.er_ids <- function(x, name, unique = TRUE) {
    x <- as.character(x)
    if (!length(x) || anyNA(x) || any(!nzchar(x)) || any(trimws(x) != x))
        stop(name, " must contain nonempty, unpadded identifiers", call. = FALSE)
    if (any(grepl("[\t\r\n]", x))) stop(name, " contains a tab or newline", call. = FALSE)
    if (unique && anyDuplicated(x)) stop(name, " contains duplicate identifiers", call. = FALSE)
    enc2utf8(x)
}

.er_int <- function(x, name, lower = 0, upper = .Machine$integer.max) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) || x < lower || x > upper)
        stop(name, " must be a whole number between ", lower, " and ", upper, call. = FALSE)
    as.integer(x)
}

.er_input <- function(x, metadata, assay, cells = FALSE) {
    if (inherits(x, "SingleCellExperiment") && !cells)
        stop("Aggregate cells with pseudobulk_extremarank() first; deletion units must be biological replicates", call. = FALSE)
    if (inherits(x, "SummarizedExperiment")) {
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Install SummarizedExperiment")
        if (is.null(assay) || length(assay) != 1L || !assay %in% SummarizedExperiment::assayNames(x))
            stop("Declare assay explicitly, using an existing assay name", call. = FALSE)
        if (is.null(metadata)) metadata <- as.data.frame(SummarizedExperiment::colData(x))
        x <- SummarizedExperiment::assay(x, assay)
    } else if (!is.null(assay)) stop("assay applies only to a SummarizedExperiment", call. = FALSE)
    if (is.data.frame(x)) x <- as.matrix(x)
    delayed <- inherits(x, "DelayedMatrix")
    if (delayed && !requireNamespace("DelayedArray", quietly = TRUE)) stop("Install DelayedArray")
    if (!is.matrix(x) && !inherits(x, "Matrix") && !delayed) stop("x must be a numeric matrix, Matrix or SummarizedExperiment", call. = FALSE)
    if (!is.numeric(x) && !inherits(x, "dMatrix") && !(delayed && DelayedArray::type(x) %in% c("double", "integer"))) stop("x must be numeric", call. = FALSE)
    if (!nrow(x) || !ncol(x)) stop("x must have features in rows and samples or cells in columns", call. = FALSE)
    rownames(x) <- .er_ids(rownames(x), "feature IDs")
    colnames(x) <- .er_ids(colnames(x), if (cells) "cell IDs" else "sample IDs")
    if (!is.null(metadata)) {
        metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
        key <- if (cells) "cell_id" else "sample_id"
        if (!key %in% names(metadata)) {
            rn <- rownames(metadata)
            if (is.null(rn) || identical(rn, as.character(seq_len(nrow(metadata)))))
                stop("metadata needs a ", key, " column or meaningful row names", call. = FALSE)
            metadata[[key]] <- rn
        }
        metadata[[key]] <- .er_ids(metadata[[key]], key)
        if (!setequal(metadata[[key]], colnames(x))) stop("matrix and metadata ID sets must match exactly", call. = FALSE)
        metadata <- metadata[match(colnames(x), metadata[[key]]), , drop = FALSE]
        rownames(metadata) <- colnames(x)
    }
    list(x = x, metadata = metadata)
}

.er_samples <- function(x, metadata, target, reference, design, incomplete_pairs) {
    if (is.null(metadata) || !"group" %in% names(metadata)) stop("metadata needs group", call. = FALSE)
    target <- .er_ids(target, "target"); reference <- .er_ids(reference, "reference")
    if (length(target) != 1L || length(reference) != 1L || target == reference) stop("declare distinct target and reference", call. = FALSE)
    metadata$group <- .er_ids(metadata$group, "group", FALSE)
    take <- metadata$group %in% c(target, reference)
    excluded <- metadata$sample_id[!take]
    metadata <- metadata[take, , drop = FALSE]; x <- x[, take, drop = FALSE]
    if (!all(c(target, reference) %in% metadata$group)) stop("both declared groups are required", call. = FALSE)
    identity <- "declared donor_id"
    if (!"donor_id" %in% names(metadata)) {
        if (design == "paired") stop("paired metadata needs donor_id", call. = FALSE)
        metadata$donor_id <- metadata$sample_id
        identity <- "sample_id assumed independent; biological identity not evaluated"
    }
    metadata$donor_id <- .er_ids(metadata$donor_id, "donor_id", FALSE)
    incomplete <- character()
    if (design == "paired") {
        key <- paste(match(metadata$donor_id, unique(metadata$donor_id)), match(metadata$group, c(target, reference)), sep = ":")
        if (anyDuplicated(key)) stop("paired analysis permits one sample per donor and group; explicitly aggregate technical replicates first", call. = FALSE)
        donors <- unique(metadata$donor_id)
        complete <- vapply(donors, function(d) sum(metadata$donor_id == d) == 2L, logical(1))
        incomplete <- donors[!complete]
        if (length(incomplete) && incomplete_pairs == "error") stop("incomplete pairs: ", paste(incomplete, collapse = ", "), call. = FALSE)
        take <- metadata$donor_id %in% donors[complete]
        excluded <- c(excluded, metadata$sample_id[!take])
        metadata <- metadata[take, , drop = FALSE]; x <- x[, take, drop = FALSE]
        units <- donors[complete]
        if (length(units) < 2L) stop("at least two complete donors are required", call. = FALSE)
    } else {
        if (anyDuplicated(metadata$donor_id)) stop("Welch requires independent donor IDs; repeated donors need a paired design or an external model", call. = FALSE)
        units <- metadata$sample_id
        if (any(table(factor(metadata$group, levels = c(target, reference))) < 2L)) stop("Welch needs at least two samples per group", call. = FALSE)
    }
    list(x = x, metadata = metadata, units = units, target = target, reference = reference,
         provenance = list(excluded_samples = excluded, incomplete_donors = incomplete, biological_identity = identity))
}

.er_counts <- function(x) {
    if (.er_native_counts(x)) { cpp_validate_counts(x); return(invisible(TRUE)) }
    if (inherits(x, "DelayedMatrix")) {
        step <- .er_block_rows(x, 1000L)
        for (start in seq.int(1L, nrow(x), by = step)) .er_counts(as.matrix(x[seq.int(start, min(nrow(x), start + step - 1L)), , drop = FALSE]))
        return(invisible(TRUE))
    }
    v <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
    if (any(!is.finite(v)) || any(v < 0) || any(v != floor(v)) || any(v > 2^53 - 1))
        stop("raw counts must be finite, nonnegative integers no larger than 2^53 - 1", call. = FALSE)
    invisible(TRUE)
}

.er_prepare <- function(x, metadata, target, reference, design, transform,
                        missing, incomplete_pairs, min_cpm, pseudocount, assay, max_dense_bytes = 512 * 1024^2, block_rows = 1000L) {
    input <- .er_input(x, metadata, assay)
    source_hash <- digest::digest(input$x, algo = "sha256")
    samples <- .er_samples(input$x, input$metadata, target, reference, design, incomplete_pairs)
    x <- .er_realize(samples$x, max_dense_bytes, block_rows); storage.mode(x) <- "double"
    original_features <- rownames(x)
    if (any(is.infinite(x))) stop("infinite input values are not allowed", call. = FALSE)
    if (anyNA(x)) {
        if (missing == "error" || transform == "cpm-log2") stop("missing values require explicit feature removal; CPM does not accept missing counts", call. = FALSE)
        x <- x[rowSums(is.na(x)) == 0L, , drop = FALSE]
    }
    if (length(min_cpm) != 1L || !is.finite(min_cpm) || min_cpm < 0 || length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount <= 0)
        stop("min_cpm must be nonnegative and pseudocount positive", call. = FALSE)
    library_total <- NULL
    if (transform == "cpm-log2") {
        .er_counts(x); library_total <- colSums(x)
        if (any(library_total <= 0) || any(!is.finite(library_total))) stop("each included library needs a finite positive total", call. = FALSE)
        x <- sweep(x, 2L, library_total, "/") * 1e6
        x <- x[rowSums(x >= min_cpm) >= ceiling(ncol(x)/2), , drop = FALSE]
        x <- log2(x + pseudocount)
    } else if (transform == "log2") {
        if (any(x < 0)) stop("log2 requires nonnegative input", call. = FALSE)
        x <- log2(x + pseudocount)
    }
    if (!nrow(x) || any(!is.finite(x))) stop("preparation must yield a nonempty finite feature universe", call. = FALSE)
    meta <- samples$metadata
    if (design == "paired") {
        ti <- match(samples$units, meta$donor_id[meta$group == target])
        ri <- match(samples$units, meta$donor_id[meta$group == reference])
        values <- x[, which(meta$group == target)[ti], drop = FALSE] - x[, which(meta$group == reference)[ri], drop = FALSE]
        colnames(values) <- samples$units
        labels <- integer(ncol(values))
    } else { values <- x; labels <- as.integer(meta$group == reference) }
    if (any(!is.finite(values))) stop("prepared contrasts overflowed; rescale inputs explicitly", call. = FALSE)
    list(values = values, labels = labels, units = samples$units,
         provenance = c(samples$provenance, list(target = target, reference = reference,
             transform = transform, preparation_frozen = TRUE, min_cpm = min_cpm, pseudocount = pseudocount,
             missing = missing, incomplete_pairs = incomplete_pairs, source_matrix_sha256 = source_hash,
             library_total = library_total,
             removed_features = setdiff(original_features, rownames(x)), prepared_sha256 = digest::digest(values, algo = "sha256"),
             hash_scope = "R serialized prepared binary64 matrix with dimnames", metadata = meta)))
}
