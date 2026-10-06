pseudobulk_extremarank <- function(x, metadata = NULL, assay = NULL, min_cells = 10L,
                                 covariates = character(), merge_technical = FALSE,
                                 missing_cell_type = c("error", "drop")) {
    input <- .er_input(x, metadata, assay, cells = TRUE)
    x <- input$x; meta <- input$metadata
    min_cells <- .er_int(min_cells, "min_cells", 1)
    required <- c("sample_id", "donor_id", "group", "cell_type")
    if (is.null(meta) || !all(required %in% names(meta))) stop("cell metadata needs cell_id, sample_id, donor_id, group and cell_type", call. = FALSE)
    if (!all(covariates %in% names(meta)) || anyDuplicated(covariates) || any(covariates %in% c(required, "cell_id", "n_cells")))
        stop("covariates must be distinct nonreserved metadata columns", call. = FALSE)
    input_cells <- ncol(x); original_metadata_hash <- digest::digest(meta, algo = "sha256")
    .er_counts(x)
    missing_cell_type <- match.arg(missing_cell_type)
    bad <- is.na(meta$cell_type) | !nzchar(as.character(meta$cell_type))
    if (any(bad) && missing_cell_type == "error") stop("missing cell-type annotations; declare missing_cell_type='drop' to exclude and record these cells", call. = FALSE)
    excluded_cells <- data.frame(cell_id = meta$cell_id[bad], reason = rep("missing cell_type", sum(bad)), stringsAsFactors = FALSE)
    if (any(bad)) { x <- x[, !bad, drop = FALSE]; meta <- meta[!bad, , drop = FALSE] }
    if (!ncol(x)) stop("no annotated cells remain", call. = FALSE)
    for (name in required) meta[[name]] <- .er_ids(meta[[name]], name, FALSE)
    for (name in covariates) if (anyNA(meta[[name]])) stop("covariates cannot contain missing values", call. = FALSE)
    # Integer codes avoid delimiter collisions in arbitrary biological IDs.
    key <- do.call(paste, c(lapply(meta[c("cell_type", "donor_id", "group")], function(z) match(z, unique(z))), sep = ":"))
    keys <- unique(key); unit <- match(key, keys); first <- match(keys, key)
    n_cells <- tabulate(unit, nbins = length(keys))
    sample_members <- lapply(seq_along(keys), function(i) unique(meta$sample_id[unit == i]))
    if (!merge_technical && any(lengths(sample_members) > 1L))
        stop("multiple sample IDs in a donor/group/cell-type unit; declare merge_technical=TRUE only for technical replicates", call. = FALSE)
    for (name in covariates) {
        ok <- vapply(seq_along(keys), function(i) length(unique(meta[[name]][unit == i])) == 1L, logical(1))
        if (!all(ok)) stop("covariate varies within a pseudobulk unit: ", name, call. = FALSE)
    }
    membership <- Matrix::sparseMatrix(i = seq_len(ncol(x)), j = unit, x = 1,
                                        dims = c(ncol(x), length(keys)))
    # With nonnegative integer summands every intermediate sum is <= its final
    # sum. Reject totals outside the exact integer range of binary64.
    counts <- x %*% membership
    .er_counts(counts)
    counts <- as.matrix(counts)
    ids <- sprintf("pb_%06d", seq_along(keys))
    colnames(counts) <- ids
    units <- meta[first, c("donor_id", "group", "cell_type", covariates), drop = FALSE]
    units$sample_id <- ids; units$n_cells <- n_cells
    units$source_samples <- vapply(sample_members, function(z) paste(z, collapse = ";"), character(1))
    units <- units[, c("sample_id", "donor_id", "group", "cell_type", "n_cells", "source_samples", covariates), drop = FALSE]
    rownames(units) <- ids
    keep <- n_cells >= min_cells
    data <- lapply(unique(units$cell_type[keep]), function(ct) {
        take <- keep & units$cell_type == ct
        list(counts = counts[, take, drop = FALSE], metadata = units[take, , drop = FALSE])
    })
    names(data) <- unique(units$cell_type[keep])
    structure(list(data = data, excluded_units = units[!keep, , drop = FALSE], excluded_cells = excluded_cells,
        provenance = list(assay = if (is.null(assay)) "explicit input count matrix" else assay,
            input_cells = input_cells, annotated_cells = ncol(x), excluded_missing_cell_type = nrow(excluded_cells),
            input_features = nrow(x), retained_units = sum(keep), missing_cell_type = missing_cell_type,
            min_cells = min_cells, merge_technical = merge_technical,
            source_samples = stats::setNames(sample_members, ids),
            counts_sha256 = digest::digest(counts, algo = "sha256"),
            metadata_sha256 = original_metadata_hash,
            arithmetic = "exact nonnegative integer sums within binary64 range",
            deletion_unit = "donor in paired analysis; independent donor/sample in Welch")),
        class = "extremarank_pseudobulk")
}

pseudobulk_experiment <- function(x, cell_type) {
    if (!inherits(x, "extremarank_pseudobulk") || length(cell_type) != 1L || !cell_type %in% names(x$data)) stop("select an existing pseudobulk cell type", call. = FALSE)
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Install SummarizedExperiment to create this object", call. = FALSE)
    z <- x$data[[cell_type]]
    SummarizedExperiment::SummarizedExperiment(assays = list(counts = z$counts), colData = z$metadata,
                                               metadata = list(extremarank = x$provenance))
}

print.extremarank_pseudobulk <- function(x, ...) {
    cat("ExtremaRank donor-level pseudobulk\n", x$provenance$input_cells, " cells; ",
        x$provenance$retained_units, " retained units; ", length(x$data), " cell types\n", sep = "")
    cat(nrow(x$excluded_units), " low-cell units excluded\n")
    invisible(x)
}
