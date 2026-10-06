.er_factor <- function(x) { x <- enc2utf8(as.character(x)); factor(x, levels = sort(unique(x), method = "radix")) }

.er_fit <- function(x, meta, keep, model, paired, categorical, numeric, contrast = NULL) {
    m <- droplevels(meta[keep, , drop = FALSE])
    m$group <- factor(m$group, levels = c("reference", "target"))
    terms <- character()
    if (paired) { m$donor_id <- .er_factor(m$donor_id); terms <- "donor_id" }
    cvs <- c(categorical, numeric)
    for (i in seq_along(cvs)) {
        name <- paste0("cv", i)
        m[[name]] <- if (cvs[i] %in% categorical) .er_factor(m[[cvs[i]]]) else as.numeric(m[[cvs[i]]])
        terms <- c(terms, name)
    }
    formula <- stats::reformulate(c(terms, "group"))
    factor_terms <- c(if (paired) "donor_id", if (length(categorical)) paste0("cv", which(cvs %in% categorical)), "group")
    if (!is.null(contrast)) for (term in factor_terms) {
        touched <- any(startsWith(names(contrast)[contrast != 0], term)) ||
            ("(Intercept)" %in% names(contrast) && contrast[["(Intercept)"]] != 0)
        baseline <- if (term == "group") c("reference", "target") else if (term == "donor_id")
            levels(.er_factor(meta$donor_id)) else levels(.er_factor(meta[[cvs[as.integer(sub("cv", "", term))]]]))
        if (touched && !baseline[1L] %in% as.character(m[[term]]))
            stop("Reference level of a declared contrast disappeared after deletion")
    }
    coding <- lapply(m[factor_terms], function(z) stats::contr.treatment(levels(z)))
    for (term in factor_terms) stats::contrasts(m[[term]]) <- coding[[term]]
    design <- stats::model.matrix(formula, m, contrasts.arg = coding)
    if (!is.null(contrast) && anyDuplicated(colnames(design))) stop("Declared contrast has ambiguous model-matrix names")
    if (qr(design)$rank < ncol(design)) stop("Declared design is rank deficient after deletion; no covariate was removed")
    if (nrow(design) <= ncol(design)) stop("No residual degrees of freedom")
    coefficient <- match("grouptarget", colnames(design))
    if (is.na(coefficient)) stop("No target group coefficient")
    cv <- rep(0, ncol(design)); cv[coefficient] <- 1
    if (!is.null(contrast)) {
        absent <- setdiff(names(contrast), colnames(design))
        if (any(contrast[absent] != 0)) stop("Declared contrast coefficient disappeared after deletion")
        cv[] <- 0; present <- intersect(names(contrast), colnames(design)); cv[match(present, colnames(design))] <- contrast[present]
        if (!any(cv != 0)) stop("Declared contrast is zero")
    }
    counts <- x[, keep, drop = FALSE]
    if (model == "limma") {
        fit <- limma::lmFit(counts, design)
        if (!is.null(contrast)) { fit <- limma::contrasts.fit(fit, contrasts = matrix(cv, ncol = 1L)); coefficient <- 1L }
        fit <- limma::eBayes(fit, trend = FALSE, robust = FALSE)
        score <- fit$t[, coefficient]; effect <- fit$coefficients[, coefficient]; pvalue <- fit$p.value[, coefficient]
        raw <- score
    } else if (model %in% c("edgeR", "limma-voom")) {
        y <- edgeR::DGEList(counts = counts, lib.size = m$library_total)
        y <- edgeR::calcNormFactors(y, method = "TMM")
        if (model == "edgeR") {
            y <- edgeR::estimateDisp(y, design, robust = FALSE)
            fit <- edgeR::glmQLFit(y, design, robust = FALSE)
            test <- if (is.null(contrast)) edgeR::glmQLFTest(fit, coef = coefficient)$table else edgeR::glmQLFTest(fit, contrast = cv)$table
            score <- sign(test$logFC) * sqrt(test$F); raw <- test$F
            effect <- test$logFC; pvalue <- test$PValue
        } else {
            v <- limma::voom(y, design, plot = FALSE)
            fit <- limma::lmFit(v, design)
            if (!is.null(contrast)) { fit <- limma::contrasts.fit(fit, contrasts = matrix(cv, ncol = 1L)); coefficient <- 1L }
            fit <- limma::eBayes(fit, trend = FALSE, robust = FALSE)
            score <- fit$t[, coefficient]; raw <- score
            effect <- fit$coefficients[, coefficient]; pvalue <- fit$p.value[, coefficient]
        }
    } else {
        if (any(counts > .Machine$integer.max)) stop("DESeq2 counts exceed its integer range")
        storage.mode(counts) <- "integer"
        dds <- DESeq2::DESeqDataSetFromMatrix(counts, m, formula)
        dds <- DESeq2::DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
        res <- DESeq2::results(dds, contrast = if (is.null(contrast)) c("group", "target", "reference") else cv, cooksCutoff = FALSE, independentFiltering = FALSE)
        score <- res$stat; raw <- score; effect <- res$log2FoldChange; pvalue <- res$pvalue
    }
    data.frame(feature_id = rownames(x), score = score, raw_statistic = raw,
               effect = effect, pvalue = pvalue, padj = stats::p.adjust(pvalue, "BH"), stringsAsFactors = FALSE)
}

refit_extremarank <- function(x, metadata = NULL, target, reference,
                             model = c("limma", "limma-voom", "edgeR", "DESeq2"),
                             design = c("paired", "welch"), k = NULL, budget = 1L,
                             direction = c("up", "down", "absolute"),
                             categorical = character(), numeric = character(),
                             min_total = 10, min_count = 1, min_count_samples = 0L,
                             max_refits = 1000L, assay = NULL,
                             incomplete_pairs = c("error", "drop"), features = NULL, all_features = FALSE,
                             contrast = NULL, delete_by = NULL, store_scores = c("memory", "disk", "none"),
                             score_path = NULL, checkpoint = NULL, resume = FALSE, progress = FALSE, workers = 1L,
                             max_dense_bytes = 512 * 1024^2, block_rows = 1000L, seed = 51006L) {
    model <- match.arg(model); design <- match.arg(design); direction <- match.arg(direction)
    packages <- switch(model, limma = "limma", `limma-voom` = c("limma", "edgeR"), edgeR = "edgeR", DESeq2 = "DESeq2")
    for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install the optional model package ", pkg, call. = FALSE)
    input <- .er_input(x, metadata, assay)
    samples <- .er_samples(input$x, input$metadata, target, reference, design, match.arg(incomplete_pairs))
    x <- .er_realize(samples$x, max_dense_bytes, block_rows); storage.mode(x) <- "double"; meta <- samples$metadata
    if (!is.null(contrast)) {
        if (!is.numeric(contrast) || any(!is.finite(contrast)) || is.null(names(contrast)) || anyDuplicated(names(contrast)) || !any(contrast != 0)) stop("contrast must be finite named, nonzero coefficients")
        names(contrast) <- .er_ids(names(contrast), "contrast coefficients")
    }
    if (any(!is.finite(x))) stop("model inputs must be finite", call. = FALSE)
    cvs <- c(categorical, numeric)
    if (anyDuplicated(cvs) || !all(cvs %in% names(meta)) || any(cvs %in% c("sample_id", "donor_id", "group", "library_total"))) stop("declare distinct nonreserved covariate columns", call. = FALSE)
    for (name in categorical) meta[[name]] <- .er_ids(meta[[name]], name, FALSE)
    for (name in numeric) {
        if (!is.numeric(meta[[name]]) || any(!is.finite(meta[[name]]))) stop("numeric covariates must be finite numeric columns", call. = FALSE)
    }
    # Covariates are read from their original columns before safe formula aliases
    # are installed, so an original column named cv1 cannot be overwritten.
    original_covariates <- meta[, cvs, drop = FALSE]
    if (length(cvs)) {
        aliases <- paste0("er_cov_", seq_along(cvs))
        for (i in seq_along(cvs)) meta[[aliases[i]]] <- original_covariates[[i]]
        categorical <- aliases[seq_along(categorical)]
        numeric <- if (length(numeric)) aliases[length(categorical) + seq_along(numeric)] else character()
    }
    original_features <- rownames(x)
    if (model != "limma") {
        .er_counts(x); meta$library_total <- colSums(x)
        if (any(meta$library_total <= 0) || any(!is.finite(meta$library_total))) stop("each included library needs a finite positive total", call. = FALSE)
        min_count_samples <- .er_int(min_count_samples, "min_count_samples", 0, ncol(x))
        if (length(min_total) != 1L || !is.finite(min_total) || min_total < 0 || length(min_count) != 1L || !is.finite(min_count) || min_count < 0) stop("count thresholds must be nonnegative", call. = FALSE)
        take <- rowSums(x) >= min_total & rowSums(x >= min_count) >= min_count_samples
        x <- x[take, , drop = FALSE]
    }
    if (!nrow(x)) stop("the frozen filter removed every feature", call. = FALSE)
    k <- if (is.null(k)) min(20L, nrow(x)) else .er_int(k, "k", 1, nrow(x))
    paired <- design == "paired"; original_units <- samples$units
    sample_units <- if (paired) meta$donor_id else meta$sample_id
    blocks <- original_units
    if (!is.null(delete_by)) {
        if (length(delete_by) != 1L || !delete_by %in% names(meta)) stop("delete_by must name a metadata column")
        block_values <- .er_ids(meta[[delete_by]], delete_by, FALSE)
        blocks <- vapply(original_units, function(u) {
            b <- unique(block_values[sample_units == u]); if (length(b) != 1L) stop("deletion block varies within a donor")
            b
        }, character(1))
    }
    units <- unique(blocks); sample_units <- blocks[match(sample_units, original_units)]
    weights <- cbind(vapply(units, function(u) if (paired) sum(blocks == u) else sum(sample_units == u & meta$group == target), integer(1)),
        vapply(units, function(u) if (paired) 0L else sum(sample_units == u & meta$group == reference), integer(1)))
    budget <- .er_int(budget, "budget", 0, max(0L, length(units) - 1L)); max_refits <- .er_int(max_refits, "max_refits", 0)
    plan <- cpp_block_plan(weights[,1L], weights[,2L], budget, paired)
    plan$deletion_unit <- if (is.null(delete_by)) if (paired) "whole donor" else "independent sample" else paste("whole", delete_by, "block")
    meta$group <- ifelse(meta$group == target, "target", "reference")
    fit <- function(keep) .er_fit(x, meta, keep, model, paired, categorical, numeric, contrast)
    prep <- c(samples$provenance, list(target = target, reference = reference, original_metadata_sha256 = digest::digest(input$metadata, algo = "sha256"), removed_features = setdiff(original_features, rownames(x)), filter_frozen = TRUE,
        fixed_count_filter = if (model == "limma") NULL else list(min_total = min_total, min_count = min_count, min_count_samples = min_count_samples),
        input_sha256 = digest::digest(x, algo = "sha256"), metadata = meta,
        categorical = cvs[seq_along(categorical)], numeric = if (length(numeric)) utils::tail(cvs, length(numeric)) else character(),
        library_totals_before_filter = meta$library_total, declared_contrast = contrast, design_coding = "fixed treatment; baseline reference protected for explicit contrasts", delete_by = delete_by))
    versions <- c(R = as.character(getRversion()), stats::setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), packages))
    .er_memory(nrow(x), ncol(x) * workers, max_dense_bytes)
    out <- .er_refit_engine(x, meta, units, sample_units, plan, weights, paired, model, direction, k, budget,
        max_refits, fit, prep, versions, features, all_features, store_scores, score_path, checkpoint,
        resume, progress, workers, seed)
    out$design <- design; out
}
