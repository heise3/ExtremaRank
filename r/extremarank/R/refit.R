.er_fit <- function(x, meta, keep, model, paired, categorical, numeric) {
    m <- droplevels(meta[keep, , drop = FALSE])
    m$group <- factor(m$group, levels = c("reference", "target"))
    terms <- character()
    if (paired) { m$donor_id <- factor(m$donor_id); terms <- "donor_id" }
    cvs <- c(categorical, numeric)
    for (i in seq_along(cvs)) {
        name <- paste0("cv", i)
        m[[name]] <- if (cvs[i] %in% categorical) factor(m[[cvs[i]]]) else as.numeric(m[[cvs[i]]])
        terms <- c(terms, name)
    }
    formula <- stats::reformulate(c(terms, "group"))
    design <- stats::model.matrix(formula, m)
    if (qr(design)$rank < ncol(design)) stop("Declared design is rank deficient after deletion; no covariate was removed")
    if (nrow(design) <= ncol(design)) stop("No residual degrees of freedom")
    coefficient <- match("grouptarget", colnames(design))
    if (is.na(coefficient)) stop("No target group coefficient")
    counts <- x[, keep, drop = FALSE]
    if (model == "limma") {
        fit <- limma::eBayes(limma::lmFit(counts, design), trend = FALSE, robust = FALSE)
        score <- fit$t[, coefficient]; effect <- fit$coefficients[, coefficient]; pvalue <- fit$p.value[, coefficient]
        raw <- score
    } else if (model %in% c("edgeR", "limma-voom")) {
        y <- edgeR::DGEList(counts = counts, lib.size = m$library_total)
        y <- edgeR::calcNormFactors(y, method = "TMM")
        if (model == "edgeR") {
            y <- edgeR::estimateDisp(y, design, robust = FALSE)
            fit <- edgeR::glmQLFit(y, design, robust = FALSE)
            test <- edgeR::glmQLFTest(fit, coef = coefficient)$table
            score <- sign(test$logFC) * sqrt(test$F); raw <- test$F
            effect <- test$logFC; pvalue <- test$PValue
        } else {
            v <- limma::voom(y, design, plot = FALSE)
            fit <- limma::eBayes(limma::lmFit(v, design), trend = FALSE, robust = FALSE)
            score <- fit$t[, coefficient]; raw <- score
            effect <- fit$coefficients[, coefficient]; pvalue <- fit$p.value[, coefficient]
        }
    } else {
        if (any(counts > .Machine$integer.max)) stop("DESeq2 counts exceed its integer range")
        storage.mode(counts) <- "integer"
        dds <- DESeq2::DESeqDataSetFromMatrix(counts, m, formula)
        dds <- DESeq2::DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
        res <- DESeq2::results(dds, contrast = c("group", "target", "reference"), cooksCutoff = FALSE, independentFiltering = FALSE)
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
                             incomplete_pairs = c("error", "drop")) {
    model <- match.arg(model); design <- match.arg(design); direction <- match.arg(direction)
    packages <- switch(model, limma = "limma", `limma-voom` = c("limma", "edgeR"), edgeR = "edgeR", DESeq2 = "DESeq2")
    for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install the optional model package ", pkg, call. = FALSE)
    input <- .er_input(x, metadata, assay)
    samples <- .er_samples(input$x, input$metadata, target, reference, design, match.arg(incomplete_pairs))
    x <- as.matrix(samples$x); storage.mode(x) <- "double"; meta <- samples$metadata
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
    units <- samples$units; paired <- design == "paired"
    labels <- if (paired) integer(length(units)) else as.integer(meta$group == reference)
    meta$group <- ifelse(meta$group == target, "target", "reference")
    counts <- tabulate(labels + 1L, nbins = if (paired) 1L else 2L)
    budget <- .er_int(budget, "budget", 0, sum(counts - 2L)); max_refits <- .er_int(max_refits, "max_refits", 0)
    masks <- list(integer()); sizes <- integer()
    # Stream combinations; never allocate choose(n,r) masks before applying cap.
    capped <- FALSE
    visit <- function(prefix, start, need) {
        if (capped) return(invisible(NULL))
        if (need == 0L) {
            if (any(counts - tabulate(labels[prefix] + 1L, nbins = length(counts)) < 2L)) return(invisible(NULL))
            if (length(masks) - 1L >= max_refits) { capped <<- TRUE; return(invisible(NULL)) }
            masks[[length(masks) + 1L]] <<- prefix; sizes <<- c(sizes, length(prefix))
        } else if (start <= length(units) - need + 1L) for (i in seq.int(start, length(units) - need + 1L)) {
            visit(c(prefix, i), i + 1L, need - 1L)
            if (capped) break
        }
        invisible(NULL)
    }
    if (budget) for (r in seq_len(budget)) { visit(integer(), 1L, r); if (capped) break }
    fit_status <- list(); scores <- list(); scenarios <- list(); baseline_topk <- NULL; baseline_signs <- NULL
    started <- proc.time()[[3L]]
    for (i in seq_along(masks)) {
        sid <- if (i == 1L) "baseline" else paste0("delete_", i - 1L)
        deleted <- units[masks[[i]]]
        keep <- if (paired) which(!meta$donor_id %in% deleted) else which(!meta$sample_id %in% deleted)
        warnings <- messages <- character()
        result <- tryCatch(withCallingHandlers(.er_fit(x, meta, keep, model, paired, categorical, numeric),
            warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
            message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }), error = function(e) e)
        detail <- "Full fixed-universe ranking fitted"
        if (inherits(result, "error")) { status <- "NOT_EVALUABLE"; detail <- conditionMessage(result); result <- NULL }
        else if (any(!is.finite(result$score))) { status <- "NOT_EVALUABLE"; detail <- paste(sum(!is.finite(result$score)), "nonfinite scores in the frozen universe") }
        else status <- "OK"
        fit_status[[i]] <- data.frame(scenario_id = sid, status = status, detail = detail,
            warnings = paste(unique(warnings), collapse = " | "), messages = paste(unique(messages), collapse = " | "), stringsAsFactors = FALSE)
        scores[sid] <- list(result)
        topk <- NULL; changed <- NA; sign_changes <- NA_integer_
        if (status == "OK") {
            ranking_score <- switch(direction, up = result$score, down = -result$score, absolute = abs(result$score))
            topk <- result$feature_id[order(-ranking_score, enc2utf8(result$feature_id), method = "radix")[seq_len(k)]]
            if (i == 1L) { baseline_topk <- topk; baseline_signs <- sign(result$effect) }
            else if (!is.null(baseline_topk)) { changed <- !setequal(topk, baseline_topk); sign_changes <- sum(sign(result$effect) != baseline_signs) }
        }
        scenarios[[i]] <- list(scenario_id = sid, deleted_units = deleted, status = status, topk = topk, topk_changed = changed, effect_sign_changes = sign_changes)
    }
    status_table <- do.call(rbind, fit_status)
    changed <- any(vapply(scenarios[-1L], function(s) isTRUE(s$topk_changed), logical(1)))
    status <- if (status_table$status[1L] != "OK") "NOT_EVALUABLE" else if (changed) "OBSERVED_CHANGED" else if (any(status_table$status != "OK")) "PARTIALLY_EVALUATED" else "OBSERVED_STABLE"
    witnesses <- Filter(function(s) isTRUE(s$topk_changed), scenarios)
    minimum_observed <- if (length(witnesses)) min(vapply(witnesses, function(s) length(s$deleted_units), integer(1))) else NA_integer_
    structure(list(status = status, method = paste("native R", model, "refits"), model = model, design = design,
        direction = direction, top_k = k, budget = budget, fit_status = status_table, scores = scores, scenarios = scenarios,
        baseline_topk = baseline_topk, minimum_observed_change = minimum_observed,
        max_refits = max_refits, enumeration_capped = capped,
        preparation = c(samples$provenance, list(removed_features = setdiff(original_features, rownames(x)), filter_frozen = TRUE,
            fixed_count_filter = if (model == "limma") NULL else list(min_total = min_total, min_count = min_count, min_count_samples = min_count_samples),
            input_sha256 = digest::digest(x, algo = "sha256"), metadata = meta,
            categorical = cvs[seq_along(categorical)], numeric = if (length(numeric)) utils::tail(cvs, length(numeric)) else character(),
            library_totals_before_filter = meta$library_total)),
        versions = c(R = as.character(getRversion()), stats::setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), packages)),
        elapsed_seconds = proc.time()[[3L]] - started, version = as.character(utils::packageVersion("extremarank")),
        scope = "Observed sensitivity for supplied fitted scenarios; no certificate for unseen deletions, FDR or biological validity"), class = "extremarank_result")
}
