.er_rows <- function(audit) {
    do.call(rbind, lapply(audit$features, function(z) {
        out <- data.frame(feature_id = z$feature_id, baseline_rank = z$baseline_rank,
            baseline_topk = z$baseline_topk, baseline_t = z$baseline_t,
            baseline_effect = z$baseline_effect, baseline_effect_sign = z$baseline_effect_sign,
            stringsAsFactors = FALSE)
        for (property in c("membership", "direction")) {
            p <- z[[property]]
            for (name in c("status", "minimum_change_lower_bound", "minimum_change_upper_bound", "exact_minimum_change", "certified_through"))
                out[[paste(property, name, sep = "_")]] <- if (is.null(p[[name]])) NA else p[[name]]
            out[[paste0(property, "_witness")]] <- I(list(p$witness))
        }
        out
    }))
}

.er_run <- function(values, labels, units, design, k, budget, direction, features,
                    all_features, max_nodes, max_scenarios, diagnostics, provenance, max_diagnostics = 100L, search_order = "input") {
    k <- if (is.null(k)) min(20L, nrow(values)) else .er_int(k, "k", 1, nrow(values))
    max_budget <- if (design == "paired") ncol(values) - 2L else sum(table(labels) - 2L)
    budget <- .er_int(budget, "budget", 0, max_budget)
    max_nodes <- .er_int(max_nodes, "max_nodes", 1)
    max_scenarios <- .er_int(max_scenarios, "max_scenarios", 1)
    if (all_features && !is.null(features)) stop("use features or all_features, not both", call. = FALSE)
    selected <- integer()
    if (all_features) selected <- seq_len(nrow(values))
    if (!is.null(features)) {
        features <- .er_ids(features, "features")
        selected <- match(features, rownames(values))
        if (anyNA(selected)) stop("features contains IDs outside the frozen universe", call. = FALSE)
    }
    started <- proc.time()[[3L]]
    audit <- cpp_audit(values, rownames(values), units, labels, design == "paired", k, budget,
                       direction, selected, max_nodes, max_scenarios, influence_order = search_order == "influence")
    search_seconds <- proc.time()[[3L]] - started
    dc <- .er_diagnostics(diagnostics, budget, max_diagnostics)
    ds <- proc.time()[[3L]]
    loo <- if (dc$limit) cpp_diagnostics(values, rownames(values), units, labels, design == "paired", k, direction, dc$limit) else NULL
    diagnostic_seconds <- proc.time()[[3L]] - ds
    structure(list(status = audit$status, features = .er_rows(audit), audit = audit,
        diagnostics = loo, preparation = provenance, method = "native exact deletion audit",
        elapsed_seconds = proc.time()[[3L]] - started,
        version = as.character(utils::packageVersion("extremarank")),
        plan = .er_plan(labels, design == "paired", budget, nrow(values), dc$limit),
        search_order = search_order, diagnostics_control = dc,
        timing = c(search_seconds = search_seconds, diagnostic_seconds = diagnostic_seconds),
        scope = "Whole Top-K membership status; separate per-feature membership/effect-sign certificates on frozen inputs. Does not certify FDR or biological validity"),
        class = "extremarank_result")
}

extremarank <- function(x, metadata = NULL, target, reference,
                       design = c("paired", "welch"), k = NULL, budget = 2L,
                       direction = c("up", "down", "absolute"), transform = c("none", "log2", "cpm-log2"),
                       features = NULL, all_features = FALSE, max_nodes = 10000L,
                       max_scenarios = 10000L, diagnostics = TRUE, assay = NULL,
                       missing = c("error", "drop_features"), incomplete_pairs = c("error", "drop"),
                       min_cpm = 0, pseudocount = 1, max_diagnostics = 100L,
                       search_order = c("input", "influence"), max_dense_bytes = 512 * 1024^2, block_rows = 1000L) {
    design <- match.arg(design); direction <- match.arg(direction); transform <- match.arg(transform)
    missing <- match.arg(missing); incomplete_pairs <- match.arg(incomplete_pairs)
    p <- .er_prepare(x, metadata, target, reference, design, transform, missing, incomplete_pairs, min_cpm, pseudocount, assay, max_dense_bytes, block_rows)
    .er_run(p$values, p$labels, p$units, design, k, budget, direction, features, all_features,
            max_nodes, max_scenarios, diagnostics, p$provenance, max_diagnostics, match.arg(search_order))
}

extremarank_effects <- function(effects, k = NULL, budget = 2L,
                               direction = c("up", "down", "absolute"), features = NULL,
                               all_features = FALSE, max_nodes = 10000L,
                               max_scenarios = 10000L, diagnostics = TRUE, max_diagnostics = 100L,
                               search_order = c("input", "influence")) {
    if (is.data.frame(effects)) effects <- as.matrix(effects)
    if (!is.matrix(effects) || !is.numeric(effects)) stop("effects must be a numeric matrix: donors in rows, features in columns", call. = FALSE)
    rownames(effects) <- .er_ids(rownames(effects), "donor IDs")
    colnames(effects) <- .er_ids(colnames(effects), "feature IDs")
    values <- t(effects); storage.mode(values) <- "double"
    .er_run(values, integer(ncol(values)), colnames(values), "paired", k, budget, match.arg(direction),
        features, all_features, max_nodes, max_scenarios, diagnostics,
        list(transform = "precomputed within-donor effects", preparation_frozen = TRUE,
             prepared_sha256 = digest::digest(values, algo = "sha256"), hash_scope = "R serialized prepared binary64 matrix with dimnames"), max_diagnostics, match.arg(search_order))
}

extremarank_extrema <- function(x, retain, included = numeric()) {
    if (!is.numeric(x) || !is.numeric(included) || any(!is.finite(c(x, included)))) stop("x and included must be finite numeric vectors", call. = FALSE)
    retain <- .er_int(retain, "retain", 0, length(x))
    if (retain + length(included) < 2L) stop("at least two values must be retained", call. = FALSE)
    cpp_extrema(as.numeric(x), as.numeric(included), retain)
}

print.extremarank_result <- function(x, ...) {
    cat("ExtremaRank native R ", x$version, "\n", x$method, ": ", x$status, "\n", sep = "")
    if (!is.null(x$audit)) {
        cat("Top-", x$audit$top_k, "; deletion budget ", x$audit$budget, "; ", x$audit$design, " / ", x$audit$direction, "\n", sep = "")
        cat(x$audit$scenarios_checked, " scenarios checked; ", x$audit$nodes, " search nodes\n", sep = "")
    } else if (!is.null(x$fit_status)) print(table(x$fit_status$status))
    cat(x$scope, "\n")
    invisible(x)
}
