library(extremarank)
set.seed(51006)
effects <- matrix(rnorm(100*30), 100, dimnames = list(paste0("d", 1:100), paste0("f", 1:30)))
z <- extremarank_effects(effects, k = 5, budget = 0)
stopifnot(is.null(z$diagnostics), z$diagnostics_control$limit == 0,
    z$plan$feasible_deletions == "0")
z <- extremarank_effects(effects, k = 5, budget = 0, diagnostics = "all", max_diagnostics = 3)
stopifnot(z$diagnostics$n_scenarios == 3, z$diagnostics$enumeration_capped,
    length(z$diagnostics$not_run_units) == 97)
for (design in c("paired", "welch")) for (direction in c("up", "down", "absolute")) {
    x <- matrix(sample(-20:20, 8*15, TRUE), 15, dimnames = list(paste0("g", 1:15), paste0("s", 1:8)))
    labels <- if (design == "paired") integer(8) else rep(0:1, each = 4)
    a <- extremarank:::cpp_audit(x, rownames(x), colnames(x), labels, design == "paired", 4L, 2L, direction, 1:15, 100000L, 100000L)
    b <- extremarank:::cpp_audit(x, rownames(x), colnames(x), labels, design == "paired", 4L, 2L, direction, 1:15, 100000L, 100000L, influence_order = TRUE)
    props <- function(y) lapply(y$features, function(f) lapply(f[c("membership", "direction")], function(p) p[c("status", "minimum_change_lower_bound", "minimum_change_upper_bound")]))
    stopifnot(identical(props(a), props(b)))
}
curve <- robustness_curve(effects[1:6, ], budgets = 0:2, k = 5, diagnostics = FALSE)
stopifnot(nrow(curve$data) == 3, inherits(summary(z), "summary.extremarank_result"),
    identical(as.data.frame(z), z$features))
plotfile <- tempfile(fileext = ".pdf"); grDevices::pdf(plotfile)
plot(z); plot(z, type = "influence"); plot(curve)
grDevices::dev.off(); stopifnot(file.info(plotfile)$size > 1000); unlink(plotfile)

# Independently count all eligible grouped deletion sets.
for (iteration in 1:20) {
    a <- sample(0:3, 5, TRUE); b <- sample(0:3, 5, TRUE)
    if (sum(a) < 2 || sum(b) < 2) next
    expected <- sapply(1:3, function(r) sum(apply(combn(1:5, r), 2, function(d) sum(a[-d]) >= 2 && sum(b[-d]) >= 2)))
    plan <- extremarank:::cpp_block_plan(a, b, 3L, FALSE)
    stopifnot(identical(as.integer(plan$by_cardinality), as.integer(expected)))
}
if (requireNamespace("limma", quietly = TRUE)) {
    x <- matrix(rnorm(30*8), 30, dimnames = list(paste0("g", 1:30), paste0("s", 1:8)))
    meta <- data.frame(sample_id = colnames(x), donor_id = colnames(x), group = rep(c("T", "R"), each = 4), site = rep(1:4, 2))
    run <- function(...) refit_extremarank(x, meta, "T", "R", design = "welch", model = "limma", k = 5, budget = 1, ...)
    zero <- run(max_refits = 0)
    stopifnot(zero$status == "NOT_EVALUATED", zero$coverage$not_run == "8", zero$coverage$evaluated == 0,
        all(zero$features$membership_status == "NOT_EVALUATED"))
    cap <- run(max_refits = 1)
    stopifnot(cap$coverage$evaluated == 1, cap$coverage$not_run == "7", cap$enumeration_capped,
        all(cap$features$membership_status %in% c("OBSERVED_CHANGED", "PARTIALLY_EVALUATED")))
    full <- run(all_features = TRUE)
    for (j in seq_len(nrow(full$features))) {
        f <- full$features[j, ]
        expected <- vapply(full$scenarios[-1L], function(s) (f$feature_id %in% s$topk) != f$baseline_topk, logical(1))
        stopifnot(f$membership_observed_changes == sum(expected))
        ranks <- vapply(full$scores, function(s) match(f$feature_id,
            s$feature_id[order(-s$score, s$feature_id, method = "radix")]), integer(1))
        selected <- vapply(full$scenarios[-1L], function(s) f$feature_id %in% s$topk, logical(1))
        stopifnot(f$best_observed_rank == min(ranks), f$worst_observed_rank == max(ranks),
            f$observed_selection_fraction == mean(selected))
        if (any(expected)) stopifnot(f$membership_status == "OBSERVED_CHANGED", length(f$membership_witness[[1]]) == 1)
    }
    tmp <- tempfile(); disk <- run(all_features = TRUE, store_scores = "disk", score_path = tmp)
    stopifnot(!length(disk$scores), length(disk$score_files) == 9, identical(full$features, disk$features))
    for (sid in names(disk$score_files)) stopifnot(identical(full$scores[[sid]], readRDS(disk$score_files[[sid]])))
    none <- run(all_features = TRUE, store_scores = "none")
    stopifnot(!length(none$scores), !length(none$score_files), identical(full$features, none$features))
    check <- tempfile(); first <- run(all_features = TRUE, checkpoint = check, max_refits = 1)
    resumed <- run(all_features = TRUE, checkpoint = check, resume = TRUE)
    stopifnot(identical(full$features, resumed$features), identical(full$scores, resumed$scores))
    bad <- try(run(all_features = TRUE, checkpoint = check, resume = TRUE, contrast = c(grouptarget = -1)), silent = TRUE)
    stopifnot(inherits(bad, "try-error"))
    blocks <- run(delete_by = "site")
    stopifnot(blocks$coverage$planned == "4", length(blocks$scenarios) == 5,
        all(vapply(blocks$scenarios[-1L], function(s) length(s$deleted_units) == 1, logical(1))))
    contrast <- run(contrast = c(grouptarget = -1))
    stopifnot(isTRUE(all.equal(contrast$scores$baseline$score, -full$scores$baseline$score)))
    old_coding <- getOption("contrasts"); options(contrasts = c("contr.sum", "contr.poly"))
    fixed <- run(all_features = TRUE); options(contrasts = old_coding)
    stopifnot(identical(full$scores, fixed$scores))
    export <- tempfile(); write_extremarank(disk, export)
    stopifnot(file.exists(file.path(export, "coverage.tsv")), file.exists(file.path(export, "baseline_scores.tsv")))
    moved <- tempfile(); stopifnot(file.rename(export, moved)); loaded <- read_extremarank(moved)
    stopifnot(identical(readRDS(loaded$score_files$baseline), full$scores$baseline))
    writeLines("changed", loaded$score_files$baseline)
    stopifnot(inherits(try(read_extremarank(moved), silent = TRUE), "try-error")); unlink(moved, recursive = TRUE)
    # Custom model: repeated observations; deletion removes every donor visit.
    visits <- meta; visits$donor_id <- rep(paste0("d", 1:4), 2)
    fitfun <- function(x, metadata) data.frame(feature_id = rownames(x), score = rowMeans(x), effect = rowMeans(x))
    custom <- refit_extremarank_custom(x, visits, fitfun, budget = 1, k = 5)
    stopifnot(custom$coverage$planned == "4", length(custom$scenarios) == 5)
    unlink(c(tmp, check, export), recursive = TRUE)
}
if (requireNamespace("DelayedArray", quietly = TRUE)) {
    counts <- matrix(1:48, 3, dimnames = list(c("A", "B", "C"), paste0("cell", 1:16)))
    cells <- data.frame(cell_id = colnames(counts), sample_id = rep(paste0("s", 1:8), each = 2),
        donor_id = rep(rep(paste0("d", 1:4), 2), each = 2), group = rep(c("T", "R"), each = 8), cell_type = "B")
    pb <- pseudobulk_extremarank(counts, cells, min_cells = 2, block_rows = 1)
    delayed <- pseudobulk_extremarank(DelayedArray::DelayedArray(counts), cells, min_cells = 2, block_rows = 1)
    stopifnot(identical(pb$data, delayed$data))
    if (requireNamespace("HDF5Array", quietly = TRUE)) {
        path <- tempfile(fileext = ".h5")
        h <- HDF5Array::writeHDF5Array(counts, filepath = path, name = "counts", with.dimnames = TRUE)
        from_disk <- pseudobulk_extremarank(h, cells, min_cells = 2, block_rows = 1)
        stopifnot(identical(pb$data, from_disk$data)); unlink(path)
        cat("HDF5-backed pseudobulk matched the independent dense reference\n")
    }
    b <- pb$data$B
    a <- extremarank(b$counts, b$metadata, "T", "R", k = 1, budget = 1, diagnostics = FALSE)
    d <- extremarank(DelayedArray::DelayedArray(b$counts), b$metadata, "T", "R", k = 1, budget = 1, diagnostics = FALSE, block_rows = 1)
    stopifnot(identical(a$features, d$features))
    bad <- try(pseudobulk_extremarank(counts, cells, min_cells = 2, max_dense_bytes = 1), silent = TRUE)
    stopifnot(inherits(bad, "try-error"))
    plan <- preflight_extremarank(b$counts, b$metadata, "T", "R", budget = 1)
    stopifnot(plan$feasible_deletions == "4")
}
cat("Coverage, candidates, diagnostics, grouped plans, contrast, disk storage, checkpoint and delayed workflows passed\n")
if (identical(Sys.getenv("EXTREMARANK_TEST_PARALLEL"), "true")) {
    x <- matrix(seq_len(40), 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:8)))
    meta <- data.frame(sample_id = colnames(x), donor_id = rep(1:4, 2))
    stochastic <- function(x, metadata) data.frame(feature_id = rownames(x), score = rowMeans(x) + stats::runif(nrow(x)), effect = rowMeans(x))
    set.seed(340); before <- .Random.seed
    a <- refit_extremarank_custom(x, meta, stochastic, workers = 1, seed = 99)
    stopifnot(identical(before, .Random.seed))
    b <- refit_extremarank_custom(x, meta, stochastic, workers = 2, seed = 99)
    stopifnot(identical(before, .Random.seed), identical(a$scores, b$scores), identical(a$features, b$features))
    cat("Serial and PSOCK stochastic refits matched; caller RNG preserved\n")
}
