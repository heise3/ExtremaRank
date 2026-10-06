library(extremarank)
root <- normalizePath(commandArgs(TRUE)[1])
out <- file.path(root, "results/v05"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
record <- list(native = list(), models = list())
props <- function(z) lapply(z$audit$features, function(f) lapply(f[c("membership", "direction")],
    function(p) p[c("status", "minimum_change_lower_bound", "minimum_change_upper_bound", "exact_minimum_change")]))
for (dataset in c("golub", "kang-b", "cptac", "nutrimouse")) for (direction in c("up", "down", "absolute")) {
    old <- readRDS(file.path(root, "results/v04", paste(dataset, direction, sep = "-"), "result.rds"))
    base <- file.path(root, "results/v03", paste0(dataset, "-up"))
    paired <- dataset == "kang-b"
    file <- file.path(base, if (paired) "effects.csv.gz" else if (dataset == "nutrimouse") "values.csv" else "values.csv.gz")
    con <- if (grepl("gz$", file)) gzfile(file, "rt") else file(file, "rt")
    tab <- read.csv(con, colClasses = "character", check.names = FALSE); close(con)
    values <- matrix(as.numeric(as.matrix(tab[,-1,drop=FALSE])), nrow = nrow(tab), dimnames = list(tab[[1]], names(tab)[-1]))
    run <- function(order) {
        if (paired) extremarank_effects(values, k = old$audit$top_k, budget = 2, direction = direction,
            diagnostics = FALSE, search_order = order) else {
            meta <- read.delim(file.path(base, "groups.tsv"), colClasses = "character")
            prep <- jsonlite::read_json(file.path(base, "preparation.json"))
            extremarank(t(values), meta, prep$target, prep$reference, design = "welch", k = old$audit$top_k,
                budget = 2, direction = direction, diagnostics = FALSE, search_order = order)
        }
    }
    a <- run("input"); b <- run("influence")
    stopifnot(identical(props(old), props(a)), identical(props(a), props(b)), identical(a$audit$baseline_topk, old$audit$baseline_topk))
    name <- paste(dataset, direction, sep = "-")
    record$native[[name]] <- list(properties = "matched v0.4 and both search orders", input_seconds = a$elapsed_seconds,
        influence_seconds = b$elapsed_seconds, input_nodes = a$audit$nodes, influence_nodes = b$audit$nodes,
        input_scenarios = a$audit$scenarios_checked, influence_scenarios = b$audit$scenarios_checked,
        minimum_change = a$audit$topk_minimum_change)
    cat(name, "both orders matched\n")
}
for (dataset in c("golub", "nutrimouse", "kang-b")) {
    base <- file.path(root, "data/prepared", switch(dataset, golub = "Golub", nutrimouse = "Nutrimouse", `kang-b` = "KangPB/celltype_1"))
    file <- file.path(base, if (dataset == "golub") "matrix.csv.gz" else if (dataset == "kang-b") "matrix.tsv.gz" else "matrix.tsv")
    con <- if (grepl("gz$", file)) gzfile(file, "rt") else file(file, "rt")
    tab <- if (dataset == "golub") read.csv(con, colClasses = "character", check.names = FALSE) else read.delim(con, colClasses = "character", check.names = FALSE)
    close(con)
    x <- matrix(as.numeric(as.matrix(tab[,-1,drop=FALSE])), nrow = nrow(tab), dimnames = list(tab[[1]], names(tab)[-1]))
    meta <- read.delim(file.path(base, "metadata.tsv"), colClasses = "character")
    for (model in if (dataset == "kang-b") c("limma-voom", "edgeR", "DESeq2") else "limma") {
        z <- refit_extremarank(x, meta, switch(dataset, golub = "AML", nutrimouse = "ppar", `kang-b` = "stim"),
            switch(dataset, golub = "ALL", nutrimouse = "wt", `kang-b` = "ctrl"), model = model,
            design = if (dataset == "kang-b") "paired" else "welch", k = if (dataset == "nutrimouse") 5 else 20,
            categorical = if (dataset == "nutrimouse") "diet" else character(), store_scores = "none")
        name <- paste(dataset, model, sep = "-")
        old <- readRDS(file.path(root, "results/v04", name, "result.rds"))
        stopifnot(identical(z$fit_status$status, old$fit_status$status), identical(z$baseline_topk, old$baseline_topk))
        for (i in seq_along(z$scenarios)) stopifnot(identical(z$scenarios[[i]]$topk, old$scenarios[[i]]$topk))
        dir <- file.path(out, name); write_extremarank(z, dir)
        record$models[[name]] <- list(status = z$status, coverage = z$coverage,
            failed_features = z$failed_features, result_bytes = as.numeric(object.size(z)),
            v04_memory_result_bytes = as.numeric(object.size(old)), store_scores = z$store_scores,
            comparison = "same baseline, failed/valid scenarios and every observed Top-K as v0.4", seconds = z$elapsed_seconds)
        cat(name, z$status, "same observed fits\n")
    }
}
record$scope <- "Unchanged frozen real analyses; model score storage disabled by explicit setting. Failed feature IDs retained. Timings are measured on this host, not universal speed claims."
jsonlite::write_json(record, file.path(out, "real_validation.json"), pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
