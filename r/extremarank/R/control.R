.er_realize <- function(x, max_dense_bytes = 512 * 1024^2, block_rows = 1000L) {
    .er_memory(nrow(x), ncol(x), max_dense_bytes)
    if (is.matrix(x)) return(x)
    z <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
    step <- .er_block_rows(x, block_rows)
    for (start in seq.int(1L, nrow(x), by = step)) {
        take <- seq.int(start, min(nrow(x), start + step - 1L))
        z[take, ] <- as.matrix(x[take, , drop = FALSE])
    }
    z
}

.er_memory <- function(rows, columns, limit) {
    if (length(limit) != 1L || !is.numeric(limit) || !is.finite(limit) || limit <= 0)
        stop("max_dense_bytes must be finite and positive", call. = FALSE)
    bytes <- as.double(rows) * columns * 8
    if (bytes > limit) stop("Prepared dense output needs approximately ", format(bytes, scientific = FALSE),
        " bytes, exceeding max_dense_bytes; explicitly increase the limit or use fewer biological units", call. = FALSE)
    invisible(bytes)
}

.er_block_rows <- function(x, block_rows) {
    block_rows <- .er_int(block_rows, "block_rows", 1)
    max(1L, min(block_rows, floor(64 * 1024^2 / (8 * ncol(x)))))
}

.er_diagnostics <- function(x, budget, limit) {
    limit <- .er_int(limit, "max_diagnostics", 0)
    if (is.logical(x) && length(x) == 1L && !is.na(x)) x <- if (x) "auto" else "none"
    x <- match.arg(x, c("auto", "all", "none"))
    list(mode = x, limit = if (x == "none" || (x == "auto" && budget == 0L)) 0L else limit)
}

.er_plan <- function(labels, paired, budget, features, max_diagnostics = 100L) {
    p <- cpp_plan(as.integer(labels), paired, budget)
    p$deletion_unit <- if (paired) "whole donor" else "independent sample"
    p$features <- features
    p$prepared_matrix_bytes <- as.double(features) * length(labels) * 8
    p$one_score_table_numeric_bytes <- as.double(features) * 5 * 8
    p$diagnostic_limit <- max_diagnostics
    p$cost_scope <- "Matrix and numeric score storage estimates exclude IDs, arbitrary precision integers, model objects and worker copies; not peak RAM guarantees"
    p
}

preflight_extremarank <- function(x, metadata = NULL, target, reference,
    design = c("paired", "welch"), budget = 2L, assay = NULL,
    incomplete_pairs = c("error", "drop"), max_diagnostics = 100L) {
    design <- match.arg(design)
    input <- .er_input(x, metadata, assay)
    s <- .er_samples(input$x, input$metadata, target, reference, design, match.arg(incomplete_pairs))
    labels <- if (design == "paired") integer(length(s$units)) else as.integer(s$metadata$group == reference)
    budget <- .er_int(budget, "budget", 0, length(labels) - if (design == "paired") 2L else 4L)
    .er_plan(labels, design == "paired", budget, nrow(s$x), .er_int(max_diagnostics, "max_diagnostics", 0))
}
