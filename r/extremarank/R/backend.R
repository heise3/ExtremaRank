.er_native_counts <- function(x) is.matrix(x) || inherits(x, c("dgCMatrix", "dgRMatrix", "dgTMatrix", "dgeMatrix"))

extremarank_backend_info <- function() {
    x <- cpp_backend_info()
    x$default <- "native"
    x$metal_scope <- "Exact integer pseudobulk; GPU buffer budget excludes R objects and framework overhead. GPU command time is not system GPU utilization."
    x$cuda_scope <- "Future backend; NVIDIA GPU execution has not been implemented or validated"
    x
}

.er_pseudobulk_compute <- function(x, unit, n_units, backend, block_rows, max_gpu_bytes) {
    backend <- match.arg(backend, c("auto", "native", "metal", "reference"))
    if (backend == "auto") backend <- "native"
    if (backend == "metal" && !isTRUE(cpp_backend_info()$metal_available))
        stop("Metal unavailable; use backend='native'. GPU access can be blocked by a sandbox", call. = FALSE)
    if (length(max_gpu_bytes) != 1L || !is.numeric(max_gpu_bytes) || !is.finite(max_gpu_bytes) || max_gpu_bytes <= 0)
        stop("max_gpu_bytes must be finite and positive", call. = FALSE)
    started <- proc.time()[[3L]]
    run <- function(z) if (backend == "metal") cpp_pseudobulk_metal(z, unit, n_units, max_gpu_bytes) else cpp_pseudobulk_native(z, unit, n_units)
    if (backend != "reference" && .er_native_counts(x)) {
        result <- run(x)
    } else {
        result <- list(counts = matrix(0, nrow(x), n_units), backend = backend,
            gpu_buffer_bytes = 0, gpu_command_seconds = 0, gpu_commands = 0L,
            stored_values_scanned = 0)
        if (backend == "reference") {
            .er_counts(x)
            take <- which(unit > 0L)
            membership <- Matrix::sparseMatrix(i = take, j = unit[take], x = 1, dims = c(ncol(x), n_units))
        }
        step <- .er_block_rows(x, block_rows)
        for (start in seq.int(1L, nrow(x), by = step)) {
            rows <- seq.int(start, min(nrow(x), start + step - 1L))
            block <- x[rows, , drop = FALSE]
            if (inherits(block, "DelayedMatrix")) block <- as.matrix(block)
            if (backend == "reference") {
                values <- block %*% membership; .er_counts(values)
                result$counts[rows, ] <- as.matrix(values)
            } else {
                if (!.er_native_counts(block)) block <- as.matrix(block)
                part <- run(block); result$counts[rows, ] <- part$counts
                result$gpu_buffer_bytes <- max(result$gpu_buffer_bytes, part$gpu_buffer_bytes)
                result$gpu_command_seconds <- result$gpu_command_seconds + part$gpu_command_seconds
                result$gpu_commands <- result$gpu_commands + part$gpu_commands
                result$stored_values_scanned <- result$stored_values_scanned + part$stored_values_scanned
                if (!is.null(part$device)) result$device <- part$device
                if (!is.null(part$unified_memory)) result$unified_memory <- part$unified_memory
            }
        }
    }
    dimnames(result$counts) <- list(rownames(x), NULL)
    result$elapsed_seconds <- proc.time()[[3L]] - started
    result
}
