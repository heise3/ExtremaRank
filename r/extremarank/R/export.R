write_extremarank <- function(x, path) {
    if (!inherits(x, "extremarank_result")) stop("x must be an extremarank_result", call. = FALSE)
    if (length(path) != 1L || !nzchar(path)) stop("declare one output directory", call. = FALSE)
    if (file.exists(path) && (!dir.exists(path) || length(list.files(path, all.files = TRUE, no.. = TRUE))))
        stop("output must be a new or empty directory", call. = FALSE)
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) stop("cannot create output directory", call. = FALSE)
    write_tsv <- function(z, name) utils::write.table(z, file.path(path, name), sep = "\t", quote = TRUE,
                                  row.names = FALSE, na = "NA", qmethod = "double")
    saveRDS(x, file.path(path, "result.rds"), version = 3)
    if (!is.null(x$features)) {
        flat <- !vapply(x$features, is.list, logical(1))
        write_tsv(x$features[, flat, drop = FALSE], "features.tsv")
        witnesses <- list()
        for (z in x$audit$features) for (property in c("membership", "direction")) {
            w <- z[[property]]$witness
            if (!is.null(w)) witnesses[[length(witnesses) + 1L]] <- data.frame(feature_id = z$feature_id,
                property = property, deleted_index = w$deleted_indices, deleted_unit = w$deleted_units, stringsAsFactors = FALSE)
        }
        if (length(witnesses)) write_tsv(do.call(rbind, witnesses), "witnesses.tsv")
        if (!is.null(x$diagnostics)) write_tsv(x$diagnostics$features, "leave_one_out.tsv")
    }
    if (!is.null(x$fit_status)) {
        write_tsv(x$fit_status, "fit_status.tsv")
        for (sid in names(x$scores)) if (!is.null(x$scores[[sid]])) write_tsv(x$scores[[sid]], paste0(sid, "_scores.tsv"))
    }
    writeLines(c(paste("ExtremaRank", x$version), paste(x$method, x$status, sep = ": "), x$scope,
        "result.rds preserves exact strings, one-based witness indices, metadata and all nested records."), file.path(path, "README.txt"))
    invisible(normalizePath(path))
}
