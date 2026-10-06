write_extremarank <- function(x, path) {
    if (!inherits(x, "extremarank_result")) stop("x must be an extremarank_result", call. = FALSE)
    if (length(path) != 1L || !nzchar(path)) stop("declare one output directory", call. = FALSE)
    if (file.exists(path) && (!dir.exists(path) || length(list.files(path, all.files = TRUE, no.. = TRUE))))
        stop("output must be a new or empty directory", call. = FALSE)
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) stop("cannot create output directory", call. = FALSE)
    write_tsv <- function(z, name) utils::write.table(z, file.path(path, name), sep = "\t", quote = TRUE,
                                  row.names = FALSE, na = "NA", qmethod = "double")
    exported <- x
    if (length(x$score_files)) {
        dir.create(file.path(path, "scores"))
        exported$score_files <- list(); exported$score_hashes <- list()
        for (sid in names(x$score_files)) {
            rel <- file.path("scores", paste0(sid, "_scores.rds"))
            if (!file.copy(x$score_files[[sid]], file.path(path, rel))) stop("Cannot copy disk score record: ", sid)
            exported$score_files[[sid]] <- rel
            exported$score_hashes[[sid]] <- digest::digest(file = file.path(path, rel), algo = "sha256")
        }
        exported$score_path_mode <- "relative_export"
    }
    if (!is.null(x$features)) {
        flat <- !vapply(x$features, is.list, logical(1))
        write_tsv(x$features[, flat, drop = FALSE], "features.tsv")
        witnesses <- list()
        for (z in x$audit$features) for (property in c("membership", "direction")) {
            w <- z[[property]]$witness
            if (!is.null(w)) witnesses[[length(witnesses) + 1L]] <- data.frame(feature_id = z$feature_id,
                property = property, deleted_index = w$deleted_indices, deleted_unit = w$deleted_units, stringsAsFactors = FALSE)
        }
        if (is.null(x$audit)) for (i in seq_len(nrow(x$features))) for (property in c("membership", "direction")) {
            w <- x$features[[paste0(property, "_witness")]][[i]]
            if (length(w)) witnesses[[length(witnesses) + 1L]] <- data.frame(feature_id = x$features$feature_id[i], property = property, deleted_unit = w, stringsAsFactors = FALSE)
        }
        if (length(witnesses)) write_tsv(do.call(rbind, witnesses), "witnesses.tsv")
        if (!is.null(x$diagnostics)) write_tsv(x$diagnostics$features, "leave_one_out.tsv")
    }
    if (!is.null(x$fit_status)) {
        for (sid in names(x$score_files)) {
            source <- x$score_files[[sid]]
            if (!file.exists(source)) stop("Missing disk-backed scores: ", source)
            write_tsv(readRDS(source), paste0(sid, "_scores.tsv"))
        }
        if (!is.null(x$deletion_units)) write_tsv(x$deletion_units, "deletion_units.tsv")
        if (!is.null(x$coverage)) write_tsv(as.data.frame(x$coverage, stringsAsFactors = FALSE), "coverage.tsv")
        write_tsv(x$fit_status, "fit_status.tsv")
        for (sid in names(x$scores)) if (!is.null(x$scores[[sid]])) write_tsv(x$scores[[sid]], paste0(sid, "_scores.tsv"))
    }
    saveRDS(exported, file.path(path, "result.rds"), version = 3)
    writeLines(c(paste("ExtremaRank", x$version), paste(x$method, x$status, sep = ": "), x$scope,
        "result.rds preserves exact strings, one-based witness indices, metadata and all nested records."), file.path(path, "README.txt"))
    invisible(normalizePath(path))
}

read_extremarank <- function(path) {
    root <- normalizePath(if (dir.exists(path)) path else dirname(path), winslash = "/")
    x <- readRDS(if (dir.exists(path)) file.path(path, "result.rds") else path)
    if (!inherits(x, "extremarank_result")) stop("Not an ExtremaRank result")
    if (identical(x$score_path_mode, "relative_export")) for (sid in names(x$score_files)) {
        file <- normalizePath(file.path(root, x$score_files[[sid]]), winslash = "/", mustWork = TRUE)
        if (!startsWith(file, paste0(root, "/"))) stop("Score path leaves the result directory")
        if (!identical(digest::digest(file = file, algo = "sha256"), x$score_hashes[[sid]])) stop("Disk score checksum mismatch: ", sid)
        x$score_files[[sid]] <- file
    }
    x
}
