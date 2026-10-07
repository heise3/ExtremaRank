# Run in R. Binary archives are for R 4.6; other versions use the source package.
install_extremarank <- function(version = "0.6.0", lib = .libPaths()[1L], prefer_binary = TRUE) {
    if (length(version) != 1L || !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", version)) stop("Declare a release version")
    if (isNamespaceLoaded("extremarank")) stop("Restart R before updating the loaded native package")
    if (getRversion() < "4.1") stop("R >= 4.1 is required")
    packages <- c("Rcpp", "BH", "Matrix", "digest")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (requireNamespace("BH", quietly = TRUE) && utils::packageVersion("BH") < "1.75.0") missing <- union(missing, "BH")
    if (length(missing)) utils::install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
    .libPaths(c(lib, .libPaths()))
    minor <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1L]][1L], sep = ".")
    os <- Sys.info()[["sysname"]]
    arch <- if (grepl("arm|aarch", R.version$arch)) "arm64" else "x86_64"
    binary <- isTRUE(prefer_binary) && minor == "4.6" &&
        ((os == "Darwin" && grepl("R.framework", R.home(), fixed = TRUE)) || (os == "Windows" && arch == "x86_64"))
    name <- if (binary) sprintf("extremarank_%s-%s-%s-R%s.%s", version,
        if (os == "Darwin") "macos" else "windows", arch, minor,
        if (os == "Darwin") "tgz" else "zip") else paste0("extremarank_", version, ".tar.gz")
    base <- paste0("https://github.com/heise3/ExtremaRank/releases/download/v", version, "/")
    lines <- readLines(paste0(base, "SHA256SUMS"), warn = FALSE)
    fields <- strsplit(lines, "[[:space:]]+")
    expected <- vapply(Filter(function(z) length(z) == 2L && z[2L] == name, fields), `[`, character(1), 1L)
    if (length(expected) != 1L) stop("No matching verified release asset: ", name)
    tmp <- tempfile(fileext = paste0(".", tools::file_ext(name))); on.exit(unlink(tmp), add = TRUE)
    utils::download.file(paste0(base, name), tmp, mode = "wb")
    if (!identical(digest::digest(file = tmp, algo = "sha256"), expected[[1L]])) stop("Release asset checksum mismatch")
    message("Installing ", name, if (binary) " (precompiled)" else " (source; C++17 compiler required)")
    utils::install.packages(tmp, lib = lib, repos = NULL, type = if (binary) "binary" else "source")
    installed <- as.character(utils::packageVersion("extremarank", lib.loc = lib))
    if (!identical(installed, version)) stop("Installed version mismatch")
    invisible(installed)
}
