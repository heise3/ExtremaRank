# R-only benchmark. Prepare inputs once, then run each method in a fresh process.
# macOS peak RSS: /usr/bin/time -l Rscript benchmarks/mac_pseudobulk.R run ...
# Windows: run the same R script; obtain PeakWorkingSet64 from the Rscript process.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3L)
suppressPackageStartupMessages(library(Matrix))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite for benchmark reports")
root <- normalizePath(".")
sha <- function(x) digest::digest(x, algo = "sha256")
count_hash <- function(x) sha(unname(as.matrix(x)))
if (args[1L] == "prepare") {
    name <- args[2L]; path <- args[3L]
    contract <- jsonlite::read_json("benchmarks/mac_gpu_contract.json")
    if (name == "kang-real") {
        folder <- "data/prepared/Kang10X"
        x <- as(readMM(gzfile(file.path(folder,"matrix.mtx.gz"))), "CsparseMatrix")
        rownames(x) <- read.delim(gzfile(file.path(folder,"features.tsv.gz")),header=FALSE,colClasses="character")[[1L]]
        colnames(x) <- read.delim(gzfile(file.path(folder,"barcodes.tsv.gz")),header=FALSE,colClasses="character")[[1L]]
        meta <- read.delim(gzfile(file.path(folder,"cell_metadata.tsv.gz")),colClasses="character")
    } else {
        index <- which(vapply(contract$cases, function(z) z$name == name, logical(1)))
        stopifnot(length(index) == 1L)
        spec <- contract$cases[[index]]; set.seed(contract$seed)
        x <- rsparsematrix(spec$features,spec$cells,spec$density,rand.x=function(n) sample.int(20L,n,TRUE))
        dimnames(x) <- list(paste0("g",seq_len(nrow(x))),paste0("c",seq_len(ncol(x))))
        unit <- sample.int(spec$units,ncol(x),TRUE)
        meta <- data.frame(cell_id=colnames(x),sample_id=paste0("s",unit),donor_id=paste0("d",unit),
            group=ifelse(unit%%2L,"T","R"),cell_type="B")
    }
    stopifnot(setequal(meta$cell_id,colnames(x)))
    meta <- meta[match(colnames(x),meta$cell_id),,drop=FALSE]
    keys <- paste(meta$cell_type,meta$donor_id,meta$group,sep="\034")
    valid <- !is.na(meta$cell_type) & nzchar(meta$cell_type)
    unit <- match(keys,unique(keys[valid])); unit[!valid] <- 0L
    n_units <- max(unit)
    keep <- which(unit>0L)
    membership <- sparseMatrix(i=keep,j=unit[keep],x=1,dims=c(ncol(x),n_units))
    expected <- x %*% membership
    expected_named <- as.matrix(expected)
    dimnames(expected_named) <- list(rownames(x),sprintf("pb_%06d",seq_len(n_units)))
    object <- list(name=name,x=x,metadata=meta,unit=unit,n_units=n_units,
        expected_hash=count_hash(expected),expected_pipeline_hash=sha(expected_named),input_sha256=sha(list(x,meta)),
        cells=ncol(x),features=nrow(x),stored_values=length(x@x))
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    saveRDS(object,path,compress=FALSE)
    cat(name,object$cells,object$features,object$stored_values,object$expected_hash,"\n")
} else if (args[1L] == "run") {
    stopifnot(length(args) == 5L)
    path <- args[2L]; method <- args[3L]; destination <- args[4L]
    lib <- normalizePath(args[5L]); .libPaths(c(lib,.libPaths()))
    suppressPackageStartupMessages(library(extremarank))
    data <- readRDS(path); x <- data$x; meta <- data$metadata
    version <- as.character(packageVersion("extremarank"))
    legacy <- method == "released-0.5"
    stopifnot(if(legacy) version == "0.5.0" else version == "0.6.0")
    budget <- switch(method,`metal-8MiB`=8*1024^2,`metal-64MiB`=64*1024^2,`metal-whole`=1024^3,64*1024^2)
    backend <- if(grepl("^metal",method)) "metal" else "native"
    if(backend == "metal") stopifnot(isTRUE(extremarank_backend_info()$metal_available))
    matrix_method <- method == "matrix-operator"
    if(matrix_method) {
        take <- which(data$unit > 0L)
        membership <- sparseMatrix(i=take,j=data$unit[take],x=1,dims=c(ncol(x),data$n_units))
    }
    operate <- function() {
        if(matrix_method) {
            extremarank:::cpp_validate_counts(x)
            counts <- x %*% membership
            extremarank:::cpp_validate_counts(counts)
            list(counts=as.matrix(counts),backend="Matrix",gpu_buffer_bytes=0,gpu_command_seconds=0,gpu_commands=0)
        } else {
            f <- if(backend == "metal") extremarank:::cpp_pseudobulk_metal else extremarank:::cpp_pseudobulk_native
            if(backend == "metal") f(x,data$unit,data$n_units,budget) else f(x,data$unit,data$n_units)
        }
    }
    pipeline <- function() {
        if(legacy) pseudobulk_extremarank(x,meta,min_cells=10,missing_cell_type="drop") else
            pseudobulk_extremarank(x,meta,min_cells=10,missing_cell_type="drop",backend=backend,max_gpu_bytes=budget)
    }
    # Warmup includes first Metal shader compilation; separately record it.
    operators <- numeric(); cold <- NULL; metrics <- NULL
    if(!legacy) {
        gc(); start <- proc.time()[[3L]]; z <- operate(); cold <- proc.time()[[3L]]-start
        stopifnot(identical(count_hash(z$counts),data$expected_hash)); z$counts <- NULL; metrics <- z; rm(z)
        for(i in 1:3) {
            gc(); start <- proc.time()[[3L]]; z <- operate(); operators[i] <- proc.time()[[3L]]-start
            stopifnot(identical(count_hash(z$counts),data$expected_hash)); z$counts <- NULL; metrics <- z; rm(z)
        }
    }
    pipelines <- numeric(); output_hash <- NULL
    if(!matrix_method) for(i in 1:3) {
        gc(); start <- proc.time()[[3L]]; z <- pipeline(); pipelines[i] <- proc.time()[[3L]]-start
        # Includes low-cell excluded units as well as retained cell types.
        output_hash <- z$provenance$counts_sha256
        stopifnot(identical(output_hash,data$expected_pipeline_hash))
        if(!legacy) metrics <- z$provenance$computation
        rm(z)
    }
    report <- list(case=data$name,method=method,version=version,platform=R.version$platform,R=R.version.string,
        Matrix=as.character(packageVersion("Matrix")),cells=data$cells,features=data$features,units=data$n_units,
        stored_values=data$stored_values,input_sha256=data$input_sha256,expected_hash=data$expected_hash,
        expected_pipeline_hash=data$expected_pipeline_hash,
        exact_match=TRUE,operator_cold_seconds=cold,operator_seconds=operators,operator_median_seconds=if(length(operators)) median(operators) else NULL,
        pipeline_seconds=pipelines,pipeline_median_seconds=if(length(pipelines)) median(pipelines) else NULL,
        gpu_budget_bytes=if(backend=="metal") budget else 0,computation=metrics,
        peak_rss_scope="Fresh process including R, libraries, input, three warm runs, validation hashes and output. External time/OS measurement; not sum of allocations.")
    dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
    jsonlite::write_json(report,destination,pretty=TRUE,auto_unbox=TRUE,digits=NA,null="null")
    cat(data$name,method,"exact hash passed; median pipeline",report$pipeline_median_seconds,"seconds\n")
} else stop("Use prepare <case> <input.rds>, or run <input.rds> <method> <report.json> <R-library>")
