library(extremarank)
fail <- function(code) stopifnot(inherits(tryCatch(force(code), error = identity), "error"))
internal <- getFromNamespace("cpp_evaluate", "extremarank")

# Independent small-integer oracle: pairwise squared differences compute
# centered variance without the production moment subtraction.
oracle <- function(x, labels, paired, direction, deleted = integer()) {
    x <- x[, setdiff(seq_len(ncol(x)), deleted), drop = FALSE]
    labels <- labels[setdiff(seq_along(labels), deleted)]
    stat <- function(v) {
        n <- length(v); s <- sum(v)
        variance_numerator <- sum(outer(v, v, "-")^2)/2
        list(s = s, n = n, v = variance_numerator)
    }
    a <- lapply(seq_len(nrow(x)), function(j) {
        t <- stat(x[j, labels == 0L]); r <- if (paired) NULL else stat(x[j, labels == 1L])
        if (paired) c(sign = sign(t$s), num = (t$n-1)*t$s^2, den = t$v)
        else {
            d <- t$s*r$n - r$s*t$n
            c(sign = sign(d), num = d^2*(t$n-1)*(r$n-1),
              den = t$v*r$n^2*(r$n-1) + r$v*t$n^2*(t$n-1))
        }
    })
    better <- function(i, j) {
        sa <- a[[i]]["sign"]; sb <- a[[j]]["sign"]
        if (direction == "absolute") { sa <- abs(sa); sb <- abs(sb) }
        if (sa != sb) c <- sign(sa - sb)
        else if (!sa) c <- 0
        else if (!a[[i]]["den"] || !a[[j]]["den"]) c <- sa*((a[[i]]["den"] == 0) - (a[[j]]["den"] == 0))
        else c <- sa*sign(a[[i]]["num"]*a[[j]]["den"] - a[[j]]["num"]*a[[i]]["den"])
        if (direction == "down") c <- -c
        c > 0 || (c == 0 && order(enc2utf8(rownames(x)[c(i,j)]), method = "radix")[1] == 1L)
    }
    order <- seq_len(nrow(x))
    for (i in seq_len(nrow(x))) for (j in seq_len(nrow(x))) if (i < j && better(order[j], order[i])) {
        tmp <- order[i]; order[i] <- order[j]; order[j] <- tmp
    }
    list(order = order, signs = vapply(a, function(z) unname(z["sign"]), numeric(1)))
}

set.seed(60319)
checks <- 0L
for (case in seq_len(35)) {
    n <- 8L
    x <- matrix(sample(-20:20, n*4L, TRUE), 4L, dimnames = list(c("d", "b", "a", "c"), paste0("u", seq_len(n))))
    for (paired in c(TRUE, FALSE)) for (direction in c("up", "down", "absolute")) {
        labels <- if (paired) integer(n) else rep(0:1, each = 4L)
        baseline <- oracle(x, labels, paired, direction); k <- 2L; budget <- 2L
        raw <- internal(x, rownames(x), labels, paired, direction, integer())
        stopifnot(identical(as.integer(raw$order), baseline$order))
        deletions <- unlist(lapply(seq_len(budget), function(r) combn(n, r, simplify = FALSE)), recursive = FALSE)
        deletions <- Filter(function(d) paired || all(tabulate(labels[-d] + 1L, 2) >= 2L), deletions)
        minimum <- matrix(NA_integer_, 4L, 2L)
        for (d in deletions) {
            truth <- oracle(x, labels, paired, direction, d)
            changed <- cbind(seq_len(4L) %in% truth$order[1:k] != (seq_len(4L) %in% baseline$order[1:k]), truth$signs != baseline$signs)
            minimum[changed & is.na(minimum)] <- length(d)
        }
        for (limit in c(1L, 3L, 100000L)) {
            out <- getFromNamespace("cpp_audit", "extremarank")(x, rownames(x), colnames(x), labels, paired, k, budget, direction, 1:4, limit, limit, 1L)
            for (j in 1:4) for (prop in 1:2) {
                z <- out$features[[j]][[c("membership", "direction")[prop]]]; actual <- minimum[j,prop]
                stopifnot(z$minimum_change_lower_bound <= if (is.na(actual)) budget+1L else actual)
                if (z$status == "CERTIFIED") stopifnot(is.na(actual))
                if (!is.null(z$minimum_change_upper_bound)) stopifnot(!is.na(actual), actual <= z$minimum_change_upper_bound)
                if (!is.null(z$exact_minimum_change)) stopifnot(z$exact_minimum_change == actual)
                if (!is.null(z$witness)) {
                    d <- z$witness$deleted_indices; truth <- oracle(x, labels, paired, direction, d)
                    stopifnot(identical(as.character(z$witness$deleted_units), colnames(x)[d]))
                    stopifnot(if (prop == 1L) (j %in% truth$order[1:k]) != (j %in% baseline$order[1:k]) else truth$signs[j] != baseline$signs[j])
                }
                if (limit == 100000L) stopifnot(z$status == if (is.na(actual)) "CERTIFIED" else "REFUTED")
                checks <- checks + 1L
            }
        }
    }
}

# Conditional extrema are compared to every retained subset, including forced
# values, all-zero and zero-variance features.
for (i in 1:100) {
    x <- sample(-20:20, 7, TRUE); fixed <- sample(-20:20, 2, TRUE)
    e <- extremarank_extrema(x, 3, fixed)
    subsets <- combn(7, 3, simplify = FALSE)
    score <- function(idx) { z <- c(fixed, x[idx]); s <- sum(z); v <- sum(outer(z,z,"-")^2)/2; if (!s) 0 else if (!v) sign(s)*Inf else s*sqrt((length(z)-1)/v) }
    values <- vapply(subsets, score, numeric(1))
    stopifnot(isTRUE(all.equal(e$minimum$t, min(values))), isTRUE(all.equal(e$maximum$t, max(values))),
              isTRUE(all.equal(score(e$minimum$optional_indices), min(values))), isTRUE(all.equal(score(e$maximum$optional_indices), max(values))))
}

effects <- matrix(c(5,4,6,5,1,2,1,2,-2,-1,-3,-2), 4,
                  dimnames = list(c("01","02","03","04"),c("A","B","C")))
r <- extremarank_effects(effects, k = 1, budget = 1)
stopifnot(r$status == "CERTIFIED", r$features$membership_exact_minimum_change %in% NA, r$audit$index_base == 1L)
path <- tempfile(); write_extremarank(r, path)
stopifnot(identical(readRDS(file.path(path, "result.rds")), r)); fail(write_extremarank(r, path)); unlink(path, recursive = TRUE)
fail(extremarank_effects(effects, budget = 3))
fail(extremarank_extrema(c(NA, 2), 2))
fail(extremarank_extrema(c(1, 2), 1))
fail(extremarank_effects(unname(effects), budget = 1))
zero <- matrix(0,4,4,dimnames=list(paste0("u",1:4),c("é","01","基因","a")))
tie <- extremarank_effects(zero, k=2,budget=2,all_features=TRUE)
stopifnot(identical(as.character(tie$audit$baseline_topk), sort(enc2utf8(colnames(zero)), method="radix")[1:2]))

counts <- matrix(seq_len(36), 3, dimnames=list(c("g1","g2","g3"), paste0("c",1:12)))
cells <- data.frame(cell_id=colnames(counts), sample_id=rep(paste0("s",1:6),each=2),
    donor_id=rep(c("01","01","02","02","03","03"),each=2), group=rep(c("T","R"),each=2,length.out=12), cell_type="B")
pb <- pseudobulk_extremarank(Matrix::Matrix(counts, sparse=TRUE), cells, min_cells=2)
expected <- vapply(split(seq_len(12), rep(1:6,each=2)), function(i) rowSums(counts[,i,drop=FALSE]), numeric(3))
stopifnot(identical(unname(pb$data$B$counts), unname(expected)))
fail(pseudobulk_extremarank(counts/2,cells,min_cells=1))
bad <- cells; bad$sample_id[1] <- "different"; fail(pseudobulk_extremarank(counts,bad,min_cells=1))
stopifnot(pseudobulk_extremarank(counts,bad,min_cells=1,merge_technical=TRUE)$provenance$retained_units == 6L)
bad <- cells; bad$cell_id[1] <- bad$cell_id[2]; fail(pseudobulk_extremarank(counts,bad,min_cells=1))
bad <- cells; bad$cell_type[1] <- NA; fail(pseudobulk_extremarank(counts,bad,min_cells=1))
filtered <- pseudobulk_extremarank(counts,bad,min_cells=1,missing_cell_type="drop")
stopifnot(filtered$provenance$excluded_missing_cell_type==1L,filtered$excluded_cells$cell_id=="c1")
bad <- counts; bad[1,1] <- 2^53-1; bad[1,2] <- 2; fail(pseudobulk_extremarank(bad,cells,min_cells=1))
pr <- extremarank(pb$data$B$counts, pb$data$B$metadata, "T", "R", k=1,budget=1)
stopifnot(pr$preparation$preparation_frozen)
fail(extremarank(pb$data$B$counts,pb$data$B$metadata,"T","R",design="welch",budget=1))
if (requireNamespace("SummarizedExperiment",quietly=TRUE)) {
    se <- pseudobulk_experiment(pb,"B")
    fail(extremarank(se,target="T",reference="R",budget=1))
    stopifnot(identical(extremarank(se,target="T",reference="R",assay="counts",k=1,budget=1)$features, pr$features))
    if (requireNamespace("SingleCellExperiment",quietly=TRUE)) {
        sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts),colData=cells)
        fail(extremarank(sce,target="T",reference="R",assay="counts",budget=1))
        stopifnot(identical(pseudobulk_extremarank(sce,assay="counts",min_cells=2)$data,pb$data))
    }
}
cat(checks, "independent property/cap/witness checks; 100 exhaustive forced extrema; native R interface tests passed\n")
