library(extremarank)
set.seed(77331)
x <- matrix(sample(20:180, 100*8, TRUE),100,dimnames=list(sprintf("g%03d",1:100),letters[1:8]))
x[1:10,1:4] <- x[1:10,1:4]+80
meta <- data.frame(sample_id=colnames(x),donor_id=colnames(x),group=rep(c("T","R"),each=4),batch=rep(c("A","B"),4),age=c(31,37,42,35,43,30,39,46))
models <- c("limma","limma-voom","edgeR","DESeq2")
for (model in models) {
    packages <- switch(model,limma="limma",`limma-voom`=c("limma","edgeR"),edgeR="edgeR",DESeq2="DESeq2")
    if (!all(vapply(packages,requireNamespace,logical(1),quietly=TRUE))) next
    z <- refit_extremarank(x,meta,"T","R",model=model,design="welch",k=5,budget=1,categorical="batch")
    stopifnot(z$status %in% c("OBSERVED_CHANGED","OBSERVED_STABLE"),nrow(z$fit_status)==9L,
              all(z$fit_status$status=="OK"),!z$enumeration_capped,length(z$scores)==9L)
    for (s in z$scenarios) {
        score <- z$scores[[s$scenario_id]]
        expected <- score$feature_id[order(-score$score,score$feature_id,method="radix")[1:5]]
        stopifnot(identical(as.character(s$topk),expected))
    }
    cap <- refit_extremarank(x,meta,"T","R",model=model,design="welch",k=5,budget=2,max_refits=1)
    stopifnot(cap$enumeration_capped,nrow(cap$fit_status)==2L)
    cat(model, "baseline, eight covariate-adjusted deletions and streamed cap passed\n")
}
if (requireNamespace("limma",quietly=TRUE)) {
    confounded <- meta; confounded$batch <- confounded$group
    z <- refit_extremarank(x,confounded,"T","R",design="welch",categorical="batch",k=5)
    stopifnot(z$status=="NOT_EVALUABLE",all(z$fit_status$status=="NOT_EVALUABLE"),grepl("rank deficient",z$fit_status$detail[1]))
    # Deletion can fail even when baseline is identifiable; retain the design.
    confounded <- meta; confounded$batch <- c("A", "A", "A", "B", "A", "A", "A", "A")
    z <- refit_extremarank(x,confounded,"T","R",design="welch",categorical="batch",k=5)
    stopifnot(z$fit_status$status[1]=="OK",any(z$fit_status$status[-1]=="NOT_EVALUABLE"))
    pair <- meta; pair$donor_id <- rep(paste0("d",1:4),2)
    z <- refit_extremarank(x,pair,"T","R",design="paired",k=5)
    stopifnot(nrow(z$fit_status)==5,all(z$fit_status$status=="OK"),all(lengths(lapply(z$scenarios[-1],`[[`,"deleted_units"))==1))
    # A covariate named cv1 or er_cov_1 cannot corrupt another declared one.
    named <- meta; named$cv1 <- named$age; named$er_cov_1 <- named$batch
    a <- refit_extremarank(x,meta,"T","R",design="welch",k=5,budget=0,categorical="batch",numeric="age")
    b <- refit_extremarank(x,named,"T","R",design="welch",k=5,budget=0,categorical="er_cov_1",numeric="cv1")
    stopifnot(identical(a$scores$baseline,b$scores$baseline))
    path <- tempfile(); write_extremarank(z,path); stopifnot(identical(z,readRDS(file.path(path,"result.rds")))); unlink(path,recursive=TRUE)
}
