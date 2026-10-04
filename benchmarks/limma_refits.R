# A real external-model adapter: refit limma including empirical Bayes at each
# deletion, then export the complete fixed feature universe for comparison.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3) stop("Usage: Rscript benchmarks/limma_refits.R matrix.csv[.gz] metadata.tsv output_directory")
suppressPackageStartupMessages(library(limma))
df <- read.csv(args[1], row.names=1, check.names=FALSE)
x <- as.matrix(df)
meta <- read.delim(args[2], stringsAsFactors=FALSE)
stopifnot(setequal(colnames(x),meta$sample_id), !anyDuplicated(meta$sample_id), all(is.finite(x)))
meta <- meta[match(colnames(x),meta$sample_id),]
stopifnot(identical(colnames(x),meta$sample_id), all(meta$group %in% c("ALL","AML")))
out <- args[3]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
refit <- function(keep) {
  group <- factor(meta$group[keep], levels=c("ALL","AML"))
  design <- model.matrix(~group)
  fit <- eBayes(lmFit(x[,keep,drop=FALSE], design), trend=FALSE, robust=FALSE)
  stopifnot(identical(colnames(fit$t)[2],"groupAML"))
  fit$t[,2]
}
baseline <- refit(seq_len(ncol(x)))
toprows <- function(scores, scenario) {
  ids <- rownames(x)
  orders <- list(up=order(-scores,ids,method="radix"), down=order(scores,ids,method="radix"),
                 absolute=order(-abs(scores),ids,method="radix"))
  do.call(rbind,lapply(names(orders),function(direction) {
    data.frame(scenario_id=scenario,direction=direction,rank=seq_len(20),feature_id=ids[orders[[direction]][seq_len(20)]])
  }))
}
expected <- list(toprows(baseline,"baseline"))
write.table(data.frame(feature_id=rownames(x), score=sprintf("%.17g",baseline)),
            file=file.path(out,"baseline.tsv"), sep="\t", row.names=FALSE, quote=TRUE, qmethod="double")
refits <- lapply(seq_len(ncol(x)), function(i) {
  scores <- refit(setdiff(seq_len(ncol(x)),i))
  expected[[i+1]] <<- toprows(scores,paste0("drop_",meta$sample_id[i]))
  data.frame(scenario_id=paste0("drop_",meta$sample_id[i]), feature_id=rownames(x),
             score=sprintf("%.17g",scores))
})
write.table(do.call(rbind,refits), file=file.path(out,"refits.tsv"), sep="\t", row.names=FALSE, quote=TRUE, qmethod="double")
write.table(data.frame(scenario_id=paste0("drop_",meta$sample_id),
                      deleted_units=paste0('["',meta$donor_id,'"]')),
            file=file.path(out,"manifest.tsv"), sep="\t", row.names=FALSE, quote=TRUE, qmethod="double")
writeLines(c(paste("R",getRversion()),paste("limma",packageVersion("limma")),
             "model: ~group, levels ALL then AML; coefficient groupAML",
             "empirical Bayes recomputed per refit; trend=FALSE, robust=FALSE",
             "source processed input/filtering frozen; all feasible single-donor deletions"),
           file.path(out,"versions.txt"))
write.table(do.call(rbind,expected),file=file.path(out,"expected_topk.tsv"),sep="\t",row.names=FALSE,quote=TRUE,qmethod="double")
cat("Exported", ncol(x), "full-feature delete-one limma refits.\n")
