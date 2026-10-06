# Provider percentages, all 21 lipids; no filtering or transformation.
args <- commandArgs(trailingOnly=TRUE)
source <- if(length(args)) args[1] else "data/source/v03/mixOmics/nutrimouse.rda"
out <- if(length(args)>1) args[2] else "data/prepared/Nutrimouse"
load(source)
x <- as.matrix(nutrimouse$lipid)
stopifnot(identical(dim(x),c(40L,21L)), all(is.finite(x)), !anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ids <- paste0("mouse_",rownames(x))
values <- matrix(sprintf("%.17g",t(x)),nrow=ncol(x))
tab <- data.frame(feature_id=colnames(x),values,check.names=FALSE)
names(tab) <- c("feature_id",ids)
write.table(tab,file.path(out,"matrix.tsv"),sep="\t",row.names=FALSE,quote=TRUE,qmethod="double")
meta <- data.frame(sample_id=ids,donor_id=ids,group=as.character(nutrimouse$genotype),diet=as.character(nutrimouse$diet))
write.table(meta,file.path(out,"metadata.tsv"),sep="\t",row.names=FALSE,quote=TRUE,qmethod="double")
writeLines(c("Nutrimouse from mixOmics 6.36.0; all 21 liver fatty-acid percentages measured by gas chromatography.",
"40 biological mice; wt versus ppar, four mice in each genotype/diet combination.",
"Values are supplied percentage compositions, not absolute concentrations. No transformation, filtering or imputation.",
"Source row IDs prefixed mouse_; raw identity NOT_EVALUATED. Frozen contract benchmarks/v03_lipid_contract.json."),file.path(out,"preparation.txt"))
