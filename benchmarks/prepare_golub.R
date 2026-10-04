# Export the public multtest data without decimal round-trip loss.
# R is needed only to regenerate these bundled tables, not to run ExtremaRank.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript benchmarks/prepare_golub.R golub.RData output_directory")
load(args[1])
stopifnot(identical(dim(golub), c(3051L,38L)), sum(golub.cl==0)==27, sum(golub.cl==1)==11)
ids <- as.character(golub.gnames[,3])
stopifnot(!anyDuplicated(ids), all(is.finite(golub)))
samples <- sprintf("Golub_train_%02d", seq_len(ncol(golub)))
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
encoded <- matrix(sprintf("%.17g", as.numeric(golub)), nrow=nrow(golub))
df <- data.frame(feature_id=ids, encoded, check.names=FALSE)
names(df) <- c("feature_id", samples)
write.table(df, file=file.path(out,"matrix.csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE, qmethod="double")
write.table(data.frame(sample_id=samples, group=ifelse(golub.cl==1,"AML","ALL"), donor_id=samples),
            file=file.path(out,"metadata.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(feature_id=ids, source_index=golub.gnames[,1], description=golub.gnames[,2]),
            file=file.path(out,"annotation.tsv"), sep="\t", row.names=FALSE, quote=TRUE, qmethod="double")
# Export original bits for an independent Python decimal parser check. R's
# read.table parser can differ by one ulp and is not this package's parser.
writeBin(as.double(t(golub)), file.path(out,"matrix.binary64"), size=8, endian="little")
cat("Exported", nrow(golub), "features and", ncol(golub), "samples; run verify_golub_roundtrip.py.\n")
