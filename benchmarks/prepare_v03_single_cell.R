# Export provider raw counts/annotations and independent R pseudobulk sums.
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(Matrix))
sce<-readRDS(args[1]);out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta<-as.data.frame(colData(sce)); keep<-which(meta$multiplets=="singlet")
meta<-meta[keep,,drop=FALSE];x<-assay(sce,"counts")[,keep,drop=FALSE]
stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(x@x>=0),all(x@x==round(x@x)))
writeMM(x,file.path(out,"matrix.mtx"))
writeLines(colnames(x),file.path(out,"barcodes.tsv"))
write.table(data.frame(id=rownames(x),name=rownames(x),type="Gene Expression"),file.path(out,"features.tsv"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
cells<-data.frame(cell_id=colnames(x),sample_id=paste(meta$ind,meta$stim,sep="_"),donor_id=as.character(meta$ind),group=meta$stim,cell_type=meta$cell)
write.table(cells,file.path(out,"cell_metadata.tsv"),sep="\t",quote=TRUE,qmethod="double",row.names=FALSE)
cell_types<-sort(unique(cells$cell_type),method="radix")
index<-data.frame(cell_type=character(),filename=character(),units=integer())
for(i in seq_along(cell_types)) {
 type<-cell_types[i];cols<-which(cells$cell_type==type)
 keys<-paste(cells$donor_id[cols],cells$group[cols],sep="|")
 counts<-table(keys);units<-sort(names(counts)[counts>=10],method="radix")
 if(!length(units)) next
 expected<-vapply(units,function(key) as.numeric(Matrix::rowSums(x[,cols[keys==key],drop=FALSE])),numeric(nrow(x)))
 filename<-paste0("expected_",i,".tsv")
 write.table(data.frame(feature_id=rownames(x),expected,check.names=FALSE),file.path(out,filename),sep="\t",quote=TRUE,qmethod="double",row.names=FALSE,col.names=c("feature_id",units))
 index<-rbind(index,data.frame(cell_type=type,filename=filename,units=length(units)))
}
write.table(index,file.path(out,"expected_index.tsv"),sep="\t",quote=TRUE,qmethod="double",row.names=FALSE)
writeLines(c(paste("source dimensions",nrow(sce),ncol(sce)),paste("provider singlets",ncol(x)),paste("raw donor IDs",paste(sort(unique(meta$ind)),collapse=",")),paste("R",getRversion()),paste("Matrix",packageVersion("Matrix")),"Filter and audit contract: benchmarks/v03_data_contract.json"),file.path(out,"export_provenance.txt"))
cat(nrow(x),"genes",ncol(x),"singlet cells;",length(unique(meta$ind)),"donors\n")
