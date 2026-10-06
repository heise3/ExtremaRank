# Model fitting is delegated to established packages, not a new DE method.
# Every refit repeats normalization and variance/dispersion estimation.
a <- commandArgs(trailingOnly=TRUE)
stopifnot(length(a)==6)
out<-a[1]; method<-a[2]; nc<-as.integer(a[3]); nn<-as.integer(a[4]); paired<-as.logical(as.integer(a[5])); k<-as.integer(a[6])
packages<-switch(method,limma="limma",`limma-voom`=c("limma","edgeR"),edgeR="edgeR",DESeq2="DESeq2",stop("Unknown model"))
for(pkg in packages) if(!requireNamespace(pkg,quietly=TRUE)) stop(paste("Missing R package",pkg))
write_tsv<-function(x,name) {
 for(j in seq_along(x)) if(is.numeric(x[[j]])) x[[j]]<-sprintf("%.17g",x[[j]])
 write.table(x,file.path(out,name),sep="\t",quote=TRUE,qmethod="double",row.names=FALSE,na="NA")
}
write_tsv(data.frame(package=c("R",packages),version=c(as.character(getRversion()),vapply(packages,function(x) as.character(packageVersion(x)),character(1)))),"versions.tsv")
x<-as.matrix(read.delim(file.path(out,"model_input.tsv"),row.names=1,check.names=FALSE)); storage.mode(x)<-"double"
meta<-read.delim(file.path(out,"model_metadata.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
masks<-read.delim(file.path(out,"model_masks.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(identical(colnames(x),meta$sample_id),all(is.finite(x)))
meta$group<-factor(meta$group,levels=c("reference","target"))
terms<-character()
if(paired) terms<-c(terms,"donor_id")
if(nc+nn>0) for(i in seq_len(nc+nn)-1L) {
 name<-paste0("cv",i)
 meta[[name]]<-if(i<nc) factor(meta[[name]]) else as.numeric(meta[[name]])
 terms<-c(terms,name)
}
formula<-reformulate(c(terms,"group"))
fit_one<-function(keep) {
 m<-droplevels(meta[keep,,drop=FALSE]); m$group<-factor(m$group,levels=c("reference","target"))
 if(paired) m$donor_id<-factor(m$donor_id)
 design<-model.matrix(formula,m)
 if(qr(design)$rank<ncol(design)) stop("Declared design is rank deficient after deletion; no covariate was removed")
 if(nrow(design)<=ncol(design)) stop("No residual degrees of freedom")
 coef<-match("grouptarget",colnames(design)); stopifnot(!is.na(coef))
 counts<-x[,keep,drop=FALSE]
 if(method=="limma") {
  fit<-limma::eBayes(limma::lmFit(counts,design),trend=FALSE,robust=FALSE)
  score<-fit$t[,coef]; effect<-fit$coefficients[,coef]; pvalue<-fit$p.value[,coef]
  raw_statistic<-score
 } else if(method %in% c("edgeR","limma-voom")) {
  y<-edgeR::DGEList(counts=counts,lib.size=m$library_total)
  y<-edgeR::calcNormFactors(y,method="TMM")
  if(method=="edgeR") {
   y<-edgeR::estimateDisp(y,design,robust=FALSE)
   fit<-edgeR::glmQLFit(y,design,robust=FALSE)
   test<-edgeR::glmQLFTest(fit,coef=coef)$table
   score<-sign(test$logFC)*sqrt(test$F); effect<-test$logFC; pvalue<-test$PValue
   raw_statistic<-test$F
  } else {
   v<-limma::voom(y,design,plot=FALSE)
   fit<-limma::eBayes(limma::lmFit(v,design),trend=FALSE,robust=FALSE)
   score<-fit$t[,coef]; effect<-fit$coefficients[,coef]; pvalue<-fit$p.value[,coef]
   raw_statistic<-score
  }
 } else {
  storage.mode(counts)<-"integer"
  dds<-DESeq2::DESeqDataSetFromMatrix(counts,m,formula)
  dds<-DESeq2::DESeq(dds,quiet=TRUE,minReplicatesForReplace=Inf)
  res<-DESeq2::results(dds,contrast=c("group","target","reference"),cooksCutoff=FALSE,independentFiltering=FALSE)
  score<-res$stat; effect<-res$log2FoldChange; pvalue<-res$pvalue
  raw_statistic<-score
 }
 # Nonfinite scores invalidate this scenario's fixed-universe ranking. They
 # are never replaced by zero or removed after observing a deletion outcome.
 data.frame(feature_id=rownames(x),score=score,raw_statistic=raw_statistic,effect=effect,pvalue=pvalue,padj=p.adjust(pvalue,"BH"))
}
statuses<-list(); refits<-list(); expected<-list(); invalid_scores<-list()
for(i in 0:nrow(masks)) {
 sid<-if(i==0) "baseline" else masks$scenario_id[i]
 keep<-if(i==0) seq_len(ncol(x)) else which(as.numeric(masks[i,-1])==1)
 warnings<-character()
 result<-tryCatch(withCallingHandlers(fit_one(keep),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 if(inherits(result,"error")) {
  statuses[[length(statuses)+1L]]<-data.frame(scenario_id=sid,status="NOT_EVALUABLE",detail=conditionMessage(result),warnings=paste(unique(warnings),collapse=" | "))
  next
 }
 bad<-!is.finite(result$score)
 if(any(bad)) {
  statuses[[length(statuses)+1L]]<-data.frame(scenario_id=sid,status="NOT_EVALUABLE",detail=paste(sum(bad),"nonfinite ranking scores in the fixed universe; see nonfinite_scores.tsv"),warnings=paste(unique(warnings),collapse=" | "))
  invalid_scores[[length(invalid_scores)+1L]]<-data.frame(scenario_id=sid,result[bad,,drop=FALSE])
  next
 }
 statuses[[length(statuses)+1L]]<-data.frame(scenario_id=sid,status="OK",detail="Full fixed-universe ranking fitted",warnings=paste(unique(warnings),collapse=" | "))
 if(i==0) write_tsv(result,"baseline.tsv") else refits[[length(refits)+1L]]<-data.frame(scenario_id=sid,result)
 orders<-list(up=order(-result$score,result$feature_id,method="radix"),down=order(result$score,result$feature_id,method="radix"),absolute=order(-abs(result$score),result$feature_id,method="radix"))
 for(direction in names(orders)) expected[[length(expected)+1L]]<-data.frame(scenario_id=sid,direction=direction,rank=seq_len(k),feature_id=result$feature_id[orders[[direction]][seq_len(k)]])
}
write_tsv(do.call(rbind,statuses),"fit_status.tsv")
write_tsv(if(length(invalid_scores)) do.call(rbind,invalid_scores) else data.frame(scenario_id=character(),feature_id=character(),score=numeric(),raw_statistic=numeric(),effect=numeric(),pvalue=numeric(),padj=numeric()),"nonfinite_scores.tsv")
write_tsv(if(length(refits)) do.call(rbind,refits) else data.frame(scenario_id=character(),feature_id=character(),score=numeric(),effect=numeric(),pvalue=numeric(),padj=numeric()),"refits.tsv")
write_tsv(if(length(expected)) do.call(rbind,expected) else data.frame(scenario_id=character(),direction=character(),rank=integer(),feature_id=character()),"expected_topk.tsv")
cat(method,"completed",length(refits),"valid refits; failures retained in fit_status.tsv\n")
