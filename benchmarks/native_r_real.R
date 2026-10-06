library(extremarank)
args <- commandArgs(TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath(".")
out <- file.path(root,"results/v04"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_matrix <- function(path, csv=FALSE) {
    con<-if(grepl("\\.gz$",path)) gzfile(path,"rt") else file(path,"rt"); on.exit(close(con))
    z<-if(csv) read.csv(con,colClasses="character",check.names=FALSE) else read.delim(con,colClasses="character",check.names=FALSE)
    x<-matrix(as.numeric(as.matrix(z[,-1,drop=FALSE])),nrow=nrow(z),dimnames=list(z[[1]],names(z)[-1]));x
}
summaries<-list()
for(dataset in c("golub","kang-b","cptac","nutrimouse")) {
    base<-file.path(root,"results/v03",paste0(dataset,"-up"));paired<-dataset=="kang-b"
    filename<-if(paired) "effects.csv.gz" else if(dataset=="nutrimouse") "values.csv" else "values.csv.gz"
    values<-read_matrix(file.path(base,filename),TRUE)
    if(!paired) {
        meta<-read.delim(file.path(base,"groups.tsv"),colClasses="character",check.names=FALSE)
        prep<-jsonlite::read_json(file.path(base,"preparation.json"))
    }
    for(direction in c("up","down","absolute")) {
        k<-if(dataset=="nutrimouse") 5L else 20L
        z<-if(paired) extremarank_effects(values,k=k,budget=2,direction=direction) else
            extremarank(t(values),meta,prep$target,prep$reference,design="welch",k=k,budget=2,direction=direction)
        name<-paste(dataset,direction,sep="-"); dest<-file.path(out,name)
        if(dir.exists(dest)) {
            prior<-readRDS(file.path(dest,"result.rds"))
            stopifnot(identical(prior$preparation$prepared_sha256,z$preparation$prepared_sha256),identical(prior$features,z$features))
        } else write_extremarank(z,dest)
        # Compare properties rather than arbitrary choice of equally valid witnesses.
        old<-jsonlite::read_json(file.path(root,"results/v03",name,"feature_robustness.json"),simplifyVector=FALSE)
        stopifnot(identical(as.character(z$audit$baseline_topk),unlist(old$baseline_topk)),z$status==old$status)
        old_ids<-vapply(old$features,`[[`,character(1),"feature_id")
        for(row in z$audit$features) {
            expected<-old$features[[match(row$feature_id,old_ids)]]
            for(prop in c("membership","direction")) for(field in c("status","minimum_change_lower_bound","minimum_change_upper_bound","exact_minimum_change"))
                stopifnot(identical(row[[prop]][[field]],expected[[prop]][[field]]))
        }
        # Hex inputs permit independent replay without decimal-parser rounding.
        frozen<-if(paired) t(values) else t(values)
        write.table(data.frame(feature_id=rownames(frozen),apply(frozen,2,function(v) sprintf("%a",v)),check.names=FALSE),
                    file.path(dest,"frozen_hex.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
        jsonlite::write_json(z$audit,file.path(dest,"audit.json"),auto_unbox=TRUE,digits=NA,null="null",pretty=TRUE)
        summaries[[name]]<-list(dataset=dataset,direction=direction,status=z$status,features=ncol(values),units=nrow(values),
            certified_members=sum(z$features$membership_status=="CERTIFIED"),certified_signs=sum(z$features$direction_status=="CERTIFIED"),
            minimum_change=z$audit$topk_minimum_change,scenarios=z$audit$scenarios_checked,nodes=z$audit$nodes,
            seconds=z$elapsed_seconds,v03_property_comparison="matched; frozen R inputs audited independently",biological_scope=if(dataset=="cptac") "technical measurement runs only" else "biological sample/donor")
        cat(name,z$status,z$elapsed_seconds,"seconds\n")
    }
}

directory<-file.path(root,"data/prepared/Kang10X")
x<-as(Matrix::readMM(gzfile(file.path(directory,"matrix.mtx.gz"))),"CsparseMatrix")
rownames(x)<-read.delim(gzfile(file.path(directory,"features.tsv.gz")),header=FALSE,colClasses="character")[[1]]
colnames(x)<-read.delim(gzfile(file.path(directory,"barcodes.tsv.gz")),header=FALSE,colClasses="character")[[1]]
meta<-read.delim(gzfile(file.path(directory,"cell_metadata.tsv.gz")),colClasses="character")
pb<-pseudobulk_extremarank(x,meta,min_cells=10,missing_cell_type="drop")
stopifnot(pb$provenance$excluded_missing_cell_type==6L)
index<-read.delim(file.path(directory,"expected_index.tsv"),colClasses="character")
sum_checks<-0L
for(i in seq_len(nrow(index))) {
    z<-pb$data[[index$cell_type[i]]]; expected<-read_matrix(file.path(directory,paste0(index$filename[i],".gz")))
    unit_keys<-paste(z$metadata$donor_id,z$metadata$group,sep="|")
    actual<-z$counts[,match(colnames(expected),unit_keys),drop=FALSE]
    stopifnot(identical(unname(actual),unname(expected)))
    sum_checks<-sum_checks+length(expected)
}
saveRDS(pb,file.path(out,"kang_pseudobulk.rds"))
cat(sum_checks,"independent Matrix pseudobulk sums matched\n")

model_runs<-list()
for(dataset in c("golub","nutrimouse","kang-b")) {
    directory<-file.path(root,"data/prepared",switch(dataset,golub="Golub",nutrimouse="Nutrimouse",`kang-b`="KangPB/celltype_1"))
    path<-file.path(directory,if(dataset=="golub") "matrix.csv.gz" else if(dataset=="kang-b") "matrix.tsv.gz" else "matrix.tsv")
    counts<-read_matrix(path,dataset=="golub")
    meta<-read.delim(file.path(directory,"metadata.tsv"),colClasses="character")
    models<-if(dataset=="kang-b") c("limma-voom","edgeR","DESeq2") else "limma"
    for(model in models) {
        z<-refit_extremarank(counts,meta,
            target=switch(dataset,golub="AML",nutrimouse="ppar",`kang-b`="stim"),
            reference=switch(dataset,golub="ALL",nutrimouse="wt",`kang-b`="ctrl"),
            model=model,design=if(dataset=="kang-b") "paired" else "welch",k=if(dataset=="nutrimouse") 5 else 20,
            categorical=if(dataset=="nutrimouse") "diet" else character(),budget=1)
        name<-paste(dataset,model,sep="-");destination<-file.path(out,name)
        if(!dir.exists(destination)) write_extremarank(z,destination)
        model_runs[[name]]<-list(status=z$status,valid_refits=sum(z$fit_status$status[-1]=="OK"),
             failed_refits=sum(z$fit_status$status[-1]!="OK"),baseline=z$fit_status$status[1],versions=z$versions,seconds=z$elapsed_seconds)
        if(dataset=="kang-b" && model=="edgeR") stopifnot(z$status=="NOT_EVALUABLE") else if(dataset!="kang-b") stopifnot(all(z$fit_status$status=="OK"))
        cat(name,z$status,"\n")
        if(dataset=="kang-b") {
            # This threshold was used in benchmarks/validate_v03.py before
            # the R port; it is not chosen after observing this run's failures.
            contract<-refit_extremarank(counts,meta,"stim","ctrl",model=model,design="paired",k=20,budget=1,min_count_samples=3)
            cname<-paste0(name,"-v03-contract");destination<-file.path(out,cname)
            if(!dir.exists(destination)) write_extremarank(contract,destination)
            model_runs[[cname]]<-list(status=contract$status,valid_refits=sum(contract$fit_status$status[-1]=="OK"),
                failed_refits=sum(contract$fit_status$status[-1]!="OK"),baseline=contract$fit_status$status[1],
                fixed_min_count_samples=3,versions=contract$versions,seconds=contract$elapsed_seconds)
            if(model=="edgeR") stopifnot(contract$status=="NOT_EVALUABLE") else stopifnot(all(contract$fit_status$status=="OK"))
            cat(cname,contract$status,"\n")
        }
    }
}
jsonlite::write_json(list(native_runs=summaries,pseudobulk_sum_checks=sum_checks,input_cells=ncol(x),models=model_runs,
    scope="R native ranking sensitivity; CPTAC technical runs are not biological replication; failed Kang edgeR scores retained"),
    file.path(out,"real_validation.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
