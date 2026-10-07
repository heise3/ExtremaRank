library(extremarank)
set.seed(61007)
metal <- isTRUE(extremarank_backend_info()$metal_available)
if (identical(Sys.getenv("EXTREMARANK_REQUIRE_METAL"),"true")) stopifnot(metal)
checks <- 0L
input_positions <- 0L
for (iteration in 1:30) {
    genes <- sample(3:25,1); cells <- sample(10:40,1)
    x <- matrix(sample(c(0L,0L,0L,1:100),genes*cells,TRUE),genes,
        dimnames=list(paste0("f",seq_len(genes)),paste0("c",seq_len(cells))))
    group <- sample(1:4,cells,TRUE)
    m <- data.frame(cell_id=colnames(x),donor_id=paste0("d",group),
        sample_id=paste0("s",group),group=ifelse(group%%2,"T","R"),cell_type="B")
    reference <- pseudobulk_extremarank(x,m,min_cells=1,backend="reference")
    representations <- list(x, Matrix::Matrix(x,sparse=TRUE))
    representations <- c(representations,list(as(representations[[2]],"RsparseMatrix"),as(representations[[2]],"TsparseMatrix")))
    for (z in representations) {
        a <- pseudobulk_extremarank(z,m,min_cells=1,backend="native")
        stopifnot(identical(a$data,reference$data),identical(a$provenance$counts_sha256,reference$provenance$counts_sha256))
        checks <- checks+sum(vapply(a$data,function(z) length(z$counts),integer(1)))
        input_positions <- input_positions+length(x)
        if (metal) {
            b <- pseudobulk_extremarank(z,m,min_cells=1,backend="metal",max_gpu_bytes=8*genes*length(unique(group))+4+256)
            stopifnot(identical(b$data,reference$data),b$provenance$computation$gpu_buffer_bytes<=8*genes*length(unique(group))+4+256)
            checks <- checks+sum(vapply(b$data,function(z) length(z$counts),integer(1)))
            input_positions <- input_positions+length(x)
        }
    }
    # Dropped annotation cells remain subject to input validation and retain IDs.
    m$cell_type[1]<-NA_character_
    a <- pseudobulk_extremarank(x,m,min_cells=1,missing_cell_type="drop")
    b <- pseudobulk_extremarank(x,m,min_cells=1,missing_cell_type="drop",backend="reference")
    stopifnot(identical(a$data,b$data),a$provenance$annotated_cells==cells-1L)
    x[1,1]<-NaN
    stopifnot(inherits(try(pseudobulk_extremarank(x,m,min_cells=1,missing_cell_type="drop"),silent=TRUE),"try-error"))
}
# Repeated carries, high limbs, exact maximum and rejection above the range.
x <- matrix(c(2^32-1,2^32+2,2^40+3,2^40+7,0,2^53-1),2,
    dimnames=list(c("low","high"),paste0("c",1:3)))
m <- data.frame(cell_id=colnames(x),donor_id=c("d1","d1","d2"),sample_id=c("s1","s1","s2"),group="T",cell_type="B")
a<-pseudobulk_extremarank(x,m,min_cells=1,backend="native")
b<-pseudobulk_extremarank(x,m,min_cells=1,backend="reference")
stopifnot(identical(a$data,b$data))
stopifnot(inherits(try(pseudobulk_extremarank(x,m,min_cells=1,block_rows=0),silent=TRUE),"try-error"))
if(metal) {
    b<-pseudobulk_extremarank(x,m,min_cells=1,backend="metal",max_gpu_bytes=128)
    stopifnot(identical(a$data,b$data),b$provenance$computation$gpu_commands>1L)
    # Many simultaneous low-limb overflows in one output.
    z<-matrix(rep(2^32-1,20000),1,dimnames=list("carry",paste0("c",1:20000)))
    meta<-data.frame(cell_id=colnames(z),donor_id="d",sample_id="s",group="T",cell_type="B")
    a<-pseudobulk_extremarank(z,meta,min_cells=1)
    b<-pseudobulk_extremarank(z,meta,min_cells=1,backend="metal",max_gpu_bytes=4096)
    stopifnot(identical(a$data,b$data),b$provenance$computation$gpu_commands>2)
    # Also detect an accumulator wrapping beyond 64 bits, not just 53 bits.
    z[,]<-2^53-1
    stopifnot(inherits(try(pseudobulk_extremarank(z,meta,min_cells=1,backend="metal"),silent=TRUE),"try-error"))
}
x[1,1]<-2^53-1;x[1,2]<-1
for (backend in if(metal)c("native","metal","reference") else c("native","reference"))
    stopifnot(inherits(try(pseudobulk_extremarank(x,m,min_cells=1,backend=backend),silent=TRUE),"try-error"))
if(requireNamespace("HDF5Array",quietly=TRUE)) {
    x<-matrix(1:72,6,dimnames=list(paste0("g",1:6),paste0("c",1:12)))
    m<-data.frame(cell_id=colnames(x),donor_id=rep(c("d1","d2"),each=6),sample_id=rep(c("s1","s2"),each=6),group="T",cell_type="B")
    path<-tempfile(fileext=".h5"); h<-HDF5Array::writeHDF5Array(x,filepath=path,name="counts",with.dimnames=TRUE)
    a<-pseudobulk_extremarank(x,m,min_cells=1)
    for(backend in if(metal)c("native","metal") else "native") {
        b<-pseudobulk_extremarank(h,m,min_cells=1,block_rows=1,backend=backend,max_gpu_bytes=256)
        stopifnot(identical(a$data,b$data))
    }
    unlink(path)
}
cat(checks,"aggregated count values compared across dense/CSC/CSR/triplet outputs;",input_positions,
    "input matrix positions represented; Metal executed:",metal,"\n")
