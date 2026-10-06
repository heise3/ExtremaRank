library(extremarank)
folder <- system.file("validation", package="extremarank")
source(file.path(folder,"cases.R"))
golden <- read.delim(file.path(folder,"golden.tsv"),colClasses="character",check.names=FALSE)
checks <- 0L
for (case in seq_along(cases)) {
    x <- t(cases[[case]]$rows); ids<-c("d","b","a","c"); rownames(x)<-ids; colnames(x)<-paste0("u",seq_len(ncol(x)))
    for (paired in c(TRUE,FALSE)) for (direction in c("up","down","absolute")) {
        expected <- golden[as.integer(golden$case)==case & golden$paired==if(paired) "True" else "False",,drop=FALSE]
        expected <- expected[expected$direction==direction,,drop=FALSE]
        labels <- if(paired) integer(ncol(x)) else cases[[case]]$labels
        budget<-as.integer(expected$budget[1])
        native<-extremarank:::cpp_evaluate(x,ids,labels,paired,direction,integer())
        ranks<-match(1:4,native$order)
        out<-extremarank:::cpp_audit(x,ids,colnames(x),labels,paired,2L,budget,direction,1:4,100000L,100000L)
        for(j in 1:4) {
            stopifnot(ranks[j]==as.integer(expected$baseline_rank[j]),native$scores[[j]]$sign==as.integer(expected$sign[j]),
                      native$scores[[j]]$effect==expected$effect[j])
            if(expected$denominator[j]!="0") stopifnot(native$scores[[j]]$numerator==expected$numerator[j],native$scores[[j]]$denominator==expected$denominator[j])
            for(prop in c("membership","direction")) {
                minimum<-as.integer(expected[[paste0(prop,"_minimum")]][j]); z<-out$features[[j]][[prop]]
                stopifnot(z$status==if(minimum) "REFUTED" else "CERTIFIED")
                if(minimum) stopifnot(z$exact_minimum_change==minimum) else stopifnot(is.null(z$exact_minimum_change),z$minimum_change_lower_bound==budget+1L)
                checks<-checks+1L
            }
        }
    }
}
# More than 2,000 feasible sets exercise the conditional DFS route. Shared
# additive gene offsets preserve every rank for every subset even when
# marginal feature bounds overlap. Forty queries exercise the indexed prover.
effects<-outer(seq_len(18),seq_len(40)/128,"+")
dimnames(effects)<-list(paste0("d",1:18),sprintf("g%02d",1:40))
full<-extremarank_effects(effects,k=20,budget=4,all_features=TRUE,max_nodes=100000,max_scenarios=100000,diagnostics=FALSE)
stopifnot(full$status=="CERTIFIED",all(full$features$membership_status=="CERTIFIED"),all(full$features$direction_status=="CERTIFIED"))
limited<-extremarank_effects(effects,k=20,budget=4,all_features=TRUE,max_nodes=3,max_scenarios=3,diagnostics=FALSE)
stopifnot(limited$status %in% c("CERTIFIED","UNRESOLVED"),!any(limited$features$membership_status=="REFUTED"))
x<-matrix(0,100,3,dimnames=list(paste0("s",1:100),c("signal","zero1","zero2")))
x[1:50,1]<-10+seq_len(50)/128; x[51:100,1]<-seq_len(50)/128
meta<-data.frame(sample_id=rownames(x),group=rep(c("T","R"),each=50))
stable<-extremarank(t(x),meta,"T","R",design="welch",k=1,budget=3,diagnostics=FALSE)
stopifnot(stable$status=="CERTIFIED",stable$audit$nodes<100L,stable$audit$scenarios_checked<100L)
cat(checks,"frozen independent Fraction golden property results replayed entirely in R; DFS and indexed bounds passed\n")
