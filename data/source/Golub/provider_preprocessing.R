###########################################################################
#
# Script for pre-processing the Golub et al. (1999) ALL AML training dataset
#
# Data available at: 
# 	http://www-genome.wi.mit.edu/mpr
#
###########################################################################

# Get data from Whitehead Institute website
URL<-"http://www-genome.wi.mit.edu/mpr/publications/projects/Leukemia/data_set_ALL_AML_train.txt"	 
golub.all<-read.table(URL,sep="\t",quote="",header=T,row.names=NULL,comment.char="")

# Gene names and tumor class labels
golub.gnames<-cbind(dimnames(golub.all)[[1]],as.character(golub.all[,1]),as.character(golub.all[,2]))
golub.cl<-c(rep(1,27),rep(2,11))

# Re-order columns
golub<-golub.all
golub<-golub[,1+2*(1:38)]
golub<-golub[,c(1:27,33:38,28:32)]
golub<-as.matrix(golub)

# Floor & ceiling
golub[golub<100]<-100
golub[golub>16000]<-16000

# Preliminary selection of genes
tmp1<-apply(golub,1,max)
tmp2<-apply(golub,1,min)
which1<-(1:7129)[(tmp1/tmp2)>5]
which2<-(1:7129)[(tmp1-tmp2)>500]
golub.sub<-intersect(which1,which2)
golub<-golub[golub.sub,]

# Log_10 transformation
golub<-log(golub,10)
	
# Normalization	
golub.expr<-scale(golub,T,T)
dimnames(golub.expr)<-list(NULL,NULL)

#export to multtest
golub<-golub.expr
golub.cl<-c(rep(0,27),rep(1,11))
golub.gnames<-golub.gnames[golub.sub,]
###########################################################################

