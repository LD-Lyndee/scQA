library(SingleCellExperiment)
library(Seurat)
library(scater)
#without counts
a<-readRDS("test.rds")
if(dim(a)[2]<500){loss=0.1}else
{loss=0.3}
l<-ceiling(dim(a)[2]*0.025)
h<-dim(a)[2]-l
numcells <- nexprs(logcounts(a), byrow=TRUE)
numcells2 <- (numcells >= l&numcells<=h)
a3<-a[numcells2,]
s1<-as.Seurat(a3,counts=NULL)
s1<-FindVariableFeatures(s1,loess.span = loss)
rownames(a3)<-rownames(s1)
m<-VariableFeatures(s1)
m3<-logcounts(a3)[m,]
m3<-as.matrix(m3)
m3<-t(m3)
write.table(m3,file = "matrix.txt", sep = "\t", quote = F, row.names = F,col.names = F)
#with counts
a<-readRDS("test.rds")
if(dim(a)[2]<500)
{loss=0.1}else
{loss=0.3}
l<-ceiling(dim(a)[2]*0.025)
h<-dim(a)[2]-l

if(dim(a)[2]<=2000 ||dim(a)[2]>=20000){

lib.size=apply(counts(a), 2,sum)
lib.sf=lib.size/mean(lib.size)
nor<-t(t(counts(a))/lib.sf)
logcounts(a)<-log2(nor+1)
numcells <- nexprs(counts(a), byrow=TRUE)
numcells2 <- (numcells >= l&numcells<=h)
a3<-a[numcells2,]

}else 
{

numcells <- nexprs(counts(a), byrow=TRUE)
numcells2 <- (numcells >= l&numcells<=h)
a3<-a[numcells2,]
lib.size=apply(counts(a3), 2,sum)
lib.sf=lib.size/mean(lib.size)
nor<-t(t(counts(a3))/lib.sf)
logcounts(a3)<-log2(nor+1)
}

s1<-as.Seurat(a3,data=NULL)
s1<-FindVariableFeatures(s1,loess.span = loss)
rownames(a3)<-rownames(s1)
m<-VariableFeatures(s1)
m3<-logcounts(a3)[m,]
m3<-as.matrix(m3)
m3<-t(m3)
write.table(m3,file = "matrix.txt", sep = "\t", quote = F, row.names = F,col.names = F)