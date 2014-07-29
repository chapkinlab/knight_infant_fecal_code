# read command line arguments
cmdArgs = commandArgs(trailingOnly = TRUE)
#sim_number = cmdArgs[1]
#cp=cmdArgs[1]
#sample_size = cmdArgs[2]
idf=as.numeric(cmdArgs)
#idf=1
#cat("Simulation number: ", sim_number, "\n")
#cat("Sample size: ", sample_size, "\n")

#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GenomicFeatures", "AnnotationDbi","edgeR"))
#biocLite("org.Rn.eg.db")
#library("org.Rn.eg.db")
library(CCA)
library(mvtnorm)
library(lattice)
library(gplots) 
#library(rstan)
#library(far)
library(coda)
library(edgeR) ## load the edgeR package
library(MASS)
#library(compositions)
#set_cppo("fast")  # for best running speed

### import the gene data set in r
fold<-"/home/roger/"

dat<-as.data.frame(read.csv( paste(fold,"karen-clean1.csv",sep="")))
dim(dat)
head(dat[,c(1:3,(ncol(dat)-3):ncol(dat))])

cat("## Move the diet and trt info to the first columns ... \n")

datnw <- data.frame(ID=dat[,1], Diet=dat[(ncol(dat)-1)], Treatmen=dat[,ncol(dat)],dat[,-c(1,(ncol(dat)-2):ncol(dat))])
head(datnw[,1:4])
datnw=datnw[order(datnw[,3]),]
datnw[,1:4]
order(datnw[,3])
#stop(" Arrete ici \n")

cat("### split the Diet column in two: Diet and trt \n")
colx=as.character(datnw[,2])

dit_fs=NULL
for(i in 1:nrow(datnw))
{
dit_fs=rbind(dit_fs, unlist(strsplit(colx[i]," "))[-(2:3)])    
}

datnw=data.frame(ID=datnw[,1],Oil=as.character(dit_fs[,1]),Diet= as.character(dit_fs[,2]),Treat=as.character(datnw[,3]),datnw[,-c(1:3)])

cat("### sort the data by treatment groups: oil + diet + treat .... \n")

idx=datnw[,1]
grp=NULL
for(i in 1:nrow(datnw))
{
grp=c(grp, as.integer(idx[i]/100)%%10)    
}

datnw=datnw[order(grp),]
datnw[,1:4]
table(grp)
#grp
#order(grp)
#stop("")

cat("#### Determine DE genes ..... \n")
izx=0
if(izx)
{
## filter some of the genes
counts=datnw[,-c(1:4)]
#grp <- factor(rep(1:8,each=6))
grp <- factor(rep(1:8,each=6)) 
dge <- DGEList(counts=t(counts), group=grp)

## Investigate the summaries for each sample
summary(dge$counts)
grp <- factor(rep(1:2,each=24))
keep <- rowSums(cpm(dge)) > 0
dge <- DGEList(counts=t(counts)[keep,], group=grp)
#summary(dge$counts)
#dim(dge$counts)
#max(rowSums(cpm(dge)))
#min(rowSums(cpm(dge)))
#quantile(rowSums(cpm(dge)),probs=c(.25,.75))
#dim(cpm(dge))

#dge <- DGEList(counts=t(counts), group=grp)
#cat("## summaries for all 48 animals .....\n") 
#summary(rowSums(counts))

cat("## summaries for all genes .....\n") 
summary(colSums(counts[1:24]))
summary(colSums(counts[-c(1:24)]))

cat("## summaries for all genes -- counts per Million .....\n") 
summary(rowSums(cpm(dge)[,1:24]))
summary(rowSums(cpm(dge)[,-c(1:24)]))

cat("## summaries for all genes .....\n")
summary(colSums(counts)[colSums(counts) > 100000])
lw_gen<-(1:ncol(counts))[colSums(counts) > 100000]
lw_gen <- (1:ncol(counts))[rowSums(cpm(dge)) < 5]
summary(colSums(counts))

dge <- DGEList(counts=t(counts[,-lw_gen]), group=grp)
summary(rowSums(cpm(dge)))
length(lw_gen)
#stop(" ")
## Filter by overall small counts

#smal_val<-apply(counts,2,min)
#larg_val<-apply(counts,2,max)

#ncol(counts)
#counts <- counts[,smal_val==0 & larg_val < 5]
#ncol(counts)
}
cat("End of Comment ... \n \n")
cat("## Recreate the DGElist ...\n")
counts=datnw[,-c(1:4)]
#head(counts[,1:10])
desig=counts[,c(1:4)]
grp <- factor(rep(1:8,each=6)) 
dge1 <- DGEList(counts=t(counts), group=grp)
dim(cpm(dge1))
summary( rowSums(cpm(dge1)) )

keep <- rowSums(cpm(dge1)) > 0
dge <- DGEList(counts=t(counts)[keep,], group=grp)
summary( rowSums(cpm(dge)) )
dim(cpm(dge))
#stop("")
#dge <- DGEList(counts=t(counts), group=grp)
design <- model.matrix(~0 + grp,data=dge$samples)
#design
#stop("Arrete \n")
#To estimate common dispersion:
D <- estimateGLMCommonDisp(dge, design,verbose=T)

cat("#To estimate tagwise dispersions: ...\n")
df_val<- 5
D <- estimateGLMTrendedDisp(D,design=design)
D <- estimateGLMTagwiseDisp(D,prior.df=df_val)
fold<-"/home/roger/scratch/Karen_code/EdgeR_plot/"
pdf(paste(fold,"bcv_df",df_val,".pdf",sep=""))
plotBCV(D,cex=.4)
dev.off()

#stop("Arrete \n")

## Fit GLM to the sequencing data
fit <- glmFit(dge, design=design,dispersion=D$tagwise.dispersion)

### do model Checking
names(fit)
qter_qt1<-qnorm(abs(.999 - pchisq(fit$deviance,df= fit$df.residual[1])))
qter_qt<- 1- pchisq(fit$deviance,df= fit$df.residual[1],lower=F) 
#qter_qt2<-1 - pchisq(fit$deviance,df= fit$df.residual[1])
head(fit$df.residual)
pdf(paste(fold,"qqd_df",df_val,".pdf",sep=""))
hist(qter_qt,col="cyan")
qqnorm(qter_qt1,cex=.4)
abline(a=0,b=1,lwd=4,col="red",lty=2)
names(dge)
dev.off()

#stop("here \n")

cat("#### Make my contrasts .... \n")
contras<-matrix(0,ncol=8,nrow=28) ## matrix of all the contrasts
cp=0
for(i in 1:7)
{
 for(j in (i+1):8)
 {
 cp=cp+1
contras[cp,c(i,j)]<- c(1,-1) ## creat a vector of 
 }
}
cat(" End .... \n")
res <- glmLRT(fit, contrast=contras[1,])
#head(res$fitted.values) ### mean mu
#head(res$dispersion)
write.csv(cbind(res$fitted.values,res$dispersion),paste(fold,"mean_and_disp.csv",sep=""))
write.csv(dge$samples, paste(fold,"lib_size.csv",sep=""))
names(res)
stop("here \n")
#**********************************************************************
cat("Find DE genes .... \n")
outfl=NULL
cp=NULL
cmp=0
gpr=NULL
for(i in 1:1)
{
 for(j in (i+1):2)
 {
 cmp=cmp+1
ix=6*(i-1)+1
jx=6*(j-1)+1
res <- glmLRT(fit, contrast=contras[cmp,])
trm=cbind(rep(paste(desig[ix,],collapse="_"),nrow(res$table)),rep(paste(desig[jx,],collapse="_"),nrow(res$table))) ## group identifiers
temp2=cbind(res$table,p.adjust(res$table[,"PValue"]))
outfl=rbind(outfl,temp2)
cp=rbind(cp,c(sum(res$padj[!is.na(res$padj)]<.1),sum(res$padj[!is.na(res$padj)]<.05)))
gpr=rbind(gpr,c(paste(desig[ix,],collapse="_"),paste(desig[jx,],collapse="_")))
cat(paste("Doing ---i= ",i,"..and..j= ",j,"\n",sep=""))
}
}

cp_dat=data.frame(fdr_.1=cp[,1],fdr_.05=cp[,2],A=gpr[,1],B=gpr[,2])


head(outfl)
dim(cp_dat)
cp_dat
write.csv(outfl,paste(fold,"All_comparisons_treatgrp_",df_val,".csv",sep=""))
write.csv(cp_dat,paste(fold,"Summary_All_comparisons_treatgrp_",df_val,".csv"))
stop(" Stop here .... \n")




cat("#### Test genes ... \n")




#my.contrast <- makeContrasts(c(1,0,0,0,0,0,-1,0), levels=design) ###
my.contrast <- makeContrasts(c(1,-1), levels=design) ### 
lrt <- glmLRT(fit, contrast=my.contrast)

## Finds DE genes
summary(de <- decideTestsDGE(lrt,adjust.method="fdr",p.value=0.1))

## Compare other groups

grp <- factor(rep(1:8,each=6))
dge <- DGEList(counts=t(counts), group=grp)
design <- model.matrix(~0+group,data=dge$samples)
lw_gen <- (1:ncol(counts))[rowSums(cpm(dge)) < 5]
dge <- DGEList(counts=t(counts[,-lw_gen]), group=grp)
design <- model.matrix(~0+group,data=dge$samples)
D <- estimateGLMCommonDisp(dge, design,verbose=T)
D <- estimateGLMTagwiseDisp(D)
fit <- glmFit(dge, design=design,dispersion=D$tagwise.dispersion)
my.contrast <- makeContrasts(c(1,-1,0,0,0,0), levels=design) ### 
lrt <- glmLRT(fit, contrast=my.contrast)
summary(de <- decideTestsDGE(lrt,adjust.method="fdr",p.value=0.1))
top <- topTags(lrt)
top

stop(" ")
top <- topTags(lrt)
top
cat("#### Produce plots ..... \n")
path_g <- "/home/roger/scratch/Karen_code/R_plot/"
pdf(paste(path_g,"outputsummary.pdf",sep=""))
plot(log(rowSums(cpm(dge)[-lw_gen,-c(1:24)])),type="l",col="red",lwd=1,cex=.19,ylab="counts per millions(AOM vs Saline)")
points(log(rowSums(cpm(dge)[-lw_gen,1:24])),type="l",col="blue",lwd=1,cex=.19)
legend("topright",c("AOM","Saline"),col=c("red","blue"),lwd=1,lty=1)

truehist(log(colSums(counts)),main="per gene sum of reads")
truehist(log(colSums(counts[1:24,])),main="per gene sum of reads (AOM)")
truehist(log(colSums(counts[-c(1:24),])),main="per gene sum of reads (Saline)")
plotBCV(D)
detags <- rownames(D)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()

#topTags(lrt)
ncol(counts)
ncol(counts[,-lw_gen])
stop("Arrete ici \n")
cat("#### Assign gene names from Refseq identifiers ..... \n")

refsq <- colnames(datnw)[-(1:3)]
head(refsq)
as.character(refsq[1:5])

#stop(" ")
# For the reverse map:
x <- org.Rn.egREFSEQ2EG[refsq]
# Get the RefSeq identifier that are mapped to an entrez gene ID
mapped_seqs <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_seqs])
if(length(xx) > 0) {
# Get the entrez gene for the first five Refseqs
xx[1:5]
# Get the first one
xx[[1]]
}
###

stop("Arrete \n")
