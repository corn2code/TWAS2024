#!/usr/bin/env Rscript

setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/")

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

args[1]
args[2]


chooseCRANmirror(ind=72)
source("http://zzlab.net/GAPIT/gapit_functions.txt") #install this directly in R and not here, it will give error

library("data.table")

pathout <- "out.TWAS/"


#load the new phe data
phe <- read.table("pheno_693.txt", head = TRUE) # BLUES.NE2020_corrected.csv  # BLUES.NE2021_corrected.csv
trait=which(colnames(phe) == args[1]) #args[1] #Anthesis.NE
colnames(phe[trait])

#covariates
covariates <- read.csv("sampling_693.order.csv", head = TRUE)
#colnames(covariates) #DaysToPollen

myCV <- covariates
colnames(myCV)
myCV2 <- myCV[,-c(2)]
colnames(myCV2)


#load counts data
counts <- fread("counts.NE2020.693.filtered.txt", data.table = F)
row.names(counts) <- counts$taxa

NROW(merge(counts[,1], phe, by = 1))


#use quantile method to handle outliers
Quantile<- apply(counts[,-1],2,  # 2 indicates it is for column and 1 indicates it is for row
                 function(A){min_x=as.numeric(quantile(A,0.05));
                 max_x=as.numeric(quantile(A,0.95));
                 out<-(2*(A-min_x)/(max_x-min_x));
                 out[out>2]<-2;out[out< 0]<- 0;return(out)})

#Quantile.t <-  as.data.frame(t(Quantile))
Quantile.t <- as.data.frame(Quantile)
Quantile.t$taxa <- row.names(Quantile.t)
myGD <-  Quantile.t[,c(ncol(Quantile.t),1: (ncol(Quantile.t)-1))]


myY <- cbind(phe[,1],phe[trait])
colnames(myY)[1] <- "taxa"


myGM <- read.table("gene_info_693.txt", head = TRUE)
unique(myGM$chr) #only cromosomes 

myGAPIT <- GAPIT(Y=myY[myY$taxa %in% myGD$taxa,],
                 GD=myGD,
                 GM=myGM,
                 CV=myCV2,
                 PCA.total=3,
                 model= "CMLM",
                 SNP.MAF=0,
                 file.output=F
)
#warnings()

#getting the important genes and Manhattan plots
values <- data.frame(myGAPIT$GWAS)
values$FDR <- p.adjust(values$P.value,method = "BH")


write.csv(values, paste0(pathout,"TWAS.CMLM_",colnames(phe[trait]),".csv"), row.names = F)





