setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/")

library(ppcor)

phe <- read.table("pheno_693.txt", head = TRUE)
counts <- fread("counts.NE2020.693.filtered.txt", data.table = F)

#zmm4 = Zm00001eb057540
#zmm15 = Zm00001eb214750

counts2 <- counts[,c("taxa","Zm00001eb057540","Zm00001eb214750")]
colnames(counts2) <- c("taxa","zmm4","zmm15")

dta1 <- merge(counts2,phe[,c("taxa","Anthesis.sp.NE")], by = "taxa")
colnames(dta1)[4] <- "Anthesis"
row.names(dta1) <- dta1[,1]
dta1 <- dta1[-1]

cor.test(dta1$zmm4,dta1$Anthesis)$estimate
cor.test(dta1$zmm15,dta1$Anthesis)$estimate

pcor(dta1)

#zcn7 = Zm00001eb293080
#zcn8 = Zm00001eb353250
counts3 <- counts[,c("taxa","Zm00001eb293080","Zm00001eb353250")]
colnames(counts3) <- c("taxa","zcn7","zcn8")

dta3 <- merge(counts3,phe[,c("taxa","Anthesis.sp.NE")], by = "taxa")
colnames(dta3)[4] <- "Anthesis"
row.names(dta3) <- dta3[,1]
dta3 <- dta3[-1]

cor.test(dta3$zcn7,dta3$Anthesis)$estimate
cor.test(dta3$zcn8,dta3$Anthesis)$estimate
cor.test(dta3$zcn8,dta3$zcn7)$estimate

pcor(dta3)


