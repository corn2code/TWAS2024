library(data.table)
library(ggplot2)
library(dplyr)
library(rMVP)
library(mice)
library(HardyWeinberg) #Try removing ‘/home/schnablelab/vladimir/R/x86_64-rocker-linux-gnu-library/4.2/00LOCK-rlang’



#install.packages("Matrix")
library(gridExtra)
library(patchwork)
library(cowplot)
library(gridGraphics)
library(ggpmisc)
library(Matrix)
#install.packages("gcookbook")
library(gcookbook) 
library(gggenes)

setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/eGWAS")


theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))


geneName = "zag6"  #args[1] #mads1 zcn26 arf34 late

index <- fread("indexGenes.csv", data.table = F)
geneid <- index$ID[which(index$Name == geneName)]

print(geneid)


print("loading data")
gwas.dat<-fread(paste0("results/",geneName,'.csv.gz'), data.table = F)
print("loading data finished")

nCHR <- length(unique(gwas.dat$CHROM))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()


for (i in sort(unique(gwas.dat$CHROM))){
  nbp[i] <- max(gwas.dat[gwas.dat$CHROM == i,]$POS)
  gwas.dat[gwas.dat$CHROM == i,"BPcum"] <- gwas.dat[gwas.dat$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.dat %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))



colnames(gwas.dat)[8] <- "gene"
gwas.dat.RNA <- gwas.dat[order(gwas.dat$gene)[1:200000],]


snps <- nrow(gwas.dat)

# to set the limite of y axis
if (max(-log10(gwas.dat.RNA$gene)) > -log10(0.05/snps) ) {
  limite = ceiling(max(-log10(gwas.dat.RNA$gene)))
} else {
  limite = ceiling(-log10(0.05/snps))
}


unique(gwas.dat.RNA$CHROM)

#to name hits in the plot 
genes.desc <- fread("/work/schnablelab/vladimir/TWAS2023.2/B73_geneModels_v5.csv", data.table = F)
hits1 <- fread("FTHits_PlusoneMichael.csv", data.table = F)
colnames(hits1)

gwas.dat.RNA2 <- gwas.dat.RNA[,c(1,9)]

hits <- merge(genes.desc, hits1[,c(1,2)] , by = 1)


# to position hit

chr.pos <- as.numeric(hits$chr[which(hits$B73_V5 == geneid)])-1


#gene.pos.x <- sum(nbp[1:chr.pos]) + hits$start[which(hits$B73_V5 == geneid)]

# to set the position of the gene
if (sum(nbp[1:chr.pos]) == 308450369 ) {
  gene.pos.x = hits$start[which(hits$B73_V5 == geneid)]
} else {
  gene.pos.x = sum(nbp[1:chr.pos]) + hits$start[which(hits$B73_V5 == geneid)]
}

#if the position of the gene is in chr 2 the formula above will place them as 1 with this we are fixing for it.
if(chr.pos == 1){
  gene.pos.x = gene.pos.x + 308450369
} else {
  gene.pos.x = gene.pos.x
}


gene.pos.y <- 2 #limite*.9
gene.id <- hits$B73_V5[which(hits$B73_V5 == geneid)]
gene.name <- hits$symbol[which(hits$B73_V5 == geneid)]

gene.name2 <- paste0(gene.name," (",gene.id,")")


hen1 <- genes.desc[which(genes.desc$B73_V5 == "Zm00001eb300700"),1:9]
hen1pos <- sum(nbp[1:6])+(hen1$start+hen1$end)/2

hen1 <- genes.desc[which(genes.desc$B73_V5 == "Zm00001eb300700"),1:9]
hen1pos <- sum(nbp[1:6])+(hen1$start+hen1$end)/2

rgd1 <- genes.desc[which(genes.desc$B73_V5 == "Zm00001eb264310"),1:9]
rgd1pos <- sum(nbp[1:5])+(rgd1$start+rgd1$end)/2


gTWAS.fdr <- ggplot() + 
  rasterise(geom_point(data=gwas.dat.RNA, aes(BPcum, -log10(gene), colour=factor(CHROM, levels = c(1:10))), size = 2),dpi=600) + 
  geom_hline(yintercept = -log10(0.05/snps), linetype=2) + 
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  geom_point(aes(y=gene.pos.y, x=gene.pos.x), fill="red", size = 5, shape = 24) +
  geom_point(aes(y=gene.pos.y, x=hen1pos), fill="black", size = 5, shape = 24) +
  geom_point(aes(y=gene.pos.y, x=rgd1pos), fill="white", size = 5, shape = 24) +
  #annotate("text", label=gene.name2, y=gene.pos.y, x=gene.pos.x, parse=T, size=5, hjust = -0.01) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](p-value))) + 
  xlab("Chromosome") + 
  scale_y_continuous(limits = c(2, limite),
                     breaks = seq(0, limite, 5))

gTWAS.fdr


asso.hits <- data.frame(matrix(NA,    # Create empty data frame
                               nrow = 1,
                               ncol = 9))
colnames(asso.hits) <- colnames(gwas.dat.RNA)

colnames(asso.hits)

#IDs for genes with FDR less than 0.05
gwas.dat.RNA$log10 <- -log10(gwas.dat.RNA$gene)

asso.hits <- gwas.dat.RNA[which(gwas.dat.RNA$log10 > -log10(0.05/snps)),]
#write.csv(asso.hits,paste0("results/associatedHits_eGWAS/",geneName,".csv"), row.names = F)

#### do the scatter plot
##
counts <- fread("../counts.NE2020.693.filtered.txt", data.table = F)
row.names(counts) <- counts$taxa

#load the new phe data
phe <- read.table("../pheno_693.txt", head = TRUE) # BLUES.NE2020_corrected.csv  # BLUES.NE2021_corrected.csv
colnames(phe)

genes <- hits1$Gene[which(hits1$symbo == geneName)]

colnames(phe)


trait=hits1$trait2plot[which(hits1$symbo == geneName)]
colnames(phe[trait])
TPM.genes <- data.frame(counts[,c(1,which(colnames(counts) 
                                          %in% genes))])

colnames(TPM.genes)[1] <- "taxa"

phenotype <- phe[,c(1,c(which(colnames(phe) %in% trait)))]
phe.TPM.genes <- merge(phenotype, TPM.genes, by="taxa")


colname2 <- gsub(".sp", "",  colnames(phe.TPM.genes)[2])

colnames(phe.TPM.genes)[2] <- colname2


##
popdata <- read.csv("/work/schnablelab/vladimir/TWAS2023.2/MarcinTPJ2023_TableS1.csv", check.names = F)
colnames(popdata)[1] <- "taxa"
unique(popdata$Group)


popdata$Group <- gsub("Mixed", "Others" , popdata$Group)
popdata$Group <- gsub("Broad origin-public", "Others" , popdata$Group)
popdata$Group <- gsub("Landrace", "Others" , popdata$Group)



groups <- c("SS","IDT", "NSS", "Sweet corn", "Popcorn", "Tropical", "Others")
groups.col <- c("darkgoldenrod1", "indianred", "deeppink", "cornflowerblue", "cornsilk2", "darkblue", "white")

unique(popdata$Group)

phe.TPM.genes.pop <- merge(popdata[,c(1,3)], phe.TPM.genes, by="taxa")


colnames(phe.TPM.genes.pop)[3:4] <- c("phenotype","expression")


colnames(phe.TPM.genes)[2]
paste0(gene.name," (TPM)")
paste0("Days to ", colnames(phe.TPM.genes)[2])

my_y_title <- expression(paste(italic("zag6"), " (TPM)"))

relation.plot <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=phenotype, y=expression, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=phenotype, y=expression,fill=Group), size = 5, shape = 21),dpi=600) +
  xlim(55,85)+
  ylab(my_y_title) + 
  xlab("Days to Silking (MI)") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none") + 
  stat_poly_line(data= phe.TPM.genes.pop, aes(x=phenotype, y=phe.TPM.genes.pop[,4]), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop, aes(x=phenotype, y=phe.TPM.genes.pop[,4]),size = 5)
relation.plot


chrMostSig <- gwas.dat.RNA$CHROM[which(gwas.dat.RNA$log10 == max(gwas.dat.RNA$log10))]
chrMostSig <- chrMostSig[1]
genodata <- paste0("1515.SNPs/mvpData/mvp.",chrMostSig,".vcf.geno.desc")
genomap <- paste0("1515.SNPs/mvpData/mvp.",chrMostSig,".vcf.geno.map")
ind <- paste0("1515.SNPs/mvpData/mvp.",chrMostSig,".vcf.geno.ind")


individuals <- fread(ind, data.table = F)
genomap <- fread(genomap)
genodata <- attach.big.matrix(genodata)

#NROW(asso.hits) == 0 modify hre if we dont want to plot not hits

gene2 <- gwas.dat.RNA$SNP[which(gwas.dat.RNA$log10 == max(gwas.dat.RNA$log10))]
gene2 <- gene2[1] # for arf34 the second most significant variation has the same outcome but is just a snp (CC and TT)

POS <- genomap$POS[which(genomap$SNP == gene2)]
REF <- genomap$REF[which(genomap$SNP == gene2)]
ALT <- genomap$ALT[which(genomap$SNP == gene2)]

### SNP effect plot as boxplot
individuals$gene2 <- genodata[which(genomap$SNP==gene2),] # from the original file NOT the simplyfied version
individuals$gene3[individuals$gene2==0] <- REF #ref
individuals$gene3[individuals$gene2==2] <- ALT #alt
individuals$gene3[individuals$gene2==1] <- NA #het
individuals$gene3 <- factor(individuals$gene3, levels = c(REF,ALT))


individuals$gene4[individuals$gene2==0] <- paste0(REF,REF) #ref
individuals$gene4[individuals$gene2==2] <- paste0(ALT,ALT) #alt
individuals$gene4[individuals$gene2==1] <- paste0(REF,ALT) #het

colnames(individuals)
colnames(individuals)[1] <- "Taxa"
individuals2 <- individuals[c("Taxa","gene3", "gene4")]

#counts <- fread("../counts.NE2020.693.filtered.txt", data.table = F)
counts1 <- counts[c("taxa",geneid)]


dta <- merge(counts1, individuals2, by =1)

table(sort(dta$gene4))
alleles.df <- as.data.frame(table(sort(dta$gene4)))

MAF <- maf(as.vector(table(sort(dta$gene4)))) # which is the same as 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB)
MAF
dta2 <- na.omit(dta)
#dta2 <- dta[!is.na(dta$gene3),]
colnames(dta2)
sum(is.na(dta2$gene2))

nREF <- NROW(which(dta2$gene3==REF))
nALT <- NROW(which(dta2$gene3==ALT))
nREF+nALT

colnames(dta2)[2] <- "tpm" # just for ploting
dta2$gene4 <- as.character(dta2$gene4)
dta2$gene4 <- factor(dta2$gene4, levels=c(paste0(REF,REF),paste0(ALT,ALT)))

paste0(ALT,ALT)

geneid.allele.plot <- ggplot(dta2, aes(gene4, tpm, fill=gene4)) +
  geom_boxplot(width=0.5) +
  annotate("text", y=-3, x=c(1, 2), label=c(paste0("n = ",nREF), paste0("n = ",nALT)), size = 5) +
  #annotate("text", y=-1, x=c(1.5), label=paste0("maf\n",round(MAF,3)), size = 5) +
  #annotate("text", y=2, x=c(1.5), label=paste0(alleles.df$Var1)) + modify if we want to add the alleles counts
  scale_fill_manual(values = c("cornflowerblue","cornsilk2")) + 
  theme(legend.position = "none") + 
  xlab("chr1:5,046,274") +
  stat_compare_means(method = "t.test", label.x.npc = "centre", aes(label = paste0("t-test, p-value ", after_stat(p.format))), hjust = 0.5) +
  ylab(my_y_title)

geneid.allele.plot

#write.csv(dta,paste0("results/Manhathan/allele_",geneName,".csv"), row.names = F)


# origin2 <- readxl::read_excel("TPJ_16123_TableS1_Classified.xlsx") # if qe tried to input as csv we will have some errors in the file
# origin3 <- merge(individuals2, origin2, by = 1)
# 
# print(colnames(origin3)[3])
# 
# write.csv(origin3, paste0("results/alleleProportion",geneName,gene2,".csv"))
# 
# origin3.1 <- origin3[!is.na(origin3$VlaCategory),]
# 
# 
# result <- origin3.1 %>% #will remove het hosted in gene3 named column
#   group_by(VlaCategory, gene4) %>%
#   summarize(count = n()) %>%
#   group_by(VlaCategory) %>%
#   mutate(proportion = count / sum(count))
# 
# 
# # Generate all possible combinations of VlaCategory and gene4
# all_combinations <- expand.grid(
#   VlaCategory = unique(result$VlaCategory),
#   gene4 = unique(result$gene4)
# )
# 
# # Find missing combinations
# missing_combinations <- anti_join(all_combinations, result, by = c("VlaCategory", "gene4"))
# missing_combinations
# # Add the missing combinations to the original data with count and proportion set to 0
# missing_rows <- missing_combinations %>%
#   mutate(count = 0, proportion = 0)
# 
# # Combine the missing rows with the original data
# result <- bind_rows(result, missing_rows)
# 
# result$gene4 <- as.character(result$gene4)
# result$gene4 <- factor(result$gene4, levels=c(paste0(REF,REF), paste0(REF,ALT), paste0(ALT,ALT)))
# 
# levels(result$gene4)
# 
# #order by reference frequency
# # result1 <- result %>% 
# #   filter(gene4==paste0(REF,REF)) %>% 
# #   arrange(desc(proportion))
# # result$VlaCategory <- factor(result$VlaCategory, levels=result1$VlaCategory)
# 
# 
# #order by closeness closeness of relation to lines in the WiDiv panel.
# category <- c("Teosinte (M)", "Teosinte (P)", "Tropical", "Temperate (E)", "Temperate (C)", "Temperate (N)")
# result$VlaCategory <- factor(result$VlaCategory, levels=category)
# result$VlaCategory
# 
# proportion <- ggplot(result, aes(x = VlaCategory, y = proportion, fill = gene4)) +
#   geom_col(colour = "black", position = "fill") +
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_manual(values = c("cornflowerblue","black","cornsilk2")) + 
#   xlab("Origin") +
#   labs(fill = "genotype") +
#   theme(legend.position="top") + 
#   theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 14), legend.box.spacing = unit(0, "cm"),
#         plot.margin = unit(c(0,1,.5,1), "cm"))
# proportion


#### combine plots
### combine all plots together
# plots <- align_plots(gTWAS.fdr, relation.plot, align = 'v', axis = 'l')
# bottom_row <- plot_grid(plots[[2]], geneid.allele.plot, proportion, rel_widths = c(.9, .5,1),
#                         nrow = 1, labels = c('B', 'C', 'D'),hjust = 0, vjust = 1)
# 
# fancyplots <- plot_grid(plots[[1]], bottom_row, ncol = 1,  labels = c('A', ''),hjust = 0, vjust = 1)
# fancyplots

#ggsave(plot=fancyplots, paste0("results/Manhathan/main_peak_",geneName,".png"), width = 14, height = 8.5)
#ggsave(plot=fancyplots, paste0("results/Manhathan/main_peak_",geneName,".pdf"), width = 14, height = 8.5)






# svg(paste0("results/Manhathan/maf_",geneName,".svg"), width = 3, height= 5)
# gg.MAF
# 
# dev.off()


##### print the second most significant peak

# Define your condition
condition <- NROW(unique(asso.hits$CHROM)) > 1  # Change this to your actual condition


# Check if the condition is true
if (condition) {
  # Code to execute if the condition is true
  cat("Condition is true. Continuing with the next lines.\n")
  
  # Your additional code goes here
  
  
  #comexar dsd aqui
  
  maxchrhit <- asso.hits$CHROM[which(asso.hits$log10 == max(asso.hits$log10))]
  maxchrhit <- maxchrhit[1]
  second.hit <- asso.hits[which(asso.hits$CHROM != maxchrhit),]
  
  #second.hit <- second.hit[2,] # just filter the SNP we want
  
  geneid2 <- second.hit$SNP[which(second.hit$log10 == max(second.hit$log10))] 
  geneid2
  maxchrhit2 <- second.hit$CHROM[which(second.hit$log10 == max(second.hit$log10))]
  maxchrhit2
  
  genodata.2 <- paste0("1515.SNPs/mvpData/mvp.",maxchrhit2,".vcf.geno.desc")
  genomap.2 <- paste0("1515.SNPs/mvpData/mvp.",maxchrhit2,".vcf.geno.map")
  ind2 <- paste0("1515.SNPs/mvpData/mvp.",maxchrhit2,".vcf.geno.ind")
  
  individuals2.1 <- fread(ind2, data.table = F)
  genomap2 <- fread(genomap.2)
  genodata2 <- attach.big.matrix(genodata.2)
  
  
  
  second.hitgene2 <- second.hit$SNP[which(second.hit$log10 == max(second.hit$log10))]
  second.hitgene2 <- second.hitgene2[1]
  second.hitPOS <- second.hit$POS[which(second.hit$SNP == second.hitgene2)]
  second.hitREF <- second.hit$REF[which(second.hit$SNP == second.hitgene2)]
  second.hitALT <- second.hit$ALT[which(second.hit$SNP == second.hitgene2)]
  
  
  ### SNP effect plot as boxplot
  individuals2.1$gene2 <- genodata2[which(genomap2$SNP==second.hitgene2),] # from the original file NOT the simplyfied version
  individuals2.1$gene3[individuals2.1$gene2==0] <- second.hitREF #ref
  individuals2.1$gene3[individuals2.1$gene2==2] <- second.hitALT #alt
  individuals2.1$gene3[individuals2.1$gene2==1] <- NA #het
  individuals2.1$gene3 <- factor(individuals2.1$gene3, levels = c(second.hitREF,second.hitALT))
  
  
  individuals2.1$gene4[individuals2.1$gene2==0] <- paste0(second.hitREF,second.hitREF) #ref
  individuals2.1$gene4[individuals2.1$gene2==2] <- paste0(second.hitALT,second.hitALT) #alt
  individuals2.1$gene4[individuals2.1$gene2==1] <- paste0(second.hitREF,second.hitALT) #het
  
  colnames(individuals2.1)
  colnames(individuals2.1)[1] <- "Taxa"
  
  
  individuals22 <- individuals2.1[c("Taxa","gene3", "gene4")]
  
  
  second.hitdta <- merge(counts1, individuals22, by =1)
  
  table(sort(second.hitdta$gene4))
  #calculate MAF using library(HardyWeinberg)
  MAF <- maf(as.vector(table(sort(second.hitdta$gene4)))) # which is the same as 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB)
  MAF
  second.hitdta2 <- na.omit(second.hitdta)
  #dta2 <- dta[!is.na(dta$gene3),]
  colnames(second.hitdta2)
  sum(is.na(second.hitdta2$gene2))
  
  nREF <- NROW(which(second.hitdta2$gene3==second.hitREF))
  nALT <- NROW(which(second.hitdta2$gene3==second.hitALT))
  nREF+nALT
  
  colnames(second.hitdta2)[2] <- "tpm" # just for ploting
  second.hitdta2$gene4 <- as.character(second.hitdta2$gene4)
  second.hitdta2$gene4 <- factor(second.hitdta2$gene4, levels=c(paste0(second.hitREF,second.hitREF),paste0(second.hitALT,second.hitALT)))
  
  table(sort(second.hitdta2$gene4))
  #calculate MAF using library(HardyWeinberg)
  MAF2 <- maf(as.vector(table(sort(second.hitdta2$gene4)))) # which is the same as 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB)
  MAF2
  
  
  
  geneid.allele.plot2 <- ggplot(second.hitdta2, aes(gene4, tpm, fill=gene4)) +
    geom_boxplot(width=0.5) +
    annotate("text", y=-3, x=c(1, 2), label=c(paste0("n = ",nREF), paste0("n = ",nALT)), size = 5) +
    #annotate("text", y=-1, x=c(1.5), label=paste0("maf\n",round(MAF,3)), size = 5) +
    #annotate("text", y=2, x=c(1.5), label=paste0(alleles.df$Var1)) + modify if we want to add the alleles counts
    scale_fill_manual(values = c("orange2","grey50")) + 
    theme(legend.position = "none") + 
    xlab("chr6:26,934,399") +
    stat_compare_means(method = "t.test", label.x.npc = "centre", aes(label = paste0("t-test, p-value = ", after_stat(p.format))), hjust = 0.5) +
    ylab(my_y_title)
  
  geneid.allele.plot2
  
  #write.csv(second.hitdta,paste0("results/Manhathan/allele_secondHit_",geneName,".csv"), row.names = F)
  
  
  
  #origin2 <- readxl::read_excel("TPJ_16123_TableS1_Classified.xlsx") # if qe tried to input as csv we will have some errors in the file
  # origin3.2 <- merge(individuals2.1, origin2, by = 1)
  # 
  # write.csv(origin3.2, paste0("results/alleleProportion_secondHit_",geneName,gene2,".csv"))
  # 
  # origin3.3 <- origin3.2[!is.na(origin3.2$VlaCategory),]
  # 
  # 
  # result2 <- origin3.3 %>% #will remove het hosted in gene3 named column
  #   group_by(VlaCategory, gene4) %>%
  #   summarize(count = n()) %>%
  #   group_by(VlaCategory) %>%
  #   mutate(proportion = count / sum(count))
  # result2
  # 
  # # Generate all possible combinations of VlaCategory and gene4
  # all_combinations2 <- expand.grid(
  #   VlaCategory = unique(result2$VlaCategory),
  #   gene4 = unique(result2$gene4)
  # )
  # 
  # # Find missing combinations
  # missing_combinations2 <- anti_join(all_combinations2, result2, by = c("VlaCategory", "gene4"))
  # missing_combinations2
  # # Add the missing combinations to the original data with count and proportion set to 0
  # missing_rows2 <- missing_combinations2 %>%
  #   mutate(count = 0, proportion = 0)
  # 
  # # Combine the missing rows with the original data
  # result2 <- bind_rows(result2, missing_rows2)
  # 
  # result2$gene4 <- as.character(result2$gene4)
  # result2$gene4 <- factor(result2$gene4, levels=c(paste0(second.hitREF,second.hitREF),
  #                                                 paste0(second.hitREF,second.hitALT),
  #                                                 paste0(second.hitALT,second.hitALT)))
  # 
  # levels(result2$gene4)
  # 
  # 
  # write.csv(result2, paste0("results/alleleProportion_summary_secondHit_",geneName,gene2,".csv"))
  # 
  # 
  # category <- c("Teosinte (M)", "Teosinte (P)", "Tropical", "Temperate (E)", "Temperate (C)", "Temperate (N)")
  # result2$VlaCategory <- factor(result2$VlaCategory, levels=category)
  # result2$VlaCategory
  # 
  # proportion2 <- ggplot(result2, aes(x = VlaCategory, y = proportion, fill = gene4)) +
  #   geom_col(colour = "black", position = "fill") +
  #   scale_y_continuous(labels = scales::percent) +
  #   scale_fill_manual(values = c("orange2","black","grey50")) + 
  #   xlab("Origin") +
  #   labs(fill = "genotype") +
  #   theme(legend.position="top") + 
  #   theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 14), legend.box.spacing = unit(0, "cm"),
  #         plot.margin = unit(c(0,1,.5,1), "cm"))
  # proportion2
  # 
  
} else {
  # Code to execute if the condition is not true
  cat("Condition is not true. Killing the job.\n")
  stop("Job terminated due to the false condition.")
}


paste0("chr",maxchrhit2)

minW <- second.hitPOS-500000
maxw <- second.hitPOS+500000

### SNPs track
second.hitPOS.snps.window <- gwas.dat %>%
  filter(CHROM==maxchrhit2, POS > minW & POS < maxw)
head(second.hitPOS.snps.window)

### Gene track
second.hitPOS.window <- genes.desc %>%
  filter(chr==maxchrhit2, start > minW &  end < maxw)
second.hitPOS.window2 <- second.hitPOS.window[,c(1,2,3,4,6,8,9,14)] 
head(second.hitPOS.window2)



second.hitPOS.window2$nothing <- NA
second.hitPOS.window2$nothing[grepl("Zm00001eb264310", second.hitPOS.window2$B73_V5)] <- "blue4"
second.hitPOS.window2$nothing[!grepl("Zm00001eb264310", second.hitPOS.window2$B73_V5)] <- "grey"

chr6_26934399.ld <- fread("chr6_26934399.ld", data.table = F)
second.hitPOS.snps.windowld <- merge(second.hitPOS.snps.window, chr6_26934399.ld , by.x=1, by.y=6)
colnames(second.hitPOS.snps.windowld)
second.hitPOS.snps.windowld$R2comp <- NA
second.hitPOS.snps.windowld$R2comp[which(second.hitPOS.snps.windowld$R2 >= 0.8)] <- "veryHigh"
second.hitPOS.snps.windowld$R2comp[which(second.hitPOS.snps.windowld$R2 >= 0.6 & second.hitPOS.snps.windowld$R2 < 0.8)] <- "High"
second.hitPOS.snps.windowld$R2comp[which(second.hitPOS.snps.windowld$R2 >= 0.5 & second.hitPOS.snps.windowld$R2 < 0.6)] <- "medium"
second.hitPOS.snps.windowld$R2comp[which(second.hitPOS.snps.windowld$R2 < 0.5)] <- "low"

zoomIn <- ggplot() + 
  geom_point(data=second.hitPOS.snps.windowld, aes(POS/1000000, -log10(gene), fill=R2comp), shape = 21, alpha=0.5, size=4) +
  ylim(0,18) +
  scale_fill_manual(values = c("red4","orange3","yellow2","white"), 
                    breaks = c("veryHigh", "High", "medium","low")) +
  geom_point(data=second.hitPOS.window2[second.hitPOS.window2$strand=="+",],
              aes(x=((start+end)/2)/1000000, y=10, colour=nothing), shape="\u25BA", size=6) +
  geom_point(data=second.hitPOS.window2[second.hitPOS.window2$strand=="-",], 
              aes(x=((start+end)/2/1000000), y=12, colour=nothing), shape="\u25C4", size=6) +
  scale_color_manual(values = c("blue4", "grey")) +
  geom_hline(yintercept = -log10(0.05/nrow(gwas.dat)), linetype=2) + 
  #geom_vline(xintercept = second.hitPOS/1000000, linetype="dotted", color = "blue", linewidth=1.5) +
  #geom_vline(xintercept = mads69pos, linetype=2, color = "green", linewidth=1) +
  #geom_point(aes(y=0, x=rgd1pos2), fill="white", size = 15, shape = 24) +
  ylab(expression(-log[10](p-value))) + 
  xlab(paste0("Chromosome ",maxchrhit2, " (Mb)")) +
  theme(legend.position = "none")
zoomIn 


########## for the second most sig snp


# Define your condition
condition2 <- NROW(unique(asso.hits$CHROM)) > 2  # Change this to your actual condition


# Check if the condition is true
if (condition2) {
  # Code to execute if the condition is true
  cat("Condition is true. Continuing with the next lines.\n")
  
  # Your additional code goes here
  
  
  #comexar dsd aqui
  
  maxchrhit <- asso.hits$CHROM[which(asso.hits$log10 == max(asso.hits$log10))]
  maxchrhit <- maxchrhit[1]
  second.hit <- asso.hits[which(asso.hits$CHROM != maxchrhit),]
  
  second.hit <- second.hit[2,] # just filter the SNP we want
  
  geneid2 <- second.hit$SNP[which(second.hit$log10 == max(second.hit$log10))] 
  geneid2
  maxchrhit2 <- second.hit$CHROM[which(second.hit$log10 == max(second.hit$log10))]
  maxchrhit2
  
  genodata.2 <- paste0("1515.SNPs/mvpData/mvp.",maxchrhit2,".vcf.geno.desc")
  genomap.2 <- paste0("1515.SNPs/mvpData/mvp.",maxchrhit2,".vcf.geno.map")
  ind2 <- paste0("1515.SNPs/mvpData/mvp.",maxchrhit2,".vcf.geno.ind")
  
  individuals2.1 <- fread(ind2, data.table = F)
  genomap2 <- fread(genomap.2)
  genodata2 <- attach.big.matrix(genodata.2)
  
  
  
  second.hitgene2 <- second.hit$SNP[which(second.hit$log10 == max(second.hit$log10))]
  second.hitgene2 <- second.hitgene2[1]
  second.hitPOS <- second.hit$POS[which(second.hit$SNP == second.hitgene2)]
  second.hitREF <- second.hit$REF[which(second.hit$SNP == second.hitgene2)]
  second.hitALT <- second.hit$ALT[which(second.hit$SNP == second.hitgene2)]
  
  
  ### SNP effect plot as boxplot
  individuals2.1$gene2 <- genodata2[which(genomap2$SNP==second.hitgene2),] # from the original file NOT the simplyfied version
  individuals2.1$gene3[individuals2.1$gene2==0] <- second.hitREF #ref
  individuals2.1$gene3[individuals2.1$gene2==2] <- second.hitALT #alt
  individuals2.1$gene3[individuals2.1$gene2==1] <- NA #het
  individuals2.1$gene3 <- factor(individuals2.1$gene3, levels = c(second.hitREF,second.hitALT))
  
  
  individuals2.1$gene4[individuals2.1$gene2==0] <- paste0(second.hitREF,second.hitREF) #ref
  individuals2.1$gene4[individuals2.1$gene2==2] <- paste0(second.hitALT,second.hitALT) #alt
  individuals2.1$gene4[individuals2.1$gene2==1] <- paste0(second.hitREF,second.hitALT) #het
  
  colnames(individuals2.1)
  colnames(individuals2.1)[1] <- "Taxa"
  
  
  individuals22 <- individuals2.1[c("Taxa","gene3", "gene4")]
  
  
  second.hitdta <- merge(counts1, individuals22, by =1)
  
  table(sort(second.hitdta$gene4))
  #calculate MAF using library(HardyWeinberg)
  MAF <- maf(as.vector(table(sort(second.hitdta$gene4)))) # which is the same as 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB)
  MAF
  second.hitdta2 <- na.omit(second.hitdta)
  #dta2 <- dta[!is.na(dta$gene3),]
  colnames(second.hitdta2)
  sum(is.na(second.hitdta2$gene2))
  
  nREF <- NROW(which(second.hitdta2$gene3==second.hitREF))
  nALT <- NROW(which(second.hitdta2$gene3==second.hitALT))
  nREF+nALT
  
  colnames(second.hitdta2)[2] <- "tpm" # just for ploting
  second.hitdta2$gene4 <- as.character(second.hitdta2$gene4)
  second.hitdta2$gene4 <- factor(second.hitdta2$gene4, levels=c(paste0(second.hitREF,second.hitREF),paste0(second.hitALT,second.hitALT)))
  
  table(sort(second.hitdta2$gene4))
  #calculate MAF using library(HardyWeinberg)
  MAF3 <- maf(as.vector(table(sort(second.hitdta2$gene4)))) # which is the same as 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB)
  MAF3
  
  
  geneid.allele.plot3 <- ggplot(second.hitdta2, aes(gene4, tpm, fill=gene4)) +
    geom_boxplot(width=0.5) +
    annotate("text", y=-3, x=c(1, 2), label=c(paste0("n = ",nREF), paste0("n = ",nALT)), size = 5) +
    #annotate("text", y=-1, x=c(1.5), label=paste0("maf\n",round(MAF,3)), size = 5) +
    #annotate("text", y=2, x=c(1.5), label=paste0(alleles.df$Var1)) + modify if we want to add the alleles counts
    scale_fill_manual(values = c("green4","cyan")) + 
    theme(legend.position = "none") + 
    xlab("chr7:7,634,465") +
    stat_compare_means(method = "t.test", label.x.npc = "centre", aes(label = paste0("t-test, p-value = ", after_stat(p.format))), hjust = 0.5) +
    ylab(my_y_title)
  
  geneid.allele.plot3
  
  #write.csv(second.hitdta,paste0("results/Manhathan/allele_selectedHit_",geneName,".csv"), row.names = F)
  
  
  
  #origin2 <- readxl::read_excel("TPJ_16123_TableS1_Classified.xlsx") # if qe tried to input as csv we will have some errors in the file
  # origin3.2 <- merge(individuals2.1, origin2, by = 1)
  # 
  # write.csv(origin3.2, paste0("results/alleleProportion_selectedHit_",geneName,gene2,".csv"))
  # 
  # origin3.3 <- origin3.2[!is.na(origin3.2$VlaCategory),]
  # 
  # 
  # result3 <- origin3.3 %>% #will remove het hosted in gene3 named column
  #   group_by(VlaCategory, gene4) %>%
  #   summarize(count = n()) %>%
  #   group_by(VlaCategory) %>%
  #   mutate(proportion = count / sum(count))
  # result3
  # 
  # # Generate all possible combinations of VlaCategory and gene4
  # all_combinations2 <- expand.grid(
  #   VlaCategory = unique(result3$VlaCategory),
  #   gene4 = unique(result3$gene4)
  # )
  # 
  # # Find missing combinations
  # missing_combinations2 <- anti_join(all_combinations2, result3, by = c("VlaCategory", "gene4"))
  # missing_combinations2
  # # Add the missing combinations to the original data with count and proportion set to 0
  # missing_rows2 <- missing_combinations2 %>%
  #   mutate(count = 0, proportion = 0)
  # 
  # # Combine the missing rows with the original data
  # result3 <- bind_rows(result3, missing_rows2)
  # 
  # result3$gene4 <- as.character(result3$gene4)
  # result3$gene4 <- factor(result3$gene4, levels=c(paste0(second.hitREF,second.hitREF),
  #                                                 paste0(second.hitREF,second.hitALT),
  #                                                 paste0(second.hitALT,second.hitALT)))
  # 
  # levels(result3$gene4)
  # 
  # 
  # write.csv(result3, paste0("results/alleleProportion_summary_selectedHit_",geneName,gene2,".csv"))
  # 
  # 
  # category <- c("Teosinte (M)", "Teosinte (P)", "Tropical", "Temperate (E)", "Temperate (C)", "Temperate (N)")
  # result3$VlaCategory <- factor(result3$VlaCategory, levels=category)
  # result3$VlaCategory
  # 
  # proportion3 <- ggplot(result3, aes(x = VlaCategory, y = proportion, fill = gene4)) +
  #   geom_col(colour = "black", position = "fill") +
  #   scale_y_continuous(labels = scales::percent) +
  #   scale_fill_manual(values = c("green4","black","cyan")) + 
  #   xlab("Origin") +
  #   labs(fill = "genotype") +
  #   theme(legend.position="top") + 
  #   theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 14), legend.box.spacing = unit(0, "cm"),
  #         plot.margin = unit(c(0,1,.5,1), "cm"))
  # proportion3
  
  
} else {
  # Code to execute if the condition is not true
  cat("Condition is not true. Killing the job.\n")
  stop("Job terminated due to the false condition.")
}


#install.packages("gggenes")
library(gggenes)



paste0("chr",maxchrhit2)

minW <- second.hitPOS-500000
maxw <- second.hitPOS+500000

### SNPs track
third.hitPOS.snps.window <- gwas.dat %>%
  filter(CHROM==maxchrhit2, POS > minW & POS < maxw)
head(third.hitPOS.snps.window)

### Gene track
third.hitPOS.window <- genes.desc %>%
  filter(chr==maxchrhit2, start > minW &  end < maxw)
third.hitPOS.window2 <- third.hitPOS.window[,c(1,2,3,4,6,8,9,14)] 
head(third.hitPOS.window2)

third.hitPOS.window2$nothing <- NA
third.hitPOS.window2$nothing[grepl("Zm00001eb300700", third.hitPOS.window2$B73_V5)] <- "blue4"
third.hitPOS.window2$nothing[!grepl("Zm00001eb300700", third.hitPOS.window2$B73_V5)] <- "grey"


chr7_7634465.ld <- fread("chr7_7634465.ld", data.table = F)
third.hitPOS.snps.windowld <- merge(third.hitPOS.snps.window, chr7_7634465.ld , by.x=1, by.y=6)
colnames(third.hitPOS.snps.windowld)
third.hitPOS.snps.windowld$R2comp <- NA
third.hitPOS.snps.windowld$R2comp[which(third.hitPOS.snps.windowld$R2 >= 0.8)] <- "veryHigh"
third.hitPOS.snps.windowld$R2comp[which(third.hitPOS.snps.windowld$R2 >= 0.6 & third.hitPOS.snps.windowld$R2 < 0.8)] <- "High"
third.hitPOS.snps.windowld$R2comp[which(third.hitPOS.snps.windowld$R2 >= 0.5 & third.hitPOS.snps.windowld$R2 < 0.6)] <- "medium"
third.hitPOS.snps.windowld$R2comp[which(third.hitPOS.snps.windowld$R2 < 0.5)] <- "low"


zoomIn3 <- ggplot() + 
  geom_point(data=third.hitPOS.snps.windowld, aes(POS/1000000, -log10(gene), fill=R2comp), shape = 21, alpha=0.5, size=4) +
  ylim(0,18) +
  scale_fill_manual(values = c("red4","orange3","yellow2","white"), 
                    breaks = c("veryHigh", "High", "medium","low")) +
  geom_point(data=third.hitPOS.window2[third.hitPOS.window2$strand=="+",],
             aes(x=((start+end)/2)/1000000, y=10, colour=nothing), shape="\u25BA", size=6) +
  geom_point(data=third.hitPOS.window2[third.hitPOS.window2$strand=="-",], 
             aes(x=((start+end)/2/1000000), y=12, colour=nothing), shape="\u25C4", size=6) +
  scale_color_manual(values = c("blue4", "grey")) +
  geom_hline(yintercept = -log10(0.05/nrow(gwas.dat)), linetype=2) + 
  #geom_vline(xintercept = second.hitPOS/1000000, linetype="dotted", color = "blue", linewidth=1.5) +
  #geom_vline(xintercept = mads69pos, linetype=2, color = "green", linewidth=1) +
  #geom_point(aes(y=0, x=rgd1pos2), fill="white", size = 15, shape = 24) +
  ylab(expression(-log[10](p-value))) + 
  xlab(paste0("Chromosome ",maxchrhit2, " (Mb)")) +
  theme(legend.position = "none")
zoomIn3 


#mads1
minW <- POS-1000000
maxw <- POS+1000000

### SNPs track
hitPOS.snps.window <- gwas.dat %>%
  filter(CHROM==1, POS > minW & POS < maxw)
head(hitPOS.snps.window)

### Gene track
hitPOS.window <- genes.desc %>%
  filter(chr==1, start > minW &  end < maxw)
hitPOS.window2 <- hitPOS.window[,c(1,2,3,4,6,8,9,14)] 
head(hitPOS.window2)

hitPOS.window2$nothing <- NA
hitPOS.window2$nothing[grepl("Zm00001eb001670", hitPOS.window2$B73_V5)] <- "blue4"
hitPOS.window2$nothing[!grepl("Zm00001eb001670", hitPOS.window2$B73_V5)] <- "grey"


chr1_5046274.ld <- fread("chr1_5046274.ld", data.table = F)
hitPOS.snps.windowld <- merge(hitPOS.snps.window, chr1_5046274.ld , by.x=1, by.y=6)
colnames(hitPOS.snps.windowld)
hitPOS.snps.windowld$R2comp <- NA
hitPOS.snps.windowld$R2comp[which(hitPOS.snps.windowld$R2 >= 0.8)] <- "veryHigh"
hitPOS.snps.windowld$R2comp[which(hitPOS.snps.windowld$R2 >= 0.6 & hitPOS.snps.windowld$R2 < 0.8)] <- "High"
hitPOS.snps.windowld$R2comp[which(hitPOS.snps.windowld$R2 >= 0.5 & hitPOS.snps.windowld$R2 < 0.6)] <- "medium"
hitPOS.snps.windowld$R2comp[which(hitPOS.snps.windowld$R2 < 0.5)] <- "low"

zoomIncis <- ggplot() + 
  geom_point(data=hitPOS.snps.windowld, aes(POS/1000000, -log10(gene), fill=R2comp), shape = 21, alpha=0.9, size=3) +
  ylim(0,50) +
  scale_fill_manual(values = c("red4","orange3","yellow2","white"), 
                    breaks = c("veryHigh", "High", "medium","low")) +
  geom_point(data=hitPOS.window2[hitPOS.window2$strand=="+",],
             aes(x=((start+end)/2)/1000000, y=49, colour=nothing), shape="\u25BA", size=6) +
  geom_point(data=hitPOS.window2[hitPOS.window2$strand=="-",], 
             aes(x=((start+end)/2/1000000), y=45, colour=nothing), shape="\u25C4", size=6) +
  scale_color_manual(values = c("blue4", "grey")) +
  geom_hline(yintercept = -log10(0.05/nrow(gwas.dat)), linetype=2) + 
  #geom_vline(xintercept = second.hitPOS/1000000, linetype="dotted", color = "blue", linewidth=1.5) +
  #geom_vline(xintercept = mads69pos, linetype=2, color = "green", linewidth=1) +
  #geom_point(aes(y=0, x=rgd1pos2), fill="white", size = 15, shape = 24) +
  ylab(expression(-log[10](p-value))) + 
  xlab(paste0("Chromosome ","1", " (Mb)")) +
  theme(legend.position = "none")
zoomIncis 


######################################
plots <- align_plots(gTWAS.fdr, relation.plot, zoomIn3, align = 'v', axis = 'l')
bottom_row <- plot_grid(plots[[2]], zoomIncis, zoomIn, rel_widths = c(1, 1,1),
                        nrow = 1, labels = c('B', 'C', 'D'),hjust = 0, vjust = 1, label_size = 18)
bottom_row2 <- plot_grid(plots[[3]], geneid.allele.plot, geneid.allele.plot2, geneid.allele.plot3, rel_widths = c(1, .63,.63,.63),
                         nrow = 1, labels = c('E', 'F', 'G', 'H'),hjust = 0, vjust = 1, label_size = 18)
#bottom_row3 <- plot_grid(plots[[4]], geneid.allele.plot3, proportion3, rel_widths = c(.9, .5,1),
#                         nrow = 1, labels = c('H', 'I', 'J'),hjust = 0, vjust = 1, label_size = 18)
fancyplots3 <- plot_grid(plots[[1]], bottom_row, bottom_row2, ncol = 1,  labels = c('A', ''),hjust = 0, vjust = 1, label_size = 18)
#fancyplots3


#ggsave(plot=fancyplots2, paste0("results/Manhathan/Fig4_",geneName,".png"), width = 15, height = 13)
ggsave(plot=fancyplots3, paste0("results/Manhathan/Fig4_",geneName,".svg"), width = 17, height = 13)
ggsave(plot=fancyplots3, paste0("results/Manhathan/Fig4_",geneName,".png"), width = 17, height = 13)

# 
# svg(paste0("results/Manhathan/Fig4_",geneName,".svg"), width = 15, height = 13)
# fancyplots2
# dev.off()
#for PNG open in inkscape then save it as png.