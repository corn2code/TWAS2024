geneName2 = "zap1"  #args[1] #mads1 zcn26 arf34 late

index <- fread("indexGenes.csv", data.table = F)
geneid2 <- index$ID[which(index$Name == geneName2)]

print(geneid2)


print("loading data")
gwas.dat2<-fread(paste0("results/",geneName2,'.csv.gz'), data.table = F)
print("loading data finished")

nCHR <- length(unique(gwas.dat2$CHROM))
gwas.dat2$BPcum <- NA
s <- 0
nbp <- c()


for (i in sort(unique(gwas.dat2$CHROM))){
  nbp[i] <- max(gwas.dat2[gwas.dat2$CHROM == i,]$POS)
  gwas.dat2[gwas.dat2$CHROM == i,"BPcum"] <- gwas.dat2[gwas.dat2$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.dat2 %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))



#install.packages("CMplot")
# library("CMplot")
# 
# pdf(paste0("results/Manhathan/qqplot",geneName2,".pdf"), height = 4, width = 4)
# par(mar=c(5.1,5.5,4.1,2.1))
# CMplot(gwas.dat2[,c(1,2,3,8)],plot.type="q",box=FALSE,file="jpg",dpi=300,
#        conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
#        file.output=F,verbose=TRUE,width=5,height=5)
# dev.off()

colnames(gwas.dat2)[8] <- "gene"
gwas.dat2.RNA2 <- gwas.dat2[order(gwas.dat2$gene)[1:200000],]
#colnames(gwas.dat2.RNA)[8] <- trait


snps <- nrow(gwas.dat2)

# to set the limite of y axis
if (max(-log10(gwas.dat2.RNA2$gene)) > -log10(0.05/snps) ) {
  limite2 = ceiling(max(-log10(gwas.dat2.RNA2$gene)))
} else {
  limite2 = ceiling(-log10(0.05/snps))
}


unique(gwas.dat2.RNA2$CHROM)

#to name hits in the plot 
genes.desc <- fread("/work/schnablelab/vladimir/TWAS2023.2/B73_geneModels_v5.csv", data.table = F)
hits1 <- fread("FTHits_PlusoneMichael.csv", data.table = F)
colnames(hits1)

hits <- merge(genes.desc, hits1[,c(1,2)] , by = 1)


# to position hit

chr.pos2 <- as.numeric(hits$chr[which(hits$B73_V5 == geneid2)])-1


#gene.pos.x <- sum(nbp[1:chr.pos2]) + hits$start[which(hits$B73_V5 == geneid2)]

# to set the position of the gene
if (sum(nbp[1:chr.pos2]) == 308450369 ) {
  gene.pos.x2 = hits$start[which(hits$B73_V5 == geneid2)]
} else {
  gene.pos.x2 = sum(nbp[1:chr.pos2]) + hits$start[which(hits$B73_V5 == geneid2)]
}

sum(1:1)

#if the position of the gene is in chr 2 the formula above will place them as 1 with this we are fixing for it.
if(chr.pos2 == 1){
  gene.pos.x2 = gene.pos.x2 + 308450369
} else {
  gene.pos.x2 = gene.pos.x2
}


gene.pos.y <- 2 #limite*.9
gene.id2 <- hits$B73_V5[which(hits$B73_V5 == geneid2)]
gene.name2 <- hits$symbol[which(hits$B73_V5 == geneid2)]

gene.name22 <- paste0(gene.name2," (",gene.id2,")")


mads69 <- genes.desc[which(genes.desc$B73_V5 == "Zm00001eb143080"),1:9]
mads69pos <- sum(nbp[1:2])+(mads69$start+mads69$end)/2


gTWAS.fdr2 <- ggplot() + 
  rasterise(geom_point(data=gwas.dat2.RNA2, aes(BPcum, -log10(gene), colour=factor(CHROM, levels = c(1:10))), size = 2),dpi=600) + 
  geom_hline(yintercept = -log10(0.05/snps), linetype=2) + 
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  geom_point(aes(y=gene.pos.y, x=gene.pos.x2), fill="red", size = 5, shape = 24) +
  geom_point(aes(y=gene.pos.y, x=mads69pos), fill="black", size = 5, shape = 24) +
  #geom_point(aes(y=gene.pos.y, x=rgd1pos), fill="white", size = 5, shape = 24) +
  #annotate("text", label=gene.name2, y=gene.pos.y, x=gene.pos.x, parse=T, size=5, hjust = -0.01) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](p-value))) + 
  xlab("Chromosome") + 
  scale_y_continuous(limits = c(2, limite2),
                     breaks = seq(0, limite2, 5))
#+ labs(title="", subtitle=paste0(geneid2))

gTWAS.fdr2

# ggsave(plot=gTWAS.fdr, paste0("results/Manhathan/Manhathan_",geneName2,".png"), width = 7, height = 4.5)
# ggsave(plot=gTWAS.fdr, paste0("results/Manhathan/Manhathan_",geneName2,".pdf"), width = 7, height = 4.5)


asso.hits2 <- data.frame(matrix(NA,    # Create empty data frame
                               nrow = 1,
                               ncol = 9))
colnames(asso.hits2) <- colnames(gwas.dat2.RNA2)

colnames(asso.hits2)

#IDs for genes with FDR less than 0.05
gwas.dat2.RNA2$log10 <- -log10(gwas.dat2.RNA2$gene)

asso.hits2 <- gwas.dat2.RNA2[which(gwas.dat2.RNA2$log10 > -log10(0.05/snps)),]
write.csv(asso.hits2,paste0("results/associatedHits_eGWAS/",geneName2,".csv"), row.names = F)

#### do the scatter plot
##
counts <- fread("../counts.NE2020.693.filtered.txt", data.table = F)
row.names(counts) <- counts$taxa

#load the new phe data
phe <- read.table("../pheno_693.txt", head = TRUE) # BLUES.NE2020_corrected.csv  # BLUES.NE2021_corrected.csv
colnames(phe)

genes <- hits1$Gene[which(hits1$symbo == geneName2)]

colnames(phe)


trait2=hits1$trait2plot[which(hits1$symbo == geneName2)]
colnames(phe[trait2])
TPM.genes2 <- data.frame(counts[,c(1,which(colnames(counts) 
                                          %in% genes))])

colnames(TPM.genes2)[1] <- "taxa"

phenotype2 <- phe[,c(1,c(which(colnames(phe) %in% trait2)))]
phe.TPM.genes2 <- merge(phenotype2, TPM.genes2, by="taxa")


colname22 <- gsub(".sp", "",  colnames(phe.TPM.genes2)[2])

colnames(phe.TPM.genes2)[2] <- colname22


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
#unique(phe.TPM.genes.pop$Group_As_On_Fig2)

#j <- ncol(phe.TPM.genes.pop)-1


phe.TPM.genes.pop22 <- merge(popdata[,c(1,3)], phe.TPM.genes2, by="taxa")


colnames(phe.TPM.genes.pop22)[3:4] <- c("phenotype","expression")


colnames(phe.TPM.genes2)[2]
paste0(gene.name2," (TPM)")
paste0("Days to ", colnames(phe.TPM.genes2)[2])

relation.plot2 <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop22[phe.TPM.genes.pop22$Group=="Others",], aes(x=phenotype, y=expression, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop22[phe.TPM.genes.pop22$Group!="Others",], aes(x=phenotype, y=expression,fill=Group), size = 5, shape = 21),dpi=600) +
  #xlim(65,90)
  ylab(paste0(gene.name2," (TPM)")) + 
  xlab(paste0("Days to ", colnames(phe.TPM.genes2)[2])) + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none") + 
  stat_poly_line(data= phe.TPM.genes.pop22, aes(x=phenotype, y=phe.TPM.genes.pop22[,4]), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop22, aes(x=phenotype, y=phe.TPM.genes.pop22[,4]),size = 5)
relation.plot2

### check from here
# gene2 <- "chr9_160163325"
# chrMostSig <- "9"

#gwas.dat2.RNA$CHROM[which(gwas.dat2.RNA$log10 == max(gwas.dat2.RNA$log10))]

chrMostSig2 <- gwas.dat2.RNA2$CHROM[which(gwas.dat2.RNA2$log10 == max(gwas.dat2.RNA2$log10))]
chrMostSig2 <- chrMostSig2[1]
genodata2 <- paste0("1515.SNPs/mvpData/mvp.",chrMostSig2,".vcf.geno.desc")
genomap2 <- paste0("1515.SNPs/mvpData/mvp.",chrMostSig2,".vcf.geno.map")
ind2 <- paste0("1515.SNPs/mvpData/mvp.",chrMostSig2,".vcf.geno.ind")


individuals2 <- fread(ind2, data.table = F)
genomap2 <- fread(genomap2)
genodata2 <- attach.big.matrix(genodata2)

#NROW(asso.hits) == 0 modify hre if we dont want to plot not hits

gene22 <- gwas.dat2.RNA2$SNP[which(gwas.dat2.RNA2$log10 == max(gwas.dat2.RNA2$log10))]
gene22 <- gene22[1] # for arf34 the second most significant variation has the same outcome but is just a snp (CC and TT)

POS2 <- genomap2$POS[which(genomap2$SNP == gene22)]
REF2 <- genomap2$REF[which(genomap2$SNP == gene22)]
ALT2 <- genomap2$ALT[which(genomap2$SNP == gene22)]

### SNP effect plot as boxplot
individuals2$gene2 <- genodata2[which(genomap2$SNP==gene22),] # from the original file NOT the simplyfied version
individuals2$gene3[individuals2$gene2==0] <- REF2 #ref
individuals2$gene3[individuals2$gene2==2] <- ALT2 #alt
individuals2$gene3[individuals2$gene2==1] <- NA #het
individuals2$gene3 <- factor(individuals2$gene3, levels = c(REF2,ALT2))


individuals2$gene4[individuals2$gene2==0] <- paste0(REF2,REF2) #ref
individuals2$gene4[individuals2$gene2==2] <- paste0(ALT2,ALT2) #alt
individuals2$gene4[individuals2$gene2==1] <- paste0(REF2,ALT2) #het

colnames(individuals2)
colnames(individuals2)[1] <- "Taxa"
individuals22 <- individuals2[c("Taxa","gene3", "gene4")]

#counts <- fread("../counts.NE2020.693.filtered.txt", data.table = F)
counts12 <- counts[c("taxa",geneid2)]


dta.2 <- merge(counts12, individuals22, by =1)

table(sort(dta.2$gene4))
alleles.df2 <- as.data.frame(table(sort(dta.2$gene4)))
#calculate MAF using library(HardyWeinberg)


#come here


MAF2 <- maf(as.vector(table(sort(dta.2$gene4)))) # which is the same as 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB)
MAF2
dta22 <- na.omit(dta.2)
#dta2 <- dta[!is.na(dta$gene3),]
colnames(dta22)
sum(is.na(dta22$gene2))

nREF2 <- NROW(which(dta22$gene3==REF2))
nALT2 <- NROW(which(dta22$gene3==ALT2))
nREF2+nALT2

colnames(dta22)[2] <- "tpm" # just for ploting
dta22$gene4 <- as.character(dta22$gene4)
dta22$gene4 <- factor(dta22$gene4, levels=c(paste0(REF2,REF2),paste0(ALT2,ALT2)))

paste0(ALT2,ALT2)

geneid2.allele.plot2 <- ggplot(dta22, aes(gene4, tpm, fill=gene4)) +
  geom_boxplot(width=0.5) +
  annotate("text", y=-3, x=c(1, 2), label=c(paste0("n = ",nREF2), paste0("n = ",nALT2)), size = 5) +
  #annotate("text", y=-1, x=c(1.5), label=paste0("maf\n",round(MAF,3)), size = 5) +
  #annotate("text", y=2, x=c(1.5), label=paste0(alleles.df$Var1)) + modify if we want to add the alleles counts
  scale_fill_manual(values = c("cornflowerblue","cornsilk2")) + 
  theme(legend.position = "none") + 
  xlab(gene22) +
  ylab(paste0(geneName2," (TPM)"))

geneid2.allele.plot2

write.csv(dta,paste0("results/Manhathan/allele_",geneName2,".csv"), row.names = F)


origin2 <- readxl::read_excel("TPJ_16123_TableS1_Classified.xlsx") # if qe tried to input as csv we will have some errors in the file
origin32 <- merge(individuals22, origin2, by = 1)

print(colnames(origin32)[3])

write.csv(origin32, paste0("results/alleleProportion",geneName2,gene22,".csv"))

origin32.1 <- origin32[!is.na(origin32$VlaCategory),]


result2 <- origin32.1 %>% #will remove het hosted in gene3 named column
  group_by(VlaCategory, gene4) %>%
  summarize(count = n()) %>%
  group_by(VlaCategory) %>%
  mutate(proportion = count / sum(count))


# Generate all possible combinations of VlaCategory and gene4
all_combinations2 <- expand.grid(
  VlaCategory = unique(result2$VlaCategory),
  gene4 = unique(result2$gene4)
)

# Find missing combinations
missing_combinations2 <- anti_join(all_combinations2, result2, by = c("VlaCategory", "gene4"))
missing_combinations2
# Add the missing combinations to the original data with count and proportion set to 0
missing_rows2 <- missing_combinations2 %>%
  mutate(count = 0, proportion = 0)

# Combine the missing rows with the original data
result2 <- bind_rows(result2, missing_rows2)

result2$gene4 <- as.character(result2$gene4)
result2$gene4 <- factor(result2$gene4, levels=c(paste0(REF2,REF2), paste0(REF2,ALT2), paste0(ALT2,ALT2)))

levels(result2$gene4)

#order by reference frequency
# result21 <- result2 %>% 
#   filter(gene4==paste0(REF,REF)) %>% 
#   arrange(desc(proportion))
# result2$VlaCategory <- factor(result2$VlaCategory, levels=result21$VlaCategory)


#order by closeness closeness of relation to lines in the WiDiv panel.
category <- c("Teosinte (M)", "Teosinte (P)", "Tropical", "Temperate (E)", "Temperate (C)", "Temperate (N)")
result2$VlaCategory <- factor(result2$VlaCategory, levels=category)
result2$VlaCategory

proportion2 <- ggplot(result2, aes(x = VlaCategory, y = proportion, fill = gene4)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("cornflowerblue","black","cornsilk2")) + 
  xlab("Origin") +
  labs(fill = "genotype") +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 14), legend.box.spacing = unit(0, "cm"),
        plot.margin = unit(c(0,1,.5,1), "cm"))
proportion2




#install.packages("gggenes")
library(gggenes)

minW2 <- POS2-1000000
maxw2 <- POS2+1000000

### SNPs track
POS.snps.window2 <- gwas.dat2 %>%
  filter(CHROM==3, POS2 > minW2 & POS2 < maxw2)
head(POS.snps.window2)

### Gene track
POS.window2 <- genes.desc %>%
  filter(chr==3, start > minW2 &  end < maxw2)
POS.window2 <- POS.window2[,c(1,2,3,4,6,8,9,14)] 
head(POS.window2)


mads69pos2 <- (mads69$start+mads69$end)/2

zoomIn2 <- ggplot() + 
  geom_point(data=POS.snps.window2, aes(POS/1000000, -log10(gene)), shape = 21, fill="white", alpha=0.5, size=4) +
  ylim(0,18) +
  geom_jitter(data=POS.window2[POS.window2$strand=="+",],
              aes(x=((start+end)/2)/1000000, y=12, colour="red"), shape="\u25BA", size=6) +
  geom_jitter(data=POS.window2[POS.window2$strand=="-",], 
              aes(x=((start+end)/2/1000000), y=12, colour="blue4"), shape="\u25C4", size=6) +
  scale_color_manual(values = c("red3", "blue4")) +
  geom_hline(yintercept = -log10(0.05/nrow(gwas.dat2)), linetype=2) + 
  #geom_vline(xintercept = second.hitPOS/1000000, linetype="dotted", color = "blue", linewidth=1.5) +
  #geom_vline(xintercept = mads69pos, linetype=2, color = "green", linewidth=1) +
  #annotate("text", x=mads69pos, y=13.5, label="paste(italic(mads69))", size=5, parse=T) +
  geom_point(aes(y=0, x=mads69pos2/1000000), fill="black", size = 15, shape = 24) +
  ylab(expression(-log[10](p-value))) + 
  xlab(paste0("Chromosome ",3, " (Mb)")) +
  theme(legend.position = "none")
zoomIn2


dev.off()

######################################
plots <- align_plots(gTWAS.fdr, gTWAS.fdr2, relation.plot, relation.plot2, align = 'v', axis = 'l')
bottom_row <- plot_grid(plots[[3]],zoomIn , geneid2.allele.plot, proportion, rel_widths = c(.9, .8,.4,1),
                        nrow = 1, labels = c('C', 'D', 'E', 'F'),hjust = 0, vjust = 1, label_size = 18)
bottom_row2 <- plot_grid(plots[[4]], zoomIn2 , geneid2.allele.plot2, proportion2, rel_widths = c(.9, .8,.4,1),
                         nrow = 1, labels = c('G', 'H', 'I', 'J'),hjust = 0, vjust = 1, label_size = 18)

fancyplots3 <- plot_grid(plots[[1]], plots[[2]], bottom_row, bottom_row2, ncol = 1,  labels = c('A', 'B'),hjust = 0, vjust = 1, label_size = 18)
fancyplots3

#fancyplots3


#ggsave(plot=fancyplots2, paste0("results/Manhathan/Fig4_",geneName2,".png"), width = 15, height = 13)
ggsave(plot=fancyplots3, paste0("results/Manhathan/Fig6_mads60hotspot.svg"), width = 15, height = 16)
# 
# svg(paste0("results/Manhathan/Fig4_",geneName2,".svg"), width = 15, height = 13)
# fancyplots2
# dev.off()
#for PNG open in inkscape then save it as png.
