## zcn8 expression 
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(cowplot)
library(gridGraphics)
library(ggpmisc)
library(Matrix)
library(gcookbook)
library(RColorBrewer)
library(ggpubr)
#install.packages("ggrastr")
library(ggrastr)

setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/out.TWAS")
genes <- read.csv("zcnGenes.csv")

head(genes)
genes$Gene

counts <- fread("../counts.NE2020.693.filtered.txt", data.table = F)
row.names(counts) <- counts$taxa


counts1 <- counts[,c("taxa",genes$Gene)]

phe <- read.table("../pheno_693.txt", head = TRUE) # BLUES.NE2020_corrected.csv  # BLUES.NE2021_corrected.csv
colnames(phe)

dta <- merge(counts1, phe[,c("taxa","Silking.sp.NE", "Anthesis.sp.NE")], by ="taxa")



for (i in genes$Gene) {
  colnames(dta)[which(colnames(dta)==i)] <- genes$symbol[which(genes$Gene==i)]
}

#if we want to use face wrap we can do this and split dataset based on genotypes then, use rbind.
# counts1.gg <- counts1 %>%
#   pivot_longer(!taxa, names_to = "Gene", values_to = "Expression")

colnames(dta)

colnames(dta)[8:9] <- gsub(".sp.NE", "",  colnames(dta)[8:9])

colnames(dta)

##
popdata <- read.csv("/work/schnablelab/vladimir/TWAS2023.2/MarcinTPJ2023_TableS1.csv", check.names = F)
colnames(popdata)[1] <- "taxa"
unique(popdata$Group)

popdata$Group <- gsub("Mixed", "Others" , popdata$Group)
popdata$Group <- gsub("Broad origin-public", "Others" , popdata$Group)
popdata$Group <- gsub("Landrace", "Others" , popdata$Group)


unique(popdata$Group_As_On_Fig2)
phe.TPM.genes.pop <- merge(popdata[,c(1,3)], dta, by="taxa")

groups <- c("SS","IDT", "NSS", "Sweet corn", "Popcorn", "Tropical", "Others")
groups.col <- c("darkgoldenrod1", "indianred", "deeppink", "cornflowerblue", "cornsilk2", "darkblue", "white")

unique(phe.TPM.genes.pop$Group)
#unique(phe.TPM.genes.pop$Group_As_On_Fig2)

#j <- ncol(phe.TPM.genes.pop)-1

colnames(phe.TPM.genes.pop)



# zcn7 <- ggplot() +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Silking, y=zcn7, fill=Group), size = 3, shape = 21),dpi=600) +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Silking, y=zcn7, fill=Group), size = 3, shape = 21),dpi=600) +
#   xlim(60,85)+
#   ylab("zcn7 (TPM)") + 
#   xlab("Days to silking") + 
#   scale_fill_manual(labels = groups, 
#                     values = groups.col, 
#                     breaks=groups,
#                     name="sub-population") + theme(legend.position="none")+ 
#   stat_poly_line(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn7), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
#                  color="black", fill="blue") +
#   stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn7),size = 5)
# zcn7


zcn8 <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Silking, y=zcn8, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Silking, y=zcn8, fill=Group), size = 5, shape = 21),dpi=600) +
  xlim(62,85)+
  ylab("zcn8 (TPM)") + 
  xlab("Days to silking.NE") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none")+ 
  stat_poly_line(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn8), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn8),size = 5)
zcn8

############

# zcn12 <- ggplot() +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Silking, y=zcn12, fill=Group), size = 3, shape = 21),dpi=600) +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Silking, y=zcn12, fill=Group), size = 3, shape = 21),dpi=600) +
#   xlim(60,85)+
#   ylab("zcn12 (TPM)") + 
#   xlab("Days to silking") + 
#   scale_fill_manual(labels = groups, 
#                     values = groups.col, 
#                     breaks=groups,
#                     name="sub-population") + theme(legend.position="none")+ 
#   stat_poly_line(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn12), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
#                  color="black", fill="blue") +
#   stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn12),size = 5)
# zcn12

############

# zcn14 <- ggplot() +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Silking, y=zcn14, fill=Group), size = 3, shape = 21),dpi=600) +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Silking, y=zcn14, fill=Group), size = 3, shape = 21),dpi=600) +
#   xlim(60,85)+
#   ylab("zcn14 (TPM)") + 
#   xlab("Days to silking") + 
#   scale_fill_manual(labels = groups, 
#                     values = groups.col, 
#                     breaks=groups,
#                     name="sub-population") + theme(legend.position="none")+ 
#   stat_poly_line(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn14), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
#                  color="black", fill="blue") +
#   stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn14),size = 5)
# zcn14

############

# zcn15 <- ggplot() +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Silking, y=zcn15, fill=Group), size = 3, shape = 21),dpi=600) +
#   rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Silking, y=zcn15, fill=Group), size = 3, shape = 21),dpi=600) +
#   xlim(60,85)+
#   ylab("zcn15 (TPM)") + 
#   xlab("Days to silking") + 
#   scale_fill_manual(labels = groups, 
#                     values = groups.col, 
#                     breaks=groups,
#                     name="sub-population") + theme(legend.position="none")+ 
#   stat_poly_line(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn15), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
#                  color="black", fill="blue") +
#   stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Silking, y=zcn15),size = 5)
# zcn15


############

zcn26 <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Anthesis, y=zcn26, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Anthesis, y=zcn26, fill=Group), size = 5, shape = 21),dpi=600) +
  xlim(62,85)+
  ylab("zcn26 (TPM)") + 
  xlab("Days to anthesis.NE") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none")+ 
  stat_poly_line(data= phe.TPM.genes.pop, aes(x=Anthesis, y=zcn26), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Anthesis, y=zcn26),size = 5)
zcn26



######################################
# plots <- align_plots(zcn7, zcn14, align = 'v', axis = 'l')
# top <- plot_grid(plots[[1]], zcn8, zcn12, rel_widths = c(1, 1, 1),
#                  nrow = 1, labels = c('A', 'B','D'),hjust = 0, vjust = 1, label_size = 18)
# middle <- plot_grid(plots[[2]], zcn15, zcn26, rel_widths = c(1,1,1),
#                     nrow = 1, labels = c('E','F','G'),hjust = 0, vjust = 1, label_size = 18)
# # bottom <- plot_grid(plots[[3]], anthesis.scatter, silk.scatter, rel_widths = c(.9,.9,.9),
# #                     nrow = 1, labels = c('E','F','G'),hjust = 0, vjust = 1, label_size = 18)
# 
# fancyplots <- plot_grid(top,middle, ncol = 1, hjust = 0, vjust = 1, label_size = 18)
# fancyplots 
# 
# ggsave("Fig3.svg", width = 15, height = 13 ) # to save text as editable
# 
# 
# 



########################################################




genes <- read.csv("CoolGenes.csv")

counts <- fread("../counts.NE2020.693.filtered.txt", data.table = F)
row.names(counts) <- counts$taxa


counts1 <- counts[,c("taxa",genes$Gene)]

phe <- read.table("../pheno_693.txt", head = TRUE) # BLUES.NE2020_corrected.csv  # BLUES.NE2021_corrected.csv
colnames(phe)

dta <- merge(counts1, phe[,c("taxa","Silking.sp.NE", "Anthesis.sp.NE","Anthesis.sp.MI")], by ="taxa")



for (i in genes$Gene) {
  colnames(dta)[which(colnames(dta)==i)] <- genes$symbol[which(genes$Gene==i)]
}

#if we want to use face wrap we can do this and split dataset based on genotypes then, use rbind.
# counts1.gg <- counts1 %>%
#   pivot_longer(!taxa, names_to = "Gene", values_to = "Expression")

colnames(dta)

#colnames(dta)[8:9] <- gsub(".sp.NE", "",  colnames(dta)[8:9])

colnames(dta)

##
popdata <- read.csv("/work/schnablelab/vladimir/TWAS2023.2/MarcinTPJ2023_TableS1.csv", check.names = F)
colnames(popdata)[1] <- "taxa"
unique(popdata$Group)

popdata$Group <- gsub("Mixed", "Others" , popdata$Group)
popdata$Group <- gsub("Broad origin-public", "Others" , popdata$Group)
popdata$Group <- gsub("Landrace", "Others" , popdata$Group)


unique(popdata$Group_As_On_Fig2)
phe.TPM.genes.pop <- merge(popdata[,c(1,3)], dta, by="taxa")

groups <- c("SS","IDT", "NSS", "Sweet corn", "Popcorn", "Tropical", "Others")
groups.col <- c("darkgoldenrod1", "indianred", "deeppink", "cornflowerblue", "cornsilk2", "darkblue", "white")

unique(phe.TPM.genes.pop$Group)
#unique(phe.TPM.genes.pop$Group_As_On_Fig2)

#j <- ncol(phe.TPM.genes.pop)-1

colnames(phe.TPM.genes.pop)



late <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Anthesis.sp.NE, y=late, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Anthesis.sp.NE, y=late, fill=Group), size = 5, shape = 21),dpi=600) +
  xlim(62,85)+
  ylab("late (TPM)") + 
  xlab("Days to Anthesis.NE") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none")+ 
  stat_poly_line(data= phe.TPM.genes.pop, aes(x=Anthesis.sp.NE, y=late), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Anthesis.sp.NE, y=late),size = 5)
late



####################

arf34 <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Anthesis.sp.MI, y=arf34, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Anthesis.sp.MI, y=arf34, fill=Group), size = 5, shape = 21),dpi=600) +
  xlim(55,80)+
  ylab("arf34 (TPM)") + 
  xlab("Days to Anthesis.MI") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none")+ 
  stat_poly_line(data= phe.TPM.genes.pop, aes(x=Anthesis.sp.MI, y=arf34), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Anthesis.sp.MI, y=arf34),size = 5)
arf34




####################

zmm4 <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Silking.sp.NE, y=zmm4, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Silking.sp.NE, y=zmm4, fill=Group), size = 5, shape = 21),dpi=600) +
  xlim(62,85)+
  ylab("zmm4 (TPM)") + 
  xlab("Days to Silking.NE") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none")+ 
  stat_poly_line(data= phe.TPM.genes.pop, aes(x=Silking.sp.NE, y=zmm4), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Silking.sp.NE, y=zmm4),size = 5)
zmm4

####################

zfp30 <- ggplot() +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group=="Others",], aes(x=Silking.sp.NE, y=zfp30, fill=Group), size = 5, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.TPM.genes.pop[phe.TPM.genes.pop$Group!="Others",], aes(x=Silking.sp.NE, y=zfp30, fill=Group), size = 5, shape = 21),dpi=600) +
  xlim(62,85)+
  ylab("zfp30 (TPM)") + 
  xlab("Days to Silking.NE") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="sub-population") + theme(legend.position="none")+ 
  stat_poly_line(data= phe.TPM.genes.pop, aes(x=Silking.sp.NE, y=zfp30), method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(data= phe.TPM.genes.pop, aes(x=Silking.sp.NE, y=zfp30),size = 5)
zfp30


######################################
# plots <- align_plots(zmm4, zfp30, align = 'v', axis = 'l')
# top <- plot_grid(plots[[1]], zcn8, zcn26, rel_widths = c(1, 1, 1),
#                  nrow = 1, labels = c('A', 'B','D'),hjust = 0, vjust = 1, label_size = 18)
# middle <- plot_grid(plots[[2]], arf34, late, rel_widths = c(1,1,1),
#                     nrow = 1, labels = c('E','F','G'),hjust = 0, vjust = 1, label_size = 18)
# # bottom <- plot_grid(plots[[3]], anthesis.scatter, silk.scatter, rel_widths = c(.9,.9,.9),
# #                     nrow = 1, labels = c('E','F','G'),hjust = 0, vjust = 1, label_size = 18)
# 
# fancyplots <- plot_grid(top,middle, ncol = 1, hjust = 0, vjust = 1, label_size = 18)
# fancyplots 
# 
# ggsave("Fig3.coolGenes.svg", width = 14, height = 7) # to save text as editable

plots2 <- align_plots(zmm4, zcn8, zcn26, zfp30,  arf34, late, align = 'hv', axis = 'l')
top2 <- plot_grid(plots2[[1]], plots2[[2]], plots2[[3]], rel_widths = c(1, 1, 1),
                 nrow = 1, labels = c('A', 'B','C'),hjust = 0, vjust = 1, label_size = 18)
middle2 <- plot_grid(plots2[[4]], plots2[[5]], plots2[[6]], rel_widths = c(1,1,1),
                    nrow = 1, labels = c('D','E','F'),hjust = 0, vjust = 1, label_size = 18)

fancyplots2 <- plot_grid(top2, middle2, ncol = 1, hjust = 0, vjust = 1, label_size = 18, align = "v")
fancyplots2 

ggsave("Fig3.coolGenes.svg", width = 14.5, height = 7) # to save text as editable

