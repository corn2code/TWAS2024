### Figure 2

#TWAS all together

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
library(RColorBrewer)
library(ggpubr)
#install.packages("ggrastr")
library(ggrastr)


setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/out.TWAS")
#list.files(pattern = "HundredKernelMassGrams")

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

index <- fread("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/eGWAS/indexGenes.csv", data.table = F)
index <- index[which(!index$Name=="Michael"),] #remove michael's gene

trait <- "Anthesis.sp.NE"
gwas.datTWAS <- fread(paste0("TWAS.CMLM_",trait,".csv"), data.table = F)


nCHR <- length(unique(gwas.datTWAS$Chromosome))
gwas.datTWAS$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.datTWAS$Chromosome))){
  nbp[i] <- max(gwas.datTWAS[gwas.datTWAS$Chromosome == i,]$`Position.`)
  gwas.datTWAS[gwas.datTWAS$Chromosome == i,"BPcum"] <- gwas.datTWAS[gwas.datTWAS$Chromosome == i,"Position."] + s
  s <- s + nbp[i]
}

axis.set <- gwas.datTWAS %>% 
  group_by(Chromosome) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

#IDs for genes with FDR less than 0.05
fdr.treshold <- -log10(0.05)

ID.fdr <- list(gwas.datTWAS$SNP[-log10(gwas.datTWAS$FDR)>fdr.treshold])

hits.fdr <- index[which(index$ID %in% unlist(ID.fdr)),] #4 or 10


gwas.datTWAS <- merge(gwas.datTWAS, hits.fdr, by = 1, all.x = T)

gTWAS.anthesis.NE <- ggplot(data=gwas.datTWAS, aes(BPcum, -log10(FDR), label=Name, colour=factor(Chromosome, levels = c(1:10)))) + 
  rasterise(geom_point(size = 2),dpi=600) + 
  geom_hline(yintercept = -log10(0.05), linetype=2) + 
  geom_text( vjust="inward",hjust="inward") +
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  annotate("text", label="Anthesis Nebraska", y=19.5, x=100, size=5,hjust="inward") +
  scale_x_continuous(label = axis.set$Chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none", axis.title.x=element_blank()) + 
  ylab(expression(-log[10](FDR))) + 
  xlab("Chromosome") +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(1, 20, 2))
#+ labs(title="", subtitle=paste0(trait))

gTWAS.anthesis.NE


######################################################33


trait <- "Silking.sp.NE"
gwas.datTWAS <- fread(paste0("TWAS.CMLM_",trait,".csv"), data.table = F)


nCHR <- length(unique(gwas.datTWAS$Chromosome))
gwas.datTWAS$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.datTWAS$Chromosome))){
  nbp[i] <- max(gwas.datTWAS[gwas.datTWAS$Chromosome == i,]$`Position.`)
  gwas.datTWAS[gwas.datTWAS$Chromosome == i,"BPcum"] <- gwas.datTWAS[gwas.datTWAS$Chromosome == i,"Position."] + s
  s <- s + nbp[i]
}

axis.set <- gwas.datTWAS %>% 
  group_by(Chromosome) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

#IDs for genes with FDR less than 0.05
fdr.treshold <- -log10(0.05)

ID.fdr <- list(gwas.datTWAS$SNP[-log10(gwas.datTWAS$FDR)>fdr.treshold])

hits.fdr <- index[which(index$ID %in% unlist(ID.fdr)),] #4 or 10


gwas.datTWAS <- merge(gwas.datTWAS, hits.fdr, by = 1, all.x = T)

gTWAS.silking.NE <- ggplot(data=gwas.datTWAS, aes(BPcum, -log10(FDR), label=Name, colour=factor(Chromosome, levels = c(1:10)))) + 
  rasterise(geom_point(size = 2),dpi=600) + 
  geom_hline(yintercept = -log10(0.05), linetype=2) + 
  geom_text( vjust="inward",hjust="inward") +
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  annotate("text", label="Silking Nebraska", y=19.5, x=100, size=5,hjust="inward") +
  scale_x_continuous(label = axis.set$Chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none", axis.title.x=element_blank()) + 
  ylab(expression(-log[10](FDR))) + 
  xlab("Chromosome") +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(1, 20, 2))
#+ labs(title="", subtitle=paste0(trait))

gTWAS.silking.NE



#################################

trait <- "Anthesis.sp.MI"
gwas.datTWAS <- fread(paste0("TWAS.CMLM_",trait,".csv"), data.table = F)


nCHR <- length(unique(gwas.datTWAS$Chromosome))
gwas.datTWAS$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.datTWAS$Chromosome))){
  nbp[i] <- max(gwas.datTWAS[gwas.datTWAS$Chromosome == i,]$`Position.`)
  gwas.datTWAS[gwas.datTWAS$Chromosome == i,"BPcum"] <- gwas.datTWAS[gwas.datTWAS$Chromosome == i,"Position."] + s
  s <- s + nbp[i]
}

axis.set <- gwas.datTWAS %>% 
  group_by(Chromosome) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

#IDs for genes with FDR less than 0.05
fdr.treshold <- -log10(0.05)

ID.fdr <- list(gwas.datTWAS$SNP[-log10(gwas.datTWAS$FDR)>fdr.treshold])

hits.fdr <- index[which(index$ID %in% unlist(ID.fdr)),] #4 or 10


gwas.datTWAS <- merge(gwas.datTWAS, hits.fdr, by = 1, all.x = T)

gTWAS.Anthesis.MI <- ggplot(data=gwas.datTWAS, aes(BPcum, -log10(FDR), label=Name, colour=factor(Chromosome, levels = c(1:10)))) + 
  rasterise(geom_point(size = 2),dpi=600) + 
  geom_hline(yintercept = -log10(0.05), linetype=2) + 
  geom_text( vjust="inward",hjust="inward") +
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  annotate("text", label="Anthesis Michigan", y=19.5, x=100, size=5, hjust="inward") +
  scale_x_continuous(label = axis.set$Chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](FDR))) + 
  xlab("Chromosome") +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(1, 20, 2))
#+ labs(title="", subtitle=paste0(trait))

gTWAS.Anthesis.MI



#################################

trait <- "Silking.sp.MI"
gwas.datTWAS <- fread(paste0("TWAS.CMLM_",trait,".csv"), data.table = F)


nCHR <- length(unique(gwas.datTWAS$Chromosome))
gwas.datTWAS$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.datTWAS$Chromosome))){
  nbp[i] <- max(gwas.datTWAS[gwas.datTWAS$Chromosome == i,]$`Position.`)
  gwas.datTWAS[gwas.datTWAS$Chromosome == i,"BPcum"] <- gwas.datTWAS[gwas.datTWAS$Chromosome == i,"Position."] + s
  s <- s + nbp[i]
}

axis.set <- gwas.datTWAS %>% 
  group_by(Chromosome) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

#IDs for genes with FDR less than 0.05
fdr.treshold <- -log10(0.05)

ID.fdr <- list(gwas.datTWAS$SNP[-log10(gwas.datTWAS$FDR)>fdr.treshold])

hits.fdr <- index[which(index$ID %in% unlist(ID.fdr)),] #4 or 10


gwas.datTWAS <- merge(gwas.datTWAS, hits.fdr, by = 1, all.x = T)


gTWAS.Silking.MI <- ggplot(data=gwas.datTWAS, aes(BPcum, -log10(FDR), label=Name, colour=factor(Chromosome, levels = c(1:10)))) + 
  rasterise(geom_point(size = 2),dpi=600) + 
  geom_hline(yintercept = -log10(0.05), linetype=2) + 
  geom_text(vjust="inward",hjust="inward") +
  scale_color_manual(values = c('#03193F','#28708C','#BF930F','#0f3bbf','#295E52','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')) + 
  annotate("text", label="Silking Michigan", y=19.5, x=100, size=5, hjust="inward") +
  scale_x_continuous(label = axis.set$Chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](FDR))) + 
  xlab("Chromosome") +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(1, 20, 2))
#+ labs(title="", subtitle=paste0(trait))

gTWAS.Silking.MI


#############################################
#scatterplot

## correlation FDR
#silking


setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/out.TWAS")

#load data
TWAS.MI.silking <- fread("TWAS.CMLM_Silking.sp.MI.csv", data.table = F)
TWAS.NE.silking <- fread("TWAS.CMLM_Silking.sp.NE.csv", data.table = F)

TWAS.MI.silking$fdr.log10 <- -log10(TWAS.MI.silking$FDR)
TWAS.NE.silking$fdr.log10 <- -log10(TWAS.NE.silking$FDR)


#merge data
silking <- merge(TWAS.NE.silking[,c(1,11)], TWAS.MI.silking[,c(1,11)], by = 1)
colnames(silking) <- c("gene", "Nebraska", "Michigan")

treshold <- -log10(0.05)

silk.low <- silking[which(!silking$Nebraska > treshold | !silking$Michigan > treshold),]
silk.high <- silking[which(silking$Nebraska > treshold | silking$Michigan > treshold ),]



silk.scatter <- ggplot() +
  rasterise(geom_point(data= silk.low, aes(x=Nebraska, y=Michigan, fill="grey"), size = 3, shape = 21),dpi=600) +
  rasterise(geom_point(data= silk.high, aes(x=Nebraska, y=Michigan, fill="orange3"), size = 3, shape = 21),dpi=600) +
  ylim(0,20) +
  xlim(0,20) +
  xlab("Nebraska (-log10 (FDR))") + 
  ylab("Michigan (-log10 (FDR))") + 
  scale_fill_manual(values = c("grey", "orange3")) +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  geom_vline(xintercept = -log10(0.05), linetype=2) +
  annotate("text", label="Silking", y=19.5, x=1.5, size=5,hjust = 0) +
  theme(legend.position = "none")
silk.scatter

####
#Anthesis
TWAS.MI.anthesis <- fread("TWAS.CMLM_Anthesis.sp.MI.csv", data.table = F)
TWAS.NE.anthesis <- fread("TWAS.CMLM_Anthesis.sp.NE.csv", data.table = F)

TWAS.MI.anthesis$fdr.log10 <- -log10(TWAS.MI.anthesis$FDR)
TWAS.NE.anthesis$fdr.log10 <- -log10(TWAS.NE.anthesis$FDR)


#merge data
anthesis <- merge(TWAS.NE.anthesis[,c(1,11)], TWAS.MI.anthesis[,c(1,11)], by = 1)
colnames(anthesis) <- c("gene", "Nebraska", "Michigan")

treshold <- -log10(0.05)

anthesis.low <- anthesis[which(anthesis$Nebraska < treshold | anthesis$Michigan < treshold),]
anthesis.high <- anthesis[which(!anthesis$Nebraska < treshold | !anthesis$Michigan < treshold ),]



anthesis.scatter <- ggplot() +
  rasterise(geom_point(data= anthesis.low, aes(x=Nebraska, y=Michigan, fill="grey"), size = 3, shape = 21),dpi=600) +
  rasterise(geom_point(data= anthesis.high, aes(x=Nebraska, y=Michigan, fill="orange3"), size = 3, shape = 21),dpi=600) +
  ylim(0,20) +
  xlim(0,20) +
  xlab("Nebraska (-log10 (FDR))") + 
  ylab("Michigan (-log10 (FDR))") + 
  scale_fill_manual(values = c("grey", "orange3")) +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  geom_vline(xintercept = -log10(0.05), linetype=2) +
  annotate("text", label="Anthesis", y=19.5, x=1.5, size=5,hjust = 0) +
  theme(legend.position = "none")
anthesis.scatter


####
# Ven diagram
#library(VennDiagram)
#install.packages("ggvenn")
library("ggvenn")


hits <- merge(silk.high, anthesis.high, by=1, all = T, suffixes = c("_silk","_anthesis"))
hits[!hits > treshold] <- NA
table(hits$Nebraska_silk)

genes.desc <- fread("/work/schnablelab/vladimir/TWAS2023.2/B73_geneModels_v5.csv", data.table = F)
colnames(genes.desc)
hits.desc.2 <- genes.desc[c(1,2,3)]
dta <- merge(hits.desc.2, hits, by =1)

index <- fread("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/eGWAS/indexGenes.csv", data.table = F)
dta2 <- merge(index, dta, by =1)


# write.csv(dta2,
#           "hits_names_pos.csv",
#           row.names = F)


Nebraska_silk <- hits$gene[!is.na(hits$Nebraska_silk)]
Michigan_silk <- hits$gene[!is.na(hits$Michigan_silk)]
Nebraska_anthesis <- hits$gene[!is.na(hits$Nebraska_anthesis)]
Michigan_anthesis <- hits$gene[!is.na(hits$Michigan_anthesis)]


genes.list <- list(`Silking\n(NE)`=Nebraska_silk,
                   `Anthesis (NE)`= Nebraska_anthesis,
                   `Anthesis (MI)`=Michigan_anthesis,
                   `Silking\n(MI)`=Michigan_silk)

venDiagram <- ggvenn(genes.list, stroke_size = 1, set_name_size = 6, show_percentage = F, text_size = 7) +
  theme(legend.box.spacing = unit(0, "cm"),
        plot.margin = unit(c(0,0,0,0), "cm"))
venDiagram


gt <- ggplot_gtable(ggplot_build(venDiagram))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid::grid.draw(gt)


######################################
plots <- align_plots(gTWAS.anthesis.NE, gTWAS.Anthesis.MI, gt, align = 'v', axis = 'l')
top <- plot_grid(plots[[1]], gTWAS.silking.NE, rel_widths = c(1, 1),
                        nrow = 1, labels = c('A', 'B'),hjust = 0, vjust = 1, label_size = 18)
middle <- plot_grid(plots[[2]], gTWAS.Silking.MI, rel_widths = c(1,1),
                         nrow = 1, labels = c('C','D'),hjust = 0, vjust = 1, label_size = 18)
bottom <- plot_grid(plots[[3]], anthesis.scatter, silk.scatter, rel_widths = c(.9,.9,.9),
                    nrow = 1, labels = c('E','F','G'),hjust = 0, vjust = 1, label_size = 18)

fancyplots <- plot_grid(top,middle,bottom, ncol = 1, hjust = 0, vjust = 1, label_size = 18)
fancyplots 

ggsave("Fig2.svg", width = 15, height = 13 ) # to save text as editable


