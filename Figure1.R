#Figure 1

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

#setwd("/work/schnablelab/vladimir/TWAS2023.2/")

setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA")

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5),legend.title=element_text(size=10), 
             legend.text=element_text(size=9))

# #load the new phe data
# order.s <- read.table("sampling_693.order.csv", head = TRUE, sep = ",")
# coor <- read.table("coordenates.693.filtered.csv", head = TRUE, sep = ",")
# dta <- merge(order.s, coor, by =1)
# colnames(dta)
# dta <- dta[,-2] #duplicated column
# colnames(dta)
# colnames(dta)[2] <- "order"
# 
# ##
# popdata <- read.csv("/work/schnablelab/vladimir/TWAS2023.2/MarcinTPJ2023_TableS1.csv", check.names = F)
# colnames(popdata)[1] <- "taxa"
# unique(popdata$Group)
# unique(popdata$Group_As_On_Fig2)
# 
# phe.genes.pop <- merge(popdata[,c(1,3)], dta, by="taxa") # Group
# 
# phe.genes.pop$color <- "NA"
# phe.genes.pop$color[which(phe.genes.pop$Group == "SS")] <- "darkgoldenrod1"
# phe.genes.pop$color[which(phe.genes.pop$Group == "IDT")] <- "indianred"
# phe.genes.pop$color[which(phe.genes.pop$Group == "Mixed")] <- "red1"
# phe.genes.pop$color[which(phe.genes.pop$Group == "NSS")] <- "orange2"
# phe.genes.pop$color[which(phe.genes.pop$Group == "Broad origin-public")] <- "lightcyan"
# phe.genes.pop$color[which(phe.genes.pop$Group == "Sweet corn")] <- "deeppink"
# phe.genes.pop$color[which(phe.genes.pop$Group == "Popcorn")] <- "cornflowerblue"
# phe.genes.pop$color[which(phe.genes.pop$Group == "Tropical")] <- "cornsilk2"
# phe.genes.pop$color[which(phe.genes.pop$Group == "Landrace")] <- "darkblue"
# 
# 
# 
# 
# colnames(phe.genes.pop)
# colnames(phe.genes.pop)[7:16] <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
# colnames(phe.genes.pop)
# write.csv(phe.genes.pop, "PCA.data.csv", row.names = F)

####### PCA

phe.genes.pop <- read.table("PCA.data.csv", head = TRUE, sep = ",")

phe.genes.pop$Group <- gsub("Mixed", "Others" , phe.genes.pop$Group)
phe.genes.pop$Group <- gsub("Broad origin-public", "Others" , phe.genes.pop$Group)
phe.genes.pop$Group <- gsub("Landrace", "Others" , phe.genes.pop$Group)

unique(phe.genes.pop$Group)

# Sort the vector and remove NAs
sorted_vector <- na.omit(sort(phe.genes.pop$PC1))
sorted_vector2 <- na.omit(sort(phe.genes.pop$PC2))
sorted_vector3 <- na.omit(sort(phe.genes.pop$PC3))

# Apply findInterval on the sorted vector
phe.genes.pop$order.PC1 <- findInterval(phe.genes.pop$PC1, sorted_vector)
phe.genes.pop$order.PC2 <- findInterval(phe.genes.pop$PC2, sorted_vector2)
phe.genes.pop$order.PC3 <- findInterval(phe.genes.pop$PC3, sorted_vector3)


groups <- c("SS","IDT", "NSS", "Sweet corn", "Popcorn", "Tropical", "Others")
groups.col <- c("darkgoldenrod1", "indianred", "deeppink", "cornflowerblue", "cornsilk2", "darkblue", "white")


PCA <- ggplot() +
  rasterise(geom_point(data= phe.genes.pop[phe.genes.pop$Group=="Others",], aes(x=PC1, y=PC2, fill=Group), size = 3, shape = 21),dpi=600) +
  rasterise(geom_point(data= phe.genes.pop[phe.genes.pop$Group!="Others",], aes(x=PC1, y=PC2, fill=Group), size = 3, shape = 21),dpi=600) +
  ylim(-150,210) +
  xlim(-150,210) +
  ylab("Principal Component 2 (PC2)") + 
  xlab("Principal Component 1 (PC1)") + 
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups,
                    name="")+ theme(legend.position="top", legend.box.spacing = unit(0, "cm"),
                                                  plot.margin = unit(c(0,.5,0,.5), "cm"))
PCA


### boxplot by subgroups

phe.genes.pop$Group <- as.character(phe.genes.pop$Group)
phe.genes.pop$Group <- factor(phe.genes.pop$Group, levels=groups)


phe.genes.popPC <- phe.genes.pop %>%
  select(c("Group", "PC1", "PC2", "PC3")) %>%
  pivot_longer(!Group, names_to = "Component", values_to = "value")



#install.packages("ggsignif")
#library(ggsignif)
library(agricolae)

t <- phe.genes.popPC %>%
  filter(Component == "PC1")
hsdt=HSD.test(aov(value~Group,data=t), c("Group"), group=T)

hsdt$groups

hsd=HSD.test(aov(value~Group+Component,data=phe.genes.popPC), c("Group","Component"), group=T)

hsdG <- as.data.frame(hsd$groups)
hsdG
phe.genes.popPC$tukey <- NA


for (i in 1:nrow(hsdG)) {
  phe.genes.popPC$tukey[which(phe.genes.popPC$Group == strsplit(rownames(hsdG), ":")[[i]][1] &
                                phe.genes.popPC$Component == strsplit(rownames(hsdG), ":")[[i]][2])][1] <- hsdG$groups[i]
}



bp <- ggplot(phe.genes.popPC, aes(Group, value, fill=Group)) +
  geom_boxplot(width=0.7, outlier.shape = NA) +
  rasterise(geom_jitter(shape=21, color="black", size=2, alpha=0.9, width = .1),dpi=600) +
  scale_y_continuous(breaks = seq(-150, 250, by = 50)) +
  scale_fill_manual(labels = groups, 
                    values = groups.col, 
                    breaks=groups) + 
  theme(legend.position="none") + 
  theme(axis.text.x = element_text(angle = 21, hjust = 1, size = 13), strip.text = element_text(size = 19), panel.border = element_rect(colour = "black", fill=NA),
        strip.background = element_blank(), strip.placement='outside') +
  facet_grid(rows = vars(Component), scales = "free", switch="both") + labs(x = NULL, y = NULL)+
  geom_text(aes(x = Group, label = tukey, y = min(value)), size = 5)

bp







############## Scatter plot PC values versus order of collection

colnames(phe.genes.pop)

phe.genes.popPC.order <- phe.genes.pop %>%
  select(c("order", "PC1", "PC2", "PC3")) %>%
  pivot_longer(!order, names_to = "Component", values_to = "value")



bp.order <- ggplot(phe.genes.popPC.order, aes(order, value, fill=Component)) +
  geom_point(size = 3, shape = 21) +
  scale_y_continuous(breaks = seq(-150, 250, by = 50)) +
  theme(legend.position="none") + 
  theme(panel.border = element_rect(colour = "black", fill=NA), strip.text = element_text(size = 19), strip.background = element_blank(), strip.placement='outside') +
  facet_grid(rows = vars(Component), scales = "free", switch="both") + labs(x = NULL, y = NULL) + scale_fill_manual(values = c("blue1", "green4", "yellow3"))+ 
  stat_poly_line(method=lm,  linetype="dashed", se=FALSE, linewidth = 2,
                 color="black", fill="blue") +
  stat_poly_eq(size = 5,rr.digits = 3) + xlab("order of collection")
bp.order


cor(phe.genes.pop$PC1,phe.genes.pop$order)^2

summary(lm(phe.genes.pop$PC3~phe.genes.pop$order))

# PC1 p-value: < 2.2e-16
# PC2 p-value: 0.06713
# PC2 p-value: 0.1378

#### heatmap

#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/

# theme_set(theme_classic(base_size = 19))
# theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
#              plot.title = element_text(hjust = 0.5), plot.subtitle=element_text(hjust=0.5),
#              panel.background = element_rect(fill = "black"))


#to introduce empty values corresponding to column 5 which was duplicated
val5 <- data.frame(matrix(NA,    # Create empty data frame
                          nrow = 30, # rows with NA in the final plot each represented one plant
                          ncol = 20)) # just to merge with the data we are using we have 13 columns

colnames(val5) <- colnames(phe.genes.pop)
val5$Row <- c(1:30) # rows with NA in the final plot each represented one plant
val5$Column <- as.factor(5) # column not sampled


#to introduce empty values corresponding to column 7 which was duplicated
val7 <- data.frame(matrix(NA,    # Create empty data frame
                          nrow = 30,
                          ncol = 20))

colnames(val7) <- colnames(phe.genes.pop)
val7$Row <- c(1:30)
val7$Column <- as.factor(7)

colnames(phe.genes.pop)

dx <- rbind(phe.genes.pop, val5, val7)
#rm(dx)
dx <- dx[order(dx$Column), ]
str(dx)

# New sorting order
desired_order <- c("1","2","3","4","5","6","7","8","9","10",
                   "11","12","13","14","15","16","17","18",
                   "19","20","21","22","23","24","25","26","27","28")
# Re-order the levels
dx$Column <- factor( as.character(dx$Column), levels=desired_order )


str(dx)
dx$Row <- as.factor(dx$Row)
dx$Column <- as.factor(dx$Column)
str(dx)

d1 <- ggplot(dx, aes(x=Row, y=Column, fill=PC1))+
  geom_tile(aes(fill = PC1)) + #colour="black", linewidth=0
  scale_fill_gradient(low="white", high="#0000FF") +
  ylab("PC1") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), panel.background = element_rect(fill = "black"), legend.key.size = unit(.8,"line"),
        legend.box.spacing = unit(0, "cm"), plot.margin = unit(c(0,.5,0,.5), "cm"))
d1



d2 <- ggplot(dx, aes(x=Row, y=Column, fill=PC2))+
  geom_tile(colour="#006600", linewidth=0) + #colour="black", linewidth=0
  scale_fill_gradient(low="white", high="#006600") + 
  ylab("PC2") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),  panel.background = element_rect(fill = "black"), legend.key.size = unit(.8,"line"),
        legend.box.spacing = unit(0, "cm"), plot.margin = unit(c(0,.5,0,.5), "cm"))
d2

d3 <- ggplot(dx, aes(x=Row, y=Column, fill=PC3))+
  geom_tile(colour="#CC9900", linewidth=0) + #colour="black", linewidth=0
  scale_fill_gradient(low="white", high="#CC9900") +
  ylab("PC3") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),  panel.background = element_rect(fill = "black"), legend.key.size = unit(.8,"line"),
        legend.box.spacing = unit(0, "cm"), plot.margin = unit(c(0,.5,0,.5), "cm"))

d3

spatial <- ggarrange(d1, d2, d3, 
          #labels = c("A", "B", "C"),
          ncol = 1, nrow = 3, widths = c(1, 1, 1))


# Finally put all together manually and adjust 
plots <- align_plots(PCA, bp.order, align = 'v', axis = 'l')

top_row <- plot_grid(plots[[1]], bp , rel_widths = c(1.1, .7),
                        nrow = 1, labels = c('A', 'B'),hjust = 0, vjust = 1, label_size = 18)
bottom_row <- plot_grid(plots[[2]], spatial, rel_widths = c(1.1, .5),
                        nrow = 1, labels = c('C', 'D'), hjust = 0, vjust = 1, label_size = 18, align = "h", axis = "b")


fancyplots <- plot_grid(top_row, bottom_row, ncol = 1,hjust = 0, vjust = 1, label_size = 18)
fancyplots


#ggsave(plot=fancyplots2, paste0("results/Manhathan/Fig1.png"), width = 15, height = 13)
#ggsave(plot=fancyplots2, paste0("results/Manhathan/Fig1.svg"), width = 15, height = 13)

ggsave("Fig1.svg", width = 10 , height = 11 ) # change hight


#