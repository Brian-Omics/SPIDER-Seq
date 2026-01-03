# --------------------
# title: FigureS13 Code
# author: Hu Zheng
# date: 2026-01-01
# --------------------

library(Seurat)
library(tidyverse)
library(caret)
library(Matrix)
library(xgboost)
library(PRROC)
library(scCustomize)
library(cowplot)
library(ggpointdensity)

library(Biorplot)
source('bin/Palettes.R')
source('bin/includes.R')

Adult.Ex <- readRDS('../data/rds/Adult.Ex.rds')
sp.PFC <- readRDS('../data/rds/sp.PFC.rds')
PFC.MERFISH <- readRDS('../data/rds/PFC.MERFISH.rds')
seu.inte <- readRDS('../data/rds/PFC.MERFISH.inte.rds')

Adult.IT.PT.barcode <- subset(Adult.Ex, cells=colnames(Adult.Ex)[which(
  (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "IT") |
    (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "PT" & Adult.Ex$sample == "Adult1")
)])

sp.PFC.Left <- subset(
  sp.PFC,
  cells = colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere=="Left")])

sp.PFC.Left.ITPT.barcode <- subset(sp.PFC, cells = colnames(sp.PFC)[which(
  sp.PFC$ABA_hemisphere=="Left" & sp.PFC$SubType_Layer %in% c("L2/3 IT","L4/5 IT","L5 IT","L6 IT", "L5 PT") & sp.PFC$BC_num>0)])



Barcode <- c('VIS-I','SSp-I','CP-I','AUD-I','RSP-I',
             'BLA-I','ACB-I','ENTl-I','AId-I','ECT-I',
             'ACB-C','PL-C','ECT-C','ENTl-C',
             'BLA-C','CP-C','AId-C','RSP-C',
             'MD-I','RE-I','DR-I','VTA-I','LHA-I','SC-I')
result_all <- readRDS('../data/rds/ML/result_target.rds')
df <- result_all
df <- df[which(df$Experiment != "spatial"),]
df$target <- factor(df$target, levels = Barcode)
df$Experiment <- factor(
  df$Experiment,
  levels = c("transcriptom + spatial","transcriptom","shuffle"))

p1 <-
  ggplot(df, aes(x=ROC_1, y=ROC_2, color=Experiment)) +
  geom_line(linewidth=1) +
  facet_wrap(~target, nrow = 4) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.spacing = unit(1,"lines"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top", legend.title = element_blank()) +
  scale_color_manual(values = c("transcriptom + spatial"="#d73027",
                                "transcriptom"="#4575b4",
                                "spatial"="#ff7f0e",
                                "shuffle"="gray")) +
  labs(x='False positive rate', y='Sensitivity',
       title = "")

ggsave("../pdf/FigureS13/p1.pdf", plot = p1, height = 8, width = 10, units = "in")


df <- df[which(df$Experiment != "shuffle"),]
unique(paste(df$target, df$AUC,df$Experiment))



spatial_MERFISH <- seu.inte@images$slice.44@coordinates[,c("x","y")]
spatial_MERFISH$CCF_metaRegion <- seu.inte@meta.data[rownames(spatial_MERFISH),'CCF_metaRegion']
MERFISH_cells <- rownames(spatial_MERFISH)[which(spatial_MERFISH$y > (-1/5*spatial_MERFISH$x-150) & spatial_MERFISH$CCF_metaRegion %in% c("MOs","ACAd","PL","ILA","DP"))]

rotation=10
x.scale=1.1
y.scale=1.1
x.move= 0.60
y.move = -2.7
spatial_test <- data.frame(
  "x" = spatial_MERFISH[MERFISH_cells,"y"]/1112 + 0.5,
  "y" = -spatial_MERFISH[MERFISH_cells,"x"]/1043 - 2
)
x0 = (spatial_test$x - median(spatial_test$x))
y0 = (spatial_test$y - median(spatial_test$y))
a = rotation / 180 * pi
x1 = x0 * cos(a) - y0 * sin(a)
y1 = x0 * sin(a) + y0 * cos(a)
x1 = x1 * x.scale + x.move
y1 = y1 * y.scale + y.move
spatial_test$x <- x1
spatial_test$y <- y1

transcriptom_test <- seu.inte@reductions$pca@cell.embeddings[MERFISH_cells,1:30]
transcriptom_train <- seu.inte@reductions$pca@cell.embeddings[colnames(sc_seu),1:30]
spatial_train <- sc_seu@meta.data[,c("ML_new","DV_new")]
colnames(spatial_train) <- c("x","y")

df_MERFISH <- spatial_test
df_MERFISH$CCF <- PFC.MERFISH$CCF_metaRegion[which(colnames(PFC.MERFISH) %in% MERFISH_cells)]
df_MERFISH$sample <- "MERFISH"
df_MERFISH$SubType <- PFC.MERFISH$L3_cluster[which(colnames(PFC.MERFISH) %in% MERFISH_cells)]
df_true <- sp.PFC.Left@meta.data[which(sp.PFC.Left$slice=="IT_slice_18"),c('ML_new','DV_new','ABA_metaRegion','SubType_Layer')]
colnames(df_true) <- c('x','y','CCF','SubType')
df_true$sample <- "our"
df <- merge(df_MERFISH, df_true, all = T)
df$y_round <- round(df$y,1)
df$noise <- ""
for (i in unique(df$y_round)){
  L23_x <- df$x[which(df$y_round == i & df$SubType == "L2/3 IT" & df$sample=="our")]
  if (length(L23_x)>0){
    L23_minx <- min(L23_x)
    df$noise[which(df$y_round == i & df$x < L23_minx & df$sample=="our")] <- "noise"
  }
}
df <- df[df$noise != "noise",]

df$CCF <- factor(df$CCF, levels = c("MOs","ACAd","ACAv","PL","ORBm","ILA","DP"))

p2_1 <-
  ggplot(df[which(df$sample=="our"),], aes(x=x,y=y,color=CCF)) +
  geom_point() +
  coord_fixed() +
  scale_color_manual(values = c("MOs"="#2166AC","ACAd"="#67A9CF","ACAv"="#92C5DE","PL"="#D1E5F0","ORBm"="#FDDBC7","ILA"="#EF8A62","DP"="#B2182B")) +
  labs(title = "SPIDER-Seq slice") +
  theme_void() +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20))

ggsave("../pdf/FigureS13/p2_1.pdf", plot = p2_1, height = 6, width = 4, units = "in")



df$CCF <- factor(df$CCF, levels = c("MOs","ACAd","ACAv","PL","ORBm","ILA","DP"))

p2_2 <-
  ggplot(df[which(df$sample=="MERFISH"),], aes(x=x,y=y,color=CCF)) +
  geom_point() +
  coord_fixed() +
  scale_color_manual(values = c("MOs"="#2166AC","ACAd"="#67A9CF","ACAv"="#92C5DE","PL"="#D1E5F0","ORBm"="#FDDBC7","ILA"="#EF8A62","DP"="#B2182B")) +
  labs(title = "MERFISH slice") +
  theme_void() +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20))

ggsave("../pdf/Figure7/p2_2.pdf", plot = p2_2, height = 6, width = 4, units = "in")



p3 <-
  DimPlot_scCustom(seurat_object = seu.inte, group.by = "Sample", reduction = "umap", figure_plot = TRUE, colors_use = c("#f8766d", "#00bfc4"),pt.size = 0.1) +
  coord_fixed()

ggsave("../pdf/FigureS13/p3.pdf", plot = p3, height = 4, width = 6, units = "in")



seu <- subset(seu.inte, cells=colnames(seu.inte)[which(seu.inte$Sample=="scRNAseq")])
colnames(seu@meta.data)[13:36] <- c(
  "ACB-C","ACB-I","AId-C","AId-I","AUD-I","BLA-C","BLA-I","CP-C","CP-I","DR-I",
  "ECT-C","ECT-I","ENTl-C","ENTl-I","LHA-I","MD-I","PL-C","RE-I","RSP-C","RSP-I",
  "SC-I","SSp-I","VIS-I","VTA-I")
seu$first_target <- "none"
seu$first_target[which(seu$BC_num>0)] <- Barcode[apply(seu@meta.data[which(seu$BC_num>0),Barcode], 1, which.max)]
seu$first_target <- factor(seu$first_target, levels = c(Barcode,"none"))

p4 <-
  DimPlot_scCustom(seurat_object = seu, group.by = "first_target", reduction = "umap",
                   colors_use = col_Barcode,
                   figure_plot = TRUE, pt.size = 0.1) +
  coord_fixed()

ggsave("../pdf/FigureS13/p4.pdf", plot = p4, height = 4, width = 6, units = "in")



i=15
module <- sp_Barcode[i]
p5_1 <-
  ggplot() +
  geom_point(pre_Barcode[which(pre_Barcode[,module]<sp_Barcode_thread[i]),],
             mapping=aes(x=x,y=y), size=0.8, color="#00204DFF") +
  geom_point(pre_Barcode[which(pre_Barcode[,module]>=sp_Barcode_thread[i]),],
             mapping=aes(x=x,y=y), size=0.8, color="#FFEA46FF") +
  ggdark::dark_theme_void() +
  labs(title = "") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_fixed() +
  xlim(0,1.5) +
  ylim(-4.3,-1.25)

ggsave(paste("../pdf/FigureS13/p5_1",module,".png",sep=""),
       plot = p5_1,
       height = 4, width = 2, units = "in")



i=15
module <- sp_Barcode[i]
slice <- 'IT_slice_18'

df <- sp.PFC.Left@meta.data[,c("ML_new","DV_new","SubType_Layer",sp_Barcode)]
colnames(df)[1:3] <- c("x","y","SubType")
df <- df[which(sp.PFC.Left$slice==slice),]
df$y_round <- round(df$y,1)
df$noise <- ""
for (i in unique(df$y_round)){
  L23_x <- df$x[which(df$y_round == i & df$SubType == "L2/3 IT")]
  if (length(L23_x)>0){
    L23_minx <- min(L23_x)
    df$noise[which(df$y_round == i & df$x < L23_minx)] <- "noise"
  }
}
df <- df[df$noise != "noise",]

p5_2 <-
  ggplot(df, aes()) +
  geom_point(df[which(df[,module]==0),], mapping=aes(x=x,y=y),
             size=0.8, color="#00204DFF") +
  geom_point(df[which(df[,module]>0),], mapping=aes(x=x,y=y),
             size=0.8, color="#FFEA46FF") +
  ggdark::dark_theme_void() +
  labs(title = "") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_fixed() +
  xlim(0,1.5) +
  ylim(-4.3,-1.25)

ggsave(paste("../pdf/FigureS13/p5_2",module,".png",sep=""),
       plot = p5_2,
       height = 4, width = 2, units = "in")




