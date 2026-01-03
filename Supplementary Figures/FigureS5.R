# --------------------
# title: FigureS5 Code
# author: Hu Zheng
# date: 2026-01-01
# --------------------

library(Seurat)
library(tidyverse)
library(viridis)
source('bin/Palettes.R')

sp.all <- readRDS('../data/rds/sp.all.rds')
sp.PFC <- readRDS('../data/rds/sp.PFC.rds')



seu <- sp.PFC
seu$SubType <- factor(
  seu$SubType,
  levels = c("L2/3_IT_1","L2/3_IT_2","L4/5_IT_1","L4/5_IT_2","L5_IT_1","L5_IT_2",
             "L6_IT_1","L6_IT_2","L5_PT_1","L5_PT_2","L5_NP","L6_CT_1","L6_CT_2"
  ))
features <- c("Ddit4l","Calb1","Otof","Hap1","Dio3","Cux2","Rorb","Tnnc1","Rspo1","Bdnf",
              "Ptn","Sema3d","Cpne7","Deptor","Fstl5","Slc24a2","Oprk1","Penk","Nnat",
              "Pou3f1","S100b","Etv1","Scn4b","Adamts2","Cdh13","Dlk1","Efnb3",
              "Tshz2","Cbln2","Grp","Syt6","Pcp4")
p1 <- DotPlot(
  seu,
  features = features,
  group.by = 'SubType',
  col.min=0, col.max=2, dot.scale = 10) +
  scale_y_discrete(limits=rev) +
  scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4")) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1, size = 17),
        axis.text.y=element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17)
  ) +
  guides(colour = guide_colorbar(title = "Expression"),
         size = guide_legend(title = "Percent")
  )

ggsave("../pdf/FigureS5/p1.pdf", plot = p1, height = 6, width = 18, units = "in")



seu <- sp.all
slice <- 'IT_slice_15'
features <- c("Ddit4l","Calb1","Otof","Hap1","Dio3","Cux2","Rorb","Tnnc1","Rspo1","Bdnf","Ptn","Sema3d","Cpne7","Deptor","Fstl5","Slc24a2","Oprk1","Penk","Nnat","Pou3f1","S100b","Etv1","Scn4b","Adamts2","Cdh13","Dlk1","Efnb3","Tshz2","Cbln2","Grp","Syt6","Pcp4")
gene <- features[2]
limits <- c(0,3)

df <- data.frame(
  X = seu$X,
  Y = seu$Y,
  Zscore = scale(log1p(seu@assays$RNA@counts[gene,]))
)
df <- df[which(seu$slice==slice),]
df$Zscore[df$Zscore<limits[1]] <- limits[1]
df$Zscore[df$Zscore>limits[2]] <- limits[2]
df <- df[order(df$Zscore),]

p2 <-
  ggplot(df,aes(x=X,y=Y)) +
  geom_point(aes(colour=Zscore), size=0.5) +
  scale_color_gradientn(colours = viridis(n = 256, option = "D", direction = 1),
                        limits = limits) +
  ggdark::dark_theme_void() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none") +
  coord_fixed()

ggsave("../pdf/FigureS5/p2_Calb1.png", plot = p2, height = 4, width = 4, units = "in")




