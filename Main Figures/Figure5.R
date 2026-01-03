# --------------------
# title: Figure5 Code
# author: Hu Zheng
# date: 2026-01-01
# --------------------

library(Seurat)
library(tidyverse)
library(hdWGCNA)
library(cowplot)
library(patchwork)
library(enrichR)
library(GeneOverlap)
library(ggpointdensity)
library(viridis)
library(ggrastr)
library(RColorBrewer)

library(Biorplot)
source('bin/Palettes.R')
source('bin/includes.R')

all.Adult <- readRDS('../data/rds/all.Adult.rds')
Adult.Ex <- readRDS('../data/rds/Adult.Ex.rds')
sp.PFC <- readRDS('../data/rds/sp.PFC.rds')
Adult.Ex.barcode <- subset(
  Adult.Ex,
  cells=colnames(Adult.Ex)[which(Adult.Ex$BC_num>0)]
)
Adult.IT.PT.barcode <- subset(Adult.Ex, cells=colnames(Adult.Ex)[which(
  (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "IT") |
    (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "PT" & Adult.Ex$sample == "Adult1")
)])



gene_lib <- read.csv('../data/csv/transmitter_and_receptor/gene_lib.csv')
# receptor
Monoamine_R <- str_to_title(gene_lib$Monoamine.Ach.receptor)
Neuropeptides_R <- str_to_title(gene_lib$Neuropeptides.receptor)
mGluR <- str_to_title(gene_lib$mGluR.Kainate.receptor)
GABA_R <- str_to_title(gene_lib$GABA.receptor)
AMPA_NMDA_R <- str_to_title(gene_lib$AMPA.NMDA)
receptor_gene <- c(Monoamine_R, Neuropeptides_R, mGluR, GABA_R, AMPA_NMDA_R)
# NT/NP
Neurotransmitter <- str_to_title(gene_lib$Neurotransmitter)
Neuropeptides <- str_to_title(gene_lib$Neuropeptides)
NTNP_gene <- c(Neurotransmitter, Neuropeptides)

all_gene <- c(receptor_gene, NTNP_gene)

# filter
seu <- Adult.IT.PT.barcode
all_gene <- all_gene[which(all_gene %in% rownames(seu))]
# gene expression filter
all_gene_exp <- AverageExpression(
  seu, features=all_gene, assays="RNA", slot="data", group.by="Proj_subtype"
)$RNA
all_gene_exp <- as.data.frame(log1p(all_gene_exp))
all_gene_exp$max <- apply(all_gene_exp, 1, max)
# gene cell percentage filter
all_gene_pct <- as.data.frame(t(as.matrix(seu@assays$RNA@data[all_gene,])))
all_gene_pct$Proj_subtype <- as.character(seu$Proj_subtype)
all_gene_pct <-
  all_gene_pct |>
  dplyr::group_by(Proj_subtype) |>
  dplyr::summarize(across(1:length(all_gene), function(x){
    length(which(x>0))/length(x)
  })) |>
  as.data.frame()
rownames(all_gene_pct) <- all_gene_pct$Proj_cluster
all_gene_pct <- as.data.frame(t(all_gene_pct[,-1]))
colnames(all_gene_pct) <- 1:33
all_gene_pct$max <- apply(all_gene_pct, 1, max)

all_gene <- all_gene[which(all_gene_exp$max>0.1 & all_gene_pct$max>0.1)]

Monoamine_R <- Monoamine_R[which(Monoamine_R %in% all_gene)]
Neuropeptides_R <- Neuropeptides_R[which(Neuropeptides_R %in% all_gene)]
mGluR <- mGluR[which(mGluR %in% all_gene)]
GABA_R <- GABA_R[which(GABA_R %in% all_gene)]
AMPA_NMDA_R <- AMPA_NMDA_R[which(AMPA_NMDA_R %in% all_gene)]

Neurotransmitter <- Neurotransmitter[which(Neurotransmitter %in% all_gene)]
Neuropeptides <- Neuropeptides[which(Neuropeptides %in% all_gene)]

receptor_gene <- receptor_gene[which(receptor_gene %in% all_gene)]
NTNP_gene <- NTNP_gene[which(NTNP_gene %in% all_gene)]




Barcode <- c('VIS-I','SSp-I','CP-I','AUD-I','RSP-I',
             'BLA-I','ACB-I','ENTl-I','AId-I','ECT-I',
             'ACB-C','PL-C','ECT-C','ENTl-C',
             'BLA-C','CP-C','AId-C','RSP-C',
             'MD-I','RE-I','DR-I','VTA-I','LHA-I','SC-I')
Proj_subtype <- c("1","14","17","22","23","25","26","27","28","29","31",
                  "7","8","9","10","15","16","18","19","20","21",
                  "2","3","11","12","13","24","30","32","33",
                  "4","5","6")
seu <- Adult.IT.PT.barcode
pct_mat <- matrix(nrow = 33, ncol = length(Barcode))
rownames(pct_mat) <- Proj_subtype
colnames(pct_mat) <- Barcode
for (i in 1:nrow(pct_mat)){
  for (j in 1:ncol(pct_mat)){
    pct_mat[i,j] <- length(which(seu@meta.data[,Barcode[j]]>0 &
                                   seu$Proj_subtype == rownames(pct_mat)[i]))/
      length(which(seu$Proj_subtype == rownames(pct_mat)[i]))
  }
}

Proj_subtype <- c("1","14","17","22","23","25","26","27","28","29","31",
                  "7","8","9","10","15","16","18","19","20","21",
                  "2","3","11","12","13","24","30","32","33",
                  "4","5","6")
seu <- Adult.IT.PT.barcode
seu$Proj_subtype <- factor(seu$Proj_subtype, levels = Proj_subtype)
avgexp <- AverageExpression(seu, features = receptor_gene,
                            group.by = "Proj_subtype", assays = "RNA")
avgexp <- avgexp$RNA
avgexp <- scale(t(avgexp))
breaks <- seq(0,1,0.01)

p1_1 <-
  pheatmap::pheatmap(
    t(avgexp),
    breaks = breaks, border_color = NA, show_rownames=F, show_colnames = F,
    color = colorRampPalette(c("white","#D73027"))(length(breaks)),
    cluster_cols = F, cluster_rows = T,
    gaps_col = c(1:32), legend=F, treeheight_row=F
  )

ggsave("../pdf/Figure5/p1_1.png", plot = p1_1, height = 0.5, width = 15, units = "in")


seu <- Adult.IT.PT.barcode
seu$Proj_subtype <- factor(seu$Proj_subtype, levels = Proj_subtype)
avgexp <- AverageExpression(seu, features = NTNP_gene,
                            group.by = "Proj_subtype", assays = "RNA")
avgexp <- avgexp$RNA
avgexp <- scale(t(avgexp))
breaks <- seq(0,1,0.01)

p1_2 <-
  pheatmap::pheatmap(
    t(avgexp),
    breaks = breaks, border_color = NA, show_rownames=F, show_colnames = F,
    color = colorRampPalette(c("white","#4575B4"))(length(breaks)),
    cluster_cols = F, cluster_rows = T,
    gaps_col = c(1:32), legend=F, treeheight_row=F
  )

ggsave("../pdf/Figure5/p1_2.png", plot = p1_2, height = 0.5, width = 15, units = "in")



seu <- Adult.IT.PT.barcode
ITi_D_DEGs <- c("Htr2a","Gabra1","Adra1b","Gabra4","Cck","Grm2","Gabrg2","Grm3",
                "Chrm3","Gria4","Slc17a7","Grin1","Grin2a","Adrb1","Gabbr1","Gabbr2",
                "Gabrb2","Adcyap1","Gabra5","Penk")
ITi_V_DEGs <- c("Gria1","Grp","Trhr","Grik2","Grm5","Grin3a","Grm1","Ntsr1",
                "Npy2r","Adra2a","Ramp3","Pdyn","Gabra2","Rxfp1","Gabrb1",
                "Cnr1","Htr2c","Npy")
gene <- c(ITi_D_DEGs, ITi_V_DEGs)
avg_exp <- AverageExpression(seu, features = gene, group.by = 'Proj_subtype',
                             assays = "RNA", slot = "data")
avg_exp <- avg_exp$RNA
avg_exp_zscore <- as.matrix(scale(t(avg_exp)))

ITi_D <- c(1,14,17,22,23,25,26,27,28,29,31)
ITi_V <- c(7,8,9,10,15,16,18,19,20,21)
ITc <- c(2,3,11,12,13,24,30,32,33)
PTi <- c(4,5,6)
Proj_subtype_order <- c(ITi_D,ITi_V)

annotation_col = data.frame(
  Gene_type = rep("Receptor", length(c(ITi_D_DEGs, ITi_V_DEGs))),
  row.names = c(ITi_D_DEGs, ITi_V_DEGs)
)

annotation_col$Gene_type[which(rownames(annotation_col) %in% NTNP_gene)] <- "NT/NP"

annotation_row = data.frame(
  Projection_module = factor(rep(c("ITi-D", "ITi-V"), c(11, 10)),
                             levels = c("ITi-D", "ITi-V")),
  row.names = Proj_subtype_order
)

ann_colors = list(
  Projection_module = c('ITi-D'='#1f77b4','ITi-V'='#ff7f0e'),
  Gene_type = c('Receptor'="#93C8C0FF",'NT/NP'="#1C3C63FF")
)

breaks <- seq(-2,2,0.01)
col <- colorRampPalette(c("lightblue3", "lightblue", "white", "red", "red4"))(length(breaks))

p2 <- pheatmap::pheatmap(avg_exp_zscore[Proj_subtype_order, gene],
                                cluster_rows = F, cluster_cols = F,
                                breaks = breaks,
                                color = col,
                                annotation_row = annotation_row, annotation_col = annotation_col,
                                annotation_colors = ann_colors,
                                gaps_col = c(20),
                                gaps_row = c(11),
                                fontsize_col = 10,
                                annotation_names_row=F,annotation_names_col=F,
                                show_colnames = T
)

ggsave("../pdf/Figure5/p2.pdf", plot = p2, height = 5, width = 9, units = "in")



seu <- Adult.IT.PT.barcode
gene <- c("Grm2","Grm3","Grin1","Grin2a","Gria4","Gria1","Grik2","Grm5","Grin3a","Grm1","Sstr2","Sstr3","Sstr1","Ramp1","Ramp2","Ramp3","Adipor1","Adipor2")
avg_exp <- AverageExpression(seu, features = gene, group.by = 'Proj_subtype',
                             assays = "RNA", slot = "data")
avg_exp <- avg_exp$RNA
avg_exp_zscore <- as.data.frame(scale(t(avg_exp)))

avg_exp_zscore$Group <- rownames(avg_exp_zscore)
df_exp <- pivot_longer(avg_exp_zscore, !Group, names_to = "Gene", values_to = "EXP")

pct <- matrix(nrow = 33, ncol = length(gene))
rownames(pct) <- rownames(avg_exp_zscore)
colnames(pct) <- gene
for (i in 1:33){
  pct[i,] <- apply(seu@assays$RNA@data[gene,seu$Proj_subtype==rownames(pct)[i]], 1,
                   function(x){
                     length(which(x>0))/length(x)
                   })
}
pct <- as.data.frame(pct)
pct$Group <- rownames(pct)
df_pct <- pivot_longer(pct,!Group, names_to = "Gene", values_to = "PCT")

df <- data.frame(
  "Group" = df_exp$Group,
  "Gene" = df_exp$Gene,
  "EXP" = df_exp$EXP,
  "PCT" = df_pct$PCT
)
df$EXP[which(df$EXP > 2)] <- 2
df$EXP[which(df$EXP < -2)] <- -2
df$Group <- factor(df$Group, levels = rev(names(col_Proj_subtype)))

df <- df[which(df$Group %in% c(14, 29, 8, 16)),]

p3 <-
  Biorplot::Bior_DotPlot(
    data = df, x = "Gene", y = "Group", size="PCT", color = "EXP",
    x.text.col = FALSE, ggtheme = theme_bw()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x="",y="") +
  scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4"), limits=c(-2,2)) +
  scale_size_continuous(range = c(1,6))

ggsave("../pdf/Figure5/p3.pdf", plot = p3, height = 4, width = 6, units = "in")



seu <- subset(Adult.IT.PT.barcode, cells = colnames(Adult.IT.PT.barcode)[which(
  Adult.IT.PT.barcode$`BLA-I` == 1 & Adult.IT.PT.barcode$Ex_subtype == "IT"
)])
seu$BLA_Layer <- "other"
seu$BLA_Layer[which(seu$SubType_Layer=="L2/3 IT")] <- "L2/3 BLA"
seu$BLA_Layer[which(seu$SubType_Layer=="L4/5 IT")] <- "L4/5 BLA"
seu$BLA_Layer[which(seu$SubType_Layer=="L5 IT")] <- "L5 BLA"
seu$BLA_Layer[which(seu$SubType_Layer=="L6 IT")] <- "L6 BLA"

gene <- NTNP_gene
gene <- c("Pdyn","Tac2","Pomc","Nppc","Grp","Slc17a6","Adcyap1","Slc17a7",
          "Dbi","Penk","Nts","Cck","Npy")
avg_exp <- AverageExpression(seu, features = gene, group.by = 'BLA_Layer',
                             assays = "RNA", slot = "data")
avg_exp <- avg_exp$RNA
avg_exp_zscore <- as.data.frame(scale(t(avg_exp)))

breaks <- seq(-2,2,0.01)
p4 <- pheatmap::pheatmap(
  t(avg_exp_zscore[,-ncol(avg_exp_zscore)]),
  cluster_rows = F, cluster_cols = F, show_rownames = T, treeheight_row = 0,
  breaks = breaks,
  color = colorRampPalette(rev(sciRcolor::pal_scircolor(83)))(length(breaks)),
  fontsize_col = 10, fontsize_row = 10,
  annotation_names_row=F, annotation_names_col=F,
  show_colnames = T
)

ggsave("../pdf/Figure5/p4.pdf", plot = p4, height = 6, width = 8.7, units = "in")



seu <- subset(
  Adult.IT.PT.barcode,
  cells=colnames(Adult.Ex)[which(Adult.Ex$Ex_subtype=="IT" &
                                   Adult.Ex$BC_num>0)])
seu$group <- "Other"
seu$group[which(seu$`AId-I`>0 & seu$`CP-I`>0 & seu$BC_num==2)] <- "AId-I + CP-I"
seu$group[which(seu$`AId-I`>0 & seu$`CP-C`>0 & seu$`PL-C`>0 & seu$BC_num==3)] <- "AId-I + CP-C + PL-C"
seu$group[which(seu$`AId-I`>0 & seu$`CP-I`>0 & seu$`CP-C`>0 & seu$`AId-C`>0 & seu$BC_num==4)] <- "AId-I + CP-I + CP-C + AId-C"
seu$group[which(seu$`AId-I`>0 & seu$`ACB-I`>0 & seu$`CP-I`>0 & seu$`CP-C`>0 & seu$`ACB-C`>0 & seu$`AId-C`>0 & seu$BC_num==6)] <- "AId-I,ACB-I,CP-I,CP-C,ACB-C,AId-C"
seu$group[which(seu$`AId-I`>0 & seu$`ACB-I`>0 & seu$BC_num==2)] <- "AId-I + ACB-I"
seu$group[which(seu$`AId-I`>0 & seu$`BLA-I`>0 & seu$`ECT-I`>0 & seu$`ENTl-I`>0 & seu$`ACB-I`>0 & seu$`CP-I`>0 & seu$BC_num==6)] <- "AId-I + BLA-I + ECT-I + LENT-I + ACB-I + CP-I"

seu$group <- factor(seu$group, levels = c("AId-I + CP-I","AId-I + CP-C + PL-C","AId-I + CP-I + CP-C + AId-C","AId-I,ACB-I,CP-I,CP-C,ACB-C,AId-C","AId-I + ACB-I","AId-I + BLA-I + ECT-I + LENT-I + ACB-I + CP-I"))
seu <- subset(seu, cells = colnames(seu)[which(seu$group != "Other")])
gene <- c("Slc17a7","Slc17a6","Penk","Pomc","Pdyn","Cck","Npy","Nppc","Grp","Adcyap1",
          "Tac2","Nts","Dbi","Nos1")
gene <- all_gene
avg_exp <- AverageExpression(seu, features = gene, group.by = 'group',
                             assays = "RNA", slot = "data")
avg_exp <- avg_exp$RNA
avg_exp_zscore <- as.data.frame(scale(t(avg_exp)))

breaks <- seq(-2,2,0.01)
p5 <- pheatmap::pheatmap(
  t(avg_exp_zscore),
  cluster_rows = T, cluster_cols = F, show_rownames = F, treeheight_row = 0,
  breaks = breaks,
  color = colorRampPalette(rev(sciRcolor::pal_scircolor(83)))(length(breaks)),
  fontsize_col = 10, fontsize_row = 3,
  annotation_names_row=F, annotation_names_col=F,
  show_colnames = T
)

ggsave("../pdf/Figure5/p5.pdf", plot = p5, height = 8, width = 7.5, units = "in")



load("../data/csv/transmitter_and_receptor/hdwgcna.RData")

set.seed(20240703)
modules <- GetModules(seu)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

p6 <-
  HubGeneNetworkPlot(
    seu,
    n_hubs = 10, n_other=100,
    edge_prop = 0.8,
    edge.alpha = 0.8,
    mods = mods,
    return_graph=FALSE
  )



Barcode_order <- c("MD-I","RE-I","DR-I","VTA-I","LHA-I","SC-I",
                   "VIS-I","SSp-I","CP-I","AUD-I","RSP-I",
                   "BLA-I","ACB-I","ENTl-I","AId-I","ECT-I",
                   "ACB-C","PL-C","ECT-C","ENTl-C","BLA-C","CP-C","AId-C","RSP-C")
Ex_BC_mat <- all.Adult@meta.data[colnames(seu), Barcode_order]
Ex_BC_mat <- log1p(10000*Ex_BC_mat/seu$nCount_RNA)

m_obj <- GetMetacellObject(seu)
MEs <- GetMEs(seu, harmonized=FALSE)
MEs <- MEs[,-grep('grey',colnames(MEs))]
meta_MEs <- matrix(nrow = nrow(m_obj@meta.data), ncol = ncol(MEs))
rownames(meta_MEs) <- rownames(m_obj@meta.data)
colnames(meta_MEs) <- colnames(MEs)
meta_BCs <- matrix(nrow = nrow(m_obj@meta.data), ncol = length(Barcode_order))
rownames(meta_BCs) <- rownames(m_obj@meta.data)
colnames(meta_BCs) <- Barcode_order
for (i in 1:nrow(meta_MEs)){
  meta_MEs[i,] <- colMeans(MEs[strsplit(m_obj$cells_merged[[i]], ',')[[1]],], na.rm = T)
  meta_BCs[i,] <- colMeans(Ex_BC_mat[strsplit(m_obj$cells_merged[[i]], ',')[[1]],], na.rm = T)
}
meta_BCs[is.nan(meta_BCs)] <- NA
meta_R <- matrix(nrow = 4, ncol = length(Barcode_order))
rownames(meta_R) <- c("M1","M2","M3","M4")
colnames(meta_R) <- Barcode_order
meta_P <- meta_R
for (i in 1:4){
  for (j in 1:24){
    cor <- cor.test(meta_MEs[,rownames(meta_R)[i]], meta_BCs[,colnames(meta_R)[j]],
                    "two.sided", "pearson")
    meta_R[i,j] <- round(cor$estimate,2)
    meta_P[i,j] <- format(cor$p.value, digits = 2)
  }
}
meta_R <- as.data.frame(meta_R)
meta_R$Module <- rownames(meta_R)
df_meta_R <- pivot_longer(meta_R, !Module, names_to = "Target", values_to = "R")
meta_P <- as.data.frame(meta_P)
meta_P$Module <- rownames(meta_P)
df_meta_P <- pivot_longer(meta_P, !Module, names_to = "Target", values_to = "Pvalue")
df <- data.frame(
  "Module" = df_meta_R$Module,
  "Target" = df_meta_R$Target,
  "R" = as.numeric(df_meta_R$R),
  "Pvalue" = as.numeric(df_meta_P$Pvalue)
)
df$Log_Pvalue <- -log10(df$Pvalue)

Barcode <- c("VIS-I","SSp-I","CP-I","AUD-I","RSP-I",
             "BLA-I","ACB-I","ENTl-I","AId-I","ECT-I",
             "ACB-C","PL-C","ECT-C","ENTl-C","BLA-C","CP-C","AId-C","RSP-C",
             "MD-I","RE-I","DR-I","VTA-I","LHA-I","SC-I")
df$Target <- factor(df$Target, levels = Barcode)
df$Module <- factor(df$Module, levels = rev(c("M1","M2","M3","M4")))
df$label <- ""
df$label[which(df$R>0 & df$Pvalue < 0.05)] <- "*"
df$label[which(df$R>0 & df$Pvalue < 0.01)] <- "**"
df$label[which(df$R>0 & df$Pvalue < 0.001)] <- "***"

breaks <- seq(-0.5,0.5,0.01)
p7 <-
  ggplot(df, aes(x=Target, y=Module, fill=R)) +
  geom_tile(color="gray", size=0.5) +
  geom_text(aes(label=label),col ="black",size = 3) +
  scale_fill_gradientn(limits=c(-0.5,0.5), colours = colorRampPalette(c("navy","white","firebrick3"))(100), na.value="firebrick3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid = element_blank()) +
  labs(x="", y="")

ggsave("../pdf/Figure5/p7.pdf", plot = p7, height = 4, width = 6.5, units = "in")



df <- as.data.frame(seu@reductions$umap@cell.embeddings)
df$M1 <- seu@misc$Adult.IT.PT.barcode$MEs$M1
df$M2 <- seu@misc$Adult.IT.PT.barcode$MEs$M2
df$`AId-I` <- seu$`AId-I`
df$`AId-I`[is.na(df$`AId-I`)] <- 0
df$`CP-I` <- seu$`CP-I`
df$`CP-I`[is.na(df$`CP-I`)] <- 0
df <- df[which(seu$Ex_subtype=="IT"),]
df <- df[df$UMAP_1 < 2 & df$UMAP_2 > -3,]

p8_1 <-
  ggplot() +
  geom_point(df, mapping = aes(x = UMAP_1, y = UMAP_2), color="lightgray", size=1) +
  geom_pointdensity(df[which(df$M2>0),], mapping = aes(x = UMAP_1, y = UMAP_2), size=1) +
  scale_color_gradientn(colours = c("lightgray","white","yellow"),
                        na.value = "yellow",
                        limits = c(0,50),
                        breaks = c(0,50)
  ) +
  coord_fixed() +
  theme_bw() +
  labs(title = "M2", x="", y="", colour = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.width  = unit(0.5, "lines"),
        legend.key.height = unit(2, "lines"))

ggsave("../pdf/Figure5/p8_1.pdf", plot = p8_1, height = 3, width = 6, units = "in")


p8_2 <-
  ggplot() +
  geom_point(df, mapping = aes(x = UMAP_1, y = UMAP_2), color="lightgray", size=1) +
  geom_pointdensity(df[which(df$`CP-I`>0),], mapping = aes(x = UMAP_1, y = UMAP_2), size=1) +
  scale_color_gradientn(colours = c("lightgray","white","red"),
                        na.value = "red",
                        limits = c(0,50),
                        breaks = c(0,50)
  ) +
  coord_fixed() +
  theme_bw() +
  labs(title = "CP-I", x="", y="", colour = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.width  = unit(0.5, "lines"),
        legend.key.height = unit(2, "lines"))

ggsave("../pdf/Figure5/p8_2.pdf", plot = p8_2, height = 3, width = 6, units = "in")



df_cor <- df[which(df$`CP-I`>0),]
df_cor[,'exp_order'] <- rank(df_cor$M2)
df_cor[,'bin'] <- cut(
  df_cor[,'exp_order'],
  seq(1, max(df_cor[,'exp_order']), length.out=11),
  labels = c(1:10)
)
df_cor <-
  df_cor |>
  group_by(bin) |>
  dplyr::summarize(across(3:6, ~ mean(.x, na.rm = TRUE)))

df_cor_norm <- as.data.frame(apply(df_cor[1:10,2:ncol(df_cor)],2,function(x){
  (x-min(x))/(max(x)-min(x))
}))
df_cor_norm$x <- df_cor$bin[1:10]

cor <- cor.test(df_cor_norm$M2, df_cor_norm$`CP-I`, "two.sided", "pearson")
R <- round(cor$estimate,2)
P <- format(cor$p.value, digits = 2)

p8_3 <-
  ggplot(df_cor_norm, aes(x=x, y=`CP-I`, group=1)) +
  geom_line(linewidth=1, color="yellow") +
  geom_point(fill="yellow",shape=21, color="black" , size=3) +
  geom_text(x=1, y=0.95, label = paste('R =',R,'\nP =',P,sep=' '), hjust=0,
            color="black") +
  scale_x_discrete(breaks = c(1,5,10),
                   labels = c("10%","50%","100%")) +
  labs(x=paste("M2 Eigengenes"),
       y=paste("CP-I projection strength"),
       title="") +
  theme_classic() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 20),
        text = element_text(size = 10))

ggsave("../pdf/Figure5/p8_3.pdf", plot = p8_3, height = 4, width = 4, units = "in")





