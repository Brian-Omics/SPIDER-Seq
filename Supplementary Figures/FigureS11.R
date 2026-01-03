# --------------------
# title: FigureS11 Code
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
library(umap)
library(scCustomize)
library(ggpointdensity)

library(Biorplot)
source('bin/Palettes.R')
source('bin/includes.R')

Adult.Ex <- readRDS('../data/rds/Adult.Ex.rds')

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



df <- data.frame(
  "gene" = all_gene,
  "max_exp" = as.numeric(all_gene_exp$max[match(all_gene, rownames(all_gene_exp))]),
  "max_pct" = as.numeric(all_gene_pct$max[match(all_gene, rownames(all_gene_exp))])
)
df$gene_type <- ""
df$gene_type[which(df$gene %in% receptor_gene)] <- "Receptor"
df$gene_type[which(df$gene %in% NTNP_gene)] <- "NT/NP"
df$gene <- factor(df$gene, levels = df$gene[order(df$max_exp)])

P1 <-
  ggplot(df, aes(x=max_exp, y=gene, color=gene_type)) +
  geom_point() +
  scale_color_manual(values = c('Receptor'="#93C8C0FF",'NT/NP'="#1C3C63FF")) +
  scale_x_continuous(breaks = c(0.1,1,2,3), limits = c(0.1,3.5)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8))

df$gene <- factor(df$gene, levels = df$gene[order(df$max_pct)])
P2 <-
  ggplot(df, aes(x=max_pct, y=gene, color=gene_type)) +
  geom_point() +
  scale_color_manual(values = c('Receptor'="#93C8C0FF",'NT/NP'="#1C3C63FF")) +
  scale_x_continuous(breaks = c(0.1,0.5,1), limits = c(0.1,1)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8))

legend <- get_legend(P1)

p <-
  plot_grid(P1 + theme(legend.position="none"),
            P2 + theme(legend.position="none"),
            legend, nrow = 1)

ggsave("../pdf/FigureS11/p1.pdf", plot = p, height = 10, width = 6, units = "in")



seu <- Adult.IT.PT.barcode
seu$ITPT <- "IT"
seu$ITPT[which(seu$Proj_module == "PTi")] <- "PT"
Idents(seu) <- "ITPT"
DefaultAssay(seu) <- "RNA"

features <- rev(c("Gria4","Grm3","Htr2a","Htr5a","Htr1b","Drd1","Npy1r","Cckbr","Npr3","Mchr1","Mc4r"))
exp <- AverageExpression(seu, features = features, assays = "RNA", slot = "data")
exp <- exp$RNA
exp <- as.data.frame(exp)
exp <- apply(exp, 1, function(x){x/max(x)})
exp <- as.data.frame(exp)
exp$Group <- rownames(exp)
df_exp <- pivot_longer(exp, !Group, names_to = "Gene", values_to = "EXP")

pct <- matrix(nrow = 2, ncol = length(features))
pct[1,] <- apply(seu@assays$RNA@data[features,seu$ITPT=="IT"],1,function(x){
  length(which(x>0))/length(x)
})
pct[2,] <- apply(seu@assays$RNA@data[features,seu$ITPT=="PT"],1,function(x){
  length(which(x>0))/length(x)
})
rownames(pct) <- c("IT","PT")
colnames(pct) <- features
pct <- as.data.frame(pct)
pct$Group <- rownames(pct)
df_pct <- pivot_longer(pct,!Group, names_to = "Gene", values_to = "PCT")

df <- data.frame(
  "Group" = df_exp$Group,
  "Gene" = df_exp$Gene,
  "EXP" = df_exp$EXP,
  "PCT" = df_pct$PCT
)
df$Group <- factor(df$Group, levels = c("IT","PT"))

p2_1 <-
  Biorplot::Bior_DotPlot(
    data = df, x = "Gene", y = "Group", size="PCT", color = "EXP",
    x.text.col = FALSE, ggtheme = theme_bw()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  labs(x="",y="") +
  scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4")) +
  scale_size_continuous(range = c(1,6))

ggsave("../pdf/FigureS11/p2_1.pdf", plot = p2_1, height = 2, width = 3, units = "in")



features <- rev(c(c("Cck","Penk","Adcyap1","Grp","Pdyn")))
features <- rev(c(c("Cck","Penk","Npy","Sst","Nppc","Gad1","Gad2","Slc17a7","Dbi","Slc17a6","Grp","Adcyap1","Pdyn")))
features <- rev(c(c("Cck","Penk","Npy","Sst","Nppc","Slc17a7","Dbi","Slc17a6","Grp","Adcyap1","Pdyn")))
exp <- AverageExpression(seu, features = features, assays = "RNA", slot = "data")
exp <- exp$RNA
exp <- as.data.frame(exp)
exp <- apply(exp, 1, function(x){x/max(x)})
exp <- as.data.frame(exp)
exp$Group <- rownames(exp)
df_exp <- pivot_longer(exp, !Group, names_to = "Gene", values_to = "EXP")

pct <- matrix(nrow = 2, ncol = length(features))
pct[1,] <- apply(seu@assays$RNA@data[features,seu$ITPT=="IT"],1,function(x){
  length(which(x>0))/length(x)
})
pct[2,] <- apply(seu@assays$RNA@data[features,seu$ITPT=="PT"],1,function(x){
  length(which(x>0))/length(x)
})
rownames(pct) <- c("IT","PT")
colnames(pct) <- features
pct <- as.data.frame(pct)
pct$Group <- rownames(pct)
df_pct <- pivot_longer(pct,!Group, names_to = "Gene", values_to = "PCT")

df <- data.frame(
  "Group" = df_exp$Group,
  "Gene" = df_exp$Gene,
  "EXP" = df_exp$EXP,
  "PCT" = df_pct$PCT
)
df$Group <- factor(df$Group, levels = c("IT","PT"))

p2_2 <-
  Biorplot::Bior_DotPlot(
    data = df, x = "Gene", y = "Group", size="PCT", color = "EXP",
    x.text.col = FALSE, ggtheme = theme_bw()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  labs(x="",y="") +
  scale_color_gradientn(colours = c("lightblue3", "lightblue", "white", "red", "red4")) +
  scale_size_continuous(range = c(1,6))

ggsave("../pdf/FigureS11/p2_2.pdf", plot = p2_2, height = 2, width = 3, units = "in")



seu <- Adult.IT.PT.barcode
gene <- c("Htr1a","Htr5a","Htr2a","Htr7","Htr1f","Htr2c","Htr1b")
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
  scale_size_continuous(range = c(1,10))

ggsave("../pdf/FigureS11/p3.pdf", plot = p3, height = 4, width = 4.5, units = "in")



seu <- subset(
  Adult.IT.PT.barcode,
  cells=colnames(Adult.Ex)[which(Adult.Ex$Ex_subtype=="IT" &
                                   Adult.Ex$BC_num>0)])
seu$group <- "Other"
seu$group[which(seu$`RSP-I`>0 & seu$`CP-I`>0 & seu$BC_num==2)] <- "RSP-I + CP-I"
seu$group[which(seu$`RSP-I`>0 & seu$`CP-I`>0 & seu$`CP-C`>0 & seu$BC_num==3)] <- "RSP-I + CP-I + CP-C"
seu$group[which(seu$`RSP-I`>0 & seu$`CP-I`>0 & seu$`ACB-I`>0 & seu$BC_num==3)] <- "RSP-I + CP-I + ACB-I"
seu$group[which(seu$`RSP-I`>0 & seu$`CP-I`>0 & seu$`SSp-I`>0 & seu$BC_num==3)] <- "RSP-I + CP-I + SSp-I"
seu$group[which(seu$`RSP-I`>0 & seu$`CP-I`>0 & seu$`VIS-I`>0 & seu$BC_num==3)] <- "RSP-I + CP-I + VIS-I"
seu$group[which(seu$`RSP-I`>0 & seu$`ACB-I`>0 & seu$BC_num==2)] <- "RSP-I + ACB-I"
seu$group[which(seu$`RSP-I`>0 & seu$`CP-I`>0 & seu$`SSp-I`>0 & seu$`VIS-I`>0 & seu$BC_num==4)] <- "RSP-I + CP-I + SSp-I + VIS-I"

seu <- subset(seu, cells = colnames(seu)[which(seu$group != "Other")])
seu$group <- factor(seu$group, levels = c("RSP-I + CP-I","RSP-I + CP-I + CP-C","RSP-I + CP-I + ACB-I","RSP-I + CP-I + SSp-I","RSP-I + CP-I + VIS-I","RSP-I + ACB-I","RSP-I + CP-I + SSp-I + VIS-I"))

gene <- all_gene
avg_exp <- AverageExpression(seu, features = gene, group.by = 'group',
                             assays = "RNA", slot = "data")
avg_exp <- avg_exp$RNA
avg_exp_zscore <- as.data.frame(scale(t(avg_exp)))
avg_exp_zscore <- avg_exp_zscore[,!is.nan(colSums(avg_exp_zscore))]

breaks <- seq(-2,2,0.01)
p4 <- pheatmap::pheatmap(
  t(avg_exp_zscore),
  cluster_rows = T, cluster_cols = T, show_rownames = F, treeheight_row = 0,
  treeheight_col = 0, breaks = breaks,
  color = colorRampPalette(rev(sciRcolor::pal_scircolor(83)))(length(breaks)),
  fontsize_col = 10, fontsize_row = 3,
  annotation_names_row=F, annotation_names_col=F,
  show_colnames = T, border_color=NA
)

ggsave("../pdf/FigureS11/p4.pdf", plot = p4, height = 8, width = 10, units = "in")



load("../data/csv/transmitter_and_receptor/hdwgcna.RData")

p5 <- PlotDendrogram(seu, main='hdWGCNA Dendrogram')



df <- as.data.frame(seu@reductions$umap@cell.embeddings)
df$M1 <- seu@misc$Adult.IT.PT.barcode$MEs$M1
df$M2 <- seu@misc$Adult.IT.PT.barcode$MEs$M2
df$M3 <- seu@misc$Adult.IT.PT.barcode$MEs$M3
df$M4 <- seu@misc$Adult.IT.PT.barcode$MEs$M4
df <- df[df$UMAP_1 < 1.2 & df$UMAP_1 > -13,]

module <- c("M1","M2","M3","M4")
col <- c("blue","yellow","turquoise","brown")
plist <- list()
for (i in 1:4){
  plist[[i]] <-
    ggplot() +
    geom_point(df, mapping = aes(x = UMAP_1, y = UMAP_2), color="lightgray", size=1) +
    geom_pointdensity(df[which(df[,module[i]]>0),], mapping = aes(x = UMAP_1, y = UMAP_2), size=1) +
    scale_color_gradientn(colours = c("lightgray","white",col[i]),
                          na.value = col[i], limits = c(0,50), breaks = c(0,50)) +
    coord_fixed() +
    theme_void() +
    labs(title = module[i], x="", y="", colour = "") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(hjust = 0.5, size = 15),
          legend.position = "right",
          legend.key.width  = unit(1, "lines"),
          legend.key.height = unit(1.5, "lines"))
}
p6 <- plot_grid(plotlist = plist, ncol=4)

ggsave("../pdf/FigureS11/p6.pdf", plot = p6, height = 3, width = 12, units = "in")



ModuleNetworkPlot(
  seu,
  outdir = '../pdf/FigureS11/ModuleNetworks'
)




