# --------------------
# title: Figure2 Code
# author: Hu Zheng
# date: 2026-01-01
# --------------------

library(Seurat)
library(scCustomize)
library(tidyverse)
library(tidydr)
library(cowplot)
library(umap)
library(pheatmap)
library(ggblur)
library(ggrepel)
library(viridis)
library(networkD3)

library(Biorplot)
source('bin/Palettes.R')
source('bin/includes.R')

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



barcode <- c("SSp-I","AUD-I","RSP-I","VIS-I","CP-I",
             "ECT-I","AId-I","BLA-I","ACB-I","LHA-I")
seu <- subset(sp.PFC, cells=colnames(sp.PFC)[which(
  sp.PFC$ABA_hemisphere=="Left" & sp.PFC$BC_num>0)])
bc_slice <- seu@meta.data[,c(barcode, 'Y')]
bc_slice <-
  bc_slice |>
  mutate(bin = cut(Y, breaks = 36))
bin <- sort(unique(bc_slice$bin))
bc_slice$bin_index <- match(bc_slice$bin, bin)

bc_slice_mean <-
  bc_slice[,c(barcode,"bin_index")] |>
  group_by(bin_index) |>
  dplyr::summarize(across(1:length(barcode), ~sum(.x>0, na.rm = TRUE))) |>
  as.data.frame()
for (i in 1:36){
  bc_slice_mean[i,2:11] <- bc_slice_mean[i,2:11]/length(which(bc_slice$bin_index==i))
}
bc_slice_mean <- bc_slice_mean[,2:11]
bc_slice_mean <- t(scale(bc_slice_mean))

annotation_row = data.frame(
  Proj_module = rep(c("ITi-D","ITi-V"),c(5,5))
)
rownames(annotation_row) <- rownames(bc_slice_mean)
ann_colors = list(
  Proj_module = c('ITi-D'='#1f77b4','ITi-V'='#ff7f0e')
)

p1 <-
  pheatmap(bc_slice_mean[barcode, rev(1:36)],
           color = colorRampPalette(c("white","#1f77b4"))(200),
           breaks = seq(-1,1,0.01), show_colnames = F, treeheight_row=10,
           cluster_rows = F, cluster_cols = F, annotation_row=annotation_row,
           annotation_colors = ann_colors, annotation_names_row=F,
           annotation_legend=F
  )

ggsave("../pdf/Figure2/p1.pdf", plot = p1, height = 3, width = 7, units = "in")



barcode <- c("SSp-I","AUD-I","RSP-I","VIS-I","CP-I",
             "ECT-I","AId-I","BLA-I","ACB-I","LHA-I")
seu <- subset(sp.PFC, cells=colnames(sp.PFC)[which(
  sp.PFC$ABA_hemisphere=="Left" & sp.PFC$BC_num>0)])
bc_slice <- seu@meta.data[,c(barcode, 'slice')]
bc_slice_mean <-
  bc_slice |>
  group_by(slice) |>
  dplyr::summarize(across(1:length(barcode), ~sum(.x>0, na.rm = TRUE))) |>
  as.data.frame()
for (i in 1:36) {
  bc_slice_mean[i,2:11] <- bc_slice_mean[i,2:11]/length(which(bc_slice$slice==bc_slice_mean$slice[i]))
}
rownames(bc_slice_mean) <- bc_slice_mean$slice
bc_slice_mean <- bc_slice_mean[,2:11]
bc_slice_mean <- t(scale(bc_slice_mean))

annotation_row = data.frame(
  Proj_module = rep(c("ITi-D","ITi-V"),c(5,5))
)
rownames(annotation_row) <- rownames(bc_slice_mean)
ann_colors = list(
  Proj_module = c('ITi-D'='#1f77b4','ITi-V'='#ff7f0e')
)

p2 <-
  pheatmap(bc_slice_mean[barcode,],
           color = colorRampPalette(c("white","#d62728"))(200),
           breaks = seq(-1,1,0.01), show_colnames = F, treeheight_row=10,
           cluster_rows = F, cluster_cols = F, annotation_row=annotation_row,
           annotation_colors = ann_colors, annotation_names_row=F,
           annotation_legend=F
  )

ggsave("../pdf/Figure2/p2.pdf", plot = p2, height = 3, width = 7, units = "in")



bg3d(color = "white")
par3d(userMatrix = rotationMatrix(-pi/6, -1, 1, 0), zoom = 0.6)
acr.list <- c("MOs","PL","ORBm","ACAd","ILA","DP","ACAv")

for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  to.del <- which(mesh$vb[1,] < 0)
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  col <- "lightgray"
  wire3d(mesh, col = col, material = list(lit=FALSE), alpha = 0.2)
}

df_plot <- sp.PFC@meta.data
df_plot$AP_new <- df_plot$AP_new + runif(n= length(df_plot$AP_new), min= -1 , max= 1)*0.055555
df_plot <- df_plot[which(df_plot$`ACB-I`>0 & df_plot$ABA_hemisphere=="Left"),]
for (i in c(1:13)){
  idx_cluster = rownames(df_plot)[which(
    df_plot$SubType==names(col_SubType)[i]
  )]
  spheres3d(x = df_plot[idx_cluster,]$ML_new,
            y = df_plot[idx_cluster,]$DV_new,
            z = df_plot[idx_cluster,]$AP_new,
            col = col_SubType[i], radius=0.01, alpha=1)
  spheres3d(x = df_plot[idx_cluster,]$ML_new,
            y = df_plot[idx_cluster,]$DV_new,
            z = -1.5,
            col = col_SubType[i], radius=0.01, alpha=1)
}

hdr <- ggdensity::get_hdr(
  data.frame(x=df_plot$ML_new, y=df_plot$AP_new),
  rangex = c(0,3),
  rangey = c(-1.5,3)
)
hdr_data <- tidyr::spread(hdr$df_est[,c("x","y","fhat")], y, fhat, fill = 0)
rownames(hdr_data) <- hdr_data$x
hdr_data <- hdr_data[,-1]

z <- as.matrix(hdr_data)
z <- z - quantile(z[z>0], probs=0.90)
z[z<0] <- 0
z <- round(z*100)
x <- as.numeric(rownames(hdr_data))
y <- as.numeric(colnames(hdr_data))

zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- colorRampPalette(
  c(rep("black",1), sciRcolor::pal_scircolor(100)))(zlen)

col <- colorlut[z - zlim[1] + 1]
z_bottom <- matrix(-5, nrow=length(x), ncol=length(y))
surface3d(x, z_bottom, y, color = col)
light3d(theta=0, phi=90)

axes3d(labels = FALSE, tick=FALSE, nticks = 0, lwd=3)
grid3d(c("x","y","z"), lwd=2)

par3d(userMatrix = rotationMatrix(-pi/4, -1/10, 1, 0), zoom = 0.65, lwd = 2)
rgl.snapshot('../pdf/Figure2/p3_1.png', top = TRUE)



bg3d(color = "white")
par3d(userMatrix = rotationMatrix(-pi/6, -1, 1, 0), zoom = 0.6)
acr.list <- c("MOs","PL","ORBm","ACAd","ILA","DP","ACAv")

for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  to.del <- which(mesh$vb[1,] < 0)
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  col <- "lightgray"
  wire3d(mesh, col = col, material = list(lit=FALSE), alpha = 0.2)
}

df_plot <- sp.PFC@meta.data
df_plot$AP_new <- df_plot$AP_new + runif(n= length(df_plot$AP_new), min= -1 , max= 1)*0.055555
df_plot <- df_plot[which(df_plot$`SSp-I`>0 & df_plot$ABA_hemisphere=="Left"),]
for (i in c(1:13)){
  idx_cluster = rownames(df_plot)[which(
    df_plot$SubType==names(col_SubType)[i]
  )]
  spheres3d(x = df_plot[idx_cluster,]$ML_new,
            y = df_plot[idx_cluster,]$DV_new,
            z = df_plot[idx_cluster,]$AP_new,
            col = col_SubType[i], radius=0.01, alpha=1)
  spheres3d(x = df_plot[idx_cluster,]$ML_new,
            y = df_plot[idx_cluster,]$DV_new,
            z = -1.5,
            col = col_SubType[i], radius=0.01, alpha=1)
}

hdr <- ggdensity::get_hdr(
  data.frame(x=df_plot$ML_new, y=df_plot$AP_new),
  rangex = c(0,3),
  rangey = c(-1.5,3)
)
hdr_data <- tidyr::spread(hdr$df_est[,c("x","y","fhat")], y, fhat, fill = 0)
rownames(hdr_data) <- hdr_data$x
hdr_data <- hdr_data[,-1]

z <- as.matrix(hdr_data)
z <- z - quantile(z[z>0], probs=0.90)
z[z<0] <- 0
z <- round(z*100)
x <- as.numeric(rownames(hdr_data))
y <- as.numeric(colnames(hdr_data))

zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- colorRampPalette(
  c(rep("black",1), sciRcolor::pal_scircolor(100)))(zlen)

col <- colorlut[z - zlim[1] + 1]
z_bottom <- matrix(-5, nrow=length(x), ncol=length(y))
surface3d(x, z_bottom, y, color = col)
light3d(theta=0, phi=90)

axes3d(labels = FALSE, tick=FALSE, nticks = 0, lwd=3)
grid3d(c("x","y","z"), lwd=2)

par3d(userMatrix = rotationMatrix(-pi/4, -1/10, 1, 0), zoom = 0.65, lwd = 2)
rgl.snapshot('../pdf/Figure2/p3_2.png', top = TRUE)



IT_SubType <- c("L2/3_IT_1", "L4/5_IT_1", "L5_IT_1", "L6_IT_1",
                "L2/3_IT_2","L4/5_IT_2", "L5_IT_2", "L6_IT_2")
seu <- subset(Adult.IT.PT.barcode, cells=colnames(Adult.IT.PT.barcode)[which(Adult.IT.PT.barcode$SubType %in% IT_SubType)])
mat <- matrix(nrow = 2, ncol = length(IT_SubType))
rownames(mat) <- c("SSp-I", "ACB-I")
colnames(mat) <- IT_SubType

for (i in 1:2){
  mat[i,] <- as.numeric(table(seu$SubType[which(seu@meta.data[,rownames(mat)[i]]>0)])[IT_SubType])
}
mat <- mat/rowSums(mat)
mat <- as.data.frame(mat)
mat$Target <- rownames(mat)

links <- pivot_longer(mat, !Target, names_to = "SubType", values_to = "value")
colnames(links) <- c('target', 'source', 'value')

links <- links[which(links$value != 0),]
nodes <- c(IT_SubType,"SSp-I","ACB-I")
nodes <- data.frame(name=nodes)
nodes$index <- 0:(nrow(nodes) - 1)
links <- merge(links, nodes, by.x="source", by.y="name")
links <- merge(links, nodes, by.x="target", by.y="name")
names(links) <- c("target","source","Value","IDsource","IDtarget")

nodes.colour <- c(
  "L2/3_IT_1"='#ffd600',"L4/5_IT_1"='#ff6d00',"L5_IT_1"='#0091ea',"L6_IT_1"='#c51162',"L2/3_IT_2"='#ffff8d',"L4/5_IT_2"='#ffd180',"L5_IT_2"='#80d8ff',"L6_IT_2"='#ff80ab','SSp-I'="#1f77b4", 'ACB-I'="#ff7f0e")
pastecolor <- paste('d3.scaleOrdinal() .domain(["', nodes$name[1], sep = '')
for (i in 2:length(nodes$name)){
  pastecolor <- paste(pastecolor, '", "', nodes$name[i], sep = '')
}
pastecolor <- paste(pastecolor, '"]) .range(["', sep = '')
pastecolor <- paste(pastecolor, nodes.colour[1], sep = '')
for (i in 2:length(nodes.colour)){
  pastecolor <- paste(pastecolor,'", "', nodes.colour[i], sep = '')
}
pastecolor <- paste(pastecolor,'"])', sep = '')
colourScale <- pastecolor

links$Group <- links$target
links$Group <- as.factor(links$Group)
colnames(links) <- c("source", "target", "Value", "IDtarget", "IDsource", "Group")

p4 <-
  sankeyNetwork(Links=links, Nodes=nodes, Source="IDsource", Target="IDtarget",
                Value="Value", NodeID="name", fontSize=20,
                nodeWidth=30, nodePadding=10, margin=NULL,
                height=600, width=400, sinksRight=TRUE,
                colourScale=colourScale, LinkGroup="Group",iterations=0)

saveNetwork(p4,"../pdf/Figure2/p4.html")



seu <- Adult.IT.PT.barcode
BLA_I <- colnames(seu)[which(seu$`ACB-I`>0 & seu$`SSp-I`==0)]
RSP_I <- colnames(seu)[which(seu$`SSp-I`>0 & seu$`ACB-I`==0)]
DEGs <- FindMarkers(seu, ident.1 = BLA_I, ident.2 = RSP_I, logfc.threshold = 0.5)

DEGs$label <- ""
top5_gene <- rownames(DEGs)[which(DEGs$avg_log2FC>0.5 & DEGs$p_val_adj<1e-2)]
if(length(top5_gene)>5){
  top5_gene <- top5_gene[1:5]
}
down5_gene <- rownames(DEGs)[which(DEGs$avg_log2FC < -0.5 &
                                     DEGs$p_val_adj < 1e-2)]
if(length(down5_gene)>5){
  down5_gene <- down5_gene[1:5]
}

DEGs$label[match(top5_gene, rownames(DEGs))] <- top5_gene
DEGs$label[match(down5_gene, rownames(DEGs))] <- down5_gene
DEGs$Type <- 'not significant'
DEGs$Type[which(DEGs$avg_log2FC>0.5 & DEGs$p_val_adj<1e-2)] <- "Up"
DEGs$Type[which(DEGs$avg_log2FC < -0.5 & DEGs$p_val_adj<1e-2)] <- "Down"

p5 <-
  ggplot(DEGs, aes(x=avg_log2FC, y= -log10(p_val_adj))) +
  geom_point(aes(color=Type), size=0.5) +
  geom_vline(aes(xintercept=0.5), colour="black", linetype="dashed",
             linewidth = 0.5) +
  geom_vline(aes(xintercept = -0.5), colour="black", linetype="dashed",
             linewidth = 0.5) +
  geom_text_repel(aes(label=label, color=Type), size=3, max.overlaps=100) +
  theme_classic() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        text = element_text(size = 15), legend.position = "none") +
  scale_color_manual(values = c("Up"='#c82423', "Down"='#2878b5',
                                'not significant'='lightgray')) +
  xlim(min(DEGs$avg_log2FC), max(DEGs$avg_log2FC)) +
  ylim(0,100) +
  labs(title = "", x = 'log2FC', y = '-log10(P value)')

ggsave("../pdf/Figure2/p5.pdf", plot = p5, height = 3, width = 4, units = "in")



BC_order <- c("AId-C","ECT-C","ENTl-C","ENTl-I","DR-I","LHA-I","VTA-I","RE-I","SC-I","MD-I","AId-I","PL-C","VIS-I","BLA-I","BLA-C","ECT-I","SSp-I","ACB-I","CP-I","ACB-C","RSP-I","RSP-C","AUD-I","CP-C")

Ex_BC_mat <- Adult.Ex.barcode@meta.data[,BC_order]
Ex_BC_mat[is.na(Ex_BC_mat)] <- 0
annotation_row <- data.frame(
  Subtype = factor(Adult.Ex.barcode$Proj_subtype, levels = 1:33),
  Module = factor(Adult.Ex.barcode$Proj_module,
                  levels = c('ITi-D', 'ITi-V', 'ITc', 'PTi'))
)
rownames(annotation_row) <- rownames(Ex_BC_mat)
annotation_colors = list(
  Subtype = col_Proj_subtype,
  Module = col_Proj_module
)

p6 <-
  pheatmap(Ex_BC_mat[,BC_order],
           cluster_rows = T, cluster_cols = F,
           show_rownames=F,
           color = colorRampPalette(c("white", "red"))(50),
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           fontsize = 10)

ggsave("../pdf/Figure2/p6.pdf", plot = p6, height = 10, width = 12, units = "in")



BC_order <- c('CP-I','VIS-I','SSp-I','RSP-I','AUD-I',
              'ACB-I','AId-I','BLA-I','ECT-I','ENTl-I',
              'AId-C','ECT-C','ENTl-C','PL-C','ACB-C','CP-C','BLA-C','RSP-C',
              'DR-I','SC-I','LHA-I','RE-I','MD-I','VTA-I')

seu <- Adult.Ex.barcode
seu@meta.data[,BC_order][is.na(seu@meta.data[,BC_order])] <- 0
seu$Proj_subtype <- factor(
  seu$Proj_subtype,
  levels = c(1,14,17,22,23,25,26,27,28,29,31,
             7,8,9,10,15,16,18,19,20,21,
             2,3,11,12,13,24,30,32,33,
             4,5,6))

p7 <-
  DotPlot(seu, features = BC_order, group.by = 'Proj_subtype',
          cols = c("lightgrey", "red")) +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  labs(x='', y='')

ggsave("../pdf/Figure2/p7.pdf", plot = p7, height = 5, width = 10, units = "in")



Barcode <- c('VIS-I','SSp-I','CP-I','AUD-I','RSP-I',
             'BLA-I','ACB-I','ENTl-I','AId-I','ECT-I',
             'ACB-C','PL-C','ECT-C','ENTl-C',
             'BLA-C','CP-C','AId-C','RSP-C',
             'MD-I','RE-I','DR-I','VTA-I','LHA-I','SC-I')
cluster_order <- c(1,14,17,22,23,25,26,27,28,29,31,
                   7,8,9,10,15,16,18,19,20,21,
                   2,3,11,12,13,24,30,32,33,
                   4,5,6)
Ex_BC_mat <- Adult.Ex.barcode@meta.data[,Barcode]
Ex_BC_mat[is.na(Ex_BC_mat)] <- 0
Ex_BC_mat <- scale(Ex_BC_mat)
Ex_BC_mat <- Ex_BC_mat[Adult.Ex.barcode$BC_num > 1,]
set.seed(20240422)
umap_out <- umap(Ex_BC_mat)

umap_result <- as.data.frame(umap_out$layout)
colnames(umap_result) = c("UMAP_1","UMAP_2")
umap_result$Proj_subtype <-
  factor(Adult.Ex.barcode@meta.data[rownames(Ex_BC_mat),'Proj_subtype'],
         levels = cluster_order)
label <- umap_result |>
  group_by(Proj_subtype) |>
  dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

p8 <-
  ggplot(umap_result, aes(UMAP_1, UMAP_2, color=Proj_subtype, fill=Proj_subtype)) +
  geom_point(size=1, alpha=0.8) +
  scale_color_manual(values = col_Proj_subtype) +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) +
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  coord_fixed()

ggsave("../pdf/Figure2/p8.pdf", plot = p8, height = 5, width = 6, units = "in")



bg3d(color = "white")
par3d(userMatrix = rotationMatrix(-pi/6, -1, 1, 0), zoom = 0.6)
acr.list <- c("MOs","PL","ORBm","ACAd","ILA","DP","ACAv")

for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  to.del <- NA
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  col <- "lightgray"
  wire3d(mesh, col = col, material = list(lit=FALSE), alpha = 0.2)
}

df_plot <- sp.PFC@meta.data
for (i in c(4)){
  idx_cluster = rownames(df_plot)[which(
    df_plot$Proj_module==names(col_Proj_module)[i]
  )]
  spheres3d(x = df_plot[idx_cluster,]$ML_new,
            y = df_plot[idx_cluster,]$DV_new,
            z = df_plot[idx_cluster,]$AP_new,
            col = col_Proj_module[i], radius=0.01, alpha=1)
}

rgl.snapshot('../pdf/Figure2/p9_PTi_white.png', top = TRUE)



seu <- Adult.IT.PT.barcode
cluster_order <- c('L2/3_IT_1', 'L4/5_IT_1', 'L5_IT_1', 'L6_IT_1',
                   'L2/3_IT_2', 'L4/5_IT_2','L5_IT_2', 'L6_IT_2',
                   'L5_PT_1',  'L5_PT_2')
seu$Ex_finetype <- factor(seu$SubType, levels = cluster_order)
seu$Proj_module <- factor(seu$Proj_module,
                          levels = c("ITi-M1","ITi-M2","ITc-M3","PTi"))
df <- table(seu$Proj_module, seu$Ex_finetype)
df_norm <- as.data.frame(apply(df,2,function(x){x/sum(x)}))
df_norm$Type <- rownames(df_norm)
df_norm_long <-
  df_norm |>
  pivot_longer(!Type, names_to = "Cluster", values_to = "Value")

df_norm_long$Cluster <- factor(
  df_norm_long$Cluster,
  levels = cluster_order)
df_norm_long$Type <- factor(df_norm_long$Type,
                            levels = c("ITi-M1","ITi-M2","ITc-M3","PTi"))

label <- round(df_norm_long$Value,2)
label[which(label<0.05)] <- ''
p10 <-
  ggplot(df_norm_long, aes(x = Cluster, y = Value, fill = Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = col_Proj_module) +
  geom_text(label = label, color = "black", size = 3,
            position = position_stack(0.5)) +
  labs(x='', y='', title='') +
  theme_classic() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"))

ggsave("../pdf/Figure2/p10.pdf", plot = p10, height = 4, width = 6, units = "in")



seu <- subset(sp.PFC, cells = colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere == "Left")])
seu$D_V_enrich <- "other"
seu$D_V_enrich[which(seu$SubType %in% c("L2/3_IT_1", "L4/5_IT_1", "L5_IT_1", "L6_IT_1"))] <- "dmPFC-enrich IT"
seu$D_V_enrich[which(seu$SubType %in% c("L2/3_IT_2", "L4/5_IT_2", "L5_IT_2", "L6_IT_2"))] <- "vmPFC-enrich IT"

df <- seu@meta.data[which(seu$slice=="IT_slice_15" & seu$D_V_enrich != "other"),
                    c("ML_new", "DV_new","D_V_enrich")]

p11 <-
  ggplot(df) +
  geom_point_blur(data = df[which(df$D_V_enrich=="vmPFC-enrich IT"),],
                  aes(x=ML_new, y=DV_new),
                  color="#fe02fd", size=2) +
  geom_point_blur(data = df[which(df$D_V_enrich=="dmPFC-enrich IT"),],
                  aes(x=ML_new, y=DV_new),
                  color="#02ffff", size=2) +
  ggdark::dark_theme_void() +
  coord_fixed() +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

ggsave("../pdf/Figure2/p11.png", plot = p11, height = 8, width = 4.58, units = "in")



seu <- subset(sp.PFC, cells = colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere == "Left")])
df <- seu@meta.data[which(seu$slice=="IT_slice_15" & seu$Proj_module %in% c("ITi-M1","ITi-M2")), c("ML_new", "DV_new","Proj_module")]
p12 <-
  ggplot(df) +
  geom_point_blur(data = df[which(df$Proj_module=="ITi-M2"),],
                  aes(x=ML_new, y=DV_new),
                  color="#ff7f0e", size=2) +
  geom_point_blur(data = df[which(df$Proj_module=="ITi-M1"),],
                  aes(x=ML_new, y=DV_new),
                  color="#1f77b4", size=2) +
  ggdark::dark_theme_void() +
  coord_fixed() +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

ggsave("../pdf/Figure2/p12.png", plot = p12, height = 8, width = 4.39, units = "in")



seu <- subset(sp.PFC, cells = colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere == "Left")])
seu$D_V_enrich <- "other"
seu$D_V_enrich[which(seu$SubType %in% c("L2/3_IT_1", "L4/5_IT_1", "L5_IT_1", "L6_IT_1"))] <- "dmPFC-enrich cluster"
seu$D_V_enrich[which(seu$SubType %in% c("L2/3_IT_2", "L4/5_IT_2", "L5_IT_2", "L6_IT_2"))] <- "vmPFC-enrich cluster"

bc_slice <- seu@meta.data[
  which(seu$D_V_enrich !="other" & seu$Proj_module %in% c("ITi-M1","ITi-M2")),
  c("D_V_enrich","Proj_module", 'Y')]
bc_slice <-
  bc_slice |>
  mutate(bin = cut(Y, breaks = 36))
bin <- sort(unique(bc_slice$bin))
bc_slice$bin_index <- match(bc_slice$bin, bin)

df <- data.frame('bin_index'=c(1:36))
for (i in 1:36){
  df$`dmPFC-enrich-cluster`[i] <- length(which(bc_slice$bin_index==i & bc_slice$D_V_enrich=="dmPFC-enrich cluster"))/
    length(which(bc_slice$bin_index==i))
  df$`vmPFC-enrich-cluster`[i] <- length(which(bc_slice$bin_index==i & bc_slice$D_V_enrich=="vmPFC-enrich cluster"))/
    length(which(bc_slice$bin_index==i))
  df$`ITi-M1`[i] <- length(which(bc_slice$bin_index==i & bc_slice$Proj_module=="ITi-M1"))/
    length(which(bc_slice$bin_index==i))
  df$`ITi-M2`[i] <- length(which(bc_slice$bin_index==i & bc_slice$Proj_module=="ITi-M2"))/
    length(which(bc_slice$bin_index==i))
}
df$x <- c(1:36)
df <- pivot_longer(df, 2:5, names_to = "type", values_to = "value")
df$type <- factor(df$type, levels = c("dmPFC-enrich-cluster","vmPFC-enrich-cluster","ITi-M1","ITi-M2"))

p13 <-
  ggplot(df, aes(x=x, y=value, color=type)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(se = F, linewidth=1.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,35,5)) +
  scale_color_manual(values = c("dmPFC-enrich-cluster"="#02ffff",
                                "vmPFC-enrich-cluster"="#fe02fd",
                                "ITi-M1"="#1f77b4",
                                "ITi-M2"="#ff7f0e")) +
  theme(text = element_text(size=15),  legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5)) +
  labs(x='D → V',y='Cell proportion') +
  xlim(36, 0)

ggsave("../pdf/Figure2/p13.pdf", plot = p13, height = 3, width = 6, units = "in")



