# --------------------
# title: FigureS8 Code
# author: Hu Zheng
# date: 2026-01-01
# --------------------

library(Seurat)
library(scCustomize)
library(tidyverse)
library(umap)
library(tidydr)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(viridis)
library(sciRcolor)
library(scRNAtoolVis)
library(networkD3)

source('bin/Palettes.R')

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



seu <- subset(sp.PFC, cells=colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere=="Left")])
slice <- unique(seu$slice)
df <- data.frame('slice'=slice)
for (i in 1:length(slice)){
  df$`ACB-I`[i] <- length(which(seu$slice==slice[i] &
                                  seu$`ACB-I`>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
  df$`SSp-I`[i] <- length(which(seu$slice==slice[i] &
                                  seu$`SSp-I`>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
}
df$x <- c(1:36)
df <- pivot_longer(df, 2:3, names_to = "target", values_to = "value")
p1_1 <-
  ggplot(df, aes(x=x, y=value, color=target)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(se = F, linewidth=1.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,35,5)) +
  scale_color_manual(values = c('SSp-I'="#1f77b4", 'ACB-I'="#ff7f0e")) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  labs(x='A → P',y='Cell proportion')

ggsave("../pdf/FigureS8/p1_1.pdf", plot = p1, height = 2.5, width = 6, units = "in")



barcode <- c("VIS-I","ACB-I","AId-I","CP-I","ECT-I",
             "BLA-I","AUD-I","SSp-I","RSP-I","LHA-I")
seu <- subset(sp.PFC, cells=colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere=="Left")])
bc_slice <- seu@meta.data[,c(barcode, 'Y','BC_num')]
bc_slice <-
  bc_slice |>
  mutate(bin = cut(Y, breaks = 36))
bin <- sort(unique(bc_slice$bin))
bc_slice$bin_index <- match(bc_slice$bin, bin)

df <- data.frame('bin_index'=c(1:36))
for (i in 1:36){
  df$`ACB-I`[i] <- length(which(bc_slice$bin_index==i & bc_slice$`ACB-I`>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
  df$`SSp-I`[i] <- length(which(bc_slice$bin_index==i &
                                  bc_slice$`SSp-I`>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
}
df$x <- c(1:36)
df <- pivot_longer(df, 2:3, names_to = "target", values_to = "value")
p1_2 <-
  ggplot(df, aes(x=x, y=value, color=target)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(se = F, linewidth=1.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,35,5)) +
  scale_color_manual(values = c('SSp-I'="#1f77b4", 'ACB-I'="#ff7f0e")) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  labs(x='V → D',y='Cell proportion')

ggsave("../pdf/FigureS8/p1_2.pdf", plot = p2, height = 2.5, width = 6, units = "in")



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
df_plot <- df_plot[which(df_plot$`CP-I`>0 & df_plot$ABA_hemisphere=="Left"),]
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
rgl.snapshot('../pdf/FigureS8/p2.png', top = TRUE)



seu <- subset(sp.PFC, cells=colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere=="Left")])
slice <- unique(seu$slice)
df <- data.frame('slice'=slice)
for (i in 1:length(slice)){
  df$`ACB-I`[i] <- length(which(seu$slice==slice[i] &
                                  seu$`ACB-I`>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
  df$`CP-I`[i] <- length(which(seu$slice==slice[i] &
                                 seu$`CP-I`>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
}
df$x <- c(1:36)
df <- pivot_longer(df, 2:3, names_to = "target", values_to = "value")

p3_1 <-
  ggplot(df, aes(x=x, y=value, color=target)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(se = F, linewidth=1.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,35,5)) +
  scale_color_manual(values = c('ACB-I'="#ff7f0e", 'CP-I'="#1f77b4")) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  labs(x='A → P',y='Cell proportion')

ggsave("../pdf/FigureS8/p3_1.pdf", plot = p3_1, height = 2.5, width = 6, units = "in")



barcode <- c("VIS-I","ACB-I","AId-I","CP-I","ECT-I",
             "BLA-I","AUD-I","SSp-I","RSP-I","LHA-I")
seu <- subset(sp.PFC, cells=colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere=="Left")])
bc_slice <- seu@meta.data[,c(barcode, 'Y','BC_num')]
bc_slice <-
  bc_slice |>
  mutate(bin = cut(Y, breaks = 36))
bin <- sort(unique(bc_slice$bin))
bc_slice$bin_index <- match(bc_slice$bin, bin)

df <- data.frame('bin_index'=c(1:36))
for (i in 1:36){
  df$`ACB-I`[i] <- length(which(bc_slice$bin_index==i & bc_slice$`ACB-I`>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
  df$`CP-I`[i] <- length(which(bc_slice$bin_index==i &
                                 bc_slice$`CP-I`>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
}
df$x <- c(1:36)
df <- pivot_longer(df, 2:3, names_to = "target", values_to = "value")
p3_2 <-
  ggplot(df, aes(x=x, y=value, color=target)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(se = F, linewidth=1.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,35,5)) +
  scale_color_manual(values = c('ACB-I'="#ff7f0e", 'CP-I'="#1f77b4")) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  labs(x='V → D',y='Cell proportion')

ggsave("../pdf/FigureS8/p3_2.pdf", plot = p3_2, height = 2.5, width = 6, units = "in")



IT_SubType <- c("L2/3_IT_1", "L4/5_IT_1", "L5_IT_1", "L6_IT_1",
                "L2/3_IT_2","L4/5_IT_2", "L5_IT_2", "L6_IT_2")
seu <- subset(Adult.IT.PT.barcode, cells=colnames(Adult.IT.PT.barcode)[which(Adult.IT.PT.barcode$SubType %in% IT_SubType)])
mat <- matrix(nrow = 2, ncol = length(IT_SubType))
rownames(mat) <- c("CP-I", "ACB-I")
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
nodes <- c(IT_SubType,"CP-I","ACB-I")
nodes <- data.frame(name=nodes)
nodes$index <- 0:(nrow(nodes) - 1)
links <- merge(links, nodes, by.x="source", by.y="name")
links <- merge(links, nodes, by.x="target", by.y="name")
names(links) <- c("target","source","Value","IDsource","IDtarget")

nodes.colour <- c(
  "L2/3_IT_1"='#ffd600',"L4/5_IT_1"='#ff6d00',"L5_IT_1"='#0091ea',"L6_IT_1"='#c51162',"L2/3_IT_2"='#ffff8d',"L4/5_IT_2"='#ffd180',"L5_IT_2"='#80d8ff',"L6_IT_2"='#ff80ab',"CP-I"="#1f77b4","ACB-I"="#ff7f0e")
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

saveNetwork(p4,"../pdf/FigureS8/p4.html")



slice <- 'IT_slice_05'
barcode <- 'PTi'
limits <- c(0,2)
seu <- subset(sp.PFC, cells=colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere=="Left")])
seu$PTi[is.na(seu$PTi)] <- 0
seu$ITi_M1[is.na(seu$ITi_M1)] <- 0
seu$ITi_M2[is.na(seu$ITi_M2)] <- 0
seu$ITC_M3[is.na(seu$ITC_M3)] <- 0

df <- data.frame(
  X = seu$X,
  Y = seu$Y,
  Zscore = scale(log1p(seu@meta.data[,barcode]))
)
df <- df[which(seu$slice==slice),]
df$Zscore[df$Zscore<limits[1]] <- limits[1]
df$Zscore[df$Zscore>limits[2]] <- limits[2]
df <- df[order(df$Zscore),]

p5 <-
  ggplot(df, aes(x=X,y=Y)) +
  geom_point(aes(colour=Zscore), size=1) +
  scale_color_gradientn(colours = viridis(n = 256, option = "E", direction = 1),
                        limits = limits) +
  ggdark::dark_theme_void() +
  theme(plot.title = element_blank(), legend.position = "none") +
  coord_fixed()

ggsave("../pdf/FigureS8/p5_PTi_slice05.png", plot = p5, height = 4, width = 3, units = "in")



seu <- sp.PFC
slice <- unique(seu$slice)
df <- data.frame('slice'=slice)
for (i in 1:length(slice)){
  df$`ITi-M1`[i] <- length(which(seu$slice==slice[i] &
                                   seu$ITi_M1>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
  df$`ITi-M2`[i] <- length(which(seu$slice==slice[i] &
                                   seu$ITi_M2>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
  df$`ITc-M3`[i] <- length(which(seu$slice==slice[i] &
                                   seu$ITC_M3>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
  df$`PTi`[i] <- length(which(seu$slice==slice[i] &
                                seu$PTi>0))/
    length(which(seu$slice==slice[i] & seu$BC_num>0))
}
df$x <- c(1:36)
df <- pivot_longer(df, 2:5, names_to = "target", values_to = "value")
p6_1 <-
  ggplot(df, aes(x=x, y=value, color=target)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(se = F, linewidth=1.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,35,5)) +
  scale_color_manual(values = col_Proj_module) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  labs(x='A → P',y='Cell proportion')

ggsave("../pdf/FigureS8/p6_1.pdf", plot = p6_1, height = 2.5, width = 6, units = "in")



bc_slice <- seu@meta.data[,c(c("ITi_M1", "ITi_M2", "ITC_M3", "PTi"), 'Y','BC_num')]
bc_slice <-
  bc_slice |>
  mutate(bin = cut(Y, breaks = 36))
bin <- sort(unique(bc_slice$bin))
bc_slice$bin_index <- match(bc_slice$bin, bin)

df <- data.frame('bin_index'=c(1:36))
for (i in 1:36){
  df$`ITi-M1`[i] <- length(which(bc_slice$bin_index==i & bc_slice$ITi_M1>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
  df$`ITi-M2`[i] <- length(which(bc_slice$bin_index==i & bc_slice$ITi_M2>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
  df$`ITc-M3`[i] <- length(which(bc_slice$bin_index==i & bc_slice$ITC_M3>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
  df$`PTi`[i] <- length(which(bc_slice$bin_index==i & bc_slice$PTi>0))/
    length(which(bc_slice$bin_index==i & bc_slice$BC_num>0))
}
df$x <- c(1:36)
df <- pivot_longer(df, 2:5, names_to = "target", values_to = "value")
p6_2 <-
  ggplot(df, aes(x=x, y=value, color=target)) +
  geom_point(alpha=0.5, size=3) +
  geom_smooth(se = F, linewidth=1.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,35,5)) +
  scale_color_manual(values = col_Proj_module) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  labs(x='D → V',y='Cell proportion') +
  xlim(36, 0)

ggsave("../pdf/FigureS8/p6_2.pdf", plot = p6_2, height = 2.5, width = 6, units = "in")



Barcode_order <- c('MD-I','RE-I','DR-I','VTA-I','LHA-I','SC-I',
                   'VIS-I','SSp-I','CP-I','AUD-I','RSP-I',
                   'BLA-I','ACB-I','ENTl-I','AId-I','ECT-I',
                   'ACB-C','PL-C','ECT-C','ENTl-C',
                   'BLA-C','CP-C','AId-C','RSP-C'
)

Ex_BC_mat <- Adult.IT.PT.barcode@meta.data[,c(Barcode_order,'Proj_module')]
BC_module_mat <-
  Ex_BC_mat |>
  group_by(Proj_module) |>
  dplyr::summarize(across(1:24, ~ mean(.x, na.rm = TRUE))) |>
  mutate_all(~replace(., is.na(.), 0))
BC_module_zscore <- scale(BC_module_mat[,2:25])
rownames(BC_module_zscore) <- BC_module_mat$Proj_module

breaks <- seq(0,1,0.01)
p7 <-
  pheatmap(BC_module_zscore[c('PTi','ITi-M1','ITi-M2','ITc-M3'),Barcode_order],
           cluster_rows = F, cluster_cols = F,
           show_colnames=T, show_rownames = T,
           breaks = breaks,
           color = colorRampPalette(c("white","firebrick3"))(length(breaks))
  )

ggsave("../pdf/FigureS8/p7.pdf", plot = p7, height = 4, width = 6, units = "in")



UMAP <- Adult.IT.PT.barcode@reductions$umap@cell.embeddings
seu <- subset(Adult.IT.PT.barcode,
              cells = colnames(Adult.IT.PT.barcode)[which(UMAP[,'UMAP_1']<1)])

p1 <- Plot_Density_Custom(
  seurat_object = seu, features = 'ITi_M1_score',
  custom_palette = colorRampPalette(c("#eeeded","#f0de36","#d71313"))(100)) +
  theme_void() +
  theme(plot.title = element_text(size = 30, hjust = 0.5)) +
  labs(title = "ITi-M1")

p2 <- Plot_Density_Custom(
  seurat_object = seu, features = 'ITi_M2_score',
  custom_palette = colorRampPalette(c("#eeeded","#f0de36","#d71313"))(100)) +
  theme_void() +
  theme(plot.title = element_text(size = 30, hjust = 0.5)) +
  labs(title = "ITi-M2")

p3 <- Plot_Density_Custom(
  seurat_object = seu, features = 'ITc_M3_score',
  custom_palette = colorRampPalette(c("#eeeded","#f0de36","#d71313"))(100)) +
  theme_void() +
  theme(plot.title = element_text(size = 30, hjust = 0.5)) +
  labs(title = "ITc-M3")

seu$PTi_score[is.nan(seu$PTi_score)] <- 0
p4 <- Plot_Density_Custom(
  seurat_object = seu, features = 'PTi_score',
  custom_palette = colorRampPalette(c("#eeeded","#f0de36","#d71313"))(100)) +
  theme_void() +
  theme(plot.title = element_text(size = 30, hjust = 0.5)) +
  labs(title = "PTi")

p <- plot_grid(p1,p2,p3,p4,ncol = 4)
ggsave("../pdf/FigureS8/p8.pdf", plot = p, height = 4, width = 16, units = "in")



#seu <- Adult.IT.PT.barcode
#seu$Proj_subtype <- factor(seu$Proj_subtype, levels = 1:33)
#Idents(seu) <- "Proj_subtype"
#DEGs <- FindAllMarkers(seu, logfc.threshold = 0.25, min.pct = 0.1)
#DEGs$p_val_adj[which(DEGs$p_val_adj==0)] <- 1e-290
#saveRDS(DEGs, '../data/rds/appeal/DEGs_Proj_subtype.rds')
DEGs <- readRDS('../data/rds/appeal/DEGs_Proj_subtype.rds')

DEGs$cluster <- factor(DEGs$cluster,
                       levels = c(1,14,17,22,23,25,26,27,28,29,31,
                                  7,8,9,10,15,16,18,19,20,21,
                                  2,3,11,12,13,24,30,32,33,
                                  4,5,6))
p9 <-
  jjVolcano(diffData = DEGs,
            topGeneN = 3,
            tile.col = rep("lightgray",33),
            aesCol = c("navy","firebrick3"),
            angle=0,
            size=4,
            pSize=1.5,
            log2FC.cutoff = 0.1,
            fontface = 'italic',
            legend.position = c(0.8,0.2),
            flip = F,
            min.segment.length = 0) +
  theme(legend.position = "none")

ggsave("../pdf/FigureS8/p9.pdf", plot = p9, height = 5, width = 20, units = "in")



#seu <- Adult.IT.PT.barcode
#Idents(seu) <- "Proj_module"
#DEGs <- FindAllMarkers(seu, logfc.threshold = 0.25, min.pct = 0.1)
#DEGs$p_val_adj[which(DEGs$p_val_adj==0)] <- 1e-290
#saveRDS(DEGs, '../data/rds/DEGs_Proj_module.rds')
DEGs <- readRDS('../data/rds/DEGs_Proj_module.rds')

proj_module <- c("ITi-D","ITi-V","ITc","PTi")
plist <- list()
for (i in 1:4){
  DEGs_module <- DEGs[DEGs$cluster == proj_module[i],]
  DEGs_module$label <- ""
  top5_gene <- DEGs_module$gene[which(DEGs_module$avg_log2FC>0.5 &
                                        DEGs_module$p_val_adj<1e-2)]
  if(length(top5_gene)>5){
    top5_gene <- top5_gene[1:5]
  }
  down5_gene <- DEGs_module$gene[which(DEGs_module$avg_log2FC< -0.5 &
                                         DEGs_module$p_val_adj<1e-2)]
  if(length(down5_gene)>5){
    down5_gene <- down5_gene[1:5]
  }

  DEGs_module$label[match(top5_gene, DEGs_module$gene)] <- top5_gene
  DEGs_module$label[match(down5_gene, DEGs_module$gene)] <- down5_gene
  DEGs_module$Type <- 'not significant'
  DEGs_module$Type[which(DEGs_module$avg_log2FC>0.5 &
                           DEGs_module$p_val_adj<1e-2)] <- "Up"
  DEGs_module$Type[which(DEGs_module$avg_log2FC < -0.5 &
                           DEGs_module$p_val_adj<1e-2)] <- "Down"

  plist[[i]] <-
    ggplot(DEGs_module, aes(x=avg_log2FC, y= -log10(p_val_adj))) +
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
    xlim(min(DEGs_module$avg_log2FC), max(DEGs_module$avg_log2FC)) +
    ylim(0,300) +
    labs(title = proj_module[i], x = 'log2FC', y = '-log10(P value)')
}
p10 <- plot_grid(plotlist = plist, ncol = 2)

ggsave("../pdf/FigureS8/p10.pdf", plot = p10, height = 6, width = 6, units = "in")



seu <- Adult.IT.PT.barcode
df <- table(seu$SubType_Layer, seu$Proj_subtype)
df_norm <- as.data.frame(apply(df, 2, function(x){x/sum(x)}))
df_norm$Layer <- rownames(df_norm)
df_plot <-
  df_norm |>
  pivot_longer(!Layer, names_to = "Proj_cluster", values_to = "value")
df_plot$Proj_cluster <- factor(
  df_plot$Proj_cluster,
  levels = names(sort(colSums(df_norm[1:2,1:(ncol(df_norm)-1)]), decreasing = T)))
df_plot$Layer <- factor(
  df_plot$Layer,
  levels = c("L2/3 IT","L4/5 IT","L5 IT","L6 IT","L5 PT"))

p11 <-
  ggplot(data=df_plot, aes(x=Proj_cluster, y=value, fill=Layer)) +
  geom_bar(stat="identity", width=0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = col_SubType_Layer) +
  labs(x='', y='')

ggsave("../pdf/FigureS8/p11.pdf", plot = p11, height = 4, width = 7, units = "in")



seu <- Adult.IT.PT.barcode
df <- table(seu$SubType, seu$Proj_subtype)
df_norm <- apply(df, 2, function(x){x/sum(x)})
df_plot <- as.data.frame(df_norm)
df_plot$Cluster <- rownames(df_plot)
df_plot <-
  df_plot |>
  pivot_longer(!Cluster, names_to = "Proj_cluster", values_to = "Value")
df_plot$Cluster <- factor(
  df_plot$Cluster,
  levels = c('L2/3_IT_1', 'L4/5_IT_1', 'L5_IT_1', 'L6_IT_1',
             'L2/3_IT_2', 'L4/5_IT_2', 'L5_IT_2', 'L6_IT_2',
             "L5_PT_1","L5_PT_2")
)
df_plot$Proj_cluster <- factor(
  df_plot$Proj_cluster,
  levels = names(sort(colSums(df_norm[c(1,3,5,9),]),decreasing = T))
)

p12 <-
  ggplot(data=df_plot, aes(x=Proj_cluster, y=Value, fill=Cluster)) +
  geom_bar(stat="identity", width=0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = col_SubType) +
  labs(x='', y='')

ggsave("../pdf/FigureS8/p12.pdf", plot = p12, height = 4, width = 7, units = "in")




