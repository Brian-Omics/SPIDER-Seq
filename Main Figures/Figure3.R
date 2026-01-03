# --------------------
# title: Figure3 Code
# author: Hu Zheng
# date: 2026-01-01
# --------------------

library(Seurat)
library(tidyverse)
library(ggsci)
library(aplot)
library(ggpointdensity)
library(scRNAtoolVis)
library(scCustomize)
library(viridis)
library(RColorBrewer)
library(cowplot)
library(ggradar)

library(Biorplot)
source('bin/Palettes.R')
source('bin/includes.R')

Adult.Ex <- readRDS('../data/rds/Adult.Ex.rds')
sp.PFC <- readRDS('../data/rds/sp.PFC.rds')
Adult.IT.PT.barcode <- subset(Adult.Ex, cells=colnames(Adult.Ex)[which(
  (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "IT") |
    (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "PT" & Adult.Ex$sample == "Adult1")
)])
Adult.IT.barcode <- subset(
  Adult.Ex,
  cells=colnames(Adult.Ex)[which(Adult.Ex$Ex_subtype=="IT" &
                                   Adult.Ex$BC_num>0)]
)
sp.PFC.Left.ITPT.barcode <- subset(sp.PFC, cells = colnames(sp.PFC)[which(
  sp.PFC$ABA_hemisphere=="Left" & sp.PFC$SubType_Layer %in% c("L2/3 IT","L4/5 IT","L5 IT","L6 IT", "L5 PT") & sp.PFC$BC_num>0)])



Barcode <- c('VIS-I','SSp-I','CP-I','AUD-I','RSP-I',
             'BLA-I','ACB-I','AId-I','ECT-I',
             'ACB-C','ECT-C',
             'CP-C','AId-C','RSP-C',
             'LHA-I')
seu <- Adult.IT.PT.barcode
seu$BC_motif <- apply(seu@meta.data[,Barcode], 1, function(x){
  paste(names(x)[which(x>0)], collapse = ',')
})

sc_motif <- table(seu$BC_motif)
sc_motif <- sc_motif[order(sc_motif, decreasing = T)]
sc_motif <- sc_motif[which(sc_motif>=10 & names(sc_motif) != "")]
sc_seu <- subset(seu, cells=colnames(seu)[which(seu$BC_motif %in% names(sc_motif))])

seu <- sp.PFC.Left.ITPT.barcode
seu$BC_motif <- apply(seu@meta.data[,Barcode], 1, function(x){
  paste(names(x)[which(x>0)], collapse = ',')
})
sp_motif <- table(seu$BC_motif)
sp_motif <- sp_motif[order(sp_motif, decreasing = T)]
sp_motif <- sp_motif[which(sp_motif>=10)]
sp_seu <- subset(seu, cells=colnames(seu)[which(seu$BC_motif %in% names(sp_motif))])

top50_motif <- names(sc_motif)
top50_motif <- top50_motif[which(!top50_motif %in% c("LHA-I","ACB-I,LHA-I"))]
top50_motif <- top50_motif[which(top50_motif %in% unique(sp_seu$BC_motif))]
top50_motif <- top50_motif[1:50]
sc_seu <- subset(sc_seu, cells=colnames(sc_seu)[which(sc_seu$BC_motif %in% top50_motif)])
sp_seu <- subset(sp_seu, cells=colnames(sp_seu)[which(sp_seu$BC_motif %in% top50_motif)])

df_transcriptom <- as.data.frame.array(table(sc_seu$BC_motif, sc_seu$SubType))
df_transcriptom <- df_transcriptom/rowSums(df_transcriptom)
df_order <- data.frame(
  "motifs" = rownames(df_transcriptom),
  "max_subtype" = apply(df_transcriptom, 1, function(x){
    colnames(df_transcriptom)[which.max(x)]
  }),
  "max_value" = apply(df_transcriptom, 1, function(x){
    max(x)
  })
)
df_order$max_subtype <- factor(
  df_order$max_subtype,
  levels = c("L2/3_IT_1","L2/3_IT_2","L4/5_IT_1","L4/5_IT_2","L5_IT_1","L5_IT_2",
             "L6_IT_1","L6_IT_2","L5_PT_1","L5_PT_2"))
df_order <- arrange(df_order, max_subtype, desc(max_value))
motifs_order <- df_order$motifs

sp_seu$BC_motif <- factor(sp_seu$BC_motif, levels = motifs_order)
df_spatial <- as.data.frame.array(table(sp_seu$BC_motif, sp_seu$ABA_metaRegion))
df_spatial <- df_spatial/rowSums(df_spatial)

Barcode_use <- c('VIS-I','SSp-I','CP-I','AUD-I','RSP-I',
                 'BLA-I','ACB-I','AId-I','ECT-I',
                 'ACB-C',
                 'CP-C','AId-C')
mat <- matrix(0, nrow = length(Barcode_use), ncol = length(motifs_order))
rownames(mat) <- Barcode_use
colnames(mat) <- motifs_order

for (j in 1:ncol(mat)){
  for (i in 1:nrow(mat)){
    if (rownames(mat)[i] %in% strsplit(colnames(mat)[j],',')[[1]]){
      mat[i,j] <- 1
    }
  }
}

df1 <- data.frame(
  'X'=rep(colnames(mat), each=nrow(mat)),
  'Y'=rep(rownames(mat), ncol(mat)),
  'value'=as.vector(mat)
)
df1$value <- factor(df1$value, levels = c(0,1))
df1$Y <- factor(df1$Y, levels = rev(Barcode_use))

df1 <- df1[which(df1$X %in% motifs_order),]
df1$X <- factor(df1$X, levels = motifs_order)

p1 <- ggplot(df1) +
  geom_point(aes(x=X, y=Y, colour=value), size=3.8) +
  theme_minimal() +
  scale_color_manual(values = c('#d7d8da','black')) +
  labs(x='', y='Projectome') +
  theme(panel.grid=element_blank(), plot.margin = margin(0,0,0,0),
        axis.ticks.x=element_blank(), axis.text.x = element_blank(),
        plot.title = element_blank(), legend.position = "none")

df_spatial$motifs <- rownames(df_spatial)
df2 <- pivot_longer(df_spatial, !motifs, names_to = "Region", values_to = "Value")
Region <- c("MOs","ILA","PL","ACAd","ORBm","ACAv", "DP")

df2$Value[is.nan(df2$Value)] <- 0
df2$motifs <- factor(df2$motifs, levels = motifs_order)
df2$Region <- factor(df2$Region, levels = c("MOs","ACAd","ACAv","PL","ORBm","ILA","DP"))

p2 <-
  ggplot(df2, aes(x = motifs, y = Value, fill = Region)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = col_subreagion) +
  labs(x='', y='', title='') +
  theme_classic() +
  theme(panel.grid=element_blank(), plot.margin = margin(0,0,0,0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        plot.title = element_blank()) +
  labs(y="Spatial region")

SubType <- c("L2/3_IT_1","L2/3_IT_2","L4/5_IT_1","L4/5_IT_2","L5_IT_1",
             "L5_IT_2","L6_IT_1","L6_IT_2")
df3 <- df_transcriptom[,SubType]
df3$motifs <- rownames(df3)
df3_long <- pivot_longer(df3, !motifs, names_to = "SubType",
                         values_to = "Value")
df3_long$motifs <- factor(df3_long$motifs, levels = motifs_order)
df3_long$SubType <- factor(df3_long$SubType, levels = rev(SubType))
df3_long$label <- ""
df3_long$label[which(df3_long$Value>0.25)] <- "*"

breaks <- seq(0,0.5,0.01)
p3 <- ggplot(df3_long, aes(x=motifs, y=SubType, fill=Value)) +
  geom_raster() +
  geom_text(aes(label=label),col ="black",size = 5) +
  scale_fill_gradientn(limits=c(0,0.5), colours = colorRampPalette(c("navy","white","firebrick3"))(100), na.value="firebrick3") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1)) +
  labs(x="", y="Transcriptome")

p <-
  p3 %>%
  insert_top(p1, height=1.3) %>%
  insert_bottom(p2, height=0.8)

ggsave("../pdf/Figure3/p1.pdf", plot = p, height = 5, width = 12, units = "in")



seu <- subset(
  Adult.Ex,
  cells=colnames(Adult.Ex)[which(Adult.Ex$Ex_subtype=="IT" &
                                   Adult.Ex$BC_num>0)])
seu$group <- "Other"
seu$group[which(seu$`AId-I`>0 & seu$`ACB-I`>0 & seu$`BLA-I`>0 & seu$BC_num==3)] <- "AId-I + ACB-I + BLA-I"
seu$group[which(seu$`AId-I`>0 & seu$`ACB-I`>0 & seu$BC_num==2)] <- "AId-I + ACB-I"
seu$group[which(seu$`AId-I`>0 & seu$`ACB-I`>0 & seu$`CP-I`>0 & seu$BC_num==3)] <- "AId-I + ACB-I + CP-I"
seu$group[which(seu$`AId-I`>0 & seu$`ACB-I`>0 & seu$`CP-I`>0 & seu$`ACB-C`>0 & seu$`CP-C`>0 & seu$`AId-C`>0 & seu$BC_num==6)] <- "AId-I,ACB-I,CP-I,AId-C,ACB-C,CP-C"
seu$group[which(seu$`AId-I`>0 & seu$`CP-I`>0 & seu$BC_num==2)] <- "AId-I + CP-I"

seu$group <- factor(seu$group, levels = c("AId-I + ACB-I + BLA-I","AId-I + ACB-I","AId-I + ACB-I + CP-I","AId-I,ACB-I,CP-I,AId-C,ACB-C,CP-C","AId-I + CP-I","Other"))

#Idents(seu) <- "group"
#DEGs <- FindAllMarkers(seu, logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)
#DEGs$Group <- "not"
#DEGs$Group[which(DEGs$avg_log2FC>0.5 & DEGs$p_val_adj < 0.01)] <- "Up"
#DEGs <- DEGs[which(DEGs$cluster != "Other"),]
#saveRDS(DEGs, '../data/rds/DEGs_AId_motifs.rds')
DEGs <- readRDS('../data/rds/DEGs_AId_motifs.rds')

p2_1 <-
  jjVolcano(diffData = DEGs,
            #myMarkers = mygene,
            topGeneN = 3,
            tile.col = rep("lightgray",6),
            aesCol = c("navy","firebrick3"),
            angle=0,
            size=3,
            pSize=1.5,
            log2FC.cutoff = 0.1,
            fontface = 'italic',
            legend.position = c(0.8,0.2),
            flip = F,
            min.segment.length = 0) +
  theme(legend.position = "none")

ggsave("../pdf/Figure3/p2_1.pdf", plot = p2_1, height = 3, width = 15, units = "in")

seu <- sp_seu
df <- seu@meta.data[,c("ML_new","DV_new","slice","BC_motif")]
# "BLA-I,ACB-I,AId-I"  "ACB-I,AId-I"  "CP-I,ACB-I,AId-I "  "CP-I,ACB-I,AId-I,ACB-C,CP-C,AId-C"  "CP-I,AId-I"
df_motif <- df[which(df$BC_motif=="BLA-I,ACB-I,AId-I"),]

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(85)

p2_2 <-
  ggplot(df, aes(ML_new, DV_new)) +
  geom_point(color="black", size=1) +
  stat_density2d(data=df_motif, geom = "raster",
                 aes(fill = ..density..),
                 contour = FALSE,na.rm = TRUE) +
  scale_fill_gradientn(colours = c(rep("transparent",15),colormap)) +
  theme_void() +
  #theme(legend.position = "none") +
  ylim(-5,0) +
  xlim(0,3) +
  coord_fixed()

ggsave("../pdf/Figure3/p2_2/BLA-I,ACB-I,AId-I.pdf", plot = p2_2,
       height = 6, width = 6, units = "in", bg="transparent")


AId_motif <- c("BLA-I,ACB-I,AId-I","ACB-I,AId-I","CP-I,ACB-I,AId-I",
               "CP-I,ACB-I,AId-I,ACB-C,CP-C,AId-C","CP-I,AId-I")
seu <- subset(sp_seu, cells = colnames(sp_seu)[which(sp_seu$BC_motif %in% AId_motif)])
seu$ABA_metaRegion <- factor(
  seu$ABA_metaRegion,
  levels = c("MOs","ACAd","ACAv","PL","ORBm","ILA","DP"))
mat <- table(seu$BC_motif, seu$ABA_metaRegion)
mat <- t(apply(mat, 1, function(x){(x-min(x))/(max(x)-min(x))}))
mat <- mat[AId_motif,]
mat <- as_tibble(mat, rownames = "group")

plist <- list()
for (i in 1:length(AId_motif)){
  plist[[i]] <-
    ggradar(mat[i,], fill = T,
            group.point.size = 2, group.line.width = 0.5,
            background.circle.colour = "white",
            group.colours = "blue") +
    theme(plot.title = element_text(hjust=0.5))
}
p2_3 <- plot_grid(plotlist=plist, ncol=5)

ggsave("../pdf/Figure3/p2_3.pdf", plot = p2_3, height = 3, width = 16, units = "in")



bg3d(color = "white")
par3d(userMatrix = rotationMatrix(0, 0, 0, 0), zoom = 0.7)
acr.list <- c("MOs","PL","ORBm","ACAd","ILA","DP","ACAv")

for(acr in acr.list){
  mesh <- mesh3d.allen.annot.from.id(id.from.acronym(acr))
  to.del <- which(mesh$vb[1,] < 0)
  mesh$it <- mesh$it[,!is.element(mesh$it[1,], to.del) & !is.element(mesh$it[2,], to.del) & !is.element(mesh$it[3,], to.del)]
  col <- "lightgray"
  wire3d(mesh, col = col, material = list(lit=FALSE), alpha = 0.2)
}

barcode <- c("VIS-I","ACB-I","CP-C","AId-I","CP-I","ECT-C","AId-C","ECT-I",
             "BLA-I","AUD-I","RSP-C","SSp-I","RSP-I","ACB-C","LHA-I")
cortex <- c("VIS-I","AId-I","ECT-C","AId-C","ECT-I",
            "BLA-I","AUD-I","RSP-C","SSp-I","RSP-I")
df_plot <- sp.PFC@meta.data
df_plot <- df_plot[df_plot$ABA_hemisphere=="Left" & df_plot$BC_num>0,]
df_plot$group <- ""
df_plot$group[which(df_plot$`ACB-I`>0 & df_plot$`CP-I`==0 &
                      rowSums(df_plot[,cortex],na.rm = T)>0)] <- "ACB-I + cortex"
df_plot$group[which(df_plot$`ACB-I`==0 & df_plot$`CP-I`>0 &
                      rowSums(df_plot[,cortex],na.rm = T)>0)] <- "CP-I + cortex"
df_plot$group[which(df_plot$`ACB-I`>0 & df_plot$`CP-I`>0 &
                      rowSums(df_plot[,cortex],na.rm = T)>0)] <- "ACB-I + CP-I + cortex"

proj_type <- c("CP-I + cortex", "ACB-I + cortex", "ACB-I + CP-I + cortex")
col <- c("red","green","blue")
for (i in c(1:3)){
  idx_cluster = rownames(df_plot)[which(
    df_plot$group==proj_type[i]
  )]
  spheres3d(x = df_plot[idx_cluster,]$ML_new,
            y = df_plot[idx_cluster,]$DV_new,
            z = df_plot[idx_cluster,]$AP_new,
            col = col[i], radius=0.01, alpha=1)
}

rgl.snapshot('../pdf/Figure3/p3_1.png', top = TRUE)


df <- df_plot[which(!df_plot$group==""),]
p3_2 <-
  ggplot(df, aes(DV_new)) +
  geom_density(aes(fill = group), alpha = 0.6) +
  scale_fill_manual(values = c("green","blue","red")) +
  labs(x="",y="") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y = element_blank(),axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("../pdf/Figure3/p3_2.pdf", plot = p3_2, height = 2, width = 8, units = "in")



sp_seu$BC_motif <- as.character(sp_seu$BC_motif)
spatial <- sp_seu@meta.data[,c("ML_new", "DV_new", "AP_new", "BC_motif")]
spatial <-
  spatial |>
  group_by(BC_motif) |>
  dplyr::summarize(across(1:3, ~ mean(.x, na.rm = TRUE)))

transcriptom <- as.data.frame(sc_seu@reductions$pca@cell.embeddings)
transcriptom$BC_motif <- sc_seu$BC_motif
transcriptom <-
  transcriptom |>
  group_by(BC_motif) |>
  dplyr::summarize(across(1:50, ~ mean(.x, na.rm = TRUE)))

projectom <- sc_seu@meta.data[,c(Barcode, "BC_motif")]
projectom <-
  projectom |>
  group_by(BC_motif) |>
  dplyr::summarize(across(1:15, ~ mean(.x, na.rm = TRUE))) |>
  mutate_all(~replace(., is.na(.), 0))

df <- data.frame(
  "spatial" = as.numeric(dist(spatial[,-1])),
  "transcriptom" = as.numeric(dist(transcriptom[,-1])),
  "projectom" = as.numeric(dist(projectom[,-1]))
)
df_norm <- as.data.frame(apply(df, 2, function(x){(x-min(x))/(max(x)-min(x))}))


cor <- cor.test(df_norm$spatial, df_norm$projectom, "two.sided", "pearson")
R <- round(cor$estimate,2)
P <- format(cor$p.value, digits = 2)

p4_1 <-
  ggplot(df_norm, aes(x=spatial, y=projectom)) +
  geom_pointdensity(adjust = 1, size=1) +
  geom_smooth(method = "lm", color='black', linewidth=0.5, se=F) +
  scale_color_distiller(palette = "Spectral", direction = -1) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 15),
        legend.position = "none") +
  labs(x="spatial distance", y="projectom distance",
       title = paste('R =',R,', P =',P,sep=' ')) +
  ylim(0,1)

ggsave("../pdf/Figure3/p4_1.pdf", plot = p4_1, height = 4, width = 4, units = "in")


cor <- cor.test(df_norm$transcriptom, df_norm$projectom, "two.sided", "pearson")
R <- round(cor$estimate,2)
P <- format(cor$p.value, digits = 2)

p4_2 <-
  ggplot(df_norm, aes(x=transcriptom, y=projectom)) +
  geom_pointdensity(adjust = 1, size=1) +
  geom_smooth(method = "lm", color='black', size=0.5, se=F) +
  scale_color_distiller(palette = "Spectral", direction = -1) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 15),
        legend.position = "none") +
  labs(x="transcriptom distance", y="projectom distance",
       title = paste('R =',R,', P =',P,sep=' ')) +
  ylim(0,1)

ggsave("../pdf/Figure3/p4_2.pdf", plot = p4_2, height = 4, width = 4, units = "in")




