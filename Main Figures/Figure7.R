# --------------------
# title: Figure7 Code
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

Adult.IT.PT.barcode <- subset(Adult.Ex, cells=colnames(Adult.Ex)[which(
  (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "IT") |
    (Adult.Ex$BC_num>0 & Adult.Ex$Ex_subtype == "PT" & Adult.Ex$sample == "Adult1")
)])

sp.PFC.Left <- subset(
  sp.PFC,
  cells = colnames(sp.PFC)[which(sp.PFC$ABA_hemisphere=="Left")])

sp.PFC.Left.ITPT.barcode <- subset(sp.PFC, cells = colnames(sp.PFC)[which(
  sp.PFC$ABA_hemisphere=="Left" & sp.PFC$SubType_Layer %in% c("L2/3 IT","L4/5 IT","L5 IT","L6 IT", "L5 PT") & sp.PFC$BC_num>0)])

seu.inte <- readRDS('../data/rds/PFC.MERFISH.inte.rds')



result_all <- readRDS('../data/rds/ML/result_all.rds')
df <- result_all
df <- df[which(df$Experiment != "spatial"),]
df$target <- factor(df$target, levels = c("PTi","ITi-D","ITi-V","ITc"))
df$Experiment <- factor(
  df$Experiment,
  levels = c("transcriptom + spatial","transcriptom","shuffle"))

p1 <-
  ggplot(df, aes(x=ROC_1, y=ROC_2, color=Experiment)) +
  geom_line(linewidth=1) +
  facet_wrap(~target, nrow = 2) +
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

unique(paste(df$target,df$AUC,df$Experiment))

ggsave("../pdf/Figure7/p1.pdf", plot = p1, height = 5.8, width = 5, units = "in")



Barcode <- c('LHA-I','CP-I','ACB-I','CP-C')
result_all <- readRDS('../data/rds/ML/result_target.rds')
df <- result_all
df <- df[df$target %in% Barcode,]
df <- df[which(df$Experiment != "spatial"),]
df$target <- factor(df$target, levels = Barcode)
df$Experiment <- factor(
  df$Experiment,
  levels = c("transcriptom + spatial","transcriptom","shuffle"))

p2 <-
  ggplot(df, aes(x=ROC_1, y=ROC_2, color=Experiment)) +
  geom_line(linewidth=1) +
  #geom_abline(slope = 1,intercept = 0,lty="dashed",color='gray') +
  facet_wrap(~target, nrow = 1) +
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

unique(paste(df$target,df$AUC))

ggsave("../pdf/Figure7/p2.pdf", plot = p2, height = 6.75, width = 5, units = "in")



result_all <- readRDS('../paper/投稿-Cell/意见/Figure/reviewer1-Q9/Q9_2/result_all.rds')
df <- result_all
df <- df[which(df$Experiment != "spatial"),]
AId_motif <- c("BLA-I,ACB-I,AId-I","ACB-I,AId-I","CP-I,ACB-I,AId-I","CP-I,ACB-I,AId-I,ACB-C,CP-C,AId-C","CP-I,AId-I")
df$target <- factor(df$target, levels = AId_motif)
df$Experiment <- factor(
  df$Experiment,
  levels = c("transcriptom + spatial","transcriptom","shuffle"))

p3 <-
  ggplot(df, aes(x=ROC_1, y=ROC_2, color=Experiment)) +
  geom_line(linewidth=1) +
  facet_wrap(~target, nrow = 1) +
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

unique(paste(df$target,df$AUC,df$Experiment))

ggsave("../pdf/Figure7/p3.pdf", plot = p3, height = 3.5, width = 12, units = "in")



col <- c("L2/3 IT"='#ffda11', "L4/5 IT"='#f68d45', "L5 IT"='#26d5ff',
         "L6 IT"= '#f05a9e', "L5 ET" = '#6d92ff', "L5/6 NP"='#8976ff',
         "L6 CT" = '#e01fff')
p4_1 <-
  ggplot(df[which(df$sample=="MERFISH"),], aes(x=x,y=y,color=SubType)) +
  geom_point() +
  coord_fixed() +
  scale_color_manual(values = col) +
  labs(title = "MERFISH slice") +
  ggdark::dark_theme_void() +
  theme(text = element_text(size = 15),
        plot.title = element_blank(), legend.position = "none")
p4_1

ggsave(paste("../pdf/Figure7/p4_1.png",sep=""), plot = p4_1, height = 6, width = 6, units = "in")



col <- c("L2/3 IT"='#ffda11', "L4/5 IT"='#f68d45', "L5 IT"='#26d5ff',
         "L6 IT"= '#f05a9e', "L5 PT" = '#6d92ff', "L5 NP"='#8976ff',
         "L6 CT" = '#e01fff')
p4_2 <-
  ggplot(df[which(df$sample=="our"),], aes(x=x,y=y,color=SubType)) +
  geom_point() +
  coord_fixed() +
  scale_color_manual(values = col) +
  labs(title = "MERFISH slice") +
  ggdark::dark_theme_void() +
  theme(text = element_text(size = 15),
        plot.title = element_blank(), legend.position = "none")

ggsave(paste("../pdf/Figure7/p4_2.png",sep=""), plot = p4_2, height = 5, width = 5, units = "in")



sp_Barcode_thread <- c("VIS-I"=0.25,"SSp-I"=0.25,"CP-I"=0.5,"AUD-I"=0.25,"RSP-I"=0.5, "BLA-I"=0.25,"ACB-I"=0.5,"AId-I"=0.75,"ECT-I"=0.5,"ACB-C"=0.25,"ECT-C"=0.5,"CP-C"=0.5,"AId-C"=0.25,"RSP-C"=0.5,"LHA-I"=0.5)
pre_Barcode <- readRDS('../data/rds/ML/pre_Barcode.rds')
df_pre <- as.data.frame(pre_Barcode)
df_pre$CCF <- PFC.MERFISH@meta.data[MERFISH_cells,"CCF_metaRegion"]
df_true <- sp.PFC.Left@meta.data[which(sp.PFC.Left$slice=="IT_slice_18"),c('ML_new','DV_new','ABA_metaRegion',sp_Barcode)]
colnames(df_true)[1:3] <- c('x','y','CCF')
for (i in 1:length(sp_Barcode)){
  df_pre[which(df_pre[,sp_Barcode[i]]>sp_Barcode_thread[i]), sp_Barcode[i]] <- 1
  df_pre[which(df_pre[,sp_Barcode[i]]<=sp_Barcode_thread[i]), sp_Barcode[i]] <- 0
  df_true[which(df_true[,sp_Barcode[i]]>0), sp_Barcode[i]] <- 1
}

df_pre_CCF <-
  df_pre[,c(sp_Barcode,"CCF")] |>
  group_by(CCF) |>
  dplyr::summarize(across(1:15, ~ sum(.x))) |>
  as.data.frame()
df_pre_CCF_norm <- apply(df_pre_CCF[,sp_Barcode], 2, function(x){x/sum(x)})
rownames(df_pre_CCF_norm) <- df_pre_CCF$CCF
df_pre_CCF_norm <- as.data.frame(df_pre_CCF_norm)
df_pre_CCF_norm$CCF <- rownames(df_pre_CCF_norm)
df_pre_CCF_norm <- pivot_longer(df_pre_CCF_norm, !CCF, names_to = "Target", values_to = "Value")

df_true_CCF <-
  df_true[,c(sp_Barcode,"CCF")] |>
  group_by(CCF) |>
  dplyr::summarize(across(1:15, ~ sum(.x))) |>
  as.data.frame()
df_true_CCF_norm <- apply(df_true_CCF[,sp_Barcode], 2, function(x){x/sum(x)})
rownames(df_true_CCF_norm) <- df_true_CCF$CCF
df_true_CCF_norm <- as.data.frame(df_true_CCF_norm)
df_true_CCF_norm$CCF <- rownames(df_true_CCF_norm)
df_true_CCF_norm <- pivot_longer(df_true_CCF_norm, !CCF, names_to = "Target", values_to = "Value")

df_true <- df_true_CCF_norm
df_true$label <- ""
df_true$label[which(df_true$Value>0.25)] <- "*"
df_true$Target <- factor(df_true$Target, levels = sp_Barcode)
df_true$CCF <- factor(df_true$CCF, levels = rev(c("MOs","ACAd","PL","ILA","DP")))

breaks <- seq(0,0.5,0.01)
p_true <-
  ggplot(df_true, aes(x=Target, y=CCF, fill=Value)) +
  geom_raster() +
  geom_text(aes(label=label),col ="black",size = 5) +
  scale_fill_gradientn(limits=c(0,0.5), colours = colorRampPalette(c("navy","white","firebrick3"))(100), na.value="firebrick3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 1)) +
  labs(x="", y="")


df_pre <- df_pre_CCF_norm
df_pre$label <- ""
df_pre$label[which(df_pre$Value>0.25)] <- "*"
df_pre$Target <- factor(df_pre$Target, levels = sp_Barcode)
df_pre$CCF <- factor(df_pre$CCF, levels = rev(c("MOs","ACAd","PL","ILA","DP")))

breaks <- seq(0,0.5,0.01)
p_pre <-
  ggplot(df_pre, aes(x=Target, y=CCF, fill=Value)) +
  geom_raster() +
  geom_text(aes(label=label),col ="black",size = 5) +
  scale_fill_gradientn(limits=c(0,0.5), colours = colorRampPalette(c("navy","white","firebrick3"))(100), na.value="firebrick3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 1)) +
  labs(x="", y="")

p5 <- plot_grid(p_true, p_pre, nrow = 2)

ggsave("../pdf/Figure7/p5.pdf", plot = p5, height = 5, width = 6, units = "in")



df <- df_true_CCF_norm
colnames(df) <- c("CCF", "Target", "True")
df$Predict <- df_pre_CCF_norm$Value
cor <- cor.test(df$True, df$Predict, "two.sided", "pearson")
R <- round(cor$estimate,2)
P <- format(cor$p.value, digits = 2)
p6 <-
  ggplot(df, aes(x=True, y=Predict)) +
  geom_pointdensity(adjust = 1, size=2) +
  geom_smooth(method = "lm", color='black', linewidth=0.5, se=F) +
  scale_color_distiller(palette = "Spectral", direction = -1) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 15),
        legend.position = "none") +
  labs(x="True", y="Predict",
       title = paste('R =',R,', P =',P,sep=' ')) +
  ylim(0,1) +
  xlim(0,1)

ggsave("../pdf/Figure7/p6.pdf", plot = p6, height = 4, width = 4, units = "in")



proj_module <- c("PTi","ITi-M1","ITi-M2","ITC-M3")
pre_module <- spatial_test
for (i in 1:length(proj_module)){
  proj_i <- proj_module[i]
  y_train <- rep(0,nrow(sc_seu@meta.data))
  y_train[which(sc_seu$Proj_module == proj_i)] <- 1

  # transcriptom + spatial
  X_train <- as.matrix(cbind(transcriptom_train, spatial_train))
  X_train <- Matrix(X_train, sparse = T)
  X_test <- as.matrix(cbind(transcriptom_test, spatial_test))
  X_test <- Matrix(X_test, sparse = T)
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  cv <- xgb.cv(data = dtrain, nrounds = 1000, nfold = 5, max_depth = 5, eta = 0.5,
               early_stopping_rounds = 5, objective = "binary:logistic",
               verbose = F)
  model_xgb <- xgboost(data=dtrain, max_depth=5, eta=0.5, nthread = 5,
                       nround = cv$best_iteration, objective = "binary:logistic",
                       verbose = F)
  pre <- predict(model_xgb, newdata = X_test)
  pre_module[,proj_i] <- pre
}

module <- proj_module[1]
pre_module$Module <- apply(pre_module[,proj_module], 1, function(x){
  proj_module[which.max(x)]
})
p7_1 <-
  ggplot() +
  geom_point(pre_Barcode[which(pre_module$Module != module),], mapping=aes(x=x,y=y),
             size=0.8, color="#00204DFF") +
  geom_point(pre_Barcode[which(pre_module$Module == module),], mapping=aes(x=x,y=y),
             size=0.8, color="#FFEA46FF") +
  ggdark::dark_theme_void() +
  labs(title = "") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_fixed() +
  xlim(0,1.5) +
  ylim(-4.3,-1.25)

ggsave(paste("../pdf/Figure7/p7_1_",module,".png",sep=""), plot = p7_1, height = 4, width = 2, units = "in")



module <- proj_module[1]
slice <- 'IT_slice_18'

df <- sp.PFC.Left@meta.data[,c("ML_new","DV_new","SubType_Layer","Proj_module")]
df$Proj_module[is.na(df$Proj_module)] <- ""
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

p7_2 <-
  ggplot(df, aes()) +
  geom_point(df[which(df$Proj_module != module),], mapping=aes(x=x,y=y),
             size=0.8, color="#00204DFF") +
  geom_point(df[which(df$Proj_module == module),], mapping=aes(x=x,y=y),
             size=0.8, color="#FFEA46FF") +
  ggdark::dark_theme_void() +
  labs(title = "") +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_fixed() +
  xlim(0,1.5) +
  ylim(-4.3,-1.25)

ggsave(paste("../pdf/Figure7/p7_2_",module,".png",sep=""), plot = p7_2, height = 4, width = 2, units = "in")



i=15
module <- sp_Barcode[i]
p8_1 <-
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

ggsave(paste("../pdf/Figure7/p8_1_",module,".png",sep=""), plot = p8_1, height = 4, width = 2, units = "in")



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

p8_2 <-
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

ggsave(paste("../pdf/Figure7/p8_2_",module,".png",sep=""), plot = p8_2, height = 4, width = 2, units = "in")




