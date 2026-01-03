# --------------------
# title: FigureS3 Code
# author: Hu Zheng
# date: 2026-01-01
# --------------------

library(Seurat)
library(tidyverse)
library(cowplot)
library(data.table)
source('bin/Palettes.R')

all.inte <- readRDS('../data/rds/all.inte.rds')
all.Adult <- readRDS('../data/rds/all.Adult.rds')
Adult.Ex <- readRDS('../data/rds/Adult.Ex.rds')

Adult.Ex.barcode <- subset(
  Adult.Ex,
  cells=colnames(Adult.Ex)[which(Adult.Ex$BC_num>0)]
)



df <- fread('../data/csv/all_3mismatch_reads.txt', header = T)
p1 <-
  ggplot(df, aes(x=mismatch, fill=Type)) +
  geom_bar(color="black",width = 0.8) +
  scale_x_continuous(breaks = seq(0,20,5)) +
  scale_fill_manual(values = c("black","white")) +
  labs(x="number of mismatch", y="UMI counts") +
  theme_classic() +
  theme(text = element_text(size = 15))

ggsave("../pdf/FigureS3/p1.pdf", plot = p1, height = 4, width = 6, units = "in")



BC_order <- c("AId-C","ECT-C","ENTl-C","ENTl-I","DR-I","LHA-I","VTA-I","RE-I","SC-I","MD-I","AId-I","PL-C","VIS-I","BLA-I","BLA-C","ECT-I","SSp-I","ACB-I","CP-I","ACB-C","RSP-I","RSP-C","AUD-I","CP-C")
BC_mat <- all.Adult@meta.data[,BC_order]
BC_mat$cellid <- rownames(BC_mat)
df <- pivot_longer(BC_mat, !cellid, names_to = "Target", values_to = "Expression")
df <- df[which(df$Expression>0),]

p2 <-
  df |>
  mutate(class = fct_reorder(Target, Expression, .fun='median')) |>
  ggplot(aes(x=reorder(Target, Expression), y=Expression, fill=Target)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = col_Barcode) +
  labs(y="UMI counts", x= "Projection targets") +
  theme_classic() +
  theme(legend.position="none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

ggsave("../pdf/FigureS3/p2.pdf", plot = p2, height = 4, width = 7, units = "in")



Barcode <- c("CP-I","SSp-I","PL-C","CP-C","ACB-C","RSP-C","AUD-I","RSP-I","LHA-I",
             "DR-I","VTA-I","RE-I","SC-I","MD-I","ECT-I","BLA-I","VIS-I","BLA-C",
             "ACB-I","AId-I","ENTl-I","ENTl-C","ECT-C","AId-C")

BC_num <- list()
for (i in 1:length(Barcode)){
  BC_num[[i]] <- length(which(Adult.Ex.barcode@meta.data[,Barcode[i]] > 0))
}
df <- data.frame('Barcode' = Barcode, 'Num' = as.numeric(BC_num))
df$Barcode <- factor(df$Barcode,
                     levels = df$Barcode[order(df$Num, decreasing = T)])

p3 <-
  ggplot(df, aes(x=Barcode, y=Num)) +
  geom_col(width=0.7, fill="#342f86") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1, size = 12),
        text = element_text(size = 15)) +
  labs(x='Projection targets', y='Cell number')

ggsave("../pdf/FigureS3/p3.pdf", plot = p3, height = 4, width = 7, units = "in")



Barcode <- c('ACB-I','ACB-C','AId-I','AId-C','AUD-I','BLA-I','BLA-C','CP-I','CP-C',
             'ECT-I','ECT-C','ENTl-I','ENTl-C','PL-C','RSP-I','RSP-C','SSp-I','VIS-I',
             'DR-I','LHA-I','MD-I','RE-I','SC-I','VTA-I')

Adult1_thread <- c('ACB-C'=3,'ACB-I'=9,'BLA-I'=4,'CP-C'=8,'CP-I'=9,'DR-I'=6,'LHA-I'=5,
                   'MD-I'=3,'RE-I'=5,'SC-I'=4,'VIS-I'=7,'VTA-I'=5)
Adult2_thread <- c('ACB-C'=5,'ACB-I'=20,'AId-C'=9,'AId-I'=11,'AUD-I'=11,'BLA-I'=11,
                   'CP-C'=30,'CP-I'=30,'ECT-I'=8,'ENTl-I'=5,'PL-C'=19,'RSP-I'=15,
                   'SSp-I'=10,'VIS-I'=3)
Adult3_thread <- c('ACB-C'=7,'ACB-I'=25,'AId-C'=8,'AId-I'=15,'BLA-C'=3,'BLA-I'=11,
                   'CP-C'=18,'CP-I'=30,'ECT-C'=9,'ECT-I'=11,'ENTl-C'=5,'ENTl-I'=3,
                   'RSP-C'=5,'RSP-I'=12,'SSp-I'=14,'VIS-I'=3)
thread <- list('Adult1'=Adult1_thread, 'Adult2'=Adult2_thread, 'Adult3'=Adult3_thread)

BC_mat <- all.inte@meta.data[,c(Barcode,'sample','IfNeuron')]
BC_mat <- BC_mat |>
  pivot_longer(!(sample:IfNeuron), names_to = "Target", values_to = "Counts") |>
  filter(Counts > 0) |>
  dplyr::count(sample, IfNeuron, Target, Counts)
BC_mat$thread <- 0
for (i in 1:nrow(BC_mat)){
  BC_mat$thread[i] <- thread[[BC_mat$sample[i]]][[BC_mat$Target[i]]]
}

plist <- list()
sample <- c('Adult1','Adult2','Adult3')
for (n in 1:3){
  sample_n <- sample[n]
  thread_n <- thread[[n]]
  for (i in 1:length(thread_n)){
    df_i <- BC_mat[BC_mat$sample==sample_n & BC_mat$Target==names(thread_n)[i],]
    thread_i <- as.numeric(thread_n[i])
    plist[[sample_n]][[i]] <-
      ggplot(df_i, aes(x=Counts, y=n, color=IfNeuron)) +
      geom_line(linewidth=1) +
      geom_vline(xintercept= thread_i, colour="#9192ab", linetype="dashed", linewidth=1) +
      labs(title = names(thread_n)[i]) +
      scale_x_continuous(name='', limits=c(0,100), breaks=c(0,thread_i,50,100)) +
      scale_y_continuous(name = '', limits = c(0,50),
                         breaks = c(0,10,20,30,40,50)) +
      theme_bw() +
      theme(plot.title = element_text(size=20, hjust=0.5),
            legend.title=element_blank(), legend.position="none") +
      scale_color_manual(values = c("#f89588","#63b2ee"))
  }
}

p4_1 <- plot_grid(plotlist = plist[['Adult1']], ncol=7)

ggsave("../pdf/FigureS3/p4_1.pdf", plot = p4_1, height = 5, width = 16, units = "in")



p4_2 <- plot_grid(plotlist = plist[['Adult2']], ncol=7)

ggsave("../pdf/FigureS3/p4_2.pdf", plot = p4_2, height = 5, width = 16, units = "in")



p4_3 <- plot_grid(plotlist = plist[['Adult3']], ncol=7)

ggsave("../pdf/FigureS3/p4_3.pdf", plot = p4_3, height = 7.5, width = 16, units = "in")




