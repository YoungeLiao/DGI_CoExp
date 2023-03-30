# === load data v1: gene interaction ===
# Hub_blue_raw <- read.csv('./03_DGI/output_data/GeneListofFilteredPathway_BvsD_HUBCluster4.csv', header = TRUE)
# Hub_blue <- Hub_blue_raw[, c(1, 13:21, 9:10, (length(Hub_blue_raw) -1) )]
# Hub_yellow_raw <- read.csv('./03_DGI/output_data/GeneListofFilteredPathway_YvsD_HUBCluster3.csv', header = TRUE)
# Hub_yellow <- Hub_yellow_raw[, c(1, 13:21, 9:10, (length(Hub_yellow_raw) -1))]
# phototransduction_blue_raw <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_BvsD_7Cluster.csv', header = T)
# phototransduction_blue<- phototransduction_blue_raw[, c(1, 13:21, 9:10, (length(phototransduction_blue_raw) -1))]
# phototransduction_yellow_raw <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_YvsD_7Cluster.csv', header = T)
# phototransduction_yellow <- phototransduction_yellow_raw[, c(1, 13:21, 9:10, (length(phototransduction_yellow_raw) -1) )]
# merge_data <- unique(rbind(Hub_blue, Hub_yellow, phototransduction_blue,phototransduction_yellow))
# === load data v2: gene interaction - split datasets, blue vs dark ===
# hub
Hub_blue_raw <- read.csv('./03_DGI/output_data/GeneListofFilteredPathway_BvsD_HUBCluster4_v2.csv', header = TRUE)
Hub_blue <- Hub_blue_raw[, c(1, 13:21, 9:10, (length(Hub_blue_raw) -1) )]
merge_data <- unique(Hub_blue)
# signaling
signal_blue_raw <- read.csv('./03_DGI/output_data/GeneListofFilteredPathway_BvsD_SignalCluster5_v2.csv', header = TRUE)
signal_blue <- signal_blue_raw[, c(1, 13:21, 9:10, (length(signal_blue_raw) -1) )]
phototransduction_blue_raw <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_BvsD_7Cluster.csv', header = T)
phototransduction_blue<- phototransduction_blue_raw[, c(1, 13:21, 9:10, (length(phototransduction_blue_raw) -1))]
merge_data <- unique(rbind(signal_blue, phototransduction_blue))

# === load data v2: gene interaction - split datasets, yellow vs dark ===
# hub
Hub_yellow_raw <- read.csv('./03_DGI/output_data/GeneListofFilteredPathway_YvsD_Cluster3.csv', header = TRUE)
Hub_yellow <- Hub_yellow_raw[, c(1, 13:21, 9:10, (length(Hub_yellow_raw) -1))]
merge_data <- unique(Hub_yellow)
# signaling
signal_yellow_raw <- read.csv('./03_DGI/output_data/GeneListofFilteredPathway_YvsD_SignalCluster1.csv', header = TRUE)
signal_yellow <- signal_yellow_raw[, c(1, 13:21, 9:10, (length(signal_yellow_raw) -1) )]
phototransduction_yellow_raw <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_YvsD_7Cluster.csv', header = T)
phototransduction_yellow <- phototransduction_yellow_raw[, c(1, 13:21, 9:10, (length(phototransduction_yellow_raw) -1) )]
merge_data <- unique(rbind(signal_yellow, phototransduction_yellow))



head(merge_data)
# gene_id     Blue1      Blue2      Blue3      Dark1      Dark2      Dark3    Yellow1
# 1 g_24930  2.412453   3.171607   3.431420   40.12658   36.59891   63.53419   49.18134
# 2 g_05066 53.302695  55.998884  38.928972 1953.27696 2101.64390 1366.70594 1096.40196
# 3 g_41985 88.519220 163.638073 225.939181 1852.31369 2007.43132 1745.86117 2833.02867
# 4 g_18966 28.812081  30.191760  33.800899  895.10850  940.92516  706.55961  845.44092
# 6 g_26961 10.731651  12.358376   9.894007  372.49442  388.76610  413.88250  261.39488
# 8 g_22023  5.337555   7.832400   7.490678  289.52171   78.02930  393.80694 1010.76070
# Yellow2    Yellow3                                                ko_des
# 1   24.29615   30.19671                             glutathione S-transferase
# 2 1958.52625 1983.24405                                             flagellin
# 3 1488.56505  440.34928                    superoxide dismutase, Fe-Mn family
# 4  801.46713  886.07421                                  glutamine synthetase
# 6  221.98411  285.76343 ubiquinol-cytochrome c reductase cytochrome b subunit
# 8  971.77936  574.33384                              molecular chaperone HtpG
gene_names <- data.frame(paste(merge_data$ko_des, ' (', merge_data$gene_id, ')', sep = ''))
colnames(gene_names) <- 'label'
corr_data <- cbind(gene_names, merge_data[,2:(length(merge_data)-1)])
rownames(corr_data) <- corr_data$label
data <- t(corr_data[, -1])
# calculate correlation matrix
library(psych)
occor = corr.test(data,use="pairwise",method="spearman",adjust="fdr",alpha=0.05)
occor.r = occor$r # 取相关性矩阵R值
occor.p = occor$p # 取相关性矩阵p值
# filtering
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0
occor.r[occor.p>0.05|abs(occor.r)<0.8] = 0
occor.r[occor.p>0.05|abs(occor.r)<0.9] = 0
# save data
write.csv(occor.r, './03_DGI/output_data/20230326_Corr_blue_cluster4_hub_r0.8.csv')
write.csv(occor.r, './03_DGI/output_data/20230326_Corr_blue_cluster4_hub_r0.9.csv')
write.csv(occor.r, './03_DGI/output_data/20230326_Corr_blue_cluster5_signal_r0.9.csv')

write.csv(occor.r, './03_DGI/output_data/20230326_Corr_yellow_cluster3_hub_r0.9.csv')
write.csv(occor.r, './03_DGI/output_data/20230326_Corr_yellow_cluster1_signal_r0.9.csv')
# write.csv(occor.r,file="FH_CK_0.05_occor.csv")

# === load data v2: clusters' pathways interaction ===
cluster4_blue <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster4.csv', header = T)
cluster5_blue <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster5.csv', header = T)
cluster7_blue <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster7.csv', header = T)
cluster1_yellow <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster1.csv', header = T)
cluster3_yellow <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster3.csv', header = T)
cluster4_yellow <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster4.csv', header = T)
pathway_names <- data.frame(unique(c(cluster4_blue[,1], cluster5_blue[,1], cluster7_blue[,1], cluster1_yellow[,1], cluster3_yellow[,1], cluster4_yellow[,1])))
colnames(pathway_names) <- 'level3_pathway_name'
# merge_data <- unique(rbind(cluster4_blue, cluster5_blue, cluster7_blue, cluster1_yellow, cluster3_yellow, cluster4_yellow))
# === data treatment ===
gene_kegg_yellow <- read.csv('./03_DGI/output_data/gene_pathway_expre_YvsD.csv', header = T)
gene_kegg_blue <- read.csv('./03_DGI/output_data/gene_pathway_expre_BvsD.csv', header = T)
DEGs_kegg <- data.frame(unique(rbind(gene_kegg_yellow,gene_kegg_blue)))

# === expression of all pathways & merge ===
# library(reshape)
# kegg_melt <- melt(DEGs_kegg[,2:11])
library(dplyr)
# library(tidyverse)
# detach('package:plyr')
# group_by(rawdata, level3_pathway_name)
rawdata <- DEGs_kegg
pathway <- rawdata %>% 
  group_by(level3_pathway_name) %>%
  summarise(B1mean = mean(Blue1), 
            B2mean = mean(Blue2),
            B3mean = mean(Blue3),
            D1mean = mean(Dark1), 
            D2mean = mean(Dark2),
            D3mean = mean(Dark3),
            Y1mean = mean(Yellow1), 
            Y2mean = mean(Yellow2),
            Y3mean = mean(Yellow3)) 

merge_data<- merge(pathway_names, pathway)
# 
head(merge_data)
# level3_pathway_name   B1mean   B2mean   B3mean    D1mean     D2mean    D3mean    Y1mean
# 1 Alanine, aspartate and glutamate metabolism 31.15849 38.90909 35.56914 17.944162 15.4665762 18.990226 18.164369
# 2             alpha-Linolenic acid metabolism 20.64777 20.70895 21.41108  1.765062  0.8976365  2.129684  2.282243
# 3                           Alzheimer disease 38.76380 23.75352 19.95845  2.047364  2.2923978  2.617527  3.699762
# 4 Amino sugar and nucleotide sugar metabolism 12.87178 12.62211 12.17253  7.888147  6.9698217  8.521604  7.000736
# 5               Amyotrophic lateral sclerosis 27.65324 13.74235 13.74070  1.644664  1.7468111  1.989740  2.997966
# 6                             Apoptosis - fly 44.74303 19.82746 18.18667  1.797269  2.0439444  2.377452  3.679412
# Y2mean    Y3mean
# 1 15.931349 19.108718
# 2  2.675228  2.403213
# 3  3.326148  1.973785
# 4  5.968388  7.405986
# 5  2.921154  1.578480
# 6  3.790025  2.285645
# gene_names <- data.frame(paste(merge_data$ko_des, ' (', merge_data$gene_id, ')', sep = ''))
# colnames(gene_names) <- 'label'
# corr_data <- cbind(gene_names, merge_data[,2:(length(merge_data)-1)])
corr_data <- merge_data[,2:10]
rownames(corr_data) <- merge_data$level3_pathway_name
# rownames(corr_data) <- corr_data$label
data <- t(corr_data[, -1])
# calculate correlation matrix
library(psych)
occor = corr.test(data,use="pairwise",method="spearman",adjust="fdr",alpha=0.05)
occor.r = occor$r # 取相关性矩阵R值
occor.p = occor$p # 取相关性矩阵p值
# filtering
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0
occor.r
# save data
write.csv(occor.r, './03_DGI/output_data/Corr_pathways.csv')
# write.csv(occor.r,file="FH_CK_0.05_occor.csv")

# === export data analysis ===
# --- v2: pathways ---
nodes <- read.csv('./03_DGI/data/Nodes_raw_pathways.csv', header = T)[,c(1, 4, 5)]
colnames(nodes) <- c('level3_pathway_name', 'degree', 'weighted degree')
nodes_cluster <- merge(nodes, merge_data, by = 'level3_pathway_name')

# --- v1: genes ---
# load data 
nodes <- read.csv('./03_DGI/data/Nodes_raw_pretreated.csv', header = T)
nodes_cluster <- merge(nodes, merge_data, by = 'gene_id')
nodes_subcellular <- nodes_cluster[,c(2:6, 16:17)]
# save data
write.csv(nodes_subcellular,'./03_DGI/output_data/Nodes_subcellular.csv', row.names = FALSE)


