## === BvsD DEGs ===
# https://mp.weixin.qq.com/s/Yj4EZTAJ7kNNUomIw_Kikg
setwd("/Users/yangliao/Documents/GitHub/OSProtein/")#设置工作目录
# --- load data ---
cluster_label <- read.csv('./03_DGI/output_data/tsne_BvsD_labels_all_7cluster.csv', header = TRUE)
fpkm <- read.table('./data/gene_fpkm_all_samples.txt', header = TRUE)
annotation <- read.csv('./data/Annotation_Expre.csv', header = TRUE)
kegg_anno <- read.csv('./data/kegg_annotation_all.csv')
seq <- read.csv('./data/sequences.csv', header = TRUE)

# --- match/merge data ---
matched <- merge(cluster_label, annotation, by.x = 'gene_id', all.x = TRUE, all.y = FALSE)
matched_kegg <- merge(matched, kegg_anno, by = 'gene_id')
table(matched_kegg$DGI_sub_label)
table(matched_kegg$Kmeans_label)
## save data
gene_pathway_expre <- matched_kegg[,c(1,13:21,29)]
write.csv(gene_pathway_expre, './03_DGI/output_data/gene_pathway_expre_BvsD.csv', row.names = FALSE)

# --- validation via KEGG ---
nitrogen <- matched_kegg[grepl("Nitrogen metabolism", matched_kegg$level3_pathway_name), ] # cluster 3
table(nitrogen$DGI_sub_label) # cluster 3
phototransduction <- matched_kegg[grepl("Phototransduction", matched_kegg$level3_pathway_name), ] # cluster 3, 5
table(phototransduction$DGI_sub_label)
# Cluster 3 Cluster 5 
# 5        13 
pathway <- data.frame(table(matched_kegg$level3_pathway_name))

# grep functions
photosynthesis <- matched_kegg[grepl("Photosynthesis", matched_kegg$level3_pathway_name), ] 
peroxisome <- matched_kegg[grepl("Peroxisome", matched_kegg$level3_pathway_name), ] 
OxidativePhosphorylation <- matched_kegg[grepl('Oxidative phosphorylation', matched_kegg$level3_pathway_name), ]  
Xebio <- matched_kegg[grepl('Metabolism of xenobiotics by cytochrome P450', matched_kegg$level3_pathway_name), ]  
# OxidativePhosphorylation possess more than 300 genes that distributed to multi-pathways, thus is not suitable for validation
Nitrogen <- matched_kegg[grepl('Nitrogen metabolism', matched_kegg$level3_pathway_name), ]  
ROS <- matched_kegg[grepl('Longevity regulating pathway - multiple species', matched_kegg$level3_pathway_name), ]  
P450 <- matched_kegg[grepl('Drug metabolism - cytochrome P450', matched_kegg$level3_pathway_name), ]  
GABA <- matched_kegg[grepl('GABAergic synapse', matched_kegg$level3_pathway_name), ]  
Exosome <- matched_kegg[grepl('Exosome', matched_kegg$level3_pathway_name), ]  
glutathione <- matched_kegg[grepl('Glutathione metabolism', matched_kegg$level3_pathway_name), ]  
glutamate <- matched_kegg[grepl('Alanine, aspartate and glutamate metabolism', matched_kegg$level3_pathway_name), ]  
Exopoly <- matched_kegg[grepl('Exopolysaccharide biosynthesis', matched_kegg$level3_pathway_name), ]  
linolenic <- matched_kegg[grepl('alpha-Linolenic acid metabolism', matched_kegg$level3_pathway_name), ]  
Galactose <- matched_kegg[grepl('Galactose metabolism', matched_kegg$level3_pathway_name), ]  
CRP <- matched_kegg[grepl('Cytochrome P450', matched_kegg$level3_pathway_name), ]  
ProDiges <- matched_kegg[grepl('Protein digestion and absorption', matched_kegg$level3_pathway_name), ]  
Insulin <- matched_kegg[grepl('Insulin secretion', matched_kegg$level3_pathway_name), ]  
Carbohydrate <- matched_kegg[grepl('Carbohydrate digestion and absorption', matched_kegg$level3_pathway_name), ]  
CO2 <- matched_kegg[grepl('Carbon fixation pathways in prokaryotes', matched_kegg$level3_pathway_name), ]  
Sulfer <- matched_kegg[grepl('Sulfur metabolism', matched_kegg$level3_pathway_name), ]  
Exopolysaccharide <- matched_kegg[grepl('Exopolysaccharide biosynthesis', matched_kegg$level3_pathway_name), ]  
Porphyrin <- matched_kegg[grepl('Porphyrin metabolism', matched_kegg$level3_pathway_name), ]  
Secretion <- matched_kegg[grepl('Secretion system', matched_kegg$level3_pathway_name), ]  
Energy <- matched_kegg[grepl('Energy metabolism', matched_kegg$level3_pathway_name), ]

neg_test <- matched_kegg[grepl('Gastric acid secretion', matched_kegg$level3_pathway_name), ]  
neg_test <- matched_kegg[grepl('Ribosome', matched_kegg$level3_pathway_name), ]  
neg_test <- matched_kegg[grepl('Huntington disease', matched_kegg$level3_pathway_name), ]  
neg_test <- matched_kegg[grepl('Chronic myeloid leukemia', matched_kegg$level3_pathway_name), ]  
neg_test <- matched_kegg[grepl('Platelet activation', matched_kegg$level3_pathway_name), ] 
neg_test <- matched_kegg[grepl('Glioma', matched_kegg$level3_pathway_name), ]  
neg_test <- matched_kegg[grepl('Leukocyte transendothelial migration', matched_kegg$level3_pathway_name), ]

neg_test1 <- matched_kegg[grepl('Transporters', matched_kegg$level3_pathway_name), ]  
neg_test2 <- matched_kegg[grepl('Amyotrophic lateral sclerosis', matched_kegg$level3_pathway_name), ]  
neg_test3 <- matched_kegg[grepl('Chromosome and associated proteins', matched_kegg$level3_pathway_name), ]  
neg_test4 <- matched_kegg[grepl('Progesterone-mediated oocyte maturation', matched_kegg$level3_pathway_name), ]  
neg_test5 <- matched_kegg[grepl('Mitophagy - animal', matched_kegg$level3_pathway_name), ]  
result_nitrate <- matched_kegg[grepl("nitrate", matched_kegg$SwissProt_Description), ] # cluster 3, the smallest clusters
result_nitrite <- matched_kegg[grepl("nitrite", matched_kegg$SwissProt_Description), ] # cluster 3, the smallest clusters
result_N <- rbind(result_nitrate, result_nitrite)

fun <- peroxisome

# reorder results
## Initialization - DGI_label 
# (-1.429994 for peroxisome) (2.272994 for phototransduction) (-1.19469) (-0.6340549)
data <-  data.frame(table(fun$DGI_label))
ref <- data.frame(table(matched_kegg$DGI_label))
## Initialization - DGI_sub_label 
# (-1.2299 for peroxisome) (2.31007 for phototransduction) (-1.289782 for OxidativePhosphorylation) (-0.5954653 for Xebio)
data <- data.frame(table(fun$DGI_sub_label))
ref <- data.frame(table(matched_kegg$DGI_sub_label))
## Initialization - Kmeans_label 
# (-1.655619 for peroxisome) (1.70732) (-1.072017) (-1.540794)
data <-  data.frame(table(fun$Kmeans_label))
ref <- data.frame(table(matched_kegg$Kmeans_label))
## Initialization - Kmeans_sub_label (-1.324743 for peroxisome) (1.521689) (0.3578049) (-1.098327)
data <-  data.frame(table(fun$Kmeans_subcel_label))
ref <- data.frame(table(matched_kegg$Kmeans_subcel_label))

# scoring 
data <-  data.frame(table(fun$DGI_label))
ref <- data.frame(table(matched_kegg$DGI_label))
data_ref <- merge(data, ref, by = 'Var1')
result <- data.frame(data_ref$Freq.x/data_ref$Freq.y) # ri
colnames(result) <- 'Freq'
result_ordered <- result[order(result$Freq, decreasing = TRUE),]
# result_ordered
# 0.002382087 0.001768242 0.001722950 0.001277427 0.001096491 0.001022495 0.000266809
if(length(result_ordered) > 2){
  score_raw <- sum(result_ordered[1:2]) - sum(result_ordered[3:length(result_ordered)])
} else{
  score_raw <- sum(result_ordered)
  }
# score_raw <- sum(result_ordered[1:2]) - sum(result_ordered[3:7])
score <- score_raw/(sum(data_ref$Freq.x)/sum(data_ref$Freq.y)) # larger is better (i.e. when value is negative, closer to zero indicates better performance)
score

# (sum(data_ref$Freq.x)/sum(data_ref$Freq.y)) # 0.001379641
# === save functions ===
fun_save <- rbind(phototransduction, peroxisome, Secretion, result_N)[,c(1, 9:21, 29:32)]
write.csv(fun_save, './03_DGI/output_data/Mechanism_genes_BvsD_20230318.csv', row.names = FALSE)

fun_save <- rbind(phototransduction, peroxisome, P450, Xebio, ROS)[,c(1, 9:21, 29:32)]
write.csv(fun_save, './03_DGI/output_data/Mechanism_genes_BvsD.csv', row.names = FALSE)

fun_save_v2 <- rbind(phototransduction, peroxisome, P450, Xebio, ROS, Secretion, Porphyrin, Exopolysaccharide, Sulfer)[,c(1, 9:21, 29:32)]
write.csv(fun_save_v2, './03_DGI/output_data/Mechanism_genes_Fig3_BvsD_v2.csv', row.names = FALSE)

fun_save_v3 <- rbind(Nitrogen, phototransduction, peroxisome, P450, Xebio, ROS, Secretion, Porphyrin, Exopolysaccharide, Sulfer)[,c(1, 9:21, 29:32)]
write.csv(fun_save_v3, './03_DGI/output_data/Mechanism_genes_Fig3Nitrogen_BvsD_v3.csv', row.names = FALSE)

fun_save_v4 <- rbind(result_N, Nitrogen, phototransduction, peroxisome, P450, Xebio, ROS, Secretion, Porphyrin, Exopolysaccharide, Sulfer)[,c(1, 9:21, 29:32)]
write.csv(fun_save_v4, './03_DGI/output_data/Mechanism_genes_NO3NO2_BvsD_v5.csv', row.names = FALSE)

# === ratio of key/functional genes ===
name_list <- list(Xebio, peroxisome, Nitrogen, phototransduction)
names <- c('Xebio', 'peroxisome', 'Nitrogen','phototransduction')
# fun <- name_list[[4]]
# data <-  data.frame(table(fun$DGI_sub_label))
# ref <- data.frame(table(matched_kegg$DGI_sub_label))
# data_ref <- merge(data, ref, by = 'Var1')
# result <- data.frame(data_ref$Freq.x/data_ref$Freq.y)
# ratio <- cbind(data_ref, result)
# colnames(ratio) <- c('Var1', 'Freq.x', 'Freq.y', 'Ratio')
# ratio_ordered <- ratio[order(ratio$Ratio, decreasing = TRUE),]
# ratio_ordered$Pathway <- names[1]
# # for i = 1
# ratio_results <- ratio_ordered
# # for i > 1
# ratio_results <- rbind(ratio_results, ratio_ordered)
for (i in 1:length(name_list)){
  fun <- name_list[[i]]
  data <-  data.frame(table(fun$DGI_sub_label))
  ref <- data.frame(table(matched_kegg$DGI_sub_label))
  data_ref <- merge(data, ref, by = 'Var1')
  result <- data.frame(data_ref$Freq.x/data_ref$Freq.y)
  ratio <- cbind(data_ref, result)
  colnames(ratio) <- c('Var1', 'Freq.x', 'Freq.y', 'Ratio')
  ratio_ordered <- ratio[order(ratio$Var1, decreasing = TRUE),]
  ratio_ordered$Pathway <- names[i]
  if (i == 1){
    ratio_results <- ratio_ordered
  }else{
    ratio_results <- rbind(ratio_results, ratio_ordered)
  }
}
# visualization
# Stacked barplot with multiple groups
pale_9_1 <- as.vector(paletteer_d('ggprism::pastels'))
pale_20_2 <- as.vector(paletteer_d('ggthemes::Tableau_20'))
pale_11_2 <- as.vector(paletteer_d('khroma::sunset'))
pale_12_2 <- as.vector(paletteer_d('RColorBrewer::Set3'))
pale_8_1 <- as.vector(paletteer_d('RColorBrewer::Pastel2'))
pale_29 <- c(pale_20_2, pale_9_1)
pale_28 <- c(pale_8_1, pale_11_2, pale_9_1)
pale_31 <- c(pale_8_1, pale_11_2, pale_12_2)
pale_32 <- c(pale_12_2, pale_20_2)
pale_25 <- c(pale_12_2, pale_11_2, '#85A4FE', 'grey')
pale_29 <- c(pale_25, pale_8_1)

top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme1 <- theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.7, angle = 15), # vjust = -0.001
        legend.text = element_text(size = 12), 
        # legend.title = element_text(size = 16),
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.title = element_blank(), 
        legend.title = element_blank(),
        legend.position = 'right', # 调整legend位置
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

# Calculate the cumulative sum of len for each pathway
library(plyr)
ratio_results$Ratio <- sprintf('%0.4f',ratio_results$Ratio*100)
ratio_results$Ratio <- as.numeric(ratio_results$Ratio)
plot_data <- ddply(ratio_results, "Pathway",
                   transform, label_ypos=cumsum(Ratio))
head(plot_data)

ggplot(data=plot_data, aes(x=Pathway, y=Ratio, fill = factor(Var1, 
                                                             levels=c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5',  'Cluster 6', 'Cluster 7')))) +
  geom_bar(stat="identity") + 
  geom_text(aes(y=label_ypos, label=Ratio), vjust=1.2, 
            color="black", size=3) +
  scale_fill_manual(values = pale_25) + 
  labs(y = 'Ratio (%)',x = 'Pathways', title = 'Blue light') +
  mytheme1


# === key/functional genes on SwissProt and save data ===
# --- co-expressed functional pathways ---
table(Xebio$DGI_sub_label)
table(peroxisome$DGI_sub_label)

write.csv(phototransduction, './03_DGI/output_data/ClusterAssign_phototransduction_BvsD_7Cluster.csv', row.names = FALSE) # cluster 5
write.csv(Xebio, './03_DGI/output_data/ClusterAssign_XebioP450_BvsD_7Cluster.csv', row.names = FALSE) # cluster 5
write.csv(peroxisome, './03_DGI/output_data/ClusterAssign_peroxisome_BvsD_7Cluster.csv', row.names = FALSE) # cluster 5

# --- nitrogen metabolism: detailed ---
# - nitrate and nitrite ---
result_nitrate <- matched[grepl("nitrate", matched$SwissProt_Description), ] # cluster 4, 7， based on narK12 and nirK
table(result_nitrate$DGI_sub_label) 
# Cluster 1 Cluster 2 Cluster 3 Cluster 4 Cluster 5 
# 2         1         8         2         1 
result_nitrite <- matched[grepl("nitrite", matched$SwissProt_Description), ] 
table(result_nitrite$DGI_sub_label) 
# Cluster 1 Cluster 2 Cluster 3 Cluster 4 Cluster 5 Cluster 6 Cluster 7 
# 3         1         7         2         3         1         5 
result_N <- unique(rbind(result_nitrate, result_nitrite))

result_N$Mean_expre <- apply(result_N[, 13:21], 1, mean)
write.csv(result_N, './03_DGI/output_data/ClusterAssign_N_BvsD_7Cluster.csv', row.names = FALSE)

# filtering for cluster identification
thre <- 3
result_N_filtered <- unique(filter(result_N, result_N$Mean_expre > thre)) # cluster 1, 5
# - nitrogen metabolism -
result_NitrogenMetabolism <- matched_kegg[grepl("Nitrogen metabolism", matched_kegg$level3_pathway_name), ] # cluster 4, 7， based on narK12 and nirK
result_NitrogenMetabolism$Mean_expre <- apply(result_NitrogenMetabolism[, 13:21], 1, mean)
# filtering for cluster identification
thre <- 3
result_NitrogenMetabolism_filtered <- unique(filter(result_NitrogenMetabolism, result_NitrogenMetabolism$Mean_expre > thre)) # cluster 1, 5
write.csv(result_NitrogenMetabolism_filtered, './03_DGI/output_data/ClusterAssign_keggNitrogenMetabo_BvsD_7Cluster.csv', row.names = FALSE)






# --- other function test ---
result_photo <- matched_kegg[grepl("photo", matched_kegg$SwissProt_Description), ] 
# result_test <- matched[grepl('Oxidative', matched$SwissProt_Description), ] 
result_dismutase <- matched[grepl("dismutase", matched$SwissProt_Description), ] # 12
# contribute to higher O2- via inhibiting the quenching of superoxide
result_oxidative <- matched[grepl("Oxidative", matched$SwissProt_Description), ]


# --- key subcellular information ---
library(dplyr)
secreted <- data.frame(filter(matched, matched$Sig0alP == 1 & matched$tmhmm == 0)) # based on the above results
secreted_kegg <- data.frame(filter(matched_kegg, matched_kegg$Sig0alP == 1 & matched_kegg$tmhmm == 0)) # based on the above results
table(cluster_label$DGI_sub_label)
# Cluster 1 Cluster 2 Cluster 3 Cluster 4 Cluster 5 Cluster 6 Cluster 7 
# 4828      1795      9533      1273      4198      2061      1589 
table(secreted$DGI_sub_label)
# Cluster 1 Cluster 2 Cluster 3 Cluster 4 Cluster 5 Cluster 6 Cluster 7 
# 219        81       403        26       204       114       198 
table(secreted_kegg$DGI_sub_label)
# Cluster 1 Cluster 2 Cluster 3 Cluster 4 Cluster 5 Cluster 6 Cluster 7 
# 148        99       289        17       210        79       110 
# --- high abundance ---
# cluster18, cluster17, etc. —— 感觉cluster分布与表达量的correlation非常大 

# === extract function cluster/genes ===
# ## V1 --- N (cluster 3): sequence 
# Cluster3 <- filter(cluster_label, cluster_label$DGI_sub_label == 'Cluster 3') 
# Cluster3_seq <- merge(Cluster3, seq, by = 'gene_id')
# Cluster3_seq <- Cluster3_seq[, -c(2:6)]
## V2 --- cluster 3, 5 (photo)
Cluster3 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 3') 
Cluster5 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 5') 
### save data
length(unique(Cluster3$gene_id)) # 4180 genes in total 
length(unique(Cluster5$gene_id)) # 2105 genes in total 
write.csv(Cluster3, './03_DGI/output_data/Cluster3_BvsD.csv', row.names = FALSE)
write.csv(Cluster5, './03_DGI/output_data/Cluster5_BvsD.csv', row.names = FALSE)


## narK12, nirK: cluster 4, 7
Cluster4 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 4') # 
Cluster7 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 7') # 
### save data
length(unique(Cluster4$gene_id)) # 394 genes in total 
length(unique(Cluster7$gene_id)) # 744 genes in total 
write.csv(Cluster4, './03_DGI/output_data/Cluster4_BvsD.csv', row.names = FALSE)
write.csv(Cluster7, './03_DGI/output_data/Cluster7_BvsD.csv', row.names = FALSE)

## subcellular - top30 un-annotated genes: cluster 2 
Cluster2 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 2')
length(unique(Cluster2$gene_id)) # 767 genes in total 
write.csv(Cluster2, './03_DGI/output_data/Cluster2_BvsD.csv', row.names = FALSE)


# cluster <- Cluster1
# write.csv(cluster, './03_DGI/output_data/YvsD/FunctionAnaly_cluster3.csv', header = TRUE)

# # Cluster2_seq <- seq
# Cluster2_seq$gene_id <- paste('>', Cluster2_seq$gene_id, sep = '')
# a <- c(Cluster2_seq$gene_id[1], Cluster2_seq$sequence[1])
# for (i in 2:length(Cluster2_seq$gene_id)){
#   a <- c(a, Cluster2_seq$gene_id[i], Cluster2_seq$sequence[i])
#   }
# a_df <- data.frame(a)

# save data
# write.table(a_df, './03_DGI/output_data/id_sequences_all.fa', col.names = FALSE, row.names = FALSE)

# --- kegg annotation ---
# cluster2_kegg <- merge(Cluster2, kegg_anno, by = 'gene_id')
# length(unique(cluster2_kegg$gene_id)) # 26
# cluster2_kegg <- merge(cluster2_kegg, fpkm, by = 'gene_id')
# pathway_table <- data.frame(table(cluster2_kegg$level3_pathway_name))
# colnames(pathway_table) <- c('level3', 'counts')
## test 
# grep <- cluster2_kegg[grepl(pathway_table$Var1[3], cluster2_kegg$level3_pathway_name), ]
# grep_mean <- apply(grep[,19:27], 1, mean)
# mean(grep_mean)

# initial abundance
# pathway_table$abundance <- 0
# pathway_table$Blue <- 0
# pathway_table$Dark <- 0
# pathway_table$Yellow <- 0
# for (i in 1:length(pathway_table$Var1)){
#   grep <- cluster2_kegg[grepl(pathway_table$Var1[i], cluster2_kegg$level3_pathway_name), ]
#   grep_mean <- sum(apply(grep[,19:27], 1, mean))
#   pathway_table$Blue[i] <- sum(apply(grep[,19:21], 1, mean))
#   pathway_table$Dark[i] <- sum(apply(grep[,22:24], 1, mean))
#   pathway_table$Yellow[i] <- sum(apply(grep[,25:27], 1, mean))
#   pathway_table$abundance[i] <- grep_mean
# }

# grep <- cluster2_kegg[grepl(pathway_table$Var1[1], cluster2_kegg$level3_pathway_name), ]
# expre_matrix <- data.frame(t(apply(grep[,19:27], 2, sum)))
# 
# for (i in 2:length(pathway_table$Var1)){
#   grep <- cluster2_kegg[grepl(pathway_table$Var1[i], cluster2_kegg$level3_pathway_name), ]
#   expre_matrix <- rbind(expre_matrix, data.frame(t(apply(grep[,19:27], 2, sum))))
# }
# 
# pathway_table <- cbind(pathway_table, expre_matrix)
