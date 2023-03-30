## === YvsD DEGs ===
# https://mp.weixin.qq.com/s/Yj4EZTAJ7kNNUomIw_Kikg
setwd("/Users/yangliao/Documents/GitHub/OSProtein")#设置工作目录
# --- load data ---
cluster_label <- read.csv('./03_DGI/output_data/tsne_YvsD_labels_all_7cluster.csv', header = TRUE)
annotation <- read.csv('./data/Annotation_Expre.csv', header = TRUE)
kegg_anno <- read.csv('./data/kegg_annotation_all.csv')
seq <- read.csv('./data/sequences.csv', header = TRUE)
fpkm <- read.table('./data/gene_fpkm_all_samples.txt', header = TRUE)

# --- match/merge data ---
matched <- merge(cluster_label, annotation, by.x = 'gene_id', all.x = TRUE, all.y = FALSE)
matched_kegg <- merge(matched, kegg_anno, by = 'gene_id')
## save data
gene_pathway_expre <- matched_kegg[,c(1,13:21,29)]
write.csv(gene_pathway_expre, './03_DGI/output_data/gene_pathway_expre_YvsD.csv', row.names = FALSE)

# --- validation via KEGG ---
nitrogen <- matched_kegg[grepl("Nitrogen metabolism", matched_kegg$level3_pathway_name), ] # cluster 3
table(nitrogen$DGI_sub_label) # cluster 3
phototransduction <- matched_kegg[grepl("Phototransduction", matched_kegg$level3_pathway_name), ] # cluster 3, 5
table(phototransduction$DGI_sub_label)
# Cluster 3 Cluster 5 
# 5        13 
## save data
# write.csv(phototransduction, './03_DGI/output_data/ClusterAssign_phototransduction_BvsD_24Cluster.csv', row.names = FALSE) # cluster 5
# photosynthesis <- matched_kegg[grepl("Photosynthesis", matched_kegg$level3_pathway_name), ] 
pathway <- data.frame(table(matched_kegg$level3_pathway_name))

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
neg_test <- matched_kegg[grepl('Gastric acid secretion', matched_kegg$level3_pathway_name), ]  
Phenylalanine <- matched_kegg[grepl('Phenylalanine, tyrosine and tryptophan biosynthesis', matched_kegg$level3_pathway_name), ]
ChemCarROS <- matched_kegg[grepl('Chemical carcinogenesis - reactive oxygen species', matched_kegg$level3_pathway_name), ]
Secretion <- matched_kegg[grepl('Secretion system', matched_kegg$level3_pathway_name), ]  
Porphyrin <- matched_kegg[grepl('Porphyrin metabolism', matched_kegg$level3_pathway_name), ] 
Sulfer <- matched_kegg[grepl('Sulfur metabolism', matched_kegg$level3_pathway_name), ]  
Exopolysaccharide <- matched_kegg[grepl('Exopolysaccharide biosynthesis', matched_kegg$level3_pathway_name), ]  
Energy <- matched_kegg[grepl('Energy metabolism', matched_kegg$level3_pathway_name), ] 
result_nitrate <- matched_kegg[grepl("nitrate", matched_kegg$SwissProt_Description), ] # cluster 3, the smallest clusters
result_nitrite <- matched_kegg[grepl("nitrite", matched_kegg$SwissProt_Description), ] # cluster 3, the smallest clusters
result_N <- rbind(result_nitrate, result_nitrite)

fun <- Energy

# pathway <- data.frame(table(matched_kegg$level3_pathway_name))
# reorder results
## Initialization - DGI_label 
# (-1.429994 for peroxisome) (2.272994 for phototransduction) (-1.19469) (-0.6340549)
data <-  data.frame(table(fun$DGI_label))
ref <- data.frame(table(matched_kegg$DGI_label))
## Initialization - DGI_sub_label 
# (-1.2299 for peroxisome) (2.31007 for phototransduction) (-1.289782 for OxidativePhosphorylation) (-0.5954653 for Xebio)
data <-  data.frame(table(fun$DGI_sub_label))
ref <- data.frame(table(matched_kegg$DGI_sub_label))
## Initialization - Kmeans_label 
# (-1.655619 for peroxisome) (1.70732) (-1.072017) (-1.540794)
data <-  data.frame(table(fun$Kmeans_label))
ref <- data.frame(table(matched_kegg$Kmeans_label))
## Initialization - Kmeans_sub_label (-1.324743 for peroxisome) (1.521689) (0.3578049) (-1.098327)
data <-  data.frame(table(fun$Kmeans_subcel_label))
ref <- data.frame(table(matched_kegg$Kmeans_subcel_label))

# # scoring (sum_n = 3)
# data_ref <- merge(data, ref, by = 'Var1')
# result <- data.frame(data_ref$Freq.x/data_ref$Freq.y)
# colnames(result) <- 'Freq'
# result_ordered <- result[order(result$Freq, decreasing = TRUE),]
# # result_ordered # 0.002506266 0.001345895 0.001231527
# if(length(result_ordered) > 3){
#   score_raw <- sum(result_ordered[1:3]) - sum(result_ordered[4:length(result_ordered)])
# } else{
#   score_raw <- sum(result_ordered)
# }
# # score_raw <- sum(result_ordered[1:2]) - sum(result_ordered[3:7])
# score <- score_raw/(sum(data_ref$Freq.x)/sum(data_ref$Freq.y)) # larger is better (i.e. when value is negative, closer to zero indicates better performance)
# score

# scoring (sum_n = 2)
data_ref <- merge(data, ref, by = 'Var1')
result <- data.frame(data_ref$Freq.x/data_ref$Freq.y)
colnames(result) <- 'Freq'
result_ordered <- result[order(result$Freq, decreasing = TRUE),]
# result_ordered # 0.002506266 0.001345895 0.001231527
if(length(result_ordered) > 2){
  score_raw <- sum(result_ordered[1:2]) - sum(result_ordered[3:length(result_ordered)])
} else{
  score_raw <- sum(result_ordered)
}
# score_raw <- sum(result_ordered[1:2]) - sum(result_ordered[3:7])
score <- score_raw/(sum(data_ref$Freq.x)/sum(data_ref$Freq.y)) # larger is better (i.e. when value is negative, closer to zero indicates better performance)
score

(sum(data_ref$Freq.x)/sum(data_ref$Freq.y)) # 0.00264109

# === save functions ===
fun_save <- rbind(phototransduction, peroxisome, Secretion, result_N)[,c(1, 9:21, 29:32)]
write.csv(fun_save, './03_DGI/output_data/Mechanism_genes_YvsD_20230318.csv', row.names = FALSE)

fun_save <- rbind(phototransduction, peroxisome, P450, Xebio, ROS)[,c(1, 9:21, 29:32)]
table(fun_save$level3_pathway_name)
write.csv(fun_save, './03_DGI/output_data/Mechanism_genes_YvsD.csv', row.names = FALSE)

fun_save_v2 <- rbind(phototransduction, peroxisome, P450, Xebio, ROS, Secretion, Porphyrin, Exopolysaccharide, Sulfer)[,c(1, 9:21, 29:32)]
write.csv(fun_save_v2, './03_DGI/output_data/Mechanism_genes_Fig3_YvsD_v2.csv', row.names = FALSE)

fun_save_v3 <- rbind(Nitrogen, phototransduction, peroxisome, P450, Xebio, ROS, Secretion, Porphyrin, Exopolysaccharide, Sulfer)[,c(1, 9:21, 29:32)]
write.csv(fun_save_v3, './03_DGI/output_data/Mechanism_genes_Fig3Nitrogen_YvsD_v3.csv', row.names = FALSE)

fun_save_v4 <- rbind(result_N, Nitrogen, phototransduction, peroxisome, P450, Xebio, ROS, Secretion, Porphyrin, Exopolysaccharide, Sulfer)[,c(1, 9:21, 29:32)]
write.csv(fun_save_v4, './03_DGI/output_data/Mechanism_genes_NO3NO2_YvsD_v5.csv', row.names = FALSE)


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

# Calculate the cumulative sum of len for each pathway
library(plyr)
ratio_results$Ratio <- sprintf('%0.4f',ratio_results$Ratio*100)
ratio_results$Ratio <- as.numeric(ratio_results$Ratio)
plot_data <- ddply(ratio_results, "Pathway",
                   transform, label_ypos=cumsum(Ratio))
head(plot_data)

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
        legend.position = c(0.8,0.8),# 'right', # 调整legend位置
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

ggplot(data=plot_data, aes(x=Pathway, y=Ratio, fill = factor(Var1, 
                                                             levels=c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5',  'Cluster 6', 'Cluster 7')))) +
  geom_bar(stat="identity") + 
  geom_text(aes(y=label_ypos, label=Ratio), vjust=1.2, 
            color="black", size=2.5) +
  scale_fill_manual(values = pale_25) + 
  labs(y = 'Ratio (%)',x = 'Pathways', title = 'Yellow light') +
  mytheme1


# --- key/functional genes ---
result_nitrate <- matched[grepl("nitrate", matched$SwissProt_Description), ] # cluster 3, the smallest clusters
result_nitrite <- matched[grepl("nitrite", matched$SwissProt_Description), ] # cluster 3, the smallest clusters
result_N <- rbind(result_nitrate, result_nitrite)
# save data
write.csv(result_N, './03_DGI/output_data/ClusterAssign_N_BvsD_7Cluster.csv', row.names = FALSE)
result_photo <- matched_kegg[grepl("photo", matched_kegg$SwissProt_Description), ] 
# result_test <- matched[grepl('Oxidative', matched$SwissProt_Description), ] 

result_dismutase <- matched[grepl("dismutase", matched$SwissProt_Description), ] # 12
# contribute to higher O2- via inhibiting the quenching of superoxide
result_oxidative <- matched[grepl("Oxidative", matched$SwissProt_Description), ]

# --- key subcellular information ---
library(dplyr)
secreted <- data.frame(filter(matched, matched$Sig0alP == 1 & matched$tmhmm == 0)) # based on the above results

# --- high abundance ---
# cluster18, cluster17, etc. —— 感觉cluster分布与表达量的correlation非常大 

# --- extract function cluster/genes ---
## V1 --- N (cluster 3): sequence 
Cluster3 <- filter(cluster_label, cluster_label$DGI_sub_label == 'Cluster 3') 
Cluster3_seq <- merge(Cluster3, seq, by = 'gene_id')
Cluster3_seq <- Cluster3_seq[, -c(2:6)]
## V2 --- N: annotation related to nitrogen
# --- cluster 3 ---
Cluster3 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 3') 
### save data
length(unique(Cluster3$gene_id)) # 26 genes in total 
write.csv(Cluster3, './03_DGI/output_data/Cluster3_YvsD.csv', row.names = FALSE)

# --- cluster 4 ---
## Phototransduction: cluster 1, 4
Cluster4 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 4') # 
### save data
length(unique(Cluster4$gene_id)) # 512 genes in total 
write.csv(Cluster4, './03_DGI/output_data/Cluster4_YvsD.csv', row.names = FALSE)

# --- cluster 1 ---
## Phototransduction: cluster 1, 4
Cluster1 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 1') # 
### save data
length(unique(Cluster1$gene_id)) # 512 genes in total 
write.csv(Cluster1, './03_DGI/output_data/Cluster1_YvsD.csv', row.names = FALSE)

# --- cluster 6 ---
## Phototransduction: cluster 1, 4
Cluster6 <- filter(matched_kegg, matched_kegg$DGI_sub_label == 'Cluster 6') # 
### save data
length(unique(Cluster6$gene_id)) # 512 genes in total 
write.csv(Cluster6, './03_DGI/output_data/Cluster6_YvsD.csv', row.names = FALSE)


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
cluster2_kegg <- merge(Cluster2, kegg_anno, by = 'gene_id')
length(unique(cluster2_kegg$gene_id)) # 26
cluster2_kegg <- merge(cluster2_kegg, fpkm, by = 'gene_id')
pathway_table <- data.frame(table(cluster2_kegg$level3_pathway_name))
colnames(pathway_table) <- c('level3', 'counts')
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

grep <- cluster2_kegg[grepl(pathway_table$Var1[1], cluster2_kegg$level3_pathway_name), ]
expre_matrix <- data.frame(t(apply(grep[,19:27], 2, sum)))

for (i in 2:length(pathway_table$Var1)){
  grep <- cluster2_kegg[grepl(pathway_table$Var1[i], cluster2_kegg$level3_pathway_name), ]
  expre_matrix <- rbind(expre_matrix, data.frame(t(apply(grep[,19:27], 2, sum))))
}

pathway_table <- cbind(pathway_table, expre_matrix)
