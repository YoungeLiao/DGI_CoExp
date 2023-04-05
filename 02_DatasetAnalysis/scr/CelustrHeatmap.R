# loading data
DEGs_exp_rpkm <- read.csv("./01_Preprocessing/output_data/DEGs_filtered.csv",header = T)
dim(DEGs_exp_rpkm)


data_raw <- DEGs_exp_rpkm[,2:length(DEGs_exp_rpkm)]
rownames(data_raw) <- DEGs_exp_rpkm$gene_id

data <- log10(data_raw[,18:26] + 1)
summary(data)


# 
library(pheatmap)

p <- pheatmap(data, cluster_row = TRUE, cluster_col = TRUE, show_rownames= FALSE,
              color=colorRampPalette(rev(c("#FE7875","white","#587DF7")))(1000),
              # color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(100),
              scale = 'row', 
              fontsize = 14,
              cutree_rows = 7)
# p
p <- pheatmap(data[c(1, 100, 1000, 1200, 1500, 2000, 3000, 4000, 5000, 6000),], cluster_row = TRUE, cluster_col = TRUE, show_rownames= FALSE,
              color=colorRampPalette(rev(c("#2979FE", "white", "#FE9F00")))(1000),
              # color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(100),
              scale = 'row', border_color = 'white',
              fontsize = 14)


# Save clusters
row_cluster <- cutree(p$tree_row, k=7)
newOrder <- data[p$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder), names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
write.csv(newOrder, "./output_data/Allgenes_HC_24cluster.csv")

# cluster statistics & save data
clus_num <- as.data.frame(table(newOrder$Cluster))
head(clus_num)
# Var1  Freq
# 1    1 14043
# 2    2  1148
# 3    3  3456
# 4    4   701
# 5    5  8650
# 6    6  1936
colnames(clus_num) <- c('HC_Cluster', 'Freq')
write.csv(clus_num, './output_data/PART3_metatranscriptomics/Fig3f_ClusterNumber_HC.csv')

# annotation 
## 1. all raw
anno <- read.csv('./data/annotation/total_annotation_v2_RemoveColumn.csv')
HCcluster_24 <- read.csv("./output_data/DEGs_HC_24cluster.csv")
loc <- match(HCcluster_24$X, anno$Gene_ID)
HCcluster_24 <- cbind(HCcluster_24, anno[loc,])
HCcluster_24 <- HCcluster_24[,-1]
## save annotated cluster
write.csv(HCcluster_24, "./output_data/DEGs_HC_24cluster_annotated_allraw.csv")

## 2. eggNOG
anno_egg <- read.csv('./data/annotation/eggNOG_annotation.csv')
table(anno_egg$Class)

HCcluster_24 <- read.csv("./output_data/DEGs_HC_24cluster.csv")
loc <- match(HCcluster_24$X, anno_egg$Gene)
HCcluster_24 <- cbind(HCcluster_24, anno_egg[loc,])
HCcluster_24 <- HCcluster_24[,-1]
## save annotated cluster
write.csv(HCcluster_24, "./output_data/DEGs_HC_24cluster_annotated_eggNOG.csv")

