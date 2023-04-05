# k-means for metatranscriptomics
library(fpc)
set.seed(0)

# 1. --- load data ---
fpkm <- read.table("./data/gene_fpkm_all_samples.txt", header = T)
DEGs <- read.csv("./01_Preprocessing/output_data/DEGs_BvsD_valid.csv", header = T)
graph_data <- read.csv('./03_DGI/output_data/graph_data_BvsD_SubCell_dgi.csv', header = T)

# v1 - basic matrix
loc <- match(DEGs$gene_id, fpkm$gene_id)
genes <- fpkm[loc,]
data <- genes[,2:7] # blue
# data <- genes[,5:10] # yellow

# v2 - graph data with subcellular information
genes <- fpkm[loc,]
data <- graph_data[,4:11]


rownames(data) <- genes$gene_id
head(data)
class(data) # data.frame


# 2. --- Normalization ---
df_numeric <- as.data.frame(apply(data[,1:length(data)], 2, as.numeric))


log_df <- log10(df_numeric+1)
summary(log_df)

# 3. --- k-means ---
n_cluster <- 24
set.seed(0)
kmean <- kmeans(log_df, n_cluster, iter.max = 500)


# results data profile
table(kmean$cluster)

# 4. --- save cluster number ---
Kmeans_clusNum <- as.data.frame(table(kmean$cluster))
colnames(Kmeans_clusNum) <- c('Kmeans_Cluster', 'Counts')
head(Kmeans_clusNum)
# Kmeans_Cluster Counts
# 1              1    245
# 2              2     55
# 3              3    267
# 4              4    129
# 5              5    282
# 6              6    446
# write.csv(Kmeans_clusNum, './output_data/ClusterNumber_Kmeans_YvsD.csv')

# 5. save genes cluster assignment 
head(kmean$cluster) # [1]  5  6  5 12 11  8
length(kmean$cluster) # 4366

KmeanClu <- cbind(rownames(data), kmean$cluster)
colnames(KmeanClu) <- c('Gene_ID', 'Kmeans_Cluster')
write.csv(KmeanClu, './03_DGI/data/10clusters/Kmeans_24cluster_BvsD_subcel.csv', row.names = FALSE)
write.csv(KmeanClu, './03_DGI/data/10clusters/Kmeans_24cluster_BvsD.csv', row.names = FALSE)
