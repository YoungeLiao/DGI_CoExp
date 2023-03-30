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
# df_scaled <- scale(df_numeric)
# head(df_scaled)
# Blue1       Blue2       Blue3        Dark1        Dark2        Dark3      Yellow1
# [1,]  0.03301536  0.05481294  0.03083476 -0.003854804 -0.002667458 -0.006791027 -0.009719643
# [2,]  0.07186596 -0.05842421 -0.05065870 -0.029161613 -0.027780067 -0.028197508 -0.028878230
# [3,]  0.04745408 -0.05882366 -0.04879494 -0.029417099 -0.029160837 -0.027951773 -0.028562388
# [4,]  0.06217157 -0.05905572 -0.05127935 -0.029123474 -0.028122684 -0.027817580 -0.029337911
# [5,]  0.02658146 -0.03836319 -0.04216777 -0.029587130 -0.029250140 -0.028356576 -0.030352154
# [6,] -0.01027961 -0.04792726 -0.04387416 -0.029632703 -0.029101252 -0.027526235 -0.030235163
# dim(df_scaled) # [1] 25886     9
log_df <- log10(df_numeric+1)
summary(log_df)
# Dark1             Dark2             Dark3             Yellow1          Yellow2          Yellow3      
# Min.   :0.00000   Min.   :0.00000   Min.   :0.000000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
# 1st Qu.:0.01265   1st Qu.:0.00000   1st Qu.:0.006166   1st Qu.:0.2425   1st Qu.:0.2616   1st Qu.:0.0345  
# Median :0.07881   Median :0.04889   Median :0.064589   Median :0.4386   Median :0.4760   Median :0.1089  
# Mean   :0.22820   Mean   :0.19753   Mean   :0.219520   Mean   :0.5330   Mean   :0.5675   Mean   :0.3799  
# 3rd Qu.:0.30781   3rd Qu.:0.23297   3rd Qu.:0.302557   3rd Qu.:0.7026   3rd Qu.:0.7552   3rd Qu.:0.6577  
# Max.   :3.26073   Max.   :3.08873   Max.   :3.148433   Max.   :3.6319   Max.   :3.8017   Max.   :3.8108  

# 3. --- k-means ---
n_cluster <- 24
set.seed(0)
kmean <- kmeans(log_df, n_cluster, iter.max = 500)


# # visualization
# plotcluster(df_scaled, kmean$cluster)

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
