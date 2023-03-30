# heatmap,URL:https://www.jianshu.com/p/ce7c8159fdbc
# annotation: URL: https://zhuanlan.zhihu.com/p/370477638
# loading data
DEGs_exp_rpkm <- read.csv("./01_Preprocessing/output_data/DEGs_BvsD_valid.csv", header = T, row.names = 1)
dim(DEGs_exp_rpkm)

data <- DEGs_exp_rpkm[,3:8]
summary(data)

library(pheatmap)
p <- pheatmap(data, cluster_row = TRUE, cluster_col = TRUE, show_rownames= FALSE,
              color=colorRampPalette(rev(c("#1D39C4", "#2E53EB", "#587DF7", "#85A4FE", "#ACC5FE","#D5E3FE", "#F8EFFE", "#B37EEB")))(1000),
              # color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(100),
              # scale = 'row', 
              fontsize = 14, # legend_breaks = c(0, 1, 2, 3, 4),
              # legend_labels=c("0","1","2","3","4"),
              cutree_rows = 24)


# Save clusters
row_cluster <- cutree(p$tree_row, k=24)
newOrder <- data[p$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder), names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
newOrder$gene_id <- rownames(newOrder)
write.csv(newOrder, "./02_DatasetAnalysis/output_data/HC_10cluster_BvsD.csv", row.names = FALSE)

# # cluster statistics & save data
# clus_num <- as.data.frame(table(newOrder$Cluster))
# head(clus_num)
# # Var1  Freq
# # 1    1 14043
# # 2    2  1148
# # 3    3  3456
# # 4    4   701
# # 5    5  8650
# # 6    6  1936
# colnames(clus_num) <- c('HC_Cluster', 'Freq')
# write.csv(clus_num, './output_data/PART3_metatranscriptomics/Fig3f_ClusterNumber_HC.csv')
# 
# # annotation 
# ## 1. all raw
# anno <- read.csv('./data/annotation/total_annotation_v2_RemoveColumn.csv')
# HCcluster_24 <- read.csv("./output_data/DEGs_HC_24cluster.csv")
# loc <- match(HCcluster_24$X, anno$Gene_ID)
# HCcluster_24 <- cbind(HCcluster_24, anno[loc,])
# HCcluster_24 <- HCcluster_24[,-1]
# ## save annotated cluster
# write.csv(HCcluster_24, "./output_data/DEGs_HC_24cluster_annotated_allraw.csv")
# 
# ## 2. eggNOG
# anno_egg <- read.csv('./data/annotation/eggNOG_annotation.csv')
# table(anno_egg$Class)
# 
# HCcluster_24 <- read.csv("./output_data/DEGs_HC_24cluster.csv")
# loc <- match(HCcluster_24$X, anno_egg$Gene)
# HCcluster_24 <- cbind(HCcluster_24, anno_egg[loc,])
# HCcluster_24 <- HCcluster_24[,-1]
# ## save annotated cluster
# write.csv(HCcluster_24, "./output_data/DEGs_HC_24cluster_annotated_eggNOG.csv")
# 
