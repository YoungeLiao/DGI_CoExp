# loading data
DEGs_exp_rpkm <- read.csv("./01_Preprocessing/output_data/DEGs_YvsD_valid.csv", header = T, row.names = 1)
dim(DEGs_exp_rpkm)

data_raw<- DEGs_exp_rpkm[,20:25]
summary(data_raw)

data <- log10(data_raw + 1)
summary(data)
data_centered <- scale(data, center = T, scale = F)
summary(data_centered)

library(pheatmap)
p <- pheatmap(data, cluster_row = TRUE, cluster_col = TRUE, show_rownames= FALSE,
              color=colorRampPalette(rev(c("#D38706", "#F9AC13", "#FEC53C", "#FED565", "#FEE48F","#FEF0B8", "#F8EFFE", "#B37EEB")))(1000),
              # color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(100),
              # scale = 'row', 
              fontsize = 14, legend_breaks = c(0, 1, 2, 3, 4),
              legend_labels=c("0","1","2","3","4"),
              cutree_rows = 10)


# Save clusters
row_cluster <- cutree(p$tree_row, k=10)
newOrder <- data[p$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder), names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
newOrder$gene_id <- rownames(newOrder)
write.csv(newOrder, "./02_DatasetAnalysis/output_data/HC_10cluster_YvsD.csv", row.names = FALSE)

