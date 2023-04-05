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

