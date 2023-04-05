
# === load data ===
genes_list <- read.csv('./02_DatasetAnalysis/output_data/top30_geneid.csv', header = TRUE)
colnames(genes_list) <- 'gene_id'

# v1: match DEGs of blue light
clusterAssaign <- read.csv('./03_DGI/output_data/tsne_BvsD_labels_all_7cluster.csv', header = TRUE)[,c(1,5)]
matched <- merge(genes_list, clusterAssaign, by = 'gene_id', all.y = FALSE) # 29 matched

# v2: match DEGs of yellow light
clusterAssaign <- read.csv('./03_DGI/output_data/tsne_YvsD_labels_all_7cluster.csv', header = TRUE)[,c(1,5)]
matched <- merge(genes_list, clusterAssaign, by = 'gene_id', all.y = FALSE) # 1 matched


# --- node data of blue light ---
node_raw <- read.csv('./03_DGI/data/topology/Nodes_raw_blue_genes.csv', header = TRUE)
kegg_anno <- read.csv('./data/kegg_annotation_all.csv')
ref <- kegg_anno[,c(1,9)]
# === preprocessing ===
labeled <- merge(node_raw, ref)

# === save data ===
write.csv(labeled, './03_DGI/output_data/nodes_LabeledGenes_blue.csv', row.names = FALSE)


