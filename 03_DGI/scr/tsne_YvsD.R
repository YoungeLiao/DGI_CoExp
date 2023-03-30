library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)


# === 1. load data & data preprocessing ===
fpkm_raw <- read.table("./data/gene_fpkm_all_samples.txt", header = T)
DEGs <- read.csv("./01_Preprocessing/output_data/DEGs_YvsD_valid.csv", header = T)
subloca_raw <- read.csv("./02_DatasetAnalysis/output_data/subloca.csv", row.names = 1)

# normalization 
fpkm <- log(fpkm_raw[2:length(fpkm_raw)] + 1, 10)
gene_id <- fpkm_raw$gene_id
fpkm <- cbind(gene_id, fpkm)

# match
loc <- match(DEGs$gene_id, fpkm_raw$gene_id)
genes <- fpkm[loc,]
data <- genes[,5:10]


# === 2. tsne ===
set.seed(0)
tsne_fpkm = Rtsne(
  data,
  dims = 2,
  pca = T,
  max_iter = 1000,
  theta = 0.4,
  perplexity = 20,
  verbose = F
)

tsne_result = as.data.frame(tsne_fpkm$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")

# === 2. visualization without label ===
library(ggplot2)
# theme
top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.title = element_blank(), 
        # legend.title = element_blank(),
        # legend.position = c(0.1, 0.8), # 调整legend位置
        legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
mytheme1 <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.title = element_blank(), 
        # legend.title = element_blank(),
        legend.position = 'left', # 调整legend位置
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
# pallete 
Pale_51 <- as.vector(paletteer_d('ggsci::default_igv'))
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

ggplot(tsne_result,aes(tSNE1,tSNE2)) + 
  geom_point(colour = '#434343', alpha = 0.8, size = 1.5) + 
  # geom_point(colour="#85A4FE", size = 1) + 
  labs(title = 'Raw') + 
  theme_bw() + 
  mytheme


# === 4. obtain graph data ===
# save data with subcellular information
position <- cbind(DEGs$gene_id, tsne_result)
colnames(position) <- c('gene_id', 'tSNE1', 'tSNE2')
graph_data <- cbind(position, data)


# === 4. Results visualization on DGI  ===
DGI <- read.table('./03_DGI/data/7clusters/YvsD_baseline/types.txt', header = F) 

# init label
tsne_plot <- cbind(position, DGI)
colnames(tsne_plot) <- c('gene_id', 'tSNE1', 'tSNE2', 'DGI_label')
tsne_plot$DGI_label <- paste('Cluster',tsne_plot$DGI_label+1)


ggplot(tsne_plot, aes(tSNE1,tSNE2,color = DGI_label)) + 
  geom_point(alpha = 0.7, size = 2) +
  # geom_density(alpha = 0.6)+
  theme_bw() + 
  scale_color_manual(values = pale_25) + 
  labs(title = 'DGI') + 
  # SCI = 0.42897538894186776 (EM)
  mytheme1
# save fig name: TSNE_YvsD_7cluster_dgi_baseline_4.8_4.5.pdf

# === 4. Results visualization on DGI with subcellular ===
DGI_sub <- read.table('./03_DGI/data/24clusters/YvsD_SubCell_dgi/types.txt', header = F) 


# init label
# tsne_plot <- cbind(tsne_plot, DGI_Sub)
# colnames(tsne_plot) <- c('otu', 'tSNE1', 'tSNE2', 'DGI_label')
tsne_plot$DGI_sub_label <- paste('Cluster',DGI_sub$V1+1)


ggplot(tsne_plot, aes(tSNE1,tSNE2, color = DGI_sub_label)) + 
  geom_point(alpha = 0.7, size = 2) +
  # geom_density(alpha = 0.6)+
  theme_bw() + 
  scale_color_manual(values = pale_25) + 
  labs(title = 'Subcellular DGI') + 
  # SCI = 0.42440431406659995 (EM)
  mytheme
# save fig name: TSNE_YvsD_7cluster_dgi_subcel_4.8_4.5.pdf

# === 5. Kmeans results visualization ===
Kmean <- read.csv('./03_DGI/data/24clusters/Kmeans_24cluster_YvsD.csv', header = TRUE) 

tsne_plot$Kmeans_label <- 'Normal'

# match cluster
for (i in 1:length(tsne_plot$gene_id)){
  if(isTRUE(match(tsne_plot[i,]$gene_id, Kmean$Gene_ID) != 'NA')){
    tsne_plot[i,]$Kmeans_label <- Kmean[match(tsne_plot[i,]$gene_id, Kmean$Gene_ID),]$Kmeans_Cluster
  }
}

tsne_plot$Kmeans_label <- paste('Cluster',tsne_plot$Kmeans_label)

ggplot(tsne_plot, aes(tSNE1,tSNE2,color = Kmeans_label)) + 
  geom_point(alpha = 0.8, size = 2) +
  # geom_density(alpha = 0.6)+
  theme_bw() + 
  scale_color_manual(values = pale_25) + 
  labs(title = 'Kmeans') + 
  # SCI = 0.42878412195481574 (Expression matrix)
  mytheme

# === 6. Subcellular Kmeans results visualization ===
Kmean_subcel <- read.csv('./03_DGI/data/24clusters/Kmeans_24cluster_YvsD_subcel.csv', header = TRUE) 

tsne_plot$Kmeans_subcel_label <- 'Normal'

# match cluster
for (i in 1:length(tsne_plot$gene_id)){
  if(isTRUE(match(tsne_plot[i,]$gene_id, Kmean_subcel$Gene_ID) != 'NA')){
    tsne_plot[i,]$Kmeans_subcel_label <- Kmean_subcel[match(tsne_plot[i,]$gene_id, Kmean$Gene_ID),]$Kmeans_Cluster
  }
}

tsne_plot$Kmeans_subcel_label <- paste('Cluster',tsne_plot$Kmeans_subcel_label)

ggplot(tsne_plot, aes(tSNE1,tSNE2,color = Kmeans_subcel_label)) + 
  geom_point(alpha = 0.8, size = 2) +
  # geom_density(alpha = 0.6)+
  theme_bw() + 
  scale_color_manual(values = pale_25) + 
  labs(title = 'Subcellular Kmeans') + # SCI = 0.12240327342437617 (Expression matrix)
  mytheme

# === 6. HC ===
HC <- read.csv('./02_DatasetAnalysis/output_data/HC_24cluster_YvsD.csv', header = TRUE) 

tsne_plot$HC_label <- 'Normal'

# match cluster
for (i in 1:length(tsne_plot$gene_id)){
  if(isTRUE(match(tsne_plot[i,]$gene_id, HC$gene_id) != 'NA')){
    tsne_plot[i,]$HC_label <- HC[match(tsne_plot[i,]$gene_id, HC$gene_id),]$Cluster
  }
}

tsne_plot$HC_label <- paste('Cluster',tsne_plot$HC_label)

ggplot(tsne_plot, aes(tSNE1,tSNE2,color = HC_label)) + 
  geom_point(alpha = 0.8, size = 2) +
  # geom_density(alpha = 0.6)+
  theme_bw() + 
  scale_color_manual(values = pale_25) + 
  labs(title = 'Hierarchical clustering') + # SCI = 0.498278177395585 (EM)
  mytheme1


# save data
write.csv(tsne_plot, './03_DGI/output_data/tsne_YvsD_labels_all_24cluster.csv', row.names = FALSE)
