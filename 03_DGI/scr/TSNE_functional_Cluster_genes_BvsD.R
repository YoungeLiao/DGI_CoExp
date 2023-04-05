# === 0. basic setting visualization ===
# packages
library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)

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
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 10, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
        # legend.title = element_text(size = 16),
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.title = element_blank(), 
        legend.title = element_blank(),
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

# === 1. load data & data preprocessing for tsne ===
fpkm_raw <- read.table("./data/gene_fpkm_all_samples.txt", header = T)
DEGs <- read.csv("./01_Preprocessing/output_data/DEGs_BvsD_valid.csv", header = T)
subloca_raw <- read.csv("./02_DatasetAnalysis/output_data/subloca.csv", row.names = 1)

# normalization 
fpkm <- log(fpkm_raw[2:length(fpkm_raw)] + 1, 10)
gene_id <- fpkm_raw$gene_id
fpkm <- cbind(gene_id, fpkm)

# match
loc <- match(DEGs$gene_id, fpkm_raw$gene_id)
genes <- fpkm[loc,]
data <- genes[,2:7]


# === 2. tsne background ===
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



# === 3. combine data ===
# save data with subcellular information
position <- cbind(DEGs$gene_id, tsne_result)
colnames(position) <- c('gene_id', 'tSNE1', 'tSNE2')
graph_data <- cbind(position, data)

# --- load DGI subcellular data---
DGI_sub <- read.table('./03_DGI/data/7clusters/BvsD_SubCell_dgi/20000epochs/types.txt', header = F) 

# labeling for plot
tsne_plot <- cbind(position, DGI_sub)
tsne_plot$DGI_sub_label <- paste('Cluster',DGI_sub$V1+1)
# tsne_plot_cluster3 <- filter(tsne_plot, tsne_plot$DGI_sub_label == 'Cluster 3') 
tsne_plot_cluster4 <- filter(tsne_plot, tsne_plot$DGI_sub_label == 'Cluster 4') 
tsne_plot_cluster5 <- filter(tsne_plot, tsne_plot$DGI_sub_label == 'Cluster 5') 
tsne_plot_cluster7 <- filter(tsne_plot, tsne_plot$DGI_sub_label == 'Cluster 7') 
tsne_result$DGI_sub_label <- 'All DEGs'

# --- load functional genes ---
phototransduction <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_BvsD_7Cluster.csv', header = T)
nitrogen <- read.csv('./03_DGI/output_data/ClusterAssign_N_BvsD_7Cluster_filtered.csv', header = T)

# labeling for plot
tsne_plot_phototrans <- phototransduction[, 2:3]
tsne_plot_phototrans$label <- 'phototransduction'
tsne_plot_nitrogen <- nitrogen[c(1:5,7), 2:3]
tsne_plot_nitrogen$label <- c('nirK (P25006)', 'narK1 (Q9RA46)', 'narG (P09152)', 'aniA (Q9JTB8)', 'nirK (Q53239)', 'narK2 (Q9RA45)')


# --- plot ---

custo_color <- c('grey', '#FF69B4FF', "#FB8072FF","#80B1D3FF", "#B3DE69FF", 
                       "#40E0D0FF", "#9370DBFF","#9370DBFF", 
                       '#1E90FFFF', '#1E90FFFF', "#FFED6FFF")

p0 <- ggplot() +  # tsne_result, aes(tSNE1,tSNE2, color = DGI_sub_label)
  # geom_point(alpha = 0.1, size = 2) + 
  geom_point(data = tsne_result, aes(tSNE1, tSNE2, color=DGI_sub_label), 
            alpha = 0.05, size = 2) + 
  geom_point(data = tsne_plot_cluster7 , aes(tSNE1, tSNE2, color=DGI_sub_label), alpha = 0.6, size = 2) + 
  geom_point(data = tsne_plot_cluster4 , aes(tSNE1, tSNE2, color=DGI_sub_label), alpha = 0.6, size = 2) + 
  geom_point(data = tsne_plot_cluster5 , aes(tSNE1, tSNE2, color=DGI_sub_label), alpha = 0.6, size = 2) + 
  geom_point(data = tsne_plot_phototrans , aes(tSNE1, tSNE2, color=label), alpha = 0.5, size = 3) + 
  geom_point(data = tsne_plot_nitrogen , aes(tSNE1, tSNE2, color=label), alpha = 0.5, size = 3) + 
  labs(title = 'DEGs of blue light') + 
  scale_color_manual(values = custo_color) + 
  theme_bw() + 
  mytheme1
p0
library(ggrepel)
# set.seed(007)
p1 <- p0 + geom_text_repel(data=tsne_plot_phototrans, aes(tSNE1, tSNE2,label=label),
                     force=20,color="grey20",
                     size=3,
                     point.padding = 0.5,hjust = 0.5,
                     # arrow = arrow(length = unit(0.01, "npc"),
                     #               type = "open", ends = "last"),
                     segment.color="grey20",
                     segment.size=0.6,
                     segment.alpha=0.8,
                     nudge_x=0,
                     nudge_y=0
                     )
p2 <- p1 + geom_text_repel(data=tsne_plot_nitrogen, aes(tSNE1, tSNE2,label=label),
                           force=20,color="grey20",
                           size=3,
                           point.padding = 0.5,hjust = 0.5,
                           # arrow = arrow(length = unit(0.01, "npc"),
                           #               type = "open", ends = "last"),
                           segment.color="grey20",
                           segment.size=0.6,
                           segment.alpha=0.8,
                           nudge_x=0,
                           nudge_y=0
)

p2
