library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)


# === 1. load data & data preprocessing ===
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

ggplot(tsne_result,aes(tSNE1,tSNE2)) + 
  geom_point(colour = '#85A4FE', alpha = 0.8, size = 1.5) + 
  # geom_point(colour="#85A4FE", size = 1) + 
  labs(title = 'Raw') + 
  theme_bw() + 
  mytheme

# === 3. construct new data ===
position <- cbind(DEGs$gene_id, tsne_result)
colnames(position) <- c('gene_id', 'tSNE1', 'tSNE2')
graph_data <- cbind(position, data)

# subcellular location
subloca <- data.frame(subloca_raw$gene_id, subloca_raw$SignalP, subloca_raw$tmhmm)
colnames(subloca) <- c('gene_id', 'Signalp', 'tmhmm')

graph_data <- merge(graph_data, subloca, by = 'gene_id', all.x = TRUE)
table(graph_data[, 10:11])
# tmhmm
# Signalp     0     1
# 0 19285  4168
# 1  1245   579

# save data with subcellular information
write.csv(graph_data, './03_DGI/output_data/graph_data_BvsD_SubCell_dgi.csv', row.names = FALSE)
# baseline dataset
baseline <- graph_data[,1:9]
write.csv(baseline, './03_DGI/output_data/graph_data_BvsD_baseline.csv', row.names = FALSE)
