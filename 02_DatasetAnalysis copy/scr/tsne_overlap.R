library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)


DEGs_all <- read.csv("01_Preprocessing/output_data/DEGs_filtered.csv", header = T) # 25886
DEGs_B <- read.csv('output_data/graph_data_BvsD_SubCel.csv', header = T)
DEGs_Y <- read.csv('output_data/graph_data_YvsD_SubCel.csv', header = T)

# === extract overlap & labeling ===
loc <- match(DEGs_Y$gene_id, DEGs_B$gene_id)
loc <- loc[!is.na(loc)]

overlap <- DEGs_B[loc, ]$gene_id
overlap <- data.frame(DEGs_B[loc,]$gene_id, rep('Overlap', length(overlap))) # 3757
colnames(overlap) <- c('gene_id', 'label')
BvsD <- DEGs_B$gene_id
BvsD <- data.frame(DEGs_B$gene_id, rep('Blue', length(BvsD))) # 25277
colnames(BvsD) <- c('gene_id', 'label')
YvsD <- DEGs_Y$gene_id
YvsD <- data.frame(DEGs_Y$gene_id, rep('Yellow', length(YvsD))) # 4366
colnames(YvsD) <- c('gene_id', 'label')

# normalization 
fpkm <- DEGs_all[,19:27]
fpkm <- log(fpkm + 1, 10)
gene_id <- DEGs_all$gene_id
fpkm <- cbind(gene_id, fpkm)

# match
# loc <- match(DEGs$gene_id, fpkm_raw$gene_id)
# genes <- fpkm[loc,]
data <- fpkm[,2:length(fpkm)]


# === 3. tsne ===

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

# theme

top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
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
        # legend.position = 'right', # 调整legend位置
        legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
mytheme2 <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 14), 
        # legend.title = element_text(size = 16),
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.title = element_blank(), 
        legend.title = element_blank(),
        legend.position = 'right', # 调整legend位置
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))


# plot
ggplot(tsne_result,aes(tSNE1,tSNE2)) + 
  geom_point(colour = 'grey', alpha = 0.5, size = 2.0) + 
  # geom_point(colour="#85A4FE", size = 1) + 
  labs(title = 'DEGs of YvsD') + 
  theme_bw() + 
  mytheme1

# === 3. construct graph data ===

## --- 2.1 tsne data ---
position <- cbind(DEGs_all$gene_id, tsne_result)
colnames(position) <- c('gene_id', 'tSNE1', 'tSNE2')
graph_data <- cbind(position, data)


## --- 2.2 YvsD SubCel ---
# --- load SignalP and tmhmm ---
signalp_raw <- read.csv("../data/SignalP.csv", header = T)
tmhmm_raw <- read.csv("../data/tmhmm.csv", header = T)

loc <- match(DEGs_all$gene_id, signalp_raw$gene_id)
signalp <- signalp_raw[loc,2:3]

loc <- match(DEGs_all$gene_id, tmhmm_raw$gene_id)
tmhmm <- tmhmm_raw[loc,2:3]

# normalization 
signalp <- log(signalp + 1, 10)
gene_id <- fpkm$gene_id
fpkm <- cbind(gene_id, fpkm)
tmhmm <- log(tmhmm + 1, 10)

graph_data <- cbind(position, data, signalp, tmhmm)
write.csv(graph_data, 'output_data/graph_data_YvsD_SubCel.csv', row.names = FALSE)

# === 3. construct label data ===
tsne_plot <- cbind(DEGs_all$gene_id, tsne_result)
colnames(tsne_plot) <- c('Gene_ID', 'tSNE1', 'tSNE2')
tsne_plot$label <- 'Normal'

for (i in 1:length(tsne_plot$Gene_ID)){
  if(isTRUE(match(tsne_plot[i,]$Gene_ID, YvsD$gene_id) != 'NA')){
    tsne_plot[i,]$label <- YvsD[match(tsne_plot[i,]$Gene_ID, YvsD$gene_id),]$label
  }
}


for (i in 1:length(tsne_plot$Gene_ID)){
  if(isTRUE(match(tsne_plot[i,]$Gene_ID, BvsD$gene_id) != 'NA')){
    tsne_plot[i,]$label <- BvsD[match(tsne_plot[i,]$Gene_ID, BvsD$gene_id),]$label
  }
}

for (i in 1:length(tsne_plot$Gene_ID)){
  if(isTRUE(match(tsne_plot[i,]$Gene_ID, overlap$gene_id) != 'NA')){
    tsne_plot[i,]$label <- overlap[match(tsne_plot[i,]$Gene_ID, overlap$gene_id),]$label
  }
}

table(tsne_plot$label)
# Blue Overlap  Yellow 
# 21520    3757     609 

# === 3. visualization ===
ggplot(tsne_plot, aes(tSNE1,tSNE2,color = label)) + 
  geom_point(alpha = 0.5, size = 2.0) +
  theme_bw() + 
  scale_color_manual(values = c("#85A4FE", '#36CFC8', '#FED565')) + 
  labs(title = 'DEGs') + 
  mytheme1

## === 4. cluster: eggNOG ===
eggNOG <- read.csv('./output_data/EggNOG_matched_singlelabel.csv', header = TRUE) 

tsne_plot$NOG_label <- 'Normal'

# match cluster
for (i in 1:length(tsne_plot$Gene_ID)){
  if(isTRUE(match(tsne_plot[i,]$Gene_ID, eggNOG$gene_id) != 'NA')){
    tsne_plot[i,]$NOG_label <- eggNOG[match(tsne_plot[i,]$Gene_ID, eggNOG$gene_id),]$eggNOG_Class
  }
}

NOG_label <- sub('Normal', 'NA', tsne_plot$NOG_label)
tsne_plot$NOG_label <- NOG_label

# theme
table(tsne_plot$NOG_label)

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
        legend.position = 'right', # 调整legend位置
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

ggplot(tsne_plot, aes(tSNE1,tSNE2,color = NOG_label)) + 
  geom_point(alpha = 0.7, size = 2.0) +
  theme_bw() + 
  scale_color_manual(values = pale_25) + 
  labs(title = 'eggNOG') + 
  mytheme2

