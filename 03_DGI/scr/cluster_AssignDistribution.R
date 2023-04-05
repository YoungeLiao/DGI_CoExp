# --- bar plot --- 
# === load packages ===
library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(paletteer)

# === load data ===
## to check the data distribution based on reads counts
cluster <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_YvsD_7Cluster.csv', header = T)

# data <- data.frame(rep(paste('cluster', 1), length(cluster$gene_id)))
data <- data.frame(c(rep('DGI', length(cluster$gene_id)),
                   rep('DGI_sub', length(cluster$gene_id)),
                   rep('Kmeans', length(cluster$gene_id)),
                   rep('Kmeans_sub', length(cluster$gene_id))))
colnames(data) <- 'methods'

clus_asign <- cluster[, 4]
for (i in 5:7){
  clus_asign <- c(clus_asign, cluster[, i])
}

data$cluster <- clus_asign
tab <- data.frame(table(data))

library(ggpubr)

pale_12 <- as.vector(paletteer_d('RColorBrewer::Set3'))
pale_20_2 <- as.vector(paletteer_d('ggthemes::Tableau_20'))
pale_51 <- as.vector(paletteer_d('ggsci::default_igv'))
pale_10 <- as.vector(paletteer_d('ggsci::default_jco'))
pale_8 <- as.vector(paletteer_d('RColorBrewer::Pastel2'))
pale_18 <- c(pale_10, pale_8)

ggbarplot(tab, x="methods", y="Freq", color = 'cluster',# add = "mean_se", color = "supp",
          palette = pale_18, position = position_dodge(0.8)) + 
  labs(x = '', y = 'Counts', title = 'Cluster number = 7 (YvsD)') + 
  mytheme1

# ====== the following are not used =====
mytheme1 <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 12), 
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


tem <- c(rep(c('Blue1'), length(fpkm[,1])),
         rep(c('Blue2'), length(fpkm[,1])),
         rep(c('Blue3'), length(fpkm[,1])),
         rep(c('Dark1'), length(fpkm[,1])),
         rep(c('Dark2'), length(fpkm[,1])),
         rep(c('Dark3'), length(fpkm[,1])),
         rep(c('Yellow1'), length(fpkm[,1])),
         rep(c('Yellow2'), length(fpkm[,1])),
         rep(c('Yellow3'), length(fpkm[,1])))
# class(tem) # check 

density_counts <- as.data.frame(c(as.data.frame(tem),
                                  as.data.frame(c(fpkm$Blue1, fpkm$Blue2, fpkm$Blue3,
                                                  fpkm$Dark1, fpkm$Dark2, fpkm$Dark3,
                                                  fpkm$Yellow1, fpkm$Yellow2, fpkm$Yellow3))))

colname <- c('samples', 'fpkm')
# class(density_counts) # check
colnames(density_counts) <- colname
head(density_counts) # check

# remove rows with counts = 0. 
temp_df <- density_counts[-which(density_counts$fpkm == 0),] 

# normalization: log 
log_density <- temp_df
log_density$fpkm <- log10(log_density$fpkm)

# plot
library(ggplot2)
p <- ggplot(log_density, aes(x=fpkm, color=samples)) +
  geom_density() + 
  scale_x_continuous("Log10(fpkm)") +
  scale_y_continuous("Density") # 设置坐标轴标签： https://blog.csdn.net/Arkardia/article/details/121851762

pale_BDY <- c('#85A4FE','#85A4FE','#85A4FE','#585858','#585858','#585858','#FEDF85','#FEDF85','#FEDF85')
p + scale_color_manual(values = pale_BDY) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face='bold'),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = c(0.9,0.7))
