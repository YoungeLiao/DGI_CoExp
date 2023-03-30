## === DEGs ===
top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme<-theme_classic()+
  theme(# plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 14),
    # legend.title = element_text(size = 16),
    legend.title = element_blank(),
    # legend.position = c(0.9, 0.9), # 调整legend位置
    legend.position = 'top',
    legend.background = element_blank(),
    plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                     units="inches"))

# === load background data ===
kegg_anno <- read.csv('./data/kegg_annotation_all.csv')
annotation <- read.csv('./data/Annotation_Expre.csv', header = TRUE)
fpkm <- read.table('./data/gene_fpkm_all_samples.txt', header = TRUE)

# === laod data and pretreatment: phototransduction ===
genes_Y <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_YvsD_7Cluster.csv', header = TRUE)
genes_B <- read.csv('./03_DGI/output_data/ClusterAssign_phototransduction_BvsD_7Cluster.csv', header = TRUE)

genes_YB <- rbind(genes_Y, genes_B)
genes_YB$Mean_expre <- apply(genes_YB[, 13:21], 1, mean)
thre <- 3
genes_YB_filtered <- unique(filter(genes_YB, genes_YB$Mean_expre > thre))
# length(unique(genes_YB_filtered$SwissProt_Description)) # 8 swissprot in total
# length(unique(genes_YB_filtered$gene_id)) # 10 genes in total
# save raw data (not summarized)
write.csv(genes_YB_filtered, './03_DGI/output_data/Phototransduction_BY_merged_raw.csv', row.names = FALSE)

# === laod data and pretreatment: phototransduction ===
# --- nitrogen metabolism ---
## v1 --- manually
# nitrogen_Y <- read.csv('./03_DGI/data/nitrogen_YvsD.csv', header = TRUE)
## v2 --- automatically
genes_Y <- read.csv('./03_DGI/output_data/ClusterAssign_N_YvsD_7Cluster.csv', header = TRUE)
genes_B <- read.csv('./03_DGI/output_data/ClusterAssign_N_BvsD_7Cluster.csv', header = TRUE)

genes_YB <- rbind(genes_Y, genes_B)
# filter data based on mean expression
genes_YB$Mean_expre <- apply(genes_YB[, 13:21], 1, mean)
thre <- 3
genes_YB_filtered <- unique(filter(genes_YB, genes_YB$Mean_expre > thre))
# length(unique(genes_YB_filtered$SwissProt_Description)) # 11 swissprot in total
# length(unique(genes_YB_filtered$gene_id)) # 19 genes in total
# save raw data (not summarized)
write.csv(genes_YB_filtered, './03_DGI/output_data/Nitrogen_BY_merged_raw.csv', row.names = FALSE)
## !!! tips: cluster is not referenable, since it includes both clusters for blue and yellow light

# merge and summary based on swissprot
library(dplyr)
## v1: by description
rawdata <- genes_YB_filtered
SwissProt <- rawdata %>% 
  group_by(SwissProt_Description) %>%
  summarise(B1mean = mean(Blue1), 
            B2mean = mean(Blue2),
            B3mean = mean(Blue3),
            D1mean = mean(Dark1), 
            D2mean = mean(Dark2),
            D3mean = mean(Dark3),
            Y1mean = mean(Yellow1), 
            Y2mean = mean(Yellow2),
            Y3mean = mean(Yellow3),
            Mean = mean(Mean_expre)) 
## v1: by ID
rawdata <- genes_YB_filtered
SwissProt_ID <- rawdata %>% 
  group_by(SwissProt_ID) %>%
  summarise(B1mean = mean(Blue1), 
            B2mean = mean(Blue2),
            B3mean = mean(Blue3),
            D1mean = mean(Dark1), 
            D2mean = mean(Dark2),
            D3mean = mean(Dark3),
            Y1mean = mean(Yellow1), 
            Y2mean = mean(Yellow2),
            Y3mean = mean(Yellow3),
            Mean = mean(Mean_expre)) 

# --- save data: phototransduction ---
top_phototran_data <- merge(SwissProt_ID, rawdata, by = 'SwissProt_ID')[, c(1:12, 16, 20:22)]
top_phototran_data <- unique(top_phototran_data)
write.csv(top_phototran_data, './03_DGI/output_data/Phtotrans_TopFiltered_7Cluster.csv', row.names = FALSE)

# --- save data: nitrogen metabolism ---
top_N_data <- merge(SwissProt_ID, rawdata, by = 'SwissProt_ID')[, c(1:12, 16, 20:22)]
top_N_data <- unique(top_N_data)
write.csv(top_N_data, './03_DGI/output_data/Nitrogen_TopFiltered_7Cluster.csv', row.names = FALSE)


# === data format ===
groupname <- c(rep('Blue', length(SwissProt_ID$SwissProt_ID)*3), 
               rep('Dark', length(SwissProt_ID$SwissProt_ID)*3),
               rep('Yellow', length(SwissProt_ID$SwissProt_ID)*3))
SwissProt_id <- SwissProt_ID$SwissProt_ID

plotdata <- data.frame(SwissProt_id,
                       groupname, 
                       as.vector(as.matrix(SwissProt_ID[,2:10])))
colnames(plotdata) <- c('group', 'SwissProt_id', 'value')
plotdata$value <- log10(plotdata$value + 1)

# save data
# write.csv(result_N_filtered, './03_DGI/output_data/ClusterAssign_N_BvsD_7Cluster.csv', row.names = FALSE)

# v1: manually nitrogen metabolism of yellow light
# light <- rep(nitrogen$light, 3)
# genes <- c(rep('NarK1', length(nitrogen$light)), rep('NarG', length(nitrogen$light)), rep('AniA', length(nitrogen$light)))
# values <- log10(c(nitrogen$NarK1, nitrogen$NarG, nitrogen$AniA)+1)
# data <- data.frame(light, genes, values)


# === 2. Visualization ===
library(ggplot2)
library(ggpubr)
Cust_palette <- c('#82B0FE', '#424242', '#FED180') 
p <- ggbarplot(plotdata, x="group", y="value", color = "SwissProt_id",
               palette = Cust_palette, 
               width = 0.7, alpha = 0.8,
               add = "mean_se", add.params = list(width = 0.4, alpha = 0.8),
               position = position_dodge(0.8)
               # add = "jitter", add.params = list(size = 0.3, alpha = 0.7)
)
p <- p + geom_jitter(data=plotdata,
                aes(x=group, y=value, color = SwissProt_id),
                size = 2, alpha = 0.7,
                height = 0.02,width = 0.2) 

p + labs(y = 'log10(fpkm)', x = '') + mytheme


# v1: manually nitrogen metabolism of yellow light
# theme(axis.text.x = element_text(angle = 45, size = 8))
# 
# Cust_palette <- c('#40A9FE', '#252525', '#FEA940') 
# p <- ggbarplot(data, x="genes", y="values", color = "light",
#                palette = Cust_palette, width = 0.7, alpha = 0.8 ,
#                add = "mean_se", add.params = list(width = 0.3, alpha = 0.7),
#                position = position_dodge(0.8)
#                # add = "jitter", add.params = list()
#                )
# p + labs(y = 'log10(fpkm)', x = '') + mytheme



