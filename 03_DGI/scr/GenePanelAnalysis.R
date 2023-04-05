# === load data ===
# --- HGP & SGP ---
# matching for topological landmark genes analysis, blue light's hub gene panel (HGP)

HGP_b <- unique(read.csv('./03_DGI/output_data/GeneListofFilteredPathway_BvsD_HUBCluster4_v2.csv', header = TRUE)[,c(1, 13:21, 30:32)])
SGP_b <- unique(read.csv('./03_DGI/output_data/GeneListofFilteredPathway_BvsD_SignalCluster5_v2.csv', header = TRUE)[,c(1, 13:21, 30:32)])
HGP_y <- unique(read.csv('./03_DGI/output_data/GeneListofFilteredPathway_YvsD_Cluster3.csv', header = TRUE)[,c(1, 13:21, 30:32)])
SGP_y <- unique(read.csv('./03_DGI/output_data/GeneListofFilteredPathway_YvsD_SignalCluster1.csv', header = TRUE)[,c(1, 13:21, 30:32)])

HGP_b$expre <- apply(HGP_b[,2:4], 1, mean)
SGP_b$expre <- apply(SGP_b[,2:4], 1, mean)
HGP_y$expre <- apply(HGP_y[,8:10], 1, mean)
SGP_y$expre <- apply(SGP_y[,8:10], 1, mean)

# --- kegg annotation ---
kegg_anno <- read.csv('./data/kegg_annotation_all.csv')[,c(1, 9, 12)]

# --- Topology ---
library(readxl)
Node_Blue_HGP <- read_excel("03_DGI/data/Data_Topology.xlsx", 
                            sheet = "Blue_HGP")
Node_Blue_SGP <- read_excel("03_DGI/data/Data_Topology.xlsx", 
                            sheet = "Blue_SGP")
Node_Yellow_HGP <- read_excel("03_DGI/data/Data_Topology.xlsx", 
                            sheet = "Yellow_HGP")
Node_Yellow_SGP <- read_excel("03_DGI/data/Data_Topology.xlsx", 
                            sheet = "Yellow_SGP")

# === match data: gene, expression, modularity class ===
module_blue_HGP <- merge(HGP_b, Node_Blue_HGP, by = 'gene_id', all.y = FALSE)[,c(1:10, 13:14, 18)] 
module_blue_SGP <- merge(SGP_b, Node_Blue_SGP, by = 'gene_id', all.y = FALSE) [,c(1:10, 13:14, 18)]
module_yellow_HGP <- merge(HGP_y, Node_Yellow_HGP, by = 'gene_id', all.y = FALSE)[,c(1:10, 13:14, 18)] 
module_yellow_SGP <- merge(SGP_y, Node_Yellow_SGP, by = 'gene_id', all.y = FALSE) [,c(1:10, 13:14, 18)]

# === landmark genes ===
library(dplyr)
GP <- module_blue_HGP
for (i in 0:(length(unique(GP$`Modularity class`))-1)){
  temp_module <- filter(GP, GP$`Modularity class` == i)
  result_ordered <- temp_module[order(temp_module$expre, decreasing = TRUE),]
  if (i == 0){
    landmark_genes <- result_ordered[1:3,]
  }else{
    landmark_i <- result_ordered[1:3,]
    landmark_genes <- rbind(landmark_genes, landmark_i)
  }
  landmark_genes <- landmark_genes[complete.cases(landmark_genes),] # remove NA row
}

gene_panels_list <- list(module_blue_HGP, module_blue_SGP, module_yellow_HGP, module_yellow_SGP)

landmark_list <- list()
for (j in 1:4){
  GP <- gene_panels_list[[j]]
  for (i in 0:(length(unique(GP$`Modularity class`))-1)){
    temp_module <- filter(GP, GP$`Modularity class` == i)
    result_ordered <- temp_module[order(temp_module$expre, decreasing = TRUE),]
    if (i == 0){
      landmark_genes <- result_ordered[1:3,]
    }else{
      landmark_i <- result_ordered[1:3,]
      landmark_genes <- rbind(landmark_genes, landmark_i)
    }
    landmark_genes <- landmark_genes[complete.cases(landmark_genes),] # remove NA row
  }
  landmark_list[[j]] <- landmark_genes
}

# save
names <- c('module_blue_HGP', 'module_blue_SGP', 'module_yellow_HGP', 'module_yellow_SGP')
xlsxfile <- list()
for (i in 1:4){
  xlsxfile[[i]] <- paste('./03_DGI/output_data/', names[i], '.csv', sep = '')
  write.csv(landmark_list[[i]], xlsxfile[[i]], row.names = FALSE)
}

# === obtain pathway table ===
test <- merge(landmark_list[[2]], kegg_anno)

# === visualization ===
# --- preprocessing ---
plotdata_raw <- landmark_list[[3]][, c(11, 2:10, 13)]
library(reshape)
plotdata <- melt(plotdata_raw[,1:10])
# merge_data <- merge(plotdata, plotdata_raw[,c(1, 11)], all = FALSE)
plotdata$`Modularity class` <- rep(plotdata_raw$`Modularity class`, 9)

plotdata$variable <- sub('Blue.', 'Blue', plotdata$variable)
plotdata$variable <- sub('Dark.', 'Dark', plotdata$variable)
plotdata$variable <- sub('Yellow.', 'Yellow', plotdata$variable)

# --- plot ---
# - default parameters -
library(ggplot2)
library(ggpubr)
Cust_palette <- c('#568cc8', '#585858', '#e8c241') 
top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.5
mytheme<-theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 10, angle = 45,  hjust = 1),
        axis.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14),
        # legend.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9), 
        # legend.position = 'right',
        legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
# - plot -
colnames(plotdata) <- c('ko_des', 'variable', 'value', 'Modularity')
p <- ggbarplot(plotdata, x="ko_des", y="value", add = c("mean_se", 'jitter'), color = "variable",
               palette = Cust_palette, position = position_dodge(0.8), width = 0.7,
               add.params = list(width = 0.5, alpha = 0.7),
               )

p + mytheme + labs(y = 'Expression level (FPKM)', x = '', title = 'HGP of yellow light')



