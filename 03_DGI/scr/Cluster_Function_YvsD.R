rm(list=ls(all=TRUE))
# === default setting - theme ===
library(ggplot2)
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

# load data
# data1: cluster 4
rawdata <- read.csv('./03_DGI/output_data/Cluster4_YvsD.csv', header = T)
# data2: cluster 6
rawdata <- read.csv('./03_DGI/output_data/Cluster6_YvsD.csv', header = T)
# data3: cluster 3
rawdata <- read.csv('./03_DGI/output_data/Cluster3_YvsD.csv', header = T)
# data4: cluster 1
rawdata <- read.csv('./03_DGI/output_data/Cluster1_YvsD.csv', header = T)

# === summary and merge based on SwissProt ===
# summary(rawdata[, 12:13])
# rawdata$level3_pathway_name
library(dplyr)
# library(tidyverse)
# detach('package:plyr')
## v1
# group_by(rawdata, level3_pathway_name)
pathway <- rawdata %>% 
  group_by(level3_pathway_name) %>%
  summarise(B1mean = mean(Blue1), 
            B2mean = mean(Blue2),
            B3mean = mean(Blue3),
            D1mean = mean(Dark1), 
            D2mean = mean(Dark2),
            D3mean = mean(Dark3),
            Y1mean = mean(Yellow1), 
            Y2mean = mean(Yellow2),
            Y3mean = mean(Yellow3)) 

# === p-value ===
temp <- data.frame(t(pathway))
colnames(temp) <- temp[1,]
temp <- temp[-1,]
name <- colnames(temp)
# colnames(temp) <- paste('Pathway', 1:length(temp), sep = '')

library(ggpubr)
temp <- data.frame(apply(temp[,1:(length(temp))], 2, as.numeric))
temp$group <- c(rep('Blue', 3), rep('Dark', 3), rep('Yellow', 3))

# initial comparison
compare_data <- temp[, c(1, length(temp))]
colnames(compare_data) <- c('Pathway', 'group')
result <- compare_means(Pathway~group, data = compare_data, ref.group = 'Dark', 
                        method = 't.test', 
                        paired = FALSE)
p_value <- data.frame(t(result$p))
colnames(p_value) <- c('Blue', 'Yellow')
# for loop
for (i in 2:(length(temp)-1)){
  compare_data <- temp[, c(i, length(temp))]
  colnames(compare_data) <- c('Pathway', 'group')
  result <- compare_means(Pathway ~ group, data = compare_data, ref.group = 'Dark', 
                          method = 't.test', 
                          paired = FALSE)
  p <- data.frame(t(result$p))
  colnames(p) <- c('Blue', 'Yellow')
  p_value <- rbind(p_value, p)
}

rownames(p_value) <- name
p_value$pathway <- rownames(p_value)
colnames(p_value) <- c('pvalue_BvsD', 'pvalue_YvsD', 'level3_pathway_name')
kegg_level3 <- merge(pathway, p_value, by = 'level3_pathway_name')

# --- statistic ---
mean <- apply(pathway[,2:10], 1, mean)

sd_B <- apply(pathway[,2:4], 1, sd)
sd_D <- apply(pathway[,5:7], 1, sd)
sd_Y <- apply(pathway[,8:10], 1, sd)

mean_B <- apply(pathway[,2:4], 1, mean)
mean_D <- apply(pathway[,5:7], 1, mean)
mean_Y <- apply(pathway[,8:10], 1, mean)

kegg_level3$YvsD <- mean_Y/mean_D
kegg_level3$BvsD <- mean_B/mean_D
kegg_level3$blue <- mean_B
kegg_level3$dark <- mean_D
kegg_level3$yellow <- mean_Y

# ===== Intra-cluster function analysis =====
# === yellow light: self-defined filtering ===
library(dplyr)
# --- cluster4: threshold: FC, p-value ---
p_thre <- 0.1
FC_up <- 2
FC_down <- 0.5
Expre <- 10
Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) & 
                      kegg_level3$pvalue_YvsD < p_thre &
                      yellow > Expre) 
# results: 11 pathways for cluster 4
# 

## v2
# p_thre <- 0.15
# FC_up <- 1.5
# FC_down <- 0.5
# Expre <- 10
# Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) & 
#                       kegg_level3$pvalue_YvsD < p_thre &
#                       blue > Expre) 
## v1
# p_thre <- 0.12
# FC_up <- 2
# FC_down <- 0.5
# Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) & kegg_level3$pvalue_YvsD < p_thre) 
write.csv(Func_YvsD, './03_DGI/output_data/YvsD/FunctionAnaly_cluster4.csv', row.names = FALSE)


# # --- cluster6: threshold: FC, p-value ---
p_thre <- 0.15
FC_up <- 2
FC_down <- 0.5
Expre <- 1.5
Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) &
                      kegg_level3$pvalue_YvsD < p_thre &
                      yellow > Expre)
# 12 for cluster 6
write.csv(Func_YvsD, './03_DGI/output_data/YvsD/FunctionAnaly_cluster6.csv', row.names = FALSE)

# --- cluster3: pvalue, FC ---
p_thre <- 0.2
FC_up <- 2
FC_down <- 0.5
Expre <- 10
Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) &
                      kegg_level3$pvalue_YvsD < p_thre &
                      yellow > Expre)
# 14 for cluster 6
# p_thre <- 0.1
# FC_up <- 4
# FC_down <- 0.3
# Func_Yvsd <- filter(kegg_level3, (kegg_level3$YvsD > 4 | kegg_level3$YvsD < 0.3) & kegg_level3$dark > 0.1) 
write.csv(Func_YvsD, './03_DGI/output_data/YvsD/FunctionAnaly_cluster3.csv', row.names = FALSE)

# --- cluster1: threshold: mean(dark) > 0.1, YvsD > 22 ---
p_thre <- 0.15
FC_up <- 2
FC_down <- 0.5
Expre <- 3
Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) &
                      kegg_level3$pvalue_YvsD < p_thre &
                      yellow > Expre)
# 14 for cluster 1
# FC_up <- 22
# dark_abun_thre <- 0.1
# Func_Yvsd <- filter(kegg_level3, kegg_level3$YvsD > FC_up & kegg_level3$dark > 0.1) 
write.csv(Func_YvsD, './03_DGI/output_data/YvsD/FunctionAnaly_cluster1.csv', row.names = FALSE)

# === Intra-cluster function plot ===
plot_data <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster1.csv', header = TRUE)
data <- plot_data[, c(1, 11:17)]
data$Expression <- apply(plot_data[,8:10], 1, mean)
# data$Expression <- apply(plot_data[,16:17], 1, mean) # mean expression of yellow light dataasets (only include dark and yellow)
head(data)

library(ggrepel)
p <- ggplot(data, aes(x = pvalue_YvsD, y = YvsD, size = Expression, color = yellow)) +
  geom_point(alpha=0.7)
p <- p + geom_text_repel(data = data, 
                    aes(pvalue_YvsD, y = YvsD,label = level3_pathway_name),
                    size=2.5,color="grey20",
                    point.padding = 0.5,hjust = 0.1,
                    segment.color="grey20",
                    segment.size=0.6,
                    segment.alpha=0.8,
                    nudge_x=0,
                    nudge_y=0, 
                    max.overlaps = 10
                    ) + 
  labs(x = 'P-value', y = 'Fold change', title = 'Yellow light (Cluster 1)')

scale_pale1 <- c('#13C2C2',"#D9D9D9",'#D2ACF7','#B37EEB', '#9253DE') #'#006C75', #36CFC8',
scale_pale2 <- c("#9253DE", "#D9D9D9", '#36CFC8', '#13C2C2', '#006C75')
scale_pale3 <- c("#9253DE", "#D9D9D9", '#36CFC8','#13C2C2')
scale_pale4 <- c("#9253DE", "#D9D9D9", '#13C2C2')
# v1
p + scale_size(range = c(1, 12), name="Expression") +
  scale_color_gradientn(colours = scale_pale3) +
  mytheme2

top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme1 <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 10, vjust = 0.5), # vjust = -0.001
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

# === Inter-cluster function plot ===
plotdata_cluster1 <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster1.csv', header = TRUE)
plotdata_cluster3 <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster3.csv', header = TRUE)
plotdata_cluster4 <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster4.csv', header = TRUE)
plotdata_cluster6 <- read.csv('./03_DGI/output_data/YvsD/FunctionAnaly_cluster6.csv', header = TRUE)
plotdata_cluster1$cluster <- 'cluster1'
plotdata_cluster3$cluster <- 'cluster3'
plotdata_cluster4$cluster <- 'cluster4'
plotdata_cluster6$cluster <- 'cluster6'
plotdata <- rbind(plotdata_cluster1, plotdata_cluster3, plotdata_cluster6, plotdata_cluster4)
plotdata$Expression <- apply(plotdata[,15:17], 1, mean)

# [1] "#8DD3C7FF" "#FFFFB3FF" "#BEBADAFF" "#FB8072FF" "#80B1D3FF" "#FDB462FF"
# [7] "#B3DE69FF" "#FCCDE5FF" "#D9D9D9FF" "#BC80BDFF" "#CCEBC5FF" "#FFED6FFF"
# custo_color <- c('grey', '#FF69B4FF', "#FB8072FF","#80B1D3FF", "#B3DE69FF", 
#                        "#40E0D0FF", "#9370DBFF", 
#                        '#1E90FFFF', '#1E90FFFF', "#FFED6FFF")
cluster_pale <- c('#8DD3C7FF','#BEBADAFF',"#FB8072FF", "#FDB462FF")
# pale_10 <- as.vector(paletteer_d('ggsci::default_jco'))
# pale_51 <- as.vector(paletteer_d('ggsci::default_igv'))

library(ggrepel)
p <- ggplot(plotdata,  
            aes(x = pvalue_YvsD, y = YvsD, size = Expression, color = cluster)) +
  geom_point(alpha=0.7)
p <- p + geom_text_repel(data = plotdata, 
                         aes(pvalue_YvsD, y = YvsD,label = level3_pathway_name),
                         size=2.5,color="grey20",
                         point.padding = 0.5,hjust = 0.1,
                         segment.color="grey20",
                         segment.size=0.6,
                         segment.alpha=0.8,
                         nudge_x=0,
                         nudge_y=0, 
                         max.overlaps = 20
) + 
  labs(x = 'P-value', y = 'Fold change', title = 'Yellow light')

# v1
p + scale_size(range = c(1, 15), name="Expression") +
  scale_color_manual(values = cluster_pale) +
  # scale_color_gradientn(colours = c("#9253DE", "#D9D9D9", '#36CFC8', '#13C2C2', '#006C75')) + 
  mytheme1

# data <- plot_data[, c(1, 11:17)]
# data$Expression <- apply(plot_data[,15:17], 1, mean)
# head(data)

# ## debug
# a <- temp[, 336:338]
# compare_means(Pathway337 ~ group, data = a, method = 't.test')
# class(a[,2]) —— 注意数据类型
