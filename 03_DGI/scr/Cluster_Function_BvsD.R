rm(list=ls(all=TRUE))
# === Default beautifying ===
# theme
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
mytheme2 <- theme_bw() +
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
        legend.position = 'right', # 调整legend位置
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

# === Prepare data, including calculating p-value === 
## ！！！！！！！！=== 注意一次只能运行一个！不然之后画点时会重合画不出来！=== ！！！！！！！
# data1: cluster 4
rawdata <- read.csv('./03_DGI/output_data/Cluster4_BvsD.csv', header = T)
# data2: cluster 5
rawdata <- read.csv('./03_DGI/output_data/Cluster5_BvsD.csv', header = T)
# data3: cluster 7
rawdata <- read.csv('./03_DGI/output_data/Cluster7_BvsD.csv', header = T)
# data4: cluster 2
rawdata <- read.csv('./03_DGI/output_data/Cluster2_BvsD.csv', header = T)


# === summary and merge based on SwissProt ===
# summary(rawdata[, 12:13])
# rawdata$level3_pathway_name
library(dplyr)
# detach('package:plyr')
## v1
# pathway <- rawdata %>% 
#   group_by(level3_pathway_name) %>%
#   summarise(B1mean = mean(Blue1), 
#             B2mean = mean(Blue2),
#             B3mean = mean(Blue3),
#             D1mean = mean(Dark1), 
#             D2mean = mean(Dark2),
#             D3mean = mean(Dark3),
#             Y1mean = mean(Yellow1), 
#             Y2mean = mean(Yellow2),
#             Y3mean = mean(Yellow3)) 
pathway <- rawdata %>% 
  group_by(level3_pathway_name) %>%
  summarise(B1mean = mean(Blue1), 
            B2mean = mean(Blue2),
            B3mean = mean(Blue3),
            D1mean = mean(Dark1), 
            D2mean = mean(Dark2),
            D3mean = mean(Dark3)) 

# === p-value ===
temp <- data.frame(t(pathway))
colnames(temp) <- temp[1,]
temp <- temp[-1,]
name <- colnames(temp)
# colnames(temp) <- paste('Pathway', 1:length(temp), sep = '')

library(ggpubr)
# temp <- data.frame(apply(temp[,1:(length(temp))], 2, as.numeric))
# temp$group <- c(rep('Blue', 3), rep('Dark', 3), rep('Yellow', 3))
temp <- data.frame(apply(temp[,1:(length(temp))], 2, as.numeric))
temp$group <- c(rep('Blue', 3), rep('Dark', 3))

# initial comparison
compare_data <- temp[, c(1, length(temp))]
colnames(compare_data) <- c('Pathway', 'group')
result <- compare_means(Pathway~group, data = compare_data, ref.group = 'Dark', 
                        method = 't.test', 
                        paired = FALSE)
p_value <- data.frame(t(result$p))
colnames(p_value) <- c('Blue')
# for loop
for (i in 2:(length(temp)-1)){
  compare_data <- temp[, c(i, length(temp))]
  colnames(compare_data) <- c('Pathway', 'group')
  result <- compare_means(Pathway ~ group, data = compare_data, ref.group = 'Dark', 
                          method = 't.test', 
                          paired = FALSE)
  p <- data.frame(t(result$p))
  colnames(p) <- c('Blue')
  p_value <- rbind(p_value, p)
}

rownames(p_value) <- name
p_value$pathway <- rownames(p_value)
colnames(p_value) <- c('pvalue_BvsD', 'level3_pathway_name')
kegg_level3 <- merge(pathway, p_value, by = 'level3_pathway_name')

# --- statistic ---
mean <- apply(pathway[,2:7], 1, mean)

sd_B <- apply(pathway[,2:4], 1, sd)
sd_D <- apply(pathway[,5:7], 1, sd)
# sd_Y <- apply(pathway[,8:10], 1, sd)

mean_B <- apply(pathway[,2:4], 1, mean)
mean_D <- apply(pathway[,5:7], 1, mean)
# mean_Y <- apply(pathway[,8:10], 1, mean)

# kegg_level3$BvsD <- mean_Y/mean_D
kegg_level3$BvsD <- mean_B/mean_D
kegg_level3$blue <- mean_B
kegg_level3$dark <- mean_D
# kegg_level3$yellow <- mean_Y

# === self-defined filtering ===
library(dplyr)
# --- cluster4: threshold: FC, p-value ---
p_thre <- 0.05
FC_up <- 2
FC_down <- 0.5
Expre <- 20
Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & 
                      kegg_level3$pvalue_BvsD < p_thre &
                      blue >Expre) 
# 26 DEGs in total

# p_thre <- 0.01
# FC_up <- 2
# FC_down <- 0.5
# Expre <- 10
# Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & 
#                       kegg_level3$pvalue_BvsD < p_thre &
#                       blue >Expre) 
# 22 DEGs in total
write.csv(Func_BvsD, './03_DGI/output_data/BvsD/FunctionAnaly_cluster4.csv', row.names = FALSE)

# --- cluster5: pvalue, FC ---
p_thre <- 0.001 
FC_up <- 10
FC_down <- 0.5
Expre <- 10
Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & 
                      kegg_level3$pvalue_BvsD < p_thre &
                      blue > Expre) 
# 15 DEGs in total 
write.csv(Func_BvsD, './03_DGI/output_data/BvsD/FunctionAnaly_cluster5.csv', row.names = FALSE)

# --- cluster7: threshold: mean(dark) > 0.1, BvsD > 22 ---
p_thre <- 0.01
FC_up <- 2
FC_down <- 0.5
Expre <- 10
# Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & kegg_level3$pvalue_BvsD < p_thre) 
Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & 
                      kegg_level3$pvalue_BvsD < p_thre &
                      blue > Expre) 
# 4 DEGs in total 
write.csv(Func_BvsD, './03_DGI/output_data/BvsD/FunctionAnaly_cluster7.csv', row.names = FALSE)

# --- cluster2: threshold: FC, p-value ---
p_thre <- 0.01
FC_up <- 2
FC_down <- 0.5
Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & kegg_level3$pvalue_BvsD < p_thre) 
write.csv(Func_BvsD, './03_DGI/output_data/BvsD/FunctionAnaly_cluster2.csv', row.names = FALSE)



# ===== Inter-cluster function plot =====
plotdata_cluster4 <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster4.csv', header = TRUE)
plotdata_cluster5 <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster5.csv', header = TRUE)
plotdata_cluster7 <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster7.csv', header = TRUE)
# plotdata_cluster2 <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster2.csv', header = TRUE)
plotdata_cluster4$cluster <- 'cluster4'
plotdata_cluster5$cluster <- 'cluster5'
plotdata_cluster7$cluster <- 'cluster7'
# plotdata_cluster2$cluster <- 'cluster2'
plotdata <- rbind(plotdata_cluster4, plotdata_cluster5, plotdata_cluster7)
# plotdata <- rbind(plotdata_cluster4, plotdata_cluster5, plotdata_cluster7, plotdata_cluster2)
plotdata$Expression <- apply(plotdata[,2:4], 1, mean)
# cluster_pale <- c("#FB8072FF","#80B1D3FF", "#B3DE69FF", '#13C2C2')
cluster_pale <- c("#FB8072FF","#80B1D3FF", "#B3DE69FF")


library(ggrepel)
# rm(list=ls(all=TRUE))
# debug start
# plotdata <- plotdata[5:35, ]
# debug end
p <- ggplot(plotdata,
            aes(x = pvalue_BvsD, y = BvsD, size = Expression, color = cluster)) +
  geom_point(alpha=0.7)
p
p <- p + geom_text_repel(data = plotdata, 
                         aes(pvalue_BvsD, y = BvsD,label = level3_pathway_name),
                         size=2.5,color="grey20",
                         point.padding = 0.5,hjust = 0.1,
                         segment.color="grey20",
                         segment.size=0.6,
                         segment.alpha=0.8,
                         nudge_x=0,
                         nudge_y=0, 
                         max.overlaps = 20
) + 
  labs(x = 'P-value', y = 'Fold change', title = 'Blue light')

# v1
library(paletteer)
pale_10 <- as.vector(paletteer_d('ggsci::default_jco'))
pale_51 <- as.vector(paletteer_d('ggsci::default_igv'))
p + scale_size(range = c(1, 15), name="Expression") +
  scale_color_manual(values = cluster_pale) +
  # scale_color_gradientn(colours = c("#9253DE", "#D9D9D9", '#36CFC8', '#13C2C2', '#006C75')) + 
  mytheme1


# ===== Intra-cluster function plot =====
plot_data <- read.csv('./03_DGI/output_data/BvsD/FunctionAnaly_cluster5.csv', header = TRUE)
# data <- plot_data[, c(1, 11:17)]
data <- plot_data
data$Expression <- apply(plot_data[,2:4], 1, mean)
head(data)

library(ggrepel)
p <- ggplot(data, aes(x = pvalue_BvsD, y = BvsD, size = Expression, color = blue)) +
  geom_point(alpha=0.7)

p <- p + geom_text_repel(data = data, 
                    aes(pvalue_BvsD, y = BvsD,label = level3_pathway_name),
                    size=2.5,color="grey20",
                    point.padding = 0.5,hjust = 0.1,
                    segment.color="grey20",
                    segment.size=0.6,
                    segment.alpha=0.8,
                    nudge_x=0,
                    nudge_y=0, 
                    max.overlaps = 10
                    ) + 
  labs(x = 'P-value', y = 'Fold change', title = 'Blue light (Cluster 5)')

scale_pale1 <- c('#13C2C2',"#D9D9D9",'#D2ACF7','#B37EEB', '#9253DE') #'#006C75', #36CFC8',
scale_pale2 <- c("#9253DE", "#D9D9D9", '#36CFC8', '#13C2C2', '#006C75')
scale_pale3 <- c("#9253DE", "#D9D9D9", '#36CFC8','#13C2C2')
scale_pale4 <- c("#9253DE", "#D9D9D9", '#13C2C2')
# scale_pale3 <- c("#9253DE", '#B5F5EB','#36CFC8','#13C2C2')
# v1
p + scale_size(range = c(1, 12), name="Expression") +
  scale_color_gradientn(colours = scale_pale4) + #'#006C75', #36CFC8',
  mytheme2






# data <- plot_data[, c(1, 11:17)]
# data$Expression <- apply(plot_data[,15:17], 1, mean)
# head(data)

# ## debug
# a <- temp[, 336:338]
# compare_means(Pathway337 ~ group, data = a, method = 't.test')
# class(a[,2]) —— 注意数据类型
