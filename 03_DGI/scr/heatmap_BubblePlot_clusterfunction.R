#!/usr/bin/env Rscript
setwd("/Users/yangliao/Documents/GitHub/OSProtein/03_DGI")#设置工作目录

#pdf(file="FunctionsBubbleChart.pdf", width=14, height=11)

# 读取文件 sep 根据文件格式确定
data <- read.table("data/bubble_FAS.txt",header = TRUE, sep = '\t')
library(ggplot2)
library(reshape)
group <- data$group
data_exp <- data.frame(group,log(1+exp(data[,2:length(data)])))
# scaletest <- data.frame(scale(data[2:length(data)]))
# data[2:length(data)] <- scaletest
# summary(scaletest)
data_melt <- melt(data)

# --- heatmap ---
library(pheatmap)
data_heat <- data.frame(t(data_exp[2:length(data)]))
colnames(data_heat) <- data_exp$group
anno_row <- data.frame(c(rep('Blue', 4), rep('Yellow', 4)))
rownames(anno_row) <- rownames(data_heat)
colnames(anno_row) <- 'Light'
# data_heat$light <- c(rep('Blue', 4), rep('Yellow', 4))
pallet1 <- c("#1D39C4", "#2E53EB", "#587DF7", "#85A4FE", "#ACC5FE","#D5E3FE", "#F8EFFE", "#B37EEB")
pallet2 <- c("#1D39C4", "#F8EFFE", "#B37EEB")
pallet3 <- c("#81C684", "#F8EFFE", "#FE7F65")

ann_colors = list(
  Light = c(Blue = "#6598FE", Yellow = "#FECC65")
  # environment = c(X1 = "slateblue3", X2 = "red2"),
  # phylum = c(phylum1 = "#7D26CD", phylum2 = "#E7298A", phylum3 = "#66A61E")
)

p <- pheatmap(data_heat[,1:5], cluster_row = TRUE, cluster_col = TRUE, show_rownames= TRUE,
              annotation_row = anno_row,
              annotation_colors = ann_colors,
              color=colorRampPalette(rev(pallet3))(1000),
              display_numbers = T,
              # color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdBu")))(100),
              # scale = 'row', 
              border_color = 'white',
              angle_col = 315,
              fontsize = 12, # legend_breaks = c(0, 1, 2, 3, 4),
              # legend_labels=c("0","1","2","3","4"),
              cutree_rows = 1)

# --- bubble plot ---
names(data_melt) = c("Functions", "Samples", "Abundances")
data_melt$Sites=substring(data_melt$Samples,1,4)# 根据第二列的样本名称提取站位信息，用于后续着色
data_melt <-as.data.frame(data_melt)
data_melt$Abundances <- log(1+exp(data_melt$Abundances))
# save data
write.csv(data_melt, './output_data/meltdata_buble.csv', row.names = FALSE)
# 做主图
bubble <- ggplot(data_melt[which(data_melt$Abundances>0),], 
                 aes(x = Samples, y = Functions, size = Abundances, color = Sites)) + geom_point(alpha = 0.8)


# 字体修饰
##windowsFonts(myFont = windowsFont("Times New Roman"))

# 修改细节 — 图注，点大小，点shape
bubble_style <- bubble + theme_classic()+
  labs(
    x = "Sampling Sites",
    y = "Functions",
    color="Sites", # 颜色图注名
    size="Abundances")+    # 大小图注名
  scale_size(range = c(0.5, 15), breaks = seq(0, 9, 3)) +  #等比修改圆圈大小
  theme(plot.title=element_text(family="Times New Roman",size=8,
                                color="red",face="italic",
                                hjust=0.5,lineheight=0.5),
        plot.subtitle = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

bubble_style
#dev.off()
