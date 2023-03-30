# === 1. load data ===
rawdata <- read.csv('./03_DGI/output_data/Cluster5_BvsD.csv', header = T)

# === 2. summary and merge based on SwissProt & filtering ===
library(dplyr)
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


# --- self-defined filtering ---
library(dplyr)
# --- 20230326_V2: cluster5 - signaling cluster --- 
p_thre <- 0.001 
FC_up <- 10
FC_down <- 0.5
Expre <- 10
Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & 
                      kegg_level3$pvalue_BvsD < p_thre &
                      blue > Expre) 

# --- 20230326_V2: cluster4 - hub cluster--- 
p_thre <- 0.05
FC_up <- 2
FC_down <- 0.5
Expre <- 20
Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) & 
                      kegg_level3$pvalue_BvsD < p_thre &
                      blue >Expre) 

# --- V1: cluster4: threshold: utlized nitrogen metabolism as threshod ---
# p_thre <- 0.025
# FC_up <- 1
# FC_down <- 0.1
# Expre <- 25
# Func_BvsD <- filter(kegg_level3, (kegg_level3$BvsD > FC_up | kegg_level3$BvsD < FC_down) &
#                       kegg_level3$pvalue_BvsD < p_thre &
#                       blue >Expre)

# --- obtain pathway list & matching & obtain gene list ---
path_list <- Func_BvsD$level3_pathway_name

# initial
path <- path_list[1]
# loc <- match(path, rawdata$level3_pathway_name) # 不能用match! match只会匹配到第一个
# loc <- loc[!is.na(loc)]
# gene_list <- rawdata[loc,]
gene_list <- rawdata[grepl(path, rawdata$level3_pathway_name), ] 
for (path in path_list[2:length(path_list)]){
  list <- rawdata[grepl(path, rawdata$level3_pathway_name), ]
  gene_list <- rbind(gene_list, list)
  }

# save data
# Hub cluster
write.csv(gene_list, './03_DGI/output_data/GeneListofFilteredPathway_BvsD_HUBCluster4_v2.csv', row.names = FALSE)
# signaling cluster
write.csv(gene_list, './03_DGI/output_data/GeneListofFilteredPathway_BvsD_SignalCluster5_v2.csv', row.names = FALSE)



