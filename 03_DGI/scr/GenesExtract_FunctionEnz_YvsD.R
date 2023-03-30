# === 1. load data ===
rawdata <- read.csv('./03_DGI/output_data/Cluster3_YvsD.csv', header = T)
rawdata <- read.csv('./03_DGI/output_data/Cluster1_YvsD.csv', header = T)
# path_filter 
# 
# result_N_filtered <- unique(filter(result_N, result_N$Mean_expre > thre))
# result_nitrate <- matched[grepl("nitrate", matched$SwissProt_Description), ]

# === 2. summary and merge based on SwissProt & filtering ===
library(dplyr)
pathway <- rawdata %>% 
  group_by(level3_pathway_name) %>%
  summarise(Y1mean = mean(Yellow1), 
            Y2mean = mean(Yellow2),
            Y3mean = mean(Yellow3),
            D1mean = mean(Dark1), 
            D2mean = mean(Dark2),
            D3mean = mean(Dark3)) 

# --- p-value ---
temp <- data.frame(t(pathway))
colnames(temp) <- temp[1,]
temp <- temp[-1,]
name <- colnames(temp)
# colnames(temp) <- paste('Pathway', 1:length(temp), sep = '')

library(ggpubr)
temp <- data.frame(apply(temp[,1:(length(temp))], 2, as.numeric))
temp$group <- c(rep('Yellow', 3), rep('Dark', 3))

# initial comparison
compare_data <- temp[, c(1, length(temp))]
colnames(compare_data) <- c('Pathway', 'group')
result <- compare_means(Pathway~group, data = compare_data, ref.group = 'Dark', 
                        method = 't.test', 
                        paired = FALSE)
p_value <- data.frame(t(result$p))
colnames(p_value) <- c('Yellow')
# for loop
for (i in 2:(length(temp)-1)){
  compare_data <- temp[, c(i, length(temp))]
  colnames(compare_data) <- c('Pathway', 'group')
  result <- compare_means(Pathway ~ group, data = compare_data, ref.group = 'Dark', 
                          method = 't.test', 
                          paired = FALSE)
  p <- data.frame(t(result$p))
  colnames(p) <- c('Yellow')
  p_value <- rbind(p_value, p)
}

rownames(p_value) <- name
p_value$pathway <- rownames(p_value)
colnames(p_value) <- c('pvalue_YvsD', 'level3_pathway_name')
kegg_level3 <- merge(pathway, p_value, by = 'level3_pathway_name')

# --- statistic ---
mean <- apply(pathway[,2:7], 1, mean)

sd_Y <- apply(pathway[,2:4], 1, sd)
sd_D <- apply(pathway[,5:7], 1, sd)


mean_Y <- apply(pathway[,2:4], 1, mean)
mean_D <- apply(pathway[,5:7], 1, mean)


kegg_level3$YvsD <- mean_Y/mean_D
kegg_level3$yellow <- mean_Y
kegg_level3$dark <- mean_D

# --- self-defined filtering ---
library(dplyr)
# --- cluster3: threshold: FC, p-value ---
p_thre <- 0.2
FC_up <- 4
FC_down <- 0.3
Expre <- 9
Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) & 
                      kegg_level3$pvalue_YvsD < p_thre &
                      yellow > Expre) 

# --- cluster1: signal ---
p_thre <- 0.15
FC_up <- 2
FC_down <- 0.5
Expre <- 3
Func_YvsD <- filter(kegg_level3, (kegg_level3$YvsD > FC_up | kegg_level3$YvsD < FC_down) &
                      kegg_level3$pvalue_YvsD < p_thre &
                      yellow > Expre)

# --- obtain pathway list & matching & obtain gene list ---
path_list <- Func_YvsD$level3_pathway_name

# initial
path <- path_list[1]
# for loop 
gene_list <- rawdata[grepl(path, rawdata$level3_pathway_name), ] 
for (path in path_list[2:length(path_list)]){
  list <- rawdata[grepl(path, rawdata$level3_pathway_name), ]
  gene_list <- rbind(gene_list, list)
}

# save data
write.csv(gene_list, './03_DGI/output_data/GeneListofFilteredPathway_YvsD_Cluster3.csv', row.names = FALSE)
write.csv(gene_list, './03_DGI/output_data/GeneListofFilteredPathway_YvsD_SignalCluster1.csv', row.names = FALSE)



