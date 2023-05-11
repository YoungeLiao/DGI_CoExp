# ===== 1. load data =====
library(readxl)
library(ggpubr)
setwd('./00_Wet-lab')
rawdata <- as.data.frame(read_excel('./data/k_wavelength.xlsx', sheet = 'Acetate'))
rawdata <- rawdata[,-1]
head(rawdata)

# ===== 2. data tCconsumpatment =====
Cconsump <- as.data.frame(t((1 - rawdata[2,]/rawdata[1,])*100))
colnames(Cconsump) <- 'AcetateUptake'
Cconsump$Group <- c(rep('Blue', 3), rep('Dark', 3), rep('Yellow', 3), rep('Solar', 3))
Cconsump

## t.test
compare_means(AcetateUptake ~ Group, data = Cconsump, method = 't.test')
# .y.           group1 group2       p  p.adj p.format p.signif method
# <chr>         <chr>  <chr>    <dbl>  <dbl> <chr>    <chr>    <chr> 
#   1 AcetateUptake Blue   Dark   0.0754  0.18   0.0754   ns       T-test
# 2 AcetateUptake Blue   Yellow 0.00128 0.0077 0.0013   **       T-test
# 3 AcetateUptake Blue   Solar  0.0426  0.17   0.0426   *        T-test
# 4 AcetateUptake Dark   Yellow 0.0207  0.1    0.0207   *        T-test
# 5 AcetateUptake Dark   Solar  0.328   0.33   0.3281   ns       T-test
# 6 AcetateUptake Yellow Solar  0.0598  0.18   0.0598   ns       T-test

Cust_palette <- c('#587DF7', '#585858', '#FED565', '#B37EEB') # blue, dark, yellow, solar
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Blue", "Dark"), c("Yellow", "Dark"), c("Yellow", "Solar"), c('Dark', 'Solar'), c("Blue", "Solar"), c('Blue', 'Yellow'))

ggboxplot(Cconsump, x="Group", y="AcetateUptake", 
          color = "Group", palette = Cust_palette, width = 0.5,
          bxp.errorbar = T, bxp.errorbar.width = 0.4,
          add = "jitter", add.params = list(size = 2.5)) +
  stat_compare_means(comparisons = my_comparisons, method = 't.test', label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.05, vjust = 0.3)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 20, method = 'anova') + # add glocal p-value
  # stat_compaCconsump_means(method = "anova", label.y = 60)+ # Add global p-value
  # stat_compaCconsump_means(label = "p.format", method = "t.test", label.y = 50) +
  labs(y = 'Acetate consumption (%)', x = '') +
  # ylim(30, 120) +
  theme(axis.title.x = element_text(size=16, face = 'bold'), axis.title.y=element_text(size=16, face = 'bold'),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        # legend.position = 'right',
        legend.position = c(0.8,0.2)
  )
