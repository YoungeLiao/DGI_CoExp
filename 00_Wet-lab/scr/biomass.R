# ===== 1. load data =====
library(readxl)
library(ggpubr)
rawdata <- as.data.frame(read_excel('./data/k_wavelength.xlsx', sheet = 'biomass'))
head(rawdata)
# Time        B1        B2        B3        D1        D2        D3        Y1        Y2        Y3       S1        S2        S3
# 1   24 0.6226006 0.7686614 0.3320137 0.4208221 0.2242275 0.3699153 0.1229206 0.1376439 0.1415068 0.175387 0.3230094 0.4727455

# ===== 2. data tbiomassatment =====
biomass <- as.data.frame(t(rawdata[,-1]))
colnames(biomass) <- 'OD600'
biomass$Group <- c(rep('Blue', 3), rep('Dark', 3), rep('Yellow', 3), rep('Solar', 3))
biomass
# OD600  Group
# B1 0.6226006   Blue
# B2 0.7686614   Blue
# B3 0.3320137   Blue
# D1 0.4208221   Dark
# D2 0.2242275   Dark
# D3 0.3699153   Dark
# Y1 0.1229206 Yellow
# Y2 0.1376439 Yellow
# Y3 0.1415068 Yellow
# S1 0.1753870  Solar
# S2 0.3230094  Solar
# S3 0.4727455  Solar

## t.test
compare_means(OD600 ~ Group, data = biomass, method = 't.test')
# .y.   group1 group2       p p.adj p.format p.signif method
# <chr> <chr>  <chr>    <dbl> <dbl> <chr>    <chr>    <chr> 
#   1 OD600   Blue   Dark   0.0784  0.16  0.0784   ns       T-test
# 2 OD600   Blue   Yellow 0.00478 0.029 0.0048   **       T-test
# 3 OD600   Blue   Solar  0.0380  0.11  0.0380   *        T-test
# 4 OD600   Dark   Yellow 0.0198  0.079 0.0198   *        T-test
# 5 OD600   Dark   Solar  0.810   0.81  0.8105   ns       T-test
# 6 OD600   Yellow Solar  0.0142  0.071 0.0142   *        T-test

Cust_palette <- c('#587DF7', '#585858', '#FED565', '#B37EEB') # blue, dark, yellow, solar
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Blue", "Dark"), c("Yellow", "Dark"), c("Yellow", "Solar"), c('Dark', 'Solar'), c("Blue", "Solar"), c('Blue', 'Yellow'))

ggboxplot(biomass, x="Group", y="OD600", 
          color = "Group", palette = Cust_palette, width = 0.5,
          bxp.errorbar = T, bxp.errorbar.width = 0.4,
          add = "jitter", add.params = list(size = 2.5)) +
  stat_compare_means(comparisons = my_comparisons, method = 't.test', label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.05, vjust = 0.3)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.3, method = 'anova') + # add glocal p-value
  # stat_compabiomass_means(method = "anova", label.y = 60)+ # Add global p-value
  # stat_compabiomass_means(label = "p.format", method = "t.test", label.y = 50) +
  labs(y = 'Biomass (OD600)', x = '') +
  # ylim(30, 120) +
  theme(axis.title.x = element_text(size=16, face = 'bold'), axis.title.y=element_text(size=16, face = 'bold'),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_blank(),
        legend.background = element_blank(),
        # legend.position = 'none',
        legend.position = c(0.8,0.2)
  )

