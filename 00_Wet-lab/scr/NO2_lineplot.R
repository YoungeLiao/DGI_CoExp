# ===== 1. load data =====
library(readxl)
df1 <- as.data.frame(read_excel('./data/fig1_lineplot.xlsx', sheet = 'NO2'))
head(df1)
# Time       B1       B2       B3       D1       D2       D3        Y1        Y2
# 1  0.0 550.0952 539.7226 546.5695 542.5035 550.5926 551.5443 544.98346 549.92562
# 2  4.0 524.3943 571.3776 537.9012 529.6360 533.8537 539.9862 525.06154 526.52868
# 3 13.0 383.6547 399.1525 401.7585 388.1502 381.4256 398.6752 326.18891 293.60073
# 4 15.0 317.7062 333.5342 350.2868 301.7553 291.8656 303.6110 223.91571  63.20087
# 5 17.0 280.5780 262.3684 282.0482 214.8330 200.2808 215.4480 114.22469  76.16548
# 6 18.5 198.8373 219.7860 222.9314 145.7617 119.4144 138.9667  42.08912  38.32639
# Y3       S1       S2       S3
# 1 538.89057 549.3925 553.2416 549.8981
# 2 526.61846 536.7917 539.2375 538.9791
# 3 326.62501 423.9881 433.3446 466.0466
# 4 227.04891 349.4712 370.1861 426.5786
# 5 115.04465 258.9608 275.7209 354.6347
# 6  40.13746 176.9294 212.4994 212.9932

# ===== 2. data tnaratment =====
df2 <- data.frame(rep(c(rep('Blue', 3), rep('Dark', 3), rep('Yellow', 3), rep('Solar', 3)), length(df1$Time)))
colnames(df2) <- 'Group'

TimeSeq <- df1$Time
for (i in TimeSeq){
  if(i == 0){
    time = rep(0, length(df1)-1)
  }
  if(i != 0){
    time <- c(time, rep(i,length(df1)-1))
  }
}
df2$time <- time


for(i in 1:length(df1$Time)){
  if(i == 1){
    nitrite <- as.numeric(df1[i,2:13])
  }
  if(i != 1){
    nitrite <- c(nitrite, as.numeric(df1[i,2:13]))
  }
}
df2$nitrite <- nitrite
head(df2)
# Group time  nitrate
# 1  Blue    0 550.0952
# 2  Blue    0 539.7226
# 3  Blue    0 546.5695
# 4  Dark    0 542.5035
# 5  Dark    0 550.5926
# 6  Dark    0 551.5443

# ===== 3. stat analysis: mean, sd —— error bar =====
attach(df2)
aov.mean<-aggregate(nitrite,by=list(Group, time),FUN=mean) # mean
aov.sd<-aggregate(nitrite,by=list(Group, time),FUN=sd) # sd
detach(df2)

plot_data <- data.frame(aov.mean, sd=aov.sd$x)
col_name <- c('Group', 'Time', 'Nitrite', 'sd')
colnames(plot_data) <- col_name
class(plot_data)
plot_data[,2:4] <- as.data.frame(apply(plot_data[,2:4], 2, as.numeric)) 

p <- ggplot(plot_data, aes(x=Time, y=Nitrite, group = Group)) + 
  geom_line(aes(color=Group)) +
  geom_point(aes(color=Group, shape = Group),size=4, alpha = 1) +
  geom_errorbar(aes(ymin=Nitrite-sd, ymax=Nitrite+sd, color = Group), width=3,
                position=position_dodge(0.05)) +
  # 手动设置line的颜色
  scale_color_manual(values=c('#587DF7', '#585858', '#B37EEB', '#FED565')) +
  scale_y_continuous(breaks = c(0, 100, 200, 300, 400, 500)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25)) +
  ylab('Nitrite (mg NO2-N·L-1)') + 
  xlab('Time (h)') 

mytheme <- theme_classic() + theme(axis.title.x = element_text(size=16, face = 'bold'), axis.title.y=element_text(size=16, face = 'bold'),
                                   axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                                   legend.text = element_text(size = 14), 
                                   legend.title = element_blank(),
                                   legend.position = c(0.15, 0.85)
)

p + mytheme
