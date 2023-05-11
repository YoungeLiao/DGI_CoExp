library(ggplot2)
rawdata <- read.csv('./data/Light.csv', header = T)
p <- ggplot(rawdata, aes(x=Wavelength, y=Intensity, group=Group)) +
  geom_line(aes(color=Group), size = 1, alpha = 0.7) +
  # geom_point(aes(color=Group))+
  scale_color_manual(values = c('#85A4FE', '#FED565')) + 
  labs(x = 'Wavelength (nm)', y = 'Intensity (a.u.)', title = 'Light spectrum')

# max(rawdata$Intensity)

# Add max line
p <- p + geom_vline(xintercept = 456.2, linetype = 2, color = 'red', linewidth = 0.7) +
  geom_vline(xintercept = 589.9, linetype = 2 , color = 'red', linewidth = 0.7)

p + mytheme
#theme(legend.position="top")

top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme <- theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.text.x = element_text(size = 14, vjust = 0.5), # vjust = -0.001
        legend.text = element_text(size = 14), 
        # legend.title = element_text(size = 16),
        # axis.text = element_blank(), 
        # axis.ticks = element_blank(), 
        # axis.title = element_blank(), 
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8), # 调整legend位置
        # legend.position = 'none', 
        # legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
