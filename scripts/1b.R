library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)
library(reshape2)
library(readxl)
library(cowplot)
library(patchwork)
library(ggsignif)

dt <- fread("data/ms.csv")

K36me <- dt %>% dplyr::filter(cond == "PA")

stats <- dt %>%
  dplyr::filter(cond == "PA") %>%
  group_by(mark) %>%
  summarise(mean.value = mean(value),
            sd.value = sd(value), count = n(),
            se.mean = sd.value/sqrt(count))

ggplot(stats, aes(x=mark, y=mean.value,fill=mark)) + 
  geom_bar(stat = "identity", position = position_dodge(),show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean.value - sd.value, ymax = mean.value + sd.value), 
                width=0.2)+
  geom_jitter(mapping=aes(x=mark,y=value,fill=mark),data = K36me,show.legend = FALSE,size=2) +
  xlab("")+ylab("% of peptides with modification")+
  scale_fill_manual(values=c("forestgreen", "blue","red")) + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y= element_text(size=9,family="Helvetica",colour = "black"),
    axis.text.x = element_text(size=9,family="Helvetica",colour = "black",hjust = 1,vjust=1),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=9,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "blue",size=9,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
      size = 9, color = "black",family = "Helvetica"),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm')) 

ggsave(filename = "1b.pdf",path="figs",device = "pdf",units = "cm",width = 10,height=7,dpi = 600,bg="white")
