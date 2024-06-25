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

K36me <- dt

stats <- dt %>%
  group_by(mark) %>%
  summarise(mean.value = mean(value),
            sd.value = sd(value), count = n(),
            se.mean = sd.value/sqrt(count))

ggplot(stats, aes(x=cond, y=mean.value,fill=cond)) + 
  geom_bar(stat = "identity", position = position_dodge(),show.legend = FALSE)+
  geom_jitter(mapping=aes(x=cond,y=value,fill=cond),data = K36me,show.legend = FALSE,size=0.7) +
  geom_errorbar(aes(ymin = mean.value - sd.value, ymax = mean.value + sd.value), 
                width=0.2)+
  xlab("")+ylab("% of peptides with modification")+ 
  facet_wrap(~mark) +
  scale_fill_manual(values=c("blue","hotpink","limegreen","cyan3","purple","yellow3","chocolate4","grey33","khaki4","orange","maroon")) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y= element_text(size=7,family="Helvetica",colour = "black"),
    axis.text.x = element_text(size=7,family="Helvetica",colour = "black",angle=45,hjust = 1,vjust=1),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=7,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "blue",size=9,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
      size = 9, color = "black",family = "Helvetica"),
    legend.text=element_text(size=8,family = "Helvetica",color = "black"),
    legend.title=element_text(size=8,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    legend.position="right",
    legend.justification="right",
    legend.box.spacing = unit(-0.001, "cm")) +
  facet_wrap(~ mark, scales = "free")

ggsave(filename = "2a.pdf",path="figs",device = "pdf",units = "cm",width = 17,height=7,dpi = 600,bg="white")