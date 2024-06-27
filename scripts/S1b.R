library(tidyverse)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggforce)
regs <- c('up' = '10kb upstream',
          'pmt' = 'Promoter',
          'e1' = 'exon',
          'i1' = 'intron',
          'e2' = 'exon',
          'i2' = 'intron',
          'e3' = 'exon',
          'i3' = 'intron',
          'eAntepenultimate' = 'exon',
          'iPenultimate' = 'intron',
          'ePenultimate' = 'exon',
          'iLast' = 'intron',
          'eLast' = 'exon',
          'dn' = '10kb downstream')

load("PA.K36me1.gene_centric.table.RData")
load("PA.K36me2.gene_centric.table.RData")
load("PA.K36me3.gene_centric.table.RData")
PA.K36me1.gene_centric.table$cond <- "H3K36me1"
PA.K36me2.gene_centric.table$cond <- "H3K36me2"
PA.K36me3.gene_centric.table$cond <- "H3K36me3"

l <- list("H3K36me1"=PA.K36me1.gene_centric.table,"H3K36me2"=PA.K36me2.gene_centric.table,"H3K36me3"=PA.K36me3.gene_centric.table)

K36me2.avg <- lapply(l,function(x){
  x %>%
    dplyr::filter(!(reg %in% c("10kb upstream","Promoter","10kb downstream")))
  x$v <- as.numeric(x$v)
  avg <- x %>%
    group_by(quartile) %>%
    summarize(avg_reg = mean(v))
  avg$quartile <- ifelse(avg$quartile == "zero","0",avg$quartile)
  avg$avg_reg <- as.numeric(avg$avg_reg)
  avg$quartile <- as.numeric(avg$quartile)
  return(avg)
})
names(K36me2.avg) <- names(l)
# Merge the data frames and add a new column with item names
c <- bind_rows(K36me2.avg, .id = "cond")
c$cond <- factor(c$cond,levels=c("H3K36me1","H3K36me2","H3K36me3"))

ggplot(c,aes(x=quartile,y=avg_reg,color=cond))+
         geom_point(show.legend = FALSE) + geom_smooth(method=lm,se=FALSE,fullrange=TRUE,
                  aes(color=cond)) +
  stat_cor(method="pearson",show.legend = FALSE,size=3.5) +
  scale_color_manual(values = c("forestgreen","blue","red"),name="Histone mark",labels=c("H3K36me1","H3K36me2","H3K36me3")) +
  labs(x="",y="") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(family="Helvetica",size=7,colour = "black"),
    axis.title.x = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=7,colour = "black"),
    axis.title.y = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=7,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=7,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 7),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    text = element_text(size = 7))

ggsave(filename = "S1b.png",path="figs",device = "png",dpi = 600,bg="white")