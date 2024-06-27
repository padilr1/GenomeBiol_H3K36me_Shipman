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

load("PA.K36me2.gene_centric.table.RData")
load("K36M.K36me2.gene_centric.table.RData")
load("SETD2KO.K36me2.gene_centric.table.RData")
load("NSD12DKO.K36me2.gene_centric.table.RData")

PA.K36me2.gene_centric.table$condition <- "PA"
SETD2KO.K36me2.gene_centric.table$condition <- "SETD2KO"
K36M.K36me2.gene_centric.table$condition <- "K36M-OE"
NSD12DKO.K36me2.gene_centric.table$condition <- "DKO"

l <- list("PA"=PA.K36me2.gene_centric.table,"SETD2KO"=SETD2KO.K36me2.gene_centric.table,"DKO"=NSD12DKO.K36me2.gene_centric.table,"K36M-OE"=K36M.K36me2.gene_centric.table)

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

c <- bind_rows(K36me2.avg, .id = "cond")
c <- c %>% dplyr::filter((cond %in% c("PA","SETD2KO","DKO","K36M-OE")))
c$cond <- factor(c$cond,levels=c("PA","SETD2KO","DKO","K36M-OE"))

dt <- fread("data/ms.csv")
K36me2 <- dt %>% dplyr::filter(mark == "H3K36me2")
K36me2$value <- K36me2$value/100
K36me2.stats <- K36me2 %>%
  group_by(cond) %>%
  summarise(mean.value = mean(value),
            sd.value = sd(value), count = n(),
            se.mean = sd.value/sqrt(count))

final <- c %>% left_join(K36me2.stats,by="cond") %>% mutate(ms_norm_avg_reg = .$avg_reg * .$mean.value)

# get slopes
summary(lm(ms_norm_avg_reg ~ quartile, data = final[final$cond=="PA",])) # 0.052612
summary(lm(ms_norm_avg_reg ~ quartile, data = final[final$cond=="SETD2KO",]))# 0.043533
summary(lm(ms_norm_avg_reg ~ quartile, data = final[final$cond=="DKO",])) #  0.04389
summary(lm(ms_norm_avg_reg ~ quartile, data = final[final$cond=="K36M-OE",])) #  0.027430

ggplot(final,aes(x=quartile,y=ms_norm_avg_reg,color=cond))+
         geom_point(show.legend = FALSE) + geom_smooth(method=lm,se=FALSE,fullrange=TRUE,
                  aes(color=cond)) +
  stat_cor(method="pearson",show.legend = FALSE,size=3.5) +
  scale_color_manual(values = c("blue","hotpink","grey33","chocolate4"),name="Condition",labels=c("PA","SETD2-KO","NSD12-DKO","K36M-OE")) +
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
ggsave(filename = "5c.png",path="figs",device = "png",dpi = 600,bg="white")