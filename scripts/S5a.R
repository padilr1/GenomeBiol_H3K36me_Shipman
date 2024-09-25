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
load("ASH1LKO.K36me2.gene_centric.table.RData")
load("NSD1KO.K36me2.gene_centric.table.RData")
load("NSD2KO.K36me2.gene_centric.table.RData")
load("NSD3KO.K36me2.gene_centric.table.RData")
load("DKO.K36me2.gene_centric.table.RData")
load("K36M.K36me2.gene_centric.table.RData")
load("SETD2KO.K36me2.gene_centric.table.RData")
load("NSD12DKO.K36me2.gene_centric.table.RData")
load("TKO.K36me2.gene_centric.table.RData")
load("QKO.K36me2.gene_centric.table.RData")
load("QuiKO.K36me2.gene_centric.table.RData")

PA.K36me2.gene_centric.table$condition <- "PA"
ASH1LKO.K36me2.gene_centric.table$condition <- "ASH1LKO"
SETD2KO.K36me2.gene_centric.table$condition <- "SETD2KO"
NSD1KO.K36me2.gene_centric.table$condition <- "NSD1KO"
NSD2KO.K36me2.gene_centric.table$condition <- "NSD2KO"
NSD3KO.K36me2.gene_centric.table$condition <- "NSD3KO"
K36M.K36me2.gene_centric.table$condition <- "K36M-OE"
NSD12DKO.K36me2.gene_centric.table$condition <- "DKO"
TKO.K36me2.gene_centric.table$condition <- "TKO"
QKO.K36me2.gene_centric.table$condition <- "QKO"
QuiKO.K36me2.gene_centric.table$condition <- "QuiKO"

l <- list("PA"=PA.K36me2.gene_centric.table,"ASH1LKO"=ASH1LKO.K36me2.gene_centric.table,"SETD2KO"=SETD2KO.K36me2.gene_centric.table,"NSD1KO"=NSD1KO.K36me2.gene_centric.table,"NSD2KO"=NSD2KO.K36me2.gene_centric.table,"NSD3KO"=NSD3KO.K36me2.gene_centric.table,"DKO"=NSD12DKO.K36me2.gene_centric.table,"K36M-OE"=K36M_OE.K36me2.gene_centric.table,"TKO"=TKO.K36me2.gene_centric.table,"QKO"=QKO.K36me2.gene_centric.table,"QuiKO"=QuiKO.K36me2.gene_centric.table)

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

c$cond <- factor(c$cond,levels=c("PA","ASH1LKO","SETD2KO","NSD1KO","NSD2KO","NSD3KO","DKO","K36M-OE","TKO","QKO","QuiKO"))

dt <- fread("data/ms.csv")
K36me2 <- dt %>% dplyr::filter(mark == "H3K36me2")
K36me2$value <- K36me2$value/100
K36me2.stats <- K36me2 %>%
  group_by(cond) %>%
  summarise(mean.value = mean(value),
            sd.value = sd(value), count = n(),
            se.mean = sd.value/sqrt(count))

final <- c %>% left_join(K36me2.stats,by="cond") %>% mutate(ms_norm_avg_reg = .$avg_reg * .$mean.value)

ggplot(c,aes(x=quartile,y=avg_reg,color=cond))+
         geom_point(show.legend = FALSE) + geom_smooth(method=lm,se=FALSE,fullrange=TRUE,
                  aes(color=cond)) +
  stat_cor(method="pearson",show.legend = FALSE,size=2.0) +
  scale_color_manual(values = c("blue", "greenyellow", "pink","cyan","purple","yellow4","grey33","chocolate4","khaki","orange","maroon3"),name="Condition") +
  labs(x="",y="") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=8,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=8,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=8,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 8, color = "black",family = "Helvetica"),
    legend.text=element_text(size=8,family = "Helvetica",color = "black"),
    legend.title=element_text(size=8,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    plot.margin = unit(c(1.2,1.2,1.2,1.2), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    legend.key.size = unit(0.2, "cm"))
ggsave(filename = "S4c.png",path="figs",device = "png",dpi = 600,bg="white")