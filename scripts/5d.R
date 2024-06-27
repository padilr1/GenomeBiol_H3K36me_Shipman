library(dplyr)
library(rtracklayer)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(data.table)
library(GenomicRanges)

s <- list.files(path = "data", pattern = '.exon_intron_ratios.',full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("samp",NA,NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path("data", f))
d <- deframe(s[,c('samp', 'f')])
raw.counts <- lapply(d,function(x){fread(x)}) %>% do.call(rbind,.)
exon <- import.bed("ref/exon.bed")
exon.sum <- sum(exon@ranges@width)
intron <- import.bed("ref/intron.bed")
intron.sum <- sum(intron@ranges@width)

exon.counts <- raw.counts %>% dplyr::filter(featureType == "exon")
colnames(exon.counts)[4] <- "exon_counts"
intron.counts <- raw.counts %>% dplyr::filter(featureType == "intron")
colnames(intron.counts)[4] <- "intron_counts"

combined.counts <- dplyr::left_join(exon.counts,intron.counts,by="file")
combined.counts$signal_exon <- (combined.counts$exon_counts/combined.counts$totalReadCount.x*1000000)/exon.sum
combined.counts$signal_intron <- (combined.counts$intron_counts/combined.counts$totalReadCount.y*1000000)/intron.sum
combined.counts$ratio <- combined.counts$signal_exon/combined.counts$signal_intron

ggplot(combined.counts[combined.counts$condition.x == "PA" | combined.counts$condition.x == "SETD2KO" | combined.counts$condition.x == "DKO" | combined.counts$condition.x == "K36M_OE"],aes(x=condition.x,y=ratio)) +
  geom_boxplot(show.legend = FALSE,aes(fill=condition.x)) +
  geom_jitter(mapping=aes(x=condition.x,y=ratio,fill=condition.x),data = combined.counts[combined.counts$condition.x == "PA" | combined.counts$condition.x == "SETD2KO" | combined.counts$condition.x == "DKO" | combined.counts$condition.x == "K36M_OE"],show.legend = FALSE,size=0.7) +
  labs(x="",y="Exon/Intron H3K36me2 signal")  +
  geom_hline(yintercept = 1,alpha=0.1) + 
  scale_fill_manual(values=c("blue","hotpink","grey33","brown")) +
  theme(axis.title.y=element_text(family="Helvetica",size=9,colour="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,family="Helvetica",size=6,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=9,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 8))

ggsave(filename = "3c.png",path="figs",device = "png",dpi = 600,bg="white")