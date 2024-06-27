library(dplyr)
library(rtracklayer)
library(tidyverse)
library(ggplot2)
library(plotly)
library(ggrepel)
library(data.table)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

s <- list.files(path = "~/Documents/10T_K36me2/work/exon.intron.ratio.analysis/data/parental", pattern = '.sorted.',full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("samp",NA,NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path("~/Documents/10T_K36me2/work/exon.intron.ratio.analysis/data/parental", f))

d <- deframe(s[,c('samp', 'f')])

raw.counts <- lapply(d,function(x){fread(x)}) %>% do.call(rbind,.)

raw.counts$mark <- raw.counts$file
raw.counts$mark <- gsub(".*_","",raw.counts$mark)
raw.counts$mark <- gsub(".sorted.bam","",raw.counts$mark)

raw.counts$condition <- "PA"

exon <- import.bed("ref/exon.bed")
exon.sum <- sum(exon@ranges@width)

intron <- import.bed("ref/intron.bed")
intron.sum <- sum(intron@ranges@width)

exon.raw.counts <- raw.counts %>% dplyr::filter(featureType == "exon")
colnames(exon.raw.counts)[4] <- "exon_counts"
intron.raw.counts <- raw.counts %>% dplyr::filter(featureType == "intron")
colnames(intron.raw.counts)[4] <- "intron_counts"

me1_me3.combined.counts <- dplyr::left_join(exon.raw.counts,intron.raw.counts,by="file")
me1_me3.combined.counts$signal_exon <- (me1_me3.combined.counts$exon_counts/me1_me3.combined.counts$totalReadCount.x*1000000)/exon.sum
me1_me3.combined.counts$signal_intron <- (me1_me3.combined.counts$intron_counts/me1_me3.combined.counts$totalReadCount.y*1000000)/intron.sum
me1_me3.combined.counts$ratio <- me1_me3.combined.counts$signal_exon/me1_me3.combined.counts$signal_intron

me1_me3.dt <- me1_me3.combined.counts %>% dplyr::select(c("file","mark.x","ratio"))
colnames(me1_me3.dt) <- c("samp","mark","ratio")

combined.counts$mark <- "K36me2"
me2.dt <- combined.counts %>% dplyr::filter(condition.x == "PA") %>% dplyr::select(c("samp","mark","ratio"))

final.dt <- rbind(me1_me3.dt,me2.dt)

final.dt$mark <- factor(final.dt$mark,levels=c("K36me1","K36me2","K36me3"))

ggplot(final.dt,aes(x=mark,y=ratio)) +
  geom_boxplot(show.legend = FALSE,aes(fill=mark)) +
  geom_jitter(mapping=aes(x=mark,y=ratio,fill=),data = final.dt,show.legend = FALSE,size=0.7) +
  labs(x="",y="")  +
  geom_hline(yintercept = 1,alpha=0.1) + 
  scale_fill_manual(values=c("forestgreen","blue","red")) +
  scale_x_discrete(labels=c("H3K36me1","H3K36me2","H3K36me3")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,family="Helvetica",size=8,colour = "black"),
    axis.title.x = element_text(angle = 45, vjust = 1, hjust=1,family="Helvetica",size=8,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=8,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 8))

ggsave(filename = "S1a.png",path="figs",device = "png",dpi = 600,bg="white")