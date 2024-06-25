library(tidyverse)
library(reactable)
library(data.table)
library(rtracklayer)
library(ggsignif)
library(ggpubr)

# sample script - DON'T RUN
bl <- import.bed("ref/blacklist.bed")
rx <- fread("data/K36me2_rx.csv")
b <- fread("10T_PA_K36me2.1kb.bed",col.names = c('chr', 'start', 'end', 'score'))
samp <- "PA_1"
librarySize <- rx$chip[rx$samp == samp]
gr <- makeGRangesFromDataFrame(b,keep.extra.columns = TRUE)
gr$score <- (gr$score/librarySize)*1000000
top1 <- gr[gr$score > (quantile(gr$score,c(0.99)))]
top1_filt <- top1[!overlapsAny(top1,bl)]
p <- data.frame(samp=samp,
                score=mean(top1_filt$score))

# load complete dataset
load("data/H3K36me2_peakiness_scores.RData")

p$condition <- factor(p$condition,levels=c("QuiKO","QKO","TKO","K36M-OE","DKO","NSD3KO","NSD2KO","NSD1KO","SETD2KO","ASH1LKO","PA"))

stats <- aggregate(score ~ condition, p, function(x) c(mean = mean(x), sd = sd(x)))

s <- stats$score %>% as.data.frame() %>% mutate(condition = stats$condition)

ggplot() + stat_summary(mapping=aes(x=condition,y=score,fill=condition),data = p,geom="col",fun=mean,show.legend = FALSE) +
  geom_jitter(mapping=aes(x=condition,y=score,fill=condition),data = p,show.legend = FALSE,size=0.7) +
  geom_errorbar(mapping = aes(y=mean,x=condition,ymin = mean - sd, ymax = mean + sd),data = s, width = 0.2) +
  coord_flip() +
  scale_fill_manual(values=c("maroon3","orange","khaki","chocolate4","grey33","yellow4","purple","cyan","hotpink","greenyellow","blue")) +
  labs(x="",y="") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.y = element_text(family="Helvetica",size=9,colour = "black"),
    axis.text.x = element_text(family="Helvetica",size=9,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=9,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 9),
    plot.margin=unit(c(0,0,0,0), "cm"))

ggsave(filename = "3b.pdf",path="figs",device = "pdf",units = "cm",width = 5.5,height=6,dpi = 600,bg="white")