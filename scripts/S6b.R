library(ggplot2)
library(dplyr)
library(DESeq2)
library(gplots)
library(hexbin)
library(GenomicFeatures)
library(rtracklayer)
library(patchwork)
library(data.table)
library(tidyverse)
library(org.Mm.eg.db)
library(ggpubr)
library(ggrepel)

# load dds object generated from featureCounts and DESeq2 from previous analysis
load("dds.RData")

# normalize raw counts by FPKM and select only the five K36MTs
norm.counts <- DESeq2::fpkm(dds,robust = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("id") %>% dplyr::filter(symbol %in% c("Nsd1","Nsd2","Nsd3","Ash1l","Setd2"))

# pivot longer and add a column grouping replicates into their respective conditions, called "cond"
mat <- norm.counts %>% pivot_longer(cols = 2:29,names_to = "samp",values_to = "FPKM") %>% mutate(cond = gsub("_.*","",.$samp))

# factorize each condition and K36MTs
mat$cond <- factor(x=mat$cond,levels = c("PA","ASH1LKO","SETD2KO","NSD2KO","NSD3KO","DKO","TKO","QKO","QuiKO"))
mat$symbol <- factor(x = mat$symbol,levels = c("Nsd1","Nsd2","Nsd3","Ash1l","Setd2"))

# summarize and generate statistics per each condition
stats <- mat %>% group_by(cond,symbol) %>%
  summarise(mean_FPKM = mean(FPKM),
            sd_FPKM = sd(FPKM), count = n(),
            se_mean = sd_FPKM/sqrt(count))
# plot
ggplot(stats, aes(x=cond, y=mean_FPKM,fill=cond)) + 
  geom_bar(stat = "identity", position = position_dodge(),show.legend = FALSE)+
  geom_jitter(mapping=aes(x=cond,y=FPKM,fill=cond),data = mat,show.legend = FALSE,size=0.7) +
  geom_errorbar(aes(ymin = mean_FPKM - sd_FPKM, ymax = mean_FPKM + sd_FPKM), 
                width=0.2)+
  xlab("")+ylab("FPKM")+
  scale_fill_manual(values=c("blue","limegreen","hotpink","cyan","purple","yellow3","grey33","darkkhaki","orange","maroon")) +
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
  facet_wrap(~ symbol, scales = "free")

# save plot
ggsave(filename = "S6b.pdf",width = 16,height=11,units="cm",dpi = 600,device = "pdf") 