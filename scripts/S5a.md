# run featureCounts on PA and QKO
```bash
featureCounts -p -s 0 -a ref-transcripts.gtf -o data/10T.RNAseq.featureCounts.txt PA_rep1.bam PA_rep2.bam PA_rep3.bam QKO_rep1.bam QKO_rep2.bam QKO_rep3.bam 
```

# run DESeq2 and plot TF families' expression
```R
df <- fread("10T.RNAseq.featureCounts.txt")
mat <- df %>% as.matrix()
metadata <- data.frame(kind=colnames(mat))
metadata$condition <- metadata$kind
metadata$condition <- gsub(pattern = "_c.*|_[1-9].*|_C.*","",metadata$condition)
metadata$condition <- factor(metadata$condition,levels=c("PA","QKO"))
metadata <- tibble::column_to_rownames(metadata,"kind")
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = metadata,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds) %>% as.data.frame() %>% tibble::rownames_to_column("gene")

norm.counts <- norm.counts %>% mutate(QKO_mean=(rowMeans(.[c(2:4)])/.$length)) %>%
  mutate(PA_mean=(rowMeans(.[c(5:7)])/.$length))
final.counts <- norm.counts %>% dplyr::select(c("symbol","PA_mean","QKO_mean")) %>% dplyr::filter(grepl("Pbx", symbol) | grepl("Zic",symbol))
final.counts$gene_family <- final.counts$symbol
final.counts$gene_family <- gsub("[1-9]","",final.counts$gene_family)
final.counts[final.counts$gene_family=="Pbxip",4] <- "Pbx"

hox.counts <- norm.counts %>% dplyr::select(c("symbol","PA_mean","QKO_mean")) %>% dplyr::filter(grepl("Hox", symbol))
hox.counts$gene_family <- "Hox"

sox.counts <- norm.counts %>% dplyr::select(c("symbol","PA_mean","QKO_mean")) %>% dplyr::filter(grepl("Sox", symbol))
sox.counts$gene_family <- "Sox"

f <- rbind(final.counts,hox.counts,sox.counts) %>% pivot_longer(c(2:3))
f$gene_family <- factor(f$gene_family,levels=c("Hox","Zic","Pbx","Sox"))
colnames(f) <- c("gene","gene_family","cond","value")
f$cond <- factor(f$cond,levels=c("PA_mean","QKO_mean"))

stats <- f %>%
  group_by(gene_family,cond) %>%
  summarise(median.value = median(value),
            sd.value = sd(value), count = n(),
            se.mean = sd.value/sqrt(count))

custom_labeller <- labeller(
  cond = c("PA_mean" = "PA", "QKO_mean" = "QKO")
)
ggplot(stats, aes(x=gene_family, y=median.value,fill=gene_family)) + 
  geom_bar(stat = "identity", position = position_dodge(),show.legend = FALSE)+
  geom_jitter(mapping=aes(x=gene_family,y=value,fill=gene_family),data = f,show.legend = FALSE,size=3,alpha=0.5) +
  geom_jitter(mapping=aes(x=gene_family,y=value,color=gene),data = subset(f, gene == "Pbx2" | gene == "Zic1" | gene == "Hoxc9"),show.legend = FALSE,size=5,alpha=30) +
  geom_text_repel(data = subset(f, gene == "Pbx2" | gene == "Zic1" | gene == "Hoxc9"),mapping=aes(x=gene_family,y=value,label=gene,color=gene),show.legend=FALSE,nudge_x = 0.3,nudge_y = 0.5,min.segment.length = 30,size=8) +
  geom_errorbar(aes(ymin = median.value - se.mean, ymax = median.value + se.mean), 
                width=0.3)+
  xlab("Gene family")+ylab("Normalized RPK") +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y= element_text(size=24,family="Helvetica",colour = "black"),
    axis.title.x= element_text(size=24,family="Helvetica",colour = "black"),
    axis.text.x = element_text(size=24,family="Helvetica",colour = "black",angle=45,hjust = 1,vjust=1),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=24,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "blue",size=14,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
      size = 24, color = "black",family = "Helvetica"),
    legend.text=element_text(size=24,family = "Helvetica",color = "black"),
    legend.title=element_text(size=24,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    legend.position="right",
    legend.justification="right",
    legend.box.spacing = unit(-0.001, "cm")) +
  facet_wrap(~ cond,labeller = custom_labeller)
ggsave(filename = "S5a.png",path = "figs",dpi = 600,device = "png") 
```