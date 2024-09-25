# run featureCounts on H3K36me2 signal within exons
```bash
featureCounts -a ref-transcripts.gtf -o 10T.K36me2.within_genes.featureCounts.txt -s 0 -p 10T.PA_1.K36me2.sorted.bam 10T.PA_2.K36me2.sorted.bam 10T.PA_3.K36me2.sorted.bam 10T.SETD2KO_1.K36me2.sorted.bam 10T.SETD2KO_2.K36me2.sorted.bam 10T.SETD2KO_3.K36me2.sorted.bam 10T.H3K36M_OE_1.K36me2.sorted.bam 10T.H3K36M_OE_2.K36me2.sorted.bam 10T.DKO_1.K36me2.sorted.bam 10T.DKO_2.K36me2.sorted.bam
```

# run featureCounts on RNAseq samples from Chao et al.
```bash
featureCounts -a ref-transcripts.gtf -o 10T.RNAseq.featureCounts.txt -s 0 -p 10T.PA_1.RNAseq.sorted.bam 10T.PA_2.RNAseq.sorted.bam 10T.PA_3.RNAseq.sorted.bam 10T.SETD2KO_1.RNAseq.sorted.bam 10T.SETD2KO_2.RNAseq.sorted.bam 10T.SETD2KO_3.RNAseq.sorted.bam 10T.H3K36M_OE_1.RNAseq.sorted.bam 10T.H3K36M_OE_2.RNAseq.sorted.bam 10T.DKO_1.RNAseq.sorted.bam 10T.DKO_2.RNAseq.sorted.bam
```

# compute H3K36me2 signal within exons using R
```R
library(ggplot2)
library(dplyr)
library(tximport)
library(DESeq2)
library(vsn)
library(gplots)
library(hexbin)
library(eulerr)
library(GenomicFeatures)
library(rtracklayer)
library(patchwork)
library(data.table)
library(plotly)
library(ggrepel)
library(stringr)
library(ggpubr)
library(ggforce)
library(grid)
library(gridExtra)

me2.counts <- fread("10T.K36me2.within_genes.featureCounts.txt")
length <- me2.counts %>% dplyr::select(c("Geneid","Length"))
colnames(length) <- c("id","length")
df <- me2.counts %>% dplyr::select(c(1,7:16)) %>% tibble::column_to_rownames("Geneid")
names(df) <- gsub("/project/6007495/shareroot/projects/cell_lines/10Ts/ChIPseq/H3K36me2/|.sorted.bam|10T_|_K36me2|_K36me2_S2.sorted.bam|.*/","",names(df))
mat <- df %>% as.matrix()
metadata <- data.frame(kind=colnames(mat))
metadata$condition <- metadata$kind
metadata$condition <- gsub(pattern = "_c.*|-[1-9].*|_C.*","",metadata$condition)
metadata[3,2] <- "Parental"
metadata$condition <- c("PA","PA","PA","SETD2KO","SETD2KO","SETD2KO","K36M_OE","K36M_OE","DKO","DKO")
metadata <- metadata %>% tibble::column_to_rownames("kind")
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = metadata,
                              design = ~condition)
dds <- DESeq(dds)

norm.counts <- DESeq2::counts(dds,normalized=TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("id") %>% left_join(length,by="id")
norm.counts <- norm.counts %>% mutate(PA_mean=(rowMeans(.[c(2:4)])/.$length)) %>%
  mutate(SETD2KO_mean=(rowMeans(.[c(5:7)])/.$length)) %>%
  mutate(K36M_mean=(rowMeans(.[c(8:9)])/.$length)) %>% mutate(DKO_mean=(rowMeans(.[c(10:11)])/.$length))
final.counts.K36me2 <- norm.counts %>% dplyr::select(c("id","PA_mean","SETD2KO_mean","K36M_mean","DKO_mean"))

# norm counts RNAseq
counts.10T <- fread("10T.RNAseq.featureCounts.txt") %>% tibble::column_to_rownames("id")
counts.10T <- counts.10T[!grepl('symbol', names(counts.10T))]
mat <- counts.10T %>% as.matrix()
mdat <- data.frame(kind=colnames(mat))
mdat$condition <- mdat$kind
mdat$condition <- gsub("_c[1-9]|_[1-9]|_c1[0-9]","",mdat$condition)
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = mdat,
                              design = ~condition)
dds <- DESeq(dds)
norm.counts <- DESeq2::counts(dds,normalized=TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("id") %>% left_join(length,by="id")
norm.counts <- norm.counts %>% mutate(K36M_mean=(rowMeans(.[c(2:4)])/.$length)) %>%
  mutate(DKO_mean=(rowMeans(.[c(5:7)])/.$length)) %>% mutate(PA_mean=(rowMeans(.[c(8,9,13)])/.$length)) %>%
  mutate(SETD2KO_mean=(rowMeans(.[c(10:12)])/.$length))
final.counts.RNAseq <- norm.counts %>% dplyr::select(c("id","PA_mean","SETD2KO_mean","K36M_mean","DKO_mean"))
colnames(final.counts.RNAseq) <- c("id","PA_gene_expression","SETD2KO_gene_expression","K36M_gene_expression","DKO_gene_expression")

# combine RNAseq with ChIP-seq H3K36me2
combined.final <- left_join(final.counts.K36me2,final.counts.RNAseq[,c(1:2)])
colnames(combined.final) <- c("id","PA_K36me2","SETD2KO_K36me2","K36M_K36me2","DKO_K36me2","PA_gene_expression")

PA <- ggplot(combined.final[combined.final$PA_K36me2 != 0,],aes(x=log2(PA_gene_expression),y=log2(PA_K36me2))) +
  labs(x="log2(Gene Expression)",y="log2(H3K36me2 Signal)",title="PA") +
  geom_bin2d(bins=100) +
  scale_fill_continuous(type = "viridis") +
  stat_cor(method="pearson",show.legend = FALSE) +
  geom_smooth(method = "lm",color="red") +
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(family="Helvetica",size=8,colour = "black"),
    axis.title.x = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.title.y = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=8,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 8),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    text = element_text(size = 8))
SETD2KO <- ggplot(combined.final[combined.final$SETD2KO_K36me2 != 0,],aes(x=log2(PA_gene_expression),y=log2(SETD2KO_K36me2))) +
  labs(x="log2(Gene Expression)",y="log2(H3K36me2 Signal)",title="SETD2KO") +
  geom_bin2d(bins=100) +
  scale_fill_continuous(type = "viridis") +
  stat_cor(method="pearson",show.legend = FALSE) +
  geom_smooth(method = "lm",color="red") +
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(family="Helvetica",size=8,colour = "black"),
    axis.title.x = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.title.y = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=8,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 8),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    text = element_text(size = 8))
K36M <- ggplot(combined.final[combined.final$K36M_K36me2 != 0,],aes(x=log2(PA_gene_expression),y=log2(K36M_K36me2))) +
  labs(x="log2(Gene Expression)",y="log2(H3K36me2 Signal)",title="K36M-OE") +
  geom_bin2d(bins=100) +
  scale_fill_continuous(type = "viridis") +
  stat_cor(method="pearson",show.legend = FALSE) +
  geom_smooth(method = "lm",color="red") +
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(family="Helvetica",size=8,colour = "black"),
    axis.title.x = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.title.y = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=8,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 8),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    text = element_text(size = 8))
DKO <- ggplot(combined.final[combined.final$DKO_K36me2 != 0,],aes(x=log2(PA_gene_expression),y=log2(DKO_K36me2))) +
  labs(x="log2(Gene Expression)",y="log2(H3K36me2 Signal)",title="DKO") +
  geom_bin2d(bins=100) +
  scale_fill_continuous(type = "viridis") +
  stat_cor(method="pearson",show.legend = FALSE) +
  geom_smooth(method = "lm",color="red") +
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(family="Helvetica",size=8,colour = "black"),
    axis.title.x = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.title.y = element_text(vjust = 0.5, hjust=0.5,family="Helvetica",size=8,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=8,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 8),
    plot.margin = unit(c(2, 2, 2, 2), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0,'cm'),
    panel.spacing.x = unit(0.1,'cm'),
    text = element_text(size = 8))
ggarrange(PA,SETD2KO,DKO,K36M)
ggsave(filename = "S4d.png",path = "figs",dpi = 600,device = "png") 
```