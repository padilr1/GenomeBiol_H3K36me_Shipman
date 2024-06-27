# run featureCounts for RNAseq wildtype cells for each cell line
```bash
featureCounts -p -s 0 -a ref-transcripts.gtf -o data/10T.RNAseq.featureCounts.txt PA_rep1.bam PA_rep2.bam PA_rep3.bam

featureCounts -p -s 0 -a ref-transcripts.gtf -o data/mESC.RNAseq.featureCounts.txt WT_rep1.bam WT_rep2.bam WT_rep3.bam

featureCounts -p -s 0 -a ref-transcripts.gtf -o data/HNSCC_Cal27.RNAseq.featureCounts.txt WT_rep1.bam WT_rep2.bam WT_rep3.bam
```

# run DESEq2 for each count matrix, then process and plot the results 
```R
### 10T ###
counts <- fread("data/10T.RNAseq.featureCounts.txt")
length <- counts$Length
counts.only <- counts %>% dplyr::select(1,7:9) %>%
  column_to_rownames("Geneid")
colnames(counts.only) <- gsub(".*/","",colnames(counts.only))
colnames(counts.only) <- gsub(".sorted.bam.*","",colnames(counts.only))
cdat <- counts.only
mdat <- data.frame(kind=colnames(cdat))
dds <- DESeqDataSetFromMatrix(countData = cdat, colData = mdat, design = ~ kind)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>%
  as.data.frame()
normalized_counts <- normalized_counts %>% mutate(unadjusted_mean=dplyr::select(.,c(1:3)) %>% rowMeans()) %>% rownames_to_column("ensembl_id")
# normalize by total length of exons combined
norm.counts <- left_join(normalized_counts,length,by="ensembl_id")
norm.counts$mean <- norm.counts$unadjusted_mean / norm.counts$Length
mMSC_norm_counts <- norm.counts %>% dplyr::filter(ensembl_id %in% c("ENSMUSG00000021488","ENSMUSG00000057406","ENSMUSG00000054823"))
### mESC ###
counts <- fread("data/mESC.RNAseq.featureCounts.txt")
length <- counts$Length
counts.only <- counts %>% dplyr::select(1,7:9) %>%
  column_to_rownames("Geneid")
colnames(counts.only) <- gsub(".*/","",colnames(counts.only))
colnames(counts.only) <- gsub(".sorted.bam.*","",colnames(counts.only))
cdat <- counts.only
mdat <- data.frame(kind=colnames(cdat))
dds <- DESeqDataSetFromMatrix(countData = cdat, colData = mdat, design = ~ kind)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>%
  as.data.frame()
normalized_counts <- normalized_counts %>% mutate(unadjusted_mean=dplyr::select(.,c(1:3)) %>% rowMeans()) %>% rownames_to_column("ensembl_id")
# normalize by total length of exons combined
norm.counts <- left_join(normalized_counts,length,by="ensembl_id")
norm.counts$mean <- norm.counts$unadjusted_mean / norm.counts$Length
mESC_norm_counts <- norm.counts %>% dplyr::filter(ensembl_id %in% c("ENSMUSG00000021488","ENSMUSG00000057406","ENSMUSG00000054823"))
### HNSCC ###
counts <- fread("data/HNSCC_Cal27.RNAseq.featureCounts.txt")
length <- counts$Length
counts.only <- counts %>% dplyr::select(1,7:9) %>%
  column_to_rownames("Geneid")
colnames(counts.only) <- gsub(".*/","",colnames(counts.only))
colnames(counts.only) <- gsub(".sorted.bam.*","",colnames(counts.only))
cdat <- counts.only
mdat <- data.frame(kind=colnames(cdat))
dds <- DESeqDataSetFromMatrix(countData = cdat, colData = mdat, design = ~ kind)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>%
  as.data.frame()
normalized_counts <- normalized_counts %>% mutate(unadjusted_mean=dplyr::select(.,c(1:3)) %>% rowMeans()) %>% rownames_to_column("ensembl_id")
# normalize by total length of exons combined
norm.counts <- left_join(normalized_counts,length,by="ensembl_id")
norm.counts$mean <- norm.counts$unadjusted_mean / norm.counts$Length
HNSCC_norm_counts <- norm.counts %>% dplyr::filter(ensembl_id %in% c("ENSG00000165671","ENSG00000109685","ENSG00000147548"))
### aggregate ###
final.combined <- rbind(mMSC_norm_counts,mESC_norm_counts,HNSCC_norm_counts) %>% dplyr::mutate(KMT = c("NSD1","NSD2","NSD3"))
### plot ###
final.combined$KMT <- factor(final.combined$KMT,levels = c("NSD3","NSD2","NSD1"))
final.combined$condition <- factor(final.combined$condition,levels=c("mESC","10T","HNSCC"))
ggplot(data=final.combined,aes(x=condition,y=adjusted_counts,fill=KMT,color=KMT)) +
  geom_bar(position="stack",stat="identity") +
  labs(y="",x="") +
  scale_fill_manual(values=c("darksalmon","steelblue","darkolivegreen")) +
  scale_color_manual(values=c("darksalmon","steelblue","darkolivegreen")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(family="Helvetica",size=8,colour = "black"),
    axis.title.x = element_text(angle=45,family="Helvetica",size=8,colour = "black"),
    axis.title.y = element_text(family="Helvetica",size=8,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=8,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 8))
ggsave(filename = "S5b.png",path = "figs",dpi = 600)
```