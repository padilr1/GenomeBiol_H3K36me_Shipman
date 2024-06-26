library(tidyverse)
library(data.table)
library(ggplot2)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(LOLA)
library(rtracklayer)
library(GenomicRanges)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# normalize each sample
## TKO
#need to set wd 
setwd(args$wd)
getwd()
# params
chrom.sizes="mm10chrom.sizes"
# sample info
cell_line="10T"
mark="H3K36me2"
samp="TKO"
input="TKO_input" #watch spaces in names
# library sizes
chipLibrary <- as.numeric(56528519)
inputLibrary <- as.numeric(28685316)
# blacklist 
blacklist=paste0("blacklist.bed")
bl <- import.bed(blacklist)
# list files
s <- list.files(path = "data/10kb.bins", pattern = '.10kb.',full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
  mutate(f = file.path('data/10kb.bins', f)) %>%
  dplyr::filter(line==args$line) %>%
  dplyr::filter(samp==args$samp | samp==args$input)
d <- deframe(s[,c('samp', 'f')])
# read in chrom sizes and chr list
keep <- fread(chrom.sizes) %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
# loop through each file
raw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  out <- dplyr::semi_join(inp,keep,by="chr")
})
pre_bl <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  imd <- dplyr::semi_join(inp,keep,by="chr")
  out <- imd[,1:3] %>% mutate(start = start + 1) 
  out <- out %>% 
    makeGRangesFromDataFrame(seqinfo = gn)
})
bw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  imd <- dplyr::semi_join(inp,keep,by="chr")
  out <- imd[,1:3] %>% mutate(start = start + 1) 
  out <- out %>% 
    makeGRangesFromDataFrame(seqinfo = gn)
  k <- !overlapsAny(out, bl)
  out <- out[!overlapsAny(out, bl)]
})
# normalize by input and scale by a factor
bw[[samp]]$score <- log2(((raw[[samp]]$score/chipLibrary) + 1e-15) / (raw[[input]]$score/inputLibrary + 1e-15))[!overlapsAny(pre_bl[[samp]], bl)]
# drop NA values
bw[[samp]] <- bw[[samp]][!is.na(bw[[samp]]$score)]
# Create the 'norm.bw' directory if it doesn't exist
if (!file.exists("data/norm.bw")) {
  dir.create("data/norm.bw")
}
# export bw
export.bw(bw[[samp]],con=sprintf("data/norm.bw/%s.%s.%s.10kb.norm.bw",cell_line,samp,mark))
## QKO
#need to set wd 
setwd(".")
getwd()
# params
chrom.sizes="mm10chrom.sizes"
# sample info
cell_line="10T"
mark="H3K36me2"
samp="QKO"
input="QKO_input" #watch spaces in names
# library sizes
chipLibrary <- as.numeric(18943359)
inputLibrary <- as.numeric(28641511)
# blacklist 
blacklist=paste0("blacklist.bed")
bl <- import.bed(blacklist)
# list files
s <- list.files(path = "data/10kb.bins", pattern = '.10kb.',full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA),'\\.',F) %>%
  mutate(f = file.path('data/10kb.bins', f)) %>%
  dplyr::filter(line==args$line) %>%
  dplyr::filter(samp==args$samp | samp==args$input)
d <- deframe(s[,c('samp', 'f')])
# read in chrom sizes and chr list
keep <- fread(chrom.sizes) %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
# loop through each file
raw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  out <- dplyr::semi_join(inp,keep,by="chr")
})
pre_bl <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  imd <- dplyr::semi_join(inp,keep,by="chr")
  out <- imd[,1:3] %>% mutate(start = start + 1) 
  out <- out %>% 
    makeGRangesFromDataFrame(seqinfo = gn)
})
bw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  imd <- dplyr::semi_join(inp,keep,by="chr")
  out <- imd[,1:3] %>% mutate(start = start + 1) 
  out <- out %>% 
    makeGRangesFromDataFrame(seqinfo = gn)
  k <- !overlapsAny(out, bl)
  out <- out[!overlapsAny(out, bl)]
})
# normalize by input and scale by a factor
bw[[samp]]$score <- log2(((raw[[samp]]$score/chipLibrary) + 1e-15) / (raw[[input]]$score/inputLibrary + 1e-15))[!overlapsAny(pre_bl[[samp]], bl)]
# drop NA values
bw[[samp]] <- bw[[samp]][!is.na(bw[[samp]]$score)]
# Create the 'norm.bw' directory if it doesn't exist
if (!file.exists("data/norm.bw")) {
  dir.create("data/norm.bw")
}
# export bw
export.bw(bw[[samp]],con=sprintf("data/norm.bw/%s.%s.%s.10kb.norm.bw",cell_line,samp,mark))

# compare the two samples
# start of code
setwd(".")
getwd()
path <- paste0("")
pattern = ".10kb."
control="TKO"
test="QKO"
cell_line="10T"
mark="H3K36me2"
s <- list.files(path = path, pattern = args$pattern,full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c("line","samp","mark",NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path(path, f)) %>%
  dplyr::filter(line == cell_line) %>%
  dplyr::filter(samp == control | samp == test)
odr <- c(test,control)
s <- s %>%
  dplyr::slice(match(odr,samp))
d <- deframe(s[,c('samp', 'f')]) %>%
  lapply(import.bw) 
lapply(d,length)
d[[control]] <- subsetByOverlaps(d[[control]],d[[test]])
d[[test]]<- subsetByOverlaps(d[[test]],d[[control]])
lapply(d,length)
r <- lapply(d, function(y) y[y$score != 0]) %>%
  Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
  granges() 
cell_line <- as.character(s$line[[1]])
mark <- as.character(s$mark[[1]])
# split
split(d, s$line) %>%
  lapply(function(x) {
    o <- lapply(x, function(y) {
      findOverlaps(r, y) %>%
        to() %>%
        {y[.]} %>%
        score()
    }) %>%
      bind_cols() %>%
      `names<-`(c('x', 'y'))
    ok <- o$x > quantile(o$x, .01) &
      o$x < quantile(o$x, .99) &
      o$y > quantile(o$y, .01) &
      o$y < quantile(o$y, .99)
    write_csv(o[ok,], sprintf('data/mat.csv/%s.%s.%s.%s.mat.csv',cell_line,control,test,mark), col_names = F) 
    export.bed(r[ok], sprintf('data/pooled.bed/%s.%s.%s.%s.pooled.bed',cell_line,control,test,mark))
  })

# read in mat.csv and identify 10kb regions enriched in TKO that is absent in QKO
lfc <- fread("data/10T.TKO.QKO.H3K36me2.mat.csv")
bed <- import.bed("10T.TKO.QKO.H3K36me2.pooled.bed")
QKO <- bed
QKO$score <- lfc$V1
QKO_filt <- QKO[QKO$score < 0,]
TKO <- bed
TKO$score <- lfc$V2
NSD3_deposited_K36me2 <- TKO[overlapsAny(TKO,QKO_filt) & TKO$score > 0,]

# run LOLA
#load gene and igr bed
g <- import.bed('ref/gene.bed')
ig <- import.bed('ref/intergenic.bed')
# background
uni <- import.bed("10T.TKO.QKO.H3K36me2.pooled.bed")
uni_igr <- uni[overlapsAny(uni, ig) & !overlapsAny(uni, g)]
uni_g <- uni[overlapsAny(uni, g) & !overlapsAny(uni, ig)]
# user set
qSet <- NSD3_deposited_K36me2
qSet.igr <- qSet[overlapsAny(qSet, ig) & !overlapsAny(qSet, g)]
qSet.g <- qSet[overlapsAny(qSet, g) & !overlapsAny(qSet, ig)]
# load functions
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}
# collection of functional annotations
ensemblDB <- loadRegionDB('regionDB/mm10',collections = "ensembl")

res <- runLOLA(qSet,uni,ensemblDB,cores = 6)

res.final <- res %>% mutate(reg = sub('.bed', '', filename),
                                         sig = as.character(signif.num(qValue)),
                                         q = case_when(qValue > 0 ~ qValue, T ~ min(qValue[qValue > 0])),
                                         q = -log10(q),
                                         info = sprintf('Odds ratio: %.3g <br>FDR: %.3g<br>Overlap: %d<br>|Annotation|: %d', oddsRatio, qValue, support, size),
                                         oddsRatio = case_when(oddsRatio > 0 ~ oddsRatio, T ~ min(oddsRatio[oddsRatio > 0])))

genome_wide_enrichment <- res.final %>%
  dplyr::select(c("oddsRatio","reg","qValue","support","sig")) %>%
  dplyr::filter(support > 1000)
genome_wide_enrichment$region <- "genome_wide_enrichment"

genome_wide_enrichment %>%
  top_n(7, -qValue) %>%
  arrange(oddsRatio) %>%
  mutate(reg = fct_inorder(reg)) %>%
  ggplot(aes(x = oddsRatio, y = reg)) +
  geom_segment(aes(x = oddsRatio, y = reg, xend = oddsRatio, yend = reg,
                   color = support)) +
  geom_point(aes(size = support),color="#ff8c00", stat = "identity") +
  geom_text(aes(label = sig), vjust = -0.4,size=4) +
  scale_size(name = "# overlaps") +
  labs(y = "Region", x = "Odds ratio") +
  scale_color_viridis_c(name = "# overlaps",
                        begin = 0.2, end = 0.8,
                        option = "A", guide = "legend",
                        direction = 1) +
  facet_wrap(. ~ region,scales = "free",labeller = region_labeller) +
  coord_cartesian(clip = "off") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.text = element_text(color = "black",size = 9),
        axis.text.x = element_text(color = "black",size = 9),
        axis.text.y = element_text(color = "black",size = 9),
        axis.title = element_text(color = "black",size = 9),
        panel.grid.major.y = element_line(color = "grey", linetype = "solid"),
        axis.ticks.y = element_line(color = "grey", linetype = "solid"),
        panel.grid.major.x = element_line(color = "grey", linetype = "dashed"),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill ="white"),
        strip.text.x = element_text(color = "white",size = 9),
        strip.background.x = element_rect(fill = "black"),
        legend.position = "none") + 
  scale_size_continuous(range = c(4, 8))
ggsave(dpi = 600,filename = "6d.png",path = "figs",units = "in",width = 6.3,height=3.2)


