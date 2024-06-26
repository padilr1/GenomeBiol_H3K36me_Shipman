# run featureCounts
```bash
featureCounts -p -s 0 -a ref-transcripts.gtf -o data/10T.PA.RNAseq.featureCounts.txt PA_rep1.bam PA_rep2.bam PA_rep3.bam
```

# compute list of coordinates for functionally annotated regions
```r
# load libraries
library(dplyr)
library(data.table)
library(tidyverse)
library(rtracklayer)
library(DESeq2)
library(biomaRt)
library(broom)

# load and process featureCounts matrix of parental RNA-seq data
counts <- fread("data/10T.PA.RNAseq.featureCounts.txt")
mat <- counts %>% dplyr::select(c(1:4)) %>% tibble::column_to_rownames("Gene_id")
mdat <- data.frame(kind=colnames(mat))
dds <- DESeqDataSetFromMatrix(countData = mat, colData = mdat, design = ~ kind)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>%
  as.data.frame()
normalized_counts <- normalized_counts %>% mutate(mean=dplyr::select(.,c(1:3)) %>% rowMeans()) %>% mutate(length=counts$Length) %>% mutate(normalized_rpk=as.numeric(.$mean)/as.numeric(.$length)) %>% rownames_to_column("Gene_id") %>%
  mutate(gene_symbol=counts$gene_symbol)

# load and process reference transcript from Encode
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ref <- fread("ref/wgEncodeGencodeCompVM25.txt")
colnames(ref) <- c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
ref$width <- abs(ref$txEnd - ref$txStart)
ref.filt <- ref %>%
  dplyr::filter(width > 50000) %>%
  dplyr::filter(exonCount >= 6)

# filter featureCounts data
colnames(normalized_counts)[1] <- "Geneid"
zero_counts <- normalized_counts %>%
  dplyr::filter(mean == 0)
normalized_counts <- normalized_counts %>%
  dplyr::filter(mean != 0)
options(scipen=999)
quantile(normalized_counts$mean)
normalized_counts <- normalized_counts %>%
  dplyr::filter(mean > 0.015477805860)

# merge featureCounts data with reference annotation and process matrix into quantiles
symbols <- fread("mouse_gene_symbols.csv") %>% as.data.frame()
normalized_counts <- left_join(normalized_counts,symbols,by="Geneid")
zero_counts <- left_join(zero_counts,symbols,by="Geneid")
stratified <- normalized_counts %>%
  mutate(quantile = ntile(mean,4))
zero <- ref.filt[ref.filt$name2 %in% zero_counts$symbol]
first <- ref.filt[ref.filt$name2 %in% stratified$symbol[stratified$quantile == 1]]
second <- ref.filt[ref.filt$name2 %in% stratified$symbol[stratified$quantile == 2]]
third <- ref.filt[ref.filt$name2 %in% stratified$symbol[stratified$quantile == 3]]
fourth <- ref.filt[ref.filt$name2 %in% stratified$symbol[stratified$quantile == 4]]

# generate bed files
regs <- c('up' = '10kb upstream',
          'pmt' = 'Promoter',
          'e1' = '1st exon',
          'i1' = '1st intron',
          'e2' = '2nd exon',
          'i2' = '2nd intron',
          'e3' = '3rd exon',
          'i3' = '3rd intron',
          'eAntepenultimate' = '3rd last exon',
          'iPenultimate' = '2nd last intron',
          'ePenultimate' = '2nd last exon',
          'iLast' = 'last intron',
          'eLast' = 'last exon',
          'dn' = '10kb downstream')
#
zero %>%
  split(., .$strand) %>%
  lapply(function(x) {
    r <- split(x, x$name) %>%
      lapply(function(y) {
        o <- list()
        ei <- strsplit(y$exonStarts, ',')[[1]] %>% as.numeric()
        ej <- strsplit(y$exonEnds, ',')[[1]] %>% as.numeric()
        exs <- rbind(ei, ej)
        ins <- rbind(ej[-length(ej)], ei[-1])
        if (x$strand[1] == '-') {
          o$pmt <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e3)
          o$up <- tibble(chr = y$chrom, start = y$txEnd + 1e3, end = start + 1e4)
          o$dn <- tibble(chr = y$chrom, start = y$txStart - 1e4, end = y$txStart)
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
        } else {
          o$pmt <- tibble(chr = y$chrom, start = y$txStart - 1e3, end = y$txStart)
          o$up <- tibble(chr = y$chrom, start = y$txStart - 1.1e4, end = y$txStart - 1e3)
          o$dn <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e4)
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate  <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
        }
        o
      })
    lapply(names(r[[1]]), function(y) {
      lapply(r, `[[`, y) %>%
        bind_rows() %>%
        fwrite(sprintf('data/%s/%s.zero.bed', ifelse(x$strand[1] == '-', 'minus', 'plus'), y),
               col.names = F, row.names = F, scipen = 999, sep = '\t')
    })
  })
#
first %>%
  split(., .$strand) %>%
  lapply(function(x) {
    r <- split(x, x$name) %>%
      lapply(function(y) {
        o <- list()
        ei <- strsplit(y$exonStarts, ',')[[1]] %>% as.numeric()
        ej <- strsplit(y$exonEnds, ',')[[1]] %>% as.numeric()
        exs <- rbind(ei, ej)
        ins <- rbind(ej[-length(ej)], ei[-1])
        if (x$strand[1] == '-') {
          o$pmt <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e3)
          o$up <- tibble(chr = y$chrom, start = y$txEnd + 1e3, end = start + 1e4)
          o$dn <- tibble(chr = y$chrom, start = y$txStart - 1e4, end = y$txStart)
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
        } else {
          o$pmt <- tibble(chr = y$chrom, start = y$txStart - 1e3, end = y$txStart)
          o$up <- tibble(chr = y$chrom, start = y$txStart - 1.1e4, end = y$txStart - 1e3)
          o$dn <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e4)
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate  <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
        }
        o
      })
    lapply(names(r[[1]]), function(y) {
      lapply(r, `[[`, y) %>%
        bind_rows() %>%
        fwrite(sprintf('data/%s/%s.first.bed', ifelse(x$strand[1] == '-', 'minus', 'plus'), y),
               col.names = F, row.names = F, scipen = 999, sep = '\t')
    })
  })
#
second %>%
  split(., .$strand) %>%
  lapply(function(x) {
    r <- split(x, x$name) %>%
      lapply(function(y) {
        o <- list()
        ei <- strsplit(y$exonStarts, ',')[[1]] %>% as.numeric()
        ej <- strsplit(y$exonEnds, ',')[[1]] %>% as.numeric()
        exs <- rbind(ei, ej)
        ins <- rbind(ej[-length(ej)], ei[-1])
        if (x$strand[1] == '-') {
          o$pmt <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e3)
          o$up <- tibble(chr = y$chrom, start = y$txEnd + 1e3, end = start + 1e4)
          o$dn <- tibble(chr = y$chrom, start = y$txStart - 1e4, end = y$txStart)
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
        } else {
          o$pmt <- tibble(chr = y$chrom, start = y$txStart - 1e3, end = y$txStart)
          o$up <- tibble(chr = y$chrom, start = y$txStart - 1.1e4, end = y$txStart - 1e3)
          o$dn <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e4)
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate  <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
        }
        o
      })
    lapply(names(r[[1]]), function(y) {
      lapply(r, `[[`, y) %>%
        bind_rows() %>%
        fwrite(sprintf('data/%s/%s.second.bed', ifelse(x$strand[1] == '-', 'minus', 'plus'), y),
               col.names = F, row.names = F, scipen = 999, sep = '\t')
    })
  })
#
third %>%
  split(., .$strand) %>%
  lapply(function(x) {
    r <- split(x, x$name) %>%
      lapply(function(y) {
        o <- list()
        ei <- strsplit(y$exonStarts, ',')[[1]] %>% as.numeric()
        ej <- strsplit(y$exonEnds, ',')[[1]] %>% as.numeric()
        exs <- rbind(ei, ej)
        ins <- rbind(ej[-length(ej)], ei[-1])
        if (x$strand[1] == '-') {
          o$pmt <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e3)
          o$up <- tibble(chr = y$chrom, start = y$txEnd + 1e3, end = start + 1e4)
          o$dn <- tibble(chr = y$chrom, start = y$txStart - 1e4, end = y$txStart)
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
        } else {
          o$pmt <- tibble(chr = y$chrom, start = y$txStart - 1e3, end = y$txStart)
          o$up <- tibble(chr = y$chrom, start = y$txStart - 1.1e4, end = y$txStart - 1e3)
          o$dn <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e4)
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate  <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
        }
        o
      })
    lapply(names(r[[1]]), function(y) {
      lapply(r, `[[`, y) %>%
        bind_rows() %>%
        fwrite(sprintf('data/%s/%s.third.bed', ifelse(x$strand[1] == '-', 'minus', 'plus'), y),
               col.names = F, row.names = F, scipen = 999, sep = '\t')
    })
  })
#
fourth %>%
  split(., .$strand) %>%
  lapply(function(x) {
    r <- split(x, x$name) %>%
      lapply(function(y) {
        o <- list()
        ei <- strsplit(y$exonStarts, ',')[[1]] %>% as.numeric()
        ej <- strsplit(y$exonEnds, ',')[[1]] %>% as.numeric()
        exs <- rbind(ei, ej)
        ins <- rbind(ej[-length(ej)], ei[-1])
        if (x$strand[1] == '-') {
          o$pmt <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e3)
          o$up <- tibble(chr = y$chrom, start = y$txEnd + 1e3, end = start + 1e4)
          o$dn <- tibble(chr = y$chrom, start = y$txStart - 1e4, end = y$txStart)
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
        } else {
          o$pmt <- tibble(chr = y$chrom, start = y$txStart - 1e3, end = y$txStart)
          o$up <- tibble(chr = y$chrom, start = y$txStart - 1.1e4, end = y$txStart - 1e3)
          o$dn <- tibble(chr = y$chrom, start = y$txEnd, end = start + 1e4)
          o$e1 <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$e2 <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$e3 <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$i1 <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$i2 <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
          o$i3 <- tibble(chr = y$chrom, start = ins[1,3], end = ins[2,3])
          exs <- exs[,ncol(exs):1]
          ins <- ins[,ncol(ins):1]
          o$eLast <- tibble(chr = y$chrom, start = exs[1,1], end = exs[2,1])
          o$ePenultimate <- tibble(chr = y$chrom, start = exs[1,2], end = exs[2,2])
          o$eAntepenultimate <- tibble(chr = y$chrom, start = exs[1,3], end = exs[2,3])
          o$iLast <- tibble(chr = y$chrom, start = ins[1,1], end = ins[2,1])
          o$iPenultimate  <- tibble(chr = y$chrom, start = ins[1,2], end = ins[2,2])
        }
        o
      })
    lapply(names(r[[1]]), function(y) {
      lapply(r, `[[`, y) %>%
        bind_rows() %>%
        fwrite(sprintf('data/%s/%s.fourth.bed', ifelse(x$strand[1] == '-', 'minus', 'plus'), y),
               col.names = F, row.names = F, scipen = 999, sep = '\t')
    })
  })
```

# use computeMatrix to count average signal in each annotated region
```bash
# process signal for each annotated bed file using deepTools computeMatrix
sample_output_dir="samp_outdir"
sample_bw="samp.cpm.bw"
sample_name="samp_name"
cd ${sample_output_dir}
parallel --jobs 6 "computeMatrix scale-regions -R {} \
-S ${sample_bw} \
-bl blacklist.bed \
-b 0 -a 0 -m 1000 -bs 5 \
-p 6 --samplesLabel ${sample_name} --verbose \
-o ${sample_name}.plus.{1/.}.mat.gz" ::: gene_body_bedfiles/plus/*.bed
# minus
parallel --jobs 6 "computeMatrix scale-regions -R {} \
-S ${sample_bw} \
-bl blacklist.bed \
-b 0 -a 0 -m 1000 -bs 5 \
-p 6 --samplesLabel ${sample_name} --verbose \
-o ${sample_name}.lfc.minus.{1/.}.mat.gz" ::: gene_body_bedfiles/minus/*.bed

# make sub-directories for each sample
for dir in samples/*; do (cd "$dir" && mkdir -p zero first second third fourth); done
# move matrix files to each subdirectory labeled zero,first,second,third,fourth
for dir in samples/*; 
do (cd "$dir" && for quartile in zero first second third fourth; do mv *.${quartile}.mat.gz ${quartile}; done); 
done
```
# aggregate the matrices and write into an R table object
```r
# aggregate matrices 
library(tidyverse)
library(ggplot2)
library(data.table)
library(ggpubr)
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
# working directories
zero_dir="samp/zero"
first_dir="samp/first"
second_dir="samp/second"
third_dir="samp/third"
fourth_dir="samp/fourth"
setwd(zero_dir)
zero <- list.files(pattern = 'mat.gz') %>%
  tibble(file = .) %>%
  separate(file, c('samp', 'norm', 'str', 'reg', NA, NA), '\\.', F) %>%
  split(., .$reg) %>%
  lapply(function(x) {
    split(x, x$str) %>%
      lapply(function(y) {
        m <- fread(y$file, skip = 1, drop = 1:6)
        if (y$str == 'minus') {
          m <- m[,ncol(m):1]
        }
        colnames(m) <- paste0('V', 1:ncol(m))
        m
      }) %>%
      bind_rows() %>%
      colMeans(na.rm = T) %>%
      tibble(v = .) %>%
      mutate(idx = 1:n())
  }) %>%
  bind_rows(.id = 'reg')
setwd(first_dir)
first <- list.files(pattern = 'mat.gz') %>%
  tibble(file = .) %>%
  separate(file, c('samp', 'norm', 'str', 'reg', NA, NA), '\\.', F) %>%
  split(., .$reg) %>%
  lapply(function(x) {
    split(x, x$str) %>%
      lapply(function(y) {
        m <- fread(y$file, skip = 1, drop = 1:6)
        if (y$str == 'minus') {
          m <- m[,ncol(m):1]
        }
        colnames(m) <- paste0('V', 1:ncol(m))
        m
      }) %>%
      bind_rows() %>%
      colMeans(na.rm = T) %>%
      tibble(v = .) %>%
      mutate(idx = 1:n())
  }) %>%
  bind_rows(.id = 'reg')
setwd(second_dir)
second <- list.files(pattern = 'mat.gz') %>%
  tibble(file = .) %>%
  separate(file, c('samp', 'norm', 'str', 'reg', NA, NA), '\\.', F) %>%
  split(., .$reg) %>%
  lapply(function(x) {
    split(x, x$str) %>%
      lapply(function(y) {
        m <- fread(y$file, skip = 1, drop = 1:6)
        if (y$str == 'minus') {
          m <- m[,ncol(m):1]
        }
        colnames(m) <- paste0('V', 1:ncol(m))
        m
      }) %>%
      bind_rows() %>%
      colMeans(na.rm = T) %>%
      tibble(v = .) %>%
      mutate(idx = 1:n())
  }) %>%
  bind_rows(.id = 'reg')
setwd(third_dir)
third <- list.files(pattern = 'mat.gz') %>%
  tibble(file = .) %>%
  separate(file, c('samp', 'norm', 'str', 'reg', NA, NA), '\\.', F) %>%
  split(., .$reg) %>%
  lapply(function(x) {
    split(x, x$str) %>%
      lapply(function(y) {
        m <- fread(y$file, skip = 1, drop = 1:6)
        if (y$str == 'minus') {
          m <- m[,ncol(m):1]
        }
        colnames(m) <- paste0('V', 1:ncol(m))
        m
      }) %>%
      bind_rows() %>%
      colMeans(na.rm = T) %>%
      tibble(v = .) %>%
      mutate(idx = 1:n())
  }) %>%
  bind_rows(.id = 'reg')
setwd(fourth_dir)
fourth <- list.files(pattern = 'mat.gz') %>%
  tibble(file = .) %>%
  separate(file, c('samp', 'norm', 'str', 'reg', NA, NA), '\\.', F) %>%
  split(., .$reg) %>%
  lapply(function(x) {
    split(x, x$str) %>%
      lapply(function(y) {
        m <- fread(y$file, skip = 1, drop = 1:6)
        if (y$str == 'minus') {
          m <- m[,ncol(m):1]
        }
        colnames(m) <- paste0('V', 1:ncol(m))
        m
      }) %>%
      bind_rows() %>%
      colMeans(na.rm = T) %>%
      tibble(v = .) %>%
      mutate(idx = 1:n())
  }) %>%
  bind_rows(.id = 'reg')
zero$quartile <- "zero"
first$quartile <- "1"
second$quartile <- "2"
third$quartile <- "3"
fourth$quartile <- "4"
zero <- zero %>%
  mutate(reg = factor(reg, names(regs))) %>%
  arrange(reg, idx) %>%
  mutate(x = 1:n(),
         reg = fct_inorder(regs[as.character(reg)]))
first <- first %>%
  mutate(reg = factor(reg, names(regs))) %>%
  arrange(reg, idx) %>%
  mutate(x = 1:n(),
         reg = fct_inorder(regs[as.character(reg)]))
second <- second %>%
  mutate(reg = factor(reg, names(regs))) %>%
  arrange(reg, idx) %>%
  mutate(x = 1:n(),
         reg = fct_inorder(regs[as.character(reg)]))
third <- third %>%
  mutate(reg = factor(reg, names(regs))) %>%
  arrange(reg, idx) %>%
  mutate(x = 1:n(),
         reg = fct_inorder(regs[as.character(reg)]))
fourth <- fourth %>%
  mutate(reg = factor(reg, names(regs))) %>%
  arrange(reg, idx) %>%
  mutate(x = 1:n(),
         reg = fct_inorder(regs[as.character(reg)]))
filt <- broom::tidy(quantile(zero$v,probs=c(0.98)))
zero <- zero %>%
  dplyr::filter(v < filt$x)
samp.gene_centric.table <- rbind(zero,first,second,third,fourth)
samp.gene_centric.table$samp <- "samp_name"
save(samp.gene_centric.table,file="samples/tables/samp.gene_centric.table.RData")
```