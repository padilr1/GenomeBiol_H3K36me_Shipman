library(data.table)
library(tidyverse)
library(rtracklayer)
library(highcharter)
library(plotly)
library(heatmaply)
library(plotly)
library(uwot)
library(smacof)
library(limma)
library(crosstalk)
library(DT)
library(future)
library(furrr)
library(future.apply)
library(ggrepel)
library(RColorBrewer)

# H3K36me1
s <- list.files(path = "data", pattern = '.10kb.',full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c(NA,"samp",NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path("data", f))
#
d <- deframe(s[,c('samp', 'f')])
# read in chrom sizes and chr list
keep <- fread("ref/mm10chrom.sizes") %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
# loop through each file and keep only signal for each 10kb bin
raw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  out <- dplyr::semi_join(inp,keep,by="chr")
  final <- out$score
})
# read in bin
bs_final <- import.bed("ref/10kb.windows.bed") %>% as.data.frame() %>% dplyr::select(1:3) %>%
  setNames(c('chr','start','end')) %>%
  dplyr::semi_join(keep,by="chr") %>%
  dplyr::filter(chr != "chrM") %>%
  makeGRangesFromDataFrame()
# get max values
mxs_final <- bind_cols(raw) %>% apply(1, max) 

bl <- import.bed('ref/blacklist.bed')
k_final <- mxs_final > 100 & !overlapsAny(bs_final,bl) 

rx <- fread("data/rx.csv")

lfc_final <- names(raw) %>%
  sub('_input', '', .) %>%
  unique() %>%
  setNames(.,.) %>%
  lapply(function(s) {
    r <- rx[rx$samp == s,]
    log2(((raw[[paste0(s) ]] / r$chip ) + 1e-15) /
      (raw[[paste0(s, '_input')]] / r$inp) + 1e-15 )[k_final] %>%
      {.[!is.finite(.)] <- NA_real_ ; .}
  }) %>%
  bind_cols() %>%
  mutate(idx = 1:n()) %>%
  na.omit()

bsk <- bs_final[k_final][lfc_final$idx]
lfc_final$idx <- NULL

rx <- fread("data/rx.csv")

md <- rx[,c(8)] 
md$cond <- md$samp


md <- md %>% column_to_rownames("samp")

res <- prcomp(t(lfc_final), scale. = T, center = T)

codr <- c("PA","SETD2KO","ASH1LKO","NSD1KO","NSD2KO","NSD3KO","K36M-OE","DKO","TKO","QKO","QuiKO")

pd <- summary(res)$importance['Proportion of Variance', ] %>%
  Map(function(x, p) sprintf('%s (%.2g%%)', p, x*100), ., names(.)) %>%
  {tmp <- res$x; colnames(tmp) <- .; tmp} %>%
  as.data.frame() %>%
  rownames_to_column('samp') %>%
  merge(rownames_to_column(md, 'samp')) %>%
  mutate(cond = factor(cond, codr)) %>%
  arrange(cond, samp)

dt <- lfc_final %>% tidyr::pivot_longer(1:29)

dt$cond <- dt$name
dt$cond <- gsub("_1.*","",dt$cond)
dt$cond <- gsub("_2.*","",dt$cond)
dt$cond <- gsub("_3.*","",dt$cond)
dt$cond <- gsub("_4.*","",dt$cond)
dt$cond <- gsub("_5.*","",dt$cond)

dt.f <- dt %>% dplyr::select(c(1,3)) %>% unique()
dt.f <- dt.f %>% column_to_rownames("name")
colnames(dt.f) <- "Condition"
dt.f$Condition <- factor(dt.f$Condition,levels=c("PA","SETD2KO","ASH1LKO","NSD1KO","NSD2KO","NSD3KO","K36M-OE","DKO","TKO","QKO","QuiKO"))

annoCol <- list("Condition"=c('PA'="blue",'SETD2KO'="hotpink",'ASH1LKO'="limegreen",'NSD1KO'="cyan",'NSD2KO'="purple",'NSD3KO'="yellow3",'K36M-OE'="chocolate4",'DKO'="grey33",'TKO'="khaki4",'QKO'="orange",'QuiKO'="maroon"))

colnames(lfc_final) <- c("ASH1LKO_1","ASH1LKO_2","ASH1LKO_3","K36M-OE_1","K36M-OE_2","NSD1KO_1","NSD1KO_2","NSD1KO_3","DKO_1","DKO_2","TKO_1","TKO_2","QuiKO_3","QuiKO_1","QuiKO_2","QKO_1","QKO_2","QKO_3","NSD2KO_1","NSD2KO_2","NSD2KO_3","NSD3KO_1","NSD3KO_2","NSD3KO_3","PA_1","PA_2","PA_3","SETD2KO_3","SETD2KO_1","SETD2KO_3")
col_heatmap <- c("limegreen","limegreen","limegreen","chocolate4","chocolate4","cyan","cyan","cyan","grey33","grey33","khaki4","khaki4","maroon","maroon","maroon","orange","orange","orange","purple","purple","purple","yellow3","yellow3","yellow3","blue","blue","blue","hotpink","hotpink","hotpink")

png(filename = "S3c_H3K36me1.png",res=600)
lfc_final %>% stats::cor(method="pearson") %>% heatmap.2(density.info = "none",key.title = "",key.xlab="",ColSideColors = col_heatmap,RowSideColors = col_heatmap,col="RdYlBu",trace="none",scale="none",cexRow=0.62,cexCol = 0.7,keysize = 2,colsep=c(1:30),sepcolor="grey",sepwidth = c(0.00000000000000000001,0.00000000000000000001),rowsep=c(1:30),offsetRow = -0.01,offsetCol = -0.01)

# H3K36me2
s <- list.files(path = "data", pattern = '.10kb.',full.names = FALSE,recursive = FALSE) %>%
  tibble(f = .) %>%
  separate(f,c(NA,"samp",NA,NA,NA),'\\.',F) %>%
  mutate(f = file.path("data", f))
d <- deframe(s[,c('samp', 'f')])
keep <- fread("ref/mm10chrom.sizes") %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
raw <- lapply(d,function(x){
  inp <- fread(x,col.names = c('chr', 'start', 'end', 'score'))
  out <- dplyr::semi_join(inp,keep,by="chr")
  final <- out$score
})

bs_final <- import.bed("ref/10kb.windows.bed") %>% as.data.frame() %>% dplyr::select(1:3) %>%
  setNames(c('chr','start','end')) %>%
  dplyr::semi_join(keep,by="chr") %>%
  dplyr::filter(chr != "chrM") %>%
  makeGRangesFromDataFrame()

mxs_final <- bind_cols(raw) %>% apply(1, max) 

bl <- import.bed('ref/blacklist.bed')
k_final <- mxs_final > 100 & !overlapsAny(bs_final,bl) 


rx <- fread("data/rx.csv")
lfc_final <- names(raw) %>%
  sub('_input', '', .) %>%
  unique() %>%
  setNames(.,.) %>%
  lapply(function(s) {
    r <- rx[rx$samp == s,]
    log2(((raw[[paste0(s) ]] / r$chip ) + 1e-15) /
      (raw[[paste0(s, '_input')]] / r$inp) + 1e-15 )[k_final] %>%
      {.[!is.finite(.)] <- NA_real_ ; .}
  }) %>%
  bind_cols() %>%
  mutate(idx = 1:n()) %>%
  na.omit()

bsk <- bs_final[k_final][lfc_final$idx]
lfc_final$idx <- NULL

rx <- fread("data/rx.csv")
md <- rx[,c(4:6)] %>% column_to_rownames('samp')
# scaling and centering
res <- prcomp(t(lfc_final), scale. = T, center = T)
codr <- c("PA","SETD2KO","ASH1LKO","NSD1KO","NSD2KO","NSD3KO","K36M-OE","DKO","TKO","QKO","QuiKO")
pd <- summary(res)$importance['Proportion of Variance', ] %>%
  Map(function(x, p) sprintf('%s (%.2g%%)', p, x*100), ., names(.)) %>%
  {tmp <- res$x; colnames(tmp) <- .; tmp} %>%
  as.data.frame() %>%
  rownames_to_column('samp') %>%
  merge(rownames_to_column(md, 'samp')) %>%
  mutate(cond = factor(cond, codr)) %>%
  arrange(cond, line, samp)

dt <- lfc_final %>% tidyr::pivot_longer(1:30)
dt$cond <- dt$name
dt$cond <- gsub("_1.*","",dt$cond)
dt$cond <- gsub("_2.*","",dt$cond)
dt$cond <- gsub("_3.*","",dt$cond)
dt$cond <- gsub("_4.*","",dt$cond)
dt$cond <- gsub("_5.*","",dt$cond)
dt.f <- dt %>% dplyr::select(c(1,3)) %>% unique() 
dt.f <- dt.f %>% column_to_rownames("name")
colnames(dt.f) <- "Condition"
dt.f$Condition <- factor(dt.f$Condition,levels=c("PA","SETD2KO","ASH1LKO","NSD1KO","NSD2KO","NSD3KO","K36M-OE","DKO","TKO","QKO","QuiKO"))

annoCol <- list("Condition"=c('PA'="blue",'SETD2KO'="hotpink",'ASH1LKO'="limegreen",'NSD1KO'="cyan",'NSD2KO'="purple",'NSD3KO'="yellow3",'K36M-OE'="chocolate4",'DKO'="grey33",'TKO'="khaki4",'QKO'="orange",'QuiKO'="maroon"))

col_heatmap <- c("limegreen","limegreen","limegreen","chocolate4","chocolate4","grey33","grey33","cyan","cyan","cyan","purple","purple","purple","blue","yellow3","yellow3","yellow3","blue","blue","orange","orange","orange","maroon","maroon","maroon","hotpink","hotpink","hotpink","khaki4","khaki4","khaki4")

png(filename = "S3c_H3K36me2.png",res=600)

lfc_final %>% stats::cor(method="pearson") %>% heatmap.2(density.info = "none",key.title = "",key.xlab="",ColSideColors = col_heatmap,RowSideColors = col_heatmap,col="RdYlBu",trace="none",scale="none",cexRow=0.62,cexCol = 0.7,keysize = 2,colsep=c(1:30),sepcolor="grey",sepwidth = c(0.00000000000000000001,0.00000000000000000001),rowsep=c(1:30),offsetRow = -0.01,offsetCol = -0.01)