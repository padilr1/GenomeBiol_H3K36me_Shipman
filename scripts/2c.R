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

pd$group <- "Multi-KO + K36M-OE"
pd[1:18,33] <- "Single-KO + PA"

pd$group <- factor(pd$group,levels=c("Single-KO + PA","Multi-KO + K36M-OE"))

pd %>%
  ggplot(aes(x = pd[,2], y = pd[,3], color = cond,label=samp,shape=group)) +
  geom_point(size = 2, show.legend = F) +
  geom_hline(yintercept = 0,alpha=0.1) +
  geom_vline(xintercept = 0,alpha=0.1) +
  xlab(label = paste(colnames(pd[2]))) +
  ylab(label = paste(colnames(pd[3]))) +
  labs(title = "",colour="") +
  scale_color_manual(values = c("blue","hotpink","limegreen","cyan3","purple","yellow3","chocolate4","grey33","khaki4","orange","maroon")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=7,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=7,family="Helvetica",colour = "black"),
    axis.text.x=element_text(size=7,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=7,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 7, color = "black",family = "Helvetica"),
    legend.text=element_text(size=7,family = "Helvetica",color = "black"),
    legend.title=element_text(size=7,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    plot.margin = unit(c(2.4,2.4,2.4,2.4), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "2c_H3K36me1.png",path="figs",device = "png",dpi = 600,bg="white")

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
pd$group <- "Multi-KO + K36M-OE"
pd[1:18,35] <- "Single-KO + PA"
pd$group <- factor(pd$group,levels=c("Single-KO + PA","Multi-KO + K36M-OE"))

pd %>%
  ggplot(aes(x = pd[,2], y = pd[,3], color = cond,label=samp,shape=group)) +
  geom_point(size = 2, show.legend = T) +
  geom_hline(yintercept = 0,alpha=0.1) +
  geom_vline(xintercept = 0,alpha=0.1) +
  xlab(label = paste(colnames(pd[2]))) +
  ylab(label = paste(colnames(pd[3]))) +
  labs(title = "",colour="") +
  scale_shape_manual(values=c(16, 17),name="Group",labels=c("Single-KO\n+ PA","Multi-KO\n + K36M-OE")) +
  scale_color_manual(values = c("blue","hotpink","limegreen","cyan3","purple","yellow3","chocolate4","grey33","khaki4","orange","maroon"),name="Condition") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=7,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=7,family="Helvetica",colour = "black"),
    axis.text.x=element_text(size=7,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=7,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 7, color = "black",family = "Helvetica"),
    legend.text=element_text(size=7,family = "Helvetica",color = "black"),
    legend.title=element_text(size=7,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    legend.key.size = unit(0.1, 'cm'),
    plot.margin = unit(c(1.2,1.2,1.2,1.2), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'))

ggsave(filename = "2c_H3K36me2.png",path="figs",device = "png",dpi = 600,bg="white")