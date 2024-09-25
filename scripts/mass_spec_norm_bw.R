library(data.table)
library(rtracklayer)
library(tidyverse)

# import bw
merged_sample <- import.bw("10T.${sample}_merged.cpm.bw")

# mass spec function
ms_norm <- function(bw, MS, exported_file_name,blacklist) {
  bw <- import.bw(bw)
  blacklist <- import.bed(blacklist)
  bw <- bw[!overlapsAny(bw,blacklist)]
  bw <- bw[!is.na(bw$score)]
  N <- as.numeric(length(bw))
  sum_bw <- as.numeric(sum(bw$score))
  MS <- as.numeric(paste0(MS))
  norm_factor <- MS * N / (sum_bw)
  bw$score <- bw$score * norm_factor
  filt_bw <- bw
  filt_bw$score[filt_bw$score > 100] <- as.numeric(100)
  print(paste0("Max value = ", max(filt_bw$score)))
  print(paste0("Post-processing quantiles:", quantile(filt_bw$score)))
  print(paste0("Norm factor = ", norm_factor))
  print(paste0("Signal mean = ", mean(filt_bw$score)))
  export.bw(filt_bw, exported_file_name)
}

# apply function to each merged sample
ms_norm(bw = "merged_sample.bw", MS = genome_wide_percentage_from_mass_spec, exported_file_name = "merged_sample_ms_scaled.bw",blacklist = "mm10_blacklist.bed")
