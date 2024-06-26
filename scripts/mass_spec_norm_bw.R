library(data.table)
library(rtracklayer)
library(tidyverse)

# import bw
merged_sample <- import.bw("10T.${sample}_merged.cpm.bw")

# multiply the bigWig scores by the mass-spec value (as a fraction)
merged_sample$score <- merged_sample$score * mass_spec_value

# export mass-spec scaled bigwig
export.bw(PA,"data/10T.${sample}_merged.${histone_mark}.ms_cpm.bw")