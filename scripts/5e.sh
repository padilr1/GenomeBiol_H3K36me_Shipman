#!/bin/bash

# promoters
computeMatrix reference-point -R data/inactive.promoter.bed data/active.promoter.bed \
-S tracks/10T.PA_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw tracks/10T.DKO_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw tracks/10T.TKO_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw tracks/10T.QKO_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw \
-bl ref/blacklist.bed \
--referencePoint center \
--binSize 50 -a 3000 -b 3000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel PA NSD12DKO NSD12-SETD2-TKO NSD123-SETD2-QKO \
-o aggr/5e.mat.gz

plotHeatmap -m aggr/5e.mat.gz \
-o figs/5e.pdf \
--dpi 600 --perGroup \
--colorMap "Reds" \
--refPointLabel "TSS" \
-T "" \
--regionsLabel Inactive_promoters Active_promoters

# enhancers
computeMatrix reference-point -R data/inactive.promoter.bed data/active.promoter.bed \
-S tracks/10T.PA_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw tracks/10T.DKO_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw tracks/10T.TKO_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw tracks/10T.QKO_merged.H3K36me2.input_norm.log2cpm.ms_scaled.bw \
-bl ref/blacklist.bed \
--referencePoint center \
--binSize 50 -a 3000 -b 3000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel PA NSD12DKO NSD12-SETD2-TKO NSD123-SETD2-QKO \
-o aggr/5e.mat.gz
