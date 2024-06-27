#!/bin/bash

# make heatmaps for NSD1-KO, NSD2-KO, DKO, K36M-OE, TKO, QKO and QuiKO at intergenic regions
computeMatrix scale-regions -R ref/igr.20000kb.bed \
-S 10T.PA_merged.H3K36me2.input_normalized.log2cpm.bw 10T.NSD1KO_merged.H3K36me2.input_normalized.log2cpm.bw 10T.NSD2KO_merged.H3K36me2.input_normalized.log2cpm.bw 10T.DKO_merged.H3K36me2.input_normalized.log2cpm.bw 10T.K36M_OE_merged.H3K36me2.input_normalized.log2cpm.bw 10T.TKO_merged.H3K36me2.input_normalized.log2cpm.bw 10T.QKO_merged.H3K36me2.input_normalized.log2cpm.bw 10T.QuiKO_merged.H3K36me2.input_normalized.log2cpm.bw \
-bl ref/blacklist.bed \
--binSize 1000 -m 20000 -a 20000 -b 20000 \
--missingDataAsZero --skipZeros \
-p 16 --samplesLabel PA NSD1KO NSD2KO NSD12DKO H3K36M_OE NSD12-SETD2-TKO NSD123-SETD2-QKO NSD123-SETD2-ASH1L-QuiKO --verbose \
-o data/4b.mat.gz

plotHeatmap -m data/4b.mat.gz -o 4b.png --dpi 600 --colorMap Reds