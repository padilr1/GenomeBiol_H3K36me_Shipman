#!/bin/bash
computeMatrix scale-regions -R igr.20000kb.bed \
-S 10T.PA_merged.H3K36me1.input_normalized.log2cpm.bw 10T.PA_merged.H3K36me2.input_normalized.log2cpm.bw 10T.PA_merged.H3K36me3.input_normalized.log2cpm.bw \
-bl ref/blacklist.bed \
--binSize 1000 -m 20000 -a 20000 -b 20000 \
--missingDataAsZero --skipZeros \
-p 16 --samplesLabel H3K36me1 H3K36me2 H3K36me3 --verbose \
-o data/1c.mat.gz

plotHeatmap -m data/1c.mat.gz -o figs/1c.pdf --dpi 600 --colorMap Reds