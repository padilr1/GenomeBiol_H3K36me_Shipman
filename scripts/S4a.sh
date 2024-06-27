#!/bin/bash

computeMatrix scale-regions -R ref/igr.20000kb.bed \
-S PA_merged.NSD3unedit.K36me2.input_normalized.log2cpm.bw ASH1LKO.K36me2.input_normalized.log2cpm.merged.final.bw SETD2KO.input_normalized.log2cpm.merged.cpm.final.bw 10T.NSD1KO.K36me2.input_norm.log2cpm.merged.final.bw 10T.NSD2KO.K36me2.input_norm.log2cpm.merged.final.bw NSD3KO.input_normalized.log2cpm.merged.cpm.final.bw NSD12DKO_merged.K36me2.input_normalized.log2cpm.bw H3K36M_OE.K36me2.input_normalized.log2cpm.merged.cpm.bw TKO.K36me2.input_normalized.log2cpm.merged.cpm.final.bw 10T.QKO.K36me2.merged.input_normalized.log2cpm.merged.bw 10T.QuiKO.K36me2.merged.input_normalized.log2cpm..bw \
-bl ref/blacklist.bed \
--binSize 1000 -m 20000 -a 20000 -b 20000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel PA ASH1LKO SETD2KO NSD1KO NSD2KO NSD3KO NSD12DKO H3K36M_OE NSD12-SETD2-TKO NSD123-SETD2-QKO NSD123-SETD2-ASH1L-QuiKO --verbose \
-o data/S4a.mat.gz

plotHeatmap -m data/S4a.mat.gz -o S4a.png --dpi 600 --colorMap Reds
