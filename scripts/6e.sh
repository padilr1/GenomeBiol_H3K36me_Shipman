#!/bin/bash

# enhancers
computeMatrix reference-point -R ref/PA.inactive_enhancer.bed ref/PA.active_enhancer.bed \
-S 10T.PA.K36me2.merged.ms_cpm.bw 10T.DKO.K36me2.merged.ms_cpm.bw 10T.TKO.K36me2.merged.ms_cpm.bw 10T.QKO.K36me2.merged.ms_cpm.bw \
-bl ref/blacklist.bed \
--referencePoint center \
--binSize 50 -a 3000 -b 3000 \
--missingDataAsZero --skipZeros \
-p 16 --samplesLabel PA DKO TKO QKO --verbose \
-o data/6e_enhancers.mat.gz

plotHeatmap -m data/6e_enhancers.mat.gz -o 6e_enhancers.png --dpi 600 --colorMap Reds -y "H3K36me2 MS norm enrichment" --regionsLabel "Inactive enhancer" "Active enhancer" --perGroup --refPointLabel ""

# promoters
computeMatrix reference-point -R ref/PA.inactive.promoter.bed ref/PA.active.promoter.bed \
-S 10T.PA.K36me2.merged.ms_cpm.bw 10T.DKO.K36me2.merged.ms_cpm.bw 10T.TKO.K36me2.merged.ms_cpm.bw 10T.QKO.K36me2.merged.ms_cpm.bw \
-bl ref/blacklist.bed \
--referencePoint center \
--binSize 50 -a 3000 -b 3000 \
--missingDataAsZero --skipZeros \
-p 16 --verbose \
-o data/6e_promoters.mat.gz

plotHeatmap -m data/6e_promoters.mat.gz -o 6e_promoters.png --dpi 600 --colorMap Reds -y "H3K36me2 MS norm enrichment" --regionsLabel "Inactive promoter" "Active promoter" --perGroup --refPointLabel ""