#!/bin/bash

plotEnrichment -b ${ChIP_sample}.sorted.bam \
--BED ref/exon.bed ref/intron.bed \
--regionLabels "exon" "intron" \
-o ${ChIP_sample}.sorted.bam.png \
--blackListFileName ref/blacklist.bed \
--centerReads \
--outRawCounts ${ChIP_sample}.exon_intron_ratios.txt