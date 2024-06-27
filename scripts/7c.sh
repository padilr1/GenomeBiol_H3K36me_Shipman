#!/bin/bash

# run IDR for ATAC-seq replicates
cat QKO_rep1.ATAC.sorted_peaks.narrowPeak QKO_rep2_ATAC.sorted_peaks.narrowPeak | sort -k8,8nr > QKO.pooled_peaks.bed

idr --samples QKO_rep1_ATAC.sorted_peaks.narrowPeak QKO_rep2_ATAC.sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file QKO.idr_peaks.bed \
--idr-threshold 0.05 \
--plot \
--log-output-file QKO.idr_peaks.log \
--peak-list QKO.pooled_peaks.bed

# find ATAC-seq peaks within remaining 119 QKO H3K36me2 peaks
bedtools intersect -wa -F 0.99 -a QKO.idr_peaks.bed -b QKO_merged.SPAN.H3K36me2_peaks.bed > 10T.QKO.ATACseq_peaks_within_K36me2.bed

findMotifsGenome.pl 10T.QKO.ATACseq_peaks_within_K36me2.bed mm10 /data -size given -mask -p 16 -bg QKO.idr_peaks.bed