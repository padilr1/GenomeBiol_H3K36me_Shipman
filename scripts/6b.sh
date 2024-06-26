#!/bin/bash

make_tracks_file --trackFiles 10T.TKO.K36me2.merged.ms_cpm.bw 10T.TKO.ATAC.merged.cpm.bw TKO_H3K27ac_1.cpm.bw enhancer.bed mm10.ncbiRefSeq.gtf -o data/6b.ini

pyGenomeTracks --tracks data/6b.ini --region chr2:71,630,381-71,667,029 --dpi 600 -o 6b.png 