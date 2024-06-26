#!/bin/bash

make_tracks_file --trackFiles 10T.TKO.K36me2.merged.ms_cpm.bw 10T.TKO.ATAC.merged.cpm.bw TKO_H3K27ac_1.cpm.bw promoter.bed mm10.ncbiRefSeq.gtf -o data/6c.ini

pyGenomeTracks --tracks data/6c.ini --region chr19:47,511,791-47,577,523 --dpi 600 -o 6c.png