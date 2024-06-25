#!/bin/bash

make_tracks_file --trackFiles 10T.PA_merged.H3K36me1.cpm.bw 10T.PA_merged.H3K36me2.cpm.bw 10T.PA_merged.H3K36me3.cpm.bw mm10.ncbiRefSeq.gtf -o data/1d.ini

pyGenomeTracks --tracks data/1d.ini --region chr10:95,930,897-96,059,439 -o figs/1d.pdf