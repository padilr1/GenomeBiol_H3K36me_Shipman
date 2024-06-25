#!/bin/bash

# make track file
make_tracks_file --trackFiles 10T.PA_merged.H3K36me1.cpm.bw 10T.PA_merged.H3K36me2.cpm.bw 10T.PA_merged.H3K36me3.cpm.bw mm10.ncbiRefSeq.gtf -o data/1a.ini

# plot track file
pyGenomeTracks --tracks data/PA_H3K36me_genome_browser_tracks.ini --region chr9:50368630-52224861 -o figs/1a.pdf
