#!/bin/bash

make_tracks_file --trackFiles 10T.PA_merged.H3K36me2.ms_cpm.bw 10T.SETD2KO_merged.H3K36me2.ms_cpm.bw 10T.DKO_merged.H3K36me2.ms_cpm.bw 10T.K36M_OE_merged.H3K36me2.ms_cpm.bw mm10.ncbiRefSeq.gtf -o data/4a.ini

pyGenomeTracks --tracks data/4a.ini --region chr6:13,620,190-13,686,722 -o figs/4a.pdf
