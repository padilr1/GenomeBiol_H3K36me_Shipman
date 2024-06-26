#!/bin/bash

# make genome browser track for PA, NSD1-KO,NSD2-KO, DKO, K36M-OE, TKO, QKO and QuiKO
make_tracks_file --trackFiles 10T.PA_merged.H3K36me2.ms_cpm.bw 10T.NSD1KO_merged.H3K36me2.ms_cpm.bw 10T.NSD2KO_merged.H3K36me2.ms_cpm.bw 10T.DKO_merged.H3K36me2.ms_cpm.bw 10T.K36M_OE_merged.H3K36me2.ms_cpm.bw 10T.TKO_merged.H3K36me2.ms_cpm.bw 10T.QKO_merged.H3K36me2.ms_cpm.bw 10T.QuiKO_merged.H3K36me2.ms_cpm.bw mm10.ncbiRefSeq.gtf -o data/3a.ini

pyGenomeTracks --tracks data/3a.ini --region chr3:27,139,141-27,204,558 -o figs/3a.pdf