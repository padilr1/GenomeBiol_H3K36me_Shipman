#!/bin/bash

make_tracks_file --trackFiles 10T.PA.K36me2.merged.ms_cpm.bw 10T.DKO.K36me2.merged.ms_cpm.bw 10T.TKO.K36me2.merged.ms_cpm.bw 10T.QKO.K36me2.merged.ms_cpm.bw mm10.ncbiRefSeq.gtf -o data/6a.ini

pyGenomeTracks --tracks data/6a.ini --region chr3:135,253,035-135,768,340 --dpi 600 -o 6a.png 