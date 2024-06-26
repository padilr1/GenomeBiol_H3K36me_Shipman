#!/bin/bash

make_tracks_file --trackFiles 10T.PA.K36me1.merged.ms_cpm.bw 10T.K36M.K36me1.merged.ms_cpm.bw 10T.DKO.K36me1.merged.ms_cpm.bw 10T.TKO.K36me1.merged.ms_cpm.bw 10T.QKO.K36me1.merged.ms_cpm.bw 10T.QuiKO.K36me1.merged.ms_cpm.bw mm10.ncbiRefSeq.gtf -o 3a.ini

pyGenomeTracks --tracks 3a.ini --region chr14:55,418,448-61,853,143 --dpi 600 -o 3a.png 