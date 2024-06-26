#!/bin/bash

# zoomed out
make_tracks_file --trackFiles 10T.QKO.K36me2.merged.ms_cpm.bw 10T.QuiKO.K36me2.merged.ms_cpm.bw mm10.ncbiRefSeq.gtf -o data/7a_1.ini

pyGenomeTracks --tracks data/7a_1.ini --region chr5:10,411,523-135,566,683 --dpi 600 -o 7a_1.png

# zoomed in
make_tracks_file --trackFiles 10T.QKO.K36me2.merged.ms_cpm.bw 10T.QuiKO.K36me2.merged.ms_cpm.bw 10T.ATAC.QKO.merged.cpm.bw mm10.ncbiRefSeq.gtf -o data/7a_2.ini

pyGenomeTracks --tracks data/7a_2.ini --region chr5:81,000,872-81,123,093 --dpi 600 -o 7a_2.png 