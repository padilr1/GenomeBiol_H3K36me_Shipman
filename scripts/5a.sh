#!/bin/bash
make_tracks_file --trackFiles 10T.PA.K36me2.merged.ms_cpm.bw 10T.SETD2KO.K36me2.merged.ms_cpm.bw 10T.DKO.K36me2.merged.ms_cpm.bw 10T.K36M_OE.K36me2.merged.ms_cpm.bw mm10.ncbiRefSeq.gtf -o data/5a.ini

pyGenomeTracks --tracks data/5a.ini --region chr6:13,620,190-13,686,722 --dpi 600 -o 5a.png 