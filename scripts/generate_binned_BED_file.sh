# convert aligned BAM file to BED file
bedtools bamtobed -i ${sample}.sorted.bam | sort -k1,1 -k2,2n > ${sample}.bed
# intersect sample BED file with binned reference files
bedtools intersect -a $WINDOW_REF -b ${sample}.bed -sorted -c -nonamecheck > $OUTPUT_NAME.binned.bed