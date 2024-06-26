#!/bin/bash

# after aligning and processing the FASTQ file, depth-normalized bigWigs (in counts per million (CPM)) for each sample was generated using:
bamCoverage -b ${sample}.sorted.bam -o ${sample}.cpm.bw --normalizeUsing CPM --centerReads -e 200

# input-normalized bigWigs were generated using:
bamCompare -b1 ${ChIP_sample}.bam -b2 ${input_sample}.bam -o ${sample}.input_normalized.log2cpm.bw -e 200 --centerReads -bs 200 --blackListFileName blacklist.bed --normalizeUsing CPM --scaleFactorsMethod None

# samples were merged using:
bigwigCompare --operation mean -bs 200 -b1 10T.${sample_rep1}.cpm.bw -b2 10T.${sample_rep2}.cpm.bw -o merged_bigwig_first_two_reps.bw
bigwigCompare -b1 merged_bigwig_first_two_reps.bw -b2 10T.${sample_rep3}.cpm.bww -o 10T.${sample}_merged.${histone_mark}.cpm.bw

