* The following programs are required:
  * [deepTools](https://github.com/deeptools/deepTools)
  * [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)
  * [SPAN](https://github.com/JetBrains-Research/span)
  * [featureCounts](https://subread.sourceforge.net/featureCounts.html)
  * [HOMER](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)
  * [R](https://cran.r-project.org/)
    * Bioconductor
      * AnnotationDbi
      * ChIPseeker
      * DESeq2
      * DiffBind
      * GenomicFeatures
      * GenomicRanges
      * SummarizedExperiment
      * TCGAbiolinks
      * TxDb.Hsapiens.UCSC.hg38.knownGene
      * fgsea
      * limma
      * org.Hs.eg.db
      * rtracklayer
      * sva
      * tximport
      * zFPKM
      * grid
      * gridExtra
    * CRAN
      * MASS
      * cowplot
      * data.table
      * forestplot
      * furrr
      * ggforce
      * ggformula
      * ggnewscale
      * ggrepel
      * ggsignif
      * hexbin
      * isoband
      * lwgeom
      * msigdbr
      * pals
      * patchwork
      * purrr
      * readxl
      * reshape2
      * RobustRankAggreg
      * scales
      * sf
      * tidyverse
      * viridis
* The scripts below should be ran to pre-process a number of files 
  * `generate_norm_bw.sh` provides examples of how to generate the normalized bigWigs used to generate the figures, including depth-normalized bigWigs <code>.cpm.bw</code> and input-normalized bigWigs <code>.input_normalized.log2cpm.bw</code>. Furthermore, it provides an example of how replicates (in bigWig format) were merged.
  * `mass_spec_norm_bw.R` demonstrates how to use mass-spectrometry values to quantitatively scale the depth-normalized bigWigs. 
  * `generate_binned_BED_file.sh` demonstrates how to generate 1kb and 10kb binned <code>BED</code> files from aligned <code>BAM</code> files.
  * `generate_tables_for_metagene_body_plots.md` demonstrates how to generate tables needed for plotting meta-gene body plots.
  * `compute_exon_intron_ratios.sh` computes exon versus intron ratios for a given aligned ChIP-seq sample file. This step will have to be run prior to generating the exon-intron ratio plots in R. 
* Figures generated using the scripts above will largely resemble the figures included in the manuscript. However, there may be minor discrepancies as a result of stochasticity in specific algorithms as well as additional manual aesthetic adjustments. Nevertheless, these minor differences do not impact the conclusions presented.
