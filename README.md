# POB_snMultiome_pipeline
This snakemake pipeline was developed for processing snMultiome data and currently under active development. It has been developed solely on Biowulf. 

## 
This pipeline includes the following steps described below.

Quality check of raw fastq reads
Trimming low quality reads
Check quality of the trimmed reads
Make a QC HTML report
Call Juicer HiC tool which involves the following steps
Generate Hi-C contact maps
Normalization of hic files
Call Arrowhead tool to identify TADs
Call Hiccups tool to identify loops
