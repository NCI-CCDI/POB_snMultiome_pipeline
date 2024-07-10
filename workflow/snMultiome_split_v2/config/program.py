##setup R version 4.3.2 and the Rlib
Rscript = "export R_LIBS_USER=/data/CCRCCDI/software/lib/R/library_4.3.2; module load R/4.3.2; echo $R_LIBS_USER; Rscript"

##setup python
python3 = "module load python/3.9;python"

##setup cellranger
cellranger_arc = "/data/CCRCCDI/software/Tools/cellranger/cellranger-arc-2.0.2/cellranger-arc"
