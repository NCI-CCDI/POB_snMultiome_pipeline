library(optparse)
# Create an option parser
option_parser <- OptionParser(usage = "Usage: %prog --arg1 <value1> --arg2 <value2>")

# Define the options
option_parser <- add_option(option_parser, "--analysis_dir", type = "character", default = NULL,
                            help = "This is the current working dir")
option_parser <- add_option(option_parser, "--sample", type = "character", default = NULL,
                            help = "This is sample name")
option_parser <- add_option(option_parser, "--nCount_RNA_LL", type = "integer", default = 500,
                            help = "Lower Limit for nCount_RNA")
option_parser <- add_option(option_parser, "--nCount_ATAC_LL", type = "integer", default = 500,
                            help = "Lower Limit for nCount_ATAC")
option_parser <- add_option(option_parser, "--nFeature_RNA_LL", type = "integer", default = 300,
                            help = "Lower Limit for nFeature_RNA")
option_parser <- add_option(option_parser, "--nFeature_ATAC_LL", type = "integer", default = 300,
                            help = "Lower Limit for nFeature_ATAC")
option_parser <- add_option(option_parser, "--nucleosome_signal_UL", type = "integer", default = 2,
                            help = "Upper Limit for nucleosome_signal")
option_parser <- add_option(option_parser, "--TSSenrichment_LL", type = "integer", default = 2,
                            help = "Lower Limit for TSSenrichment")
option_parser <- add_option(option_parser, "--pct_reads_in_peaks_LL", type = "numeric", default = 0.15,
                            help = "Lower Limit for pct_reads_in_peaks_LL")
option_parser <- add_option(option_parser, "--blacklist_fraction_UL", type = "numeric", default = 0.05,
                            help = "Upper Limit for blacklist")
option_parser <- add_option(option_parser, "--mt_UL", type = "integer", default = 10,
                            help = "Upper Limit for mt percentage")
opt <- parse_args(option_parser)


working_dir = opt$analysis_dir
sample = opt$sample
nCount_RNA_LL = opt$nCount_RNA_LL
nCount_ATAC_LL = opt$nCount_ATAC_LL 
nFeature_RNA_LL = opt$nFeature_RNA_LL
nFeature_ATAC_LL = opt$nFeature_ATAC_LL
nucleosome_signal_UL = opt$nucleosome_signal_UL
TSSenrichment_LL = opt$TSSenrichment_LL 
pct_reads_in_peaks_LL = opt$pct_reads_in_peaks_LL
blacklist_fraction_UL = opt$blacklist_fraction_UL
mt_UL = opt$mt_UL
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/Snake"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
cat(working_dir,"\n")
cat("nCount_RNA_LL:", nCount_RNA_LL, ", nFeature_RNA_LL:",nFeature_RNA_LL,", mt_UL:",mt_UL,"\n",
    "nCount_ATAC_LL:", nCount_ATAC_LL,", nFeature_ATAC_LL:", nFeature_ATAC_LL,", pct_reads_in_peaks_LL:",pct_reads_in_peaks_LL,"\n",
    "nucleosome_signal_UL:",nucleosome_signal_UL,", TSSenrichment_LL:",TSSenrichment_LL,", blacklist_fraction_UL:",blacklist_fraction_UL,"\n")
setwd(working_dir)
source('scripts/R/main/load_packages.R')

output_path = paste0(working_dir,"/result_wSoupX")
plot_path = paste0(output_path,"/qc_result_plot")
sample_plot_path =  paste0(output_path,"/sample_result_plot")

sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet))
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist

#mt_UL = 10 #outliers_mad(object$percent.mt,threshold = 5)$UL_CI_MAD
obj_name = paste0(sample,"_GEXATACobj_raw.rds")
#print(obj_name)
file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
#print(sample_outdir_list[[seob]])
object = readRDS(file_name) 
#sample = Project(object)
m = object@meta.data
nCount_RNA_UL = outliers_mad(object$nCount_RNA,threshold = 3)$UL_CI_MAD
nCount_ATAC_UL = outliers_mad(object$nCount_ATAC,threshold = 3)$UL_CI_MAD  #max(object$nCount_ATAC) + 1000  #outliers_mad(object$nCount_ATAC,threshold = 3)$UL_CI_MAD
percentile_now = 1 - (ncol(object) / 1000) * 0.008  #96% percentile was used for CCDI4
doublet_threshold <- quantile(m$doublet_scores, percentile_now)
#do filtering now
object_filter = object %>% subset(
    subset=
      nCount_RNA > nCount_RNA_LL & nCount_RNA < nCount_RNA_UL &
      nCount_ATAC > nCount_ATAC_LL &  nCount_ATAC < nCount_ATAC_UL &
      nFeature_RNA > nFeature_RNA_LL & nFeature_ATAC > nFeature_ATAC_LL &
      nucleosome_signal < nucleosome_signal_UL & TSS.enrichment > TSSenrichment_LL &
      percent.mt < mt_UL & pct_reads_in_peaks > pct_reads_in_peaks_LL & blacklist_fraction < blacklist_fraction_UL &
      doublet_scores < doublet_threshold)

obj_file = paste0(sample,"_GEXATACobj_filt.rds")
saveRDS(object_filter,paste0(output_path,"/seurat/preprocessed/",obj_file))
cat("total cell after filtering for",sample, ncol(object_filter),"\n")
