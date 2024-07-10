library(optparse)
# Create an option parser
option_parser <- OptionParser(usage = "Usage: %prog --arg1 <value1> --arg2 <value2>")

# Define the options
option_parser <- add_option(option_parser, "--analysis_dir", type = "character", default = NULL,
                            help = "This is the current working dir")
#option_parser <- add_option(option_parser, "--sample", type = "character", default = NULL,
#                            help = "This is sample name")
option_parser <- add_option(option_parser, "--nCount_RNA_LL", type = "integer", default = 500,
                            help = "Lower Limit for nCount_RNA")
option_parser <- add_option(option_parser, "--nCount_RNA_UL", type = "character", default = "3MADS",
                            help = "Lower Limit for nCount_RNA")
option_parser <- add_option(option_parser, "--nCount_ATAC_LL", type = "integer", default = 0,
                            help = "Lower Limit for nCount_ATAC for minimal ATAC filtering")
option_parser <- add_option(option_parser, "--nCount_ATAC_UL", type = "character", default = "3MADS",
                            help = "Lower Limit for nCount_ATAC for minimal ATAC filtering")
option_parser <- add_option(option_parser, "--nFeature_RNA_LL", type = "integer", default = 300,
                            help = "Lower Limit for nFeature_RNA")
option_parser <- add_option(option_parser, "--nFeature_ATAC_LL", type = "integer", default = 0,
                            help = "Lower Limit for nFeature_ATAC for minimal ATAC filtering")
option_parser <- add_option(option_parser, "--nucleosome_signal_UL", type = "integer", default = 2,
                            help = "Upper Limit for nucleosome_signal")
option_parser <- add_option(option_parser, "--TSSenrichment_LL", type = "integer", default = 2,
                            help = "Lower Limit for TSSenrichment")
option_parser <- add_option(option_parser, "--pct_reads_in_peaks_LL", type = "numeric", default = 0,
                            help = "Lower Limit for pct_reads_in_peaks_LL for minimal ATAC filtering")
option_parser <- add_option(option_parser, "--blacklist_fraction_UL", type = "numeric", default = 0.05,
                            help = "Upper Limit for blacklist")
option_parser <- add_option(option_parser, "--mt_UL", type = "integer", default = 10,
                            help = "Upper Limit for mt percentage")
opt <- parse_args(option_parser)


working_dir = opt$analysis_dir
sample = opt$sample
nCount_RNA_LL = opt$nCount_RNA_LL
nCount_ATAC_LL = opt$nCount_ATAC_LL 
nCount_RNA_UL = opt$nCount_RNA_UL
nCount_ATAC_UL = opt$nCount_ATAC_UL 
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
cat("nCount_RNA_UL:", nCount_RNA_UL, ", nCount_ATAC_UL:",nCount_ATAC_UL,"\n")


setwd(working_dir)
source('scripts/R/main/load_packages.R')

output_path = paste0(working_dir,"/result_wSoupX")
plot_path = paste0(output_path,"/qc_result_plot")
sample_plot_path =  paste0(output_path,"/sample_result_plot")

sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet))
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist

object.list=list()
for (sample in samplelist) {
  obj_name = paste0(sample,"_GEXATACobj_raw.rds")
  #print(obj_name)
  file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
  #print(sample_outdir_list[[seob]])
  object = readRDS(file_name) 
 object.list[[sample]] = object
}
 
#saveRDS(object.list,paste0(output_path,"/seurat/2.object_list_rawGEXATAC.rds"))


#mt_UL = 10 #outliers_mad(object$percent.mt,threshold = 5)$UL_CI_MAD

df_threshold = data.frame("Sample" = character(), "nCount_RNA_LL" = numeric(), "nCount_RNA_UL" = numeric(), "nFeature_RNA_LL" = numeric(), "mt_UL" = numeric(),
                         "nCount_ATAC_LL" = numeric(), "nCount_ATAC_UL" = numeric(), "nFeature_ATAC_LL"= numeric(),
                         "nucleosome_signal_UL" = numeric(),  "TSSenrichment_LL" = numeric(),
                         "pct_fragments_in_peaks_LL"= numeric(), "blacklist_fraction_LL"= numeric(),"scrublet_criteria"=character(),"doublet_score_cutoff"=numeric(),
                          stringsAsFactors = FALSE)
df_filter = data.frame("Sample" = character(), "nCount_RNA_LL"=numeric(), "nCount_RNA_UL"=numeric(),"nFeature_RNA_LL"=numeric(),"percent_mt"=numeric(),
                       "nCount_ATAC_LL"=numeric(),"nCount_ATAC_UL"=numeric(),"nFeature_ATAC_LL"=numeric(),
                       "nucleosome_signal"=numeric(), "TSSenrichment"=numeric(),
                       "pct_fragments_in_peaks"=numeric(), "blacklist_fraction"=numeric(),
                       "total_filt"=numeric(), "scrablet_doublet"=numeric(),"total_with_doublet_filter"=numeric(),"scrublet_doublet_percentile"=numeric(),"total_with_doublet_filterpercentile"=numeric(),
                       "total"=numeric(),
                       stringsAsFactors = FALSE)
nCount_RNA_UL_threshold = nCount_RNA_UL
nCount_ATAC_UL_threshold = nCount_ATAC_UL
for (sample in samplelist) {
  print(sample)
  object = object.list[[sample]]
  m = object@meta.data
  print(object)
  MADS_threshold_RNA = nCount_RNA_UL_threshold
  if (grepl("MAD", nCount_RNA_UL_threshold )) {
    MADS_threshold_RNA <- as.numeric(gsub("MADS", "", MADS_threshold_RNA))
    nCount_RNA_UL = outliers_mad(object$nCount_RNA,threshold = MADS_threshold_RNA)$UL_CI_MAD
    cat ("RNA MADS threshold:", MADS_threshold_RNA,"\n","nCount_RNA_UL:",nCount_RNA_UL,"\n")
  } else if (grepl("max", nCount_RNA_UL_threshold )) {
    nCount_RNA_UL = max(object$nCount_RNA) + 1000
    print("nCount_RNA_UL is max+1000")
  } else if (grepl("^[0-9]+$",nCount_RNA_UL_threshold)) {
    nCount_RNA_UL = as.numeric(nCount_RNA_UL_threshold)
    print ("nCount_RNA_UL is a hard cutoff")
    cat ("nCount_RNA_UL:",nCount_RNA_UL,"\n")
  } else {
    print("nCount_RNA_UL is not in correct format!!!!!!!!!!!!!!!!!!")
    print(nCount_RNA_UL)
    print("nCount_RNA_UL is not in correct format!!!!!!!!!!!!!!!!!!")
  }
  
  MADS_threshold_ATAC = nCount_ATAC_UL_threshold
  if (grepl("MAD", nCount_ATAC_UL_threshold)) {
    MADS_threshold_ATAC <- as.numeric(gsub("MADS", "", MADS_threshold_ATAC))
    nCount_ATAC_UL = outliers_mad(object$nCount_ATAC,threshold = MADS_threshold_ATAC)$UL_CI_MAD
    print("nCount_ATAC_UL is MADS method")
    cat ("nCount_ATAC_UL:",nCount_ATAC_UL,"\n")
  } else if (grepl("max", nCount_ATAC_UL_threshold)) {
    nCount_ATAC_UL = max(object$nCount_ATAC) + 1000
    print("nCount_ATAC_UL is max + 1000")
    cat ("nCount_ATAC_UL:",nCount_ATAC_UL,"\n")
  } else if (grepl("^[0-9]+$",nCount_ATAC_UL_threshold)) {
    nCount_ATAC_UL = as.numeric(nCount_ATAC_UL_threshold)
    print ("nCount_RNA_UL is a hard cutoff")
    cat ("nCount_ATAC_UL:",nCount_ATAC_UL,"\n")
  } else {
    print("nCount_ATAC_UL is not in correct format!!!!!!!!!!!!!!!!!!")
    print(nCount_ATAC_UL)
    print("nCount_ATAC_UL is not in correct format!!!!!!!!!!!!!!!!!!")
  }
  
  #nCount_ATAC_UL =  outliers_mad(object$nCount_ATAC,threshold = 3)$UL_CI_MAD
  #percentile_now = 1 - (ncol(object) / 1000) * 0.008
  #percentile_now = round(percentile_now,2)
  percentile_now = unique(object@meta.data$scrub.percentile)
  doublet_criteria = unique(object@meta.data$scrub.criteria)
  #doublet_threshold <- quantile(m$doublet_scores, percentile_now)
  doublet_threshold = unique(object@meta.data$scrub.cutoff_final)
  thresholdStat = c(sample,nCount_RNA_LL,nCount_RNA_UL,nFeature_RNA_LL,mt_UL,
                nCount_ATAC_LL,nCount_ATAC_UL,nFeature_ATAC_LL,
                nucleosome_signal_UL,TSSenrichment_LL,pct_reads_in_peaks_LL,blacklist_fraction_UL,doublet_criteria,doublet_threshold)

  df_threshold_now = rbind(thresholdStat=thresholdStat)

  colnames(df_threshold_now) = c("Sample", "nCount_RNA_LL","nCount_RNA_UL", "nFeature_RNA_LL", "mt_UL",
                       "nCount_ATAC_LL", "nCount_ATAC_UL", "nFeature_ATAC_LL",
                       "nucleosome_signal_UL",  "TSSenrichment_LL",
                       "pct_fragments_in_peaks_LL", "blacklist_fraction_LL","scrublet_criteria","doublet_score_cutoff")
  df_threshold = rbind(df_threshold, df_threshold_now)
  

  total_cell = nrow(m)
  nCountRNA_LL_filt = sum(m$nCount_RNA <= nCount_RNA_LL)
  nCountRNA_UL_filt = sum(m$nCount_RNA >= nCount_RNA_UL)
  nFeatureRNA_LL_filt = sum(m$nFeature_RNA <= nFeature_RNA_LL)
  mt_UL_filt = sum(m$percent.mt >= mt_UL)
  
  nCountATAC_LL_filt = sum(m$nCount_ATAC <= nCount_ATAC_LL)
  nCountATAC_UL_filt = sum(m$nCount_ATAC >= nCount_ATAC_UL)
  nFeatureATAC_LL_filt = sum(m$nFeature_ATAC <= nFeature_ATAC_LL)
  
  nucleasome_UL_filt = sum(m$nucleosome_signal >= nucleosome_signal_UL,na.rm = TRUE)
  TSS_UL_filt = sum(m$TSS.enrichment <= TSSenrichment_LL)
  Frip_LL_filt = sum(m$pct_reads_in_peaks <= pct_reads_in_peaks_LL)
  Blacklist_filt = sum(m$blacklist_fraction >= blacklist_fraction_UL)
  #Doublet_filt = sum(m$if.doublet.gmm == "doublet")
  #Doublet_filt2 = sum(m$doublet_scores>=doublet_threshold)
  Doublet_filt = sum(m$scrub.if_doublet ==  "Doublet")
  Doublet_filt2 = sum(m$scrub.if_doublet_percentile == "Doublet")
                     
  #do filtering now
  object_filt = object %>% subset(
    subset=
      nCount_RNA > nCount_RNA_LL & nCount_RNA < nCount_RNA_UL &
      nCount_ATAC > nCount_ATAC_LL &  nCount_ATAC < nCount_ATAC_UL &
      nFeature_RNA > nFeature_RNA_LL & nFeature_ATAC > nFeature_ATAC_LL &
      nucleosome_signal < nucleosome_signal_UL & TSS.enrichment > TSSenrichment_LL &
      percent.mt < mt_UL & pct_reads_in_peaks > pct_reads_in_peaks_LL & blacklist_fraction < blacklist_fraction_UL)
  object_filt_doublet = object %>% subset(
    subset=
      nCount_RNA > nCount_RNA_LL & nCount_RNA < nCount_RNA_UL &
      nCount_ATAC > nCount_ATAC_LL &  nCount_ATAC < nCount_ATAC_UL &
      nFeature_RNA > nFeature_RNA_LL & nFeature_ATAC > nFeature_ATAC_LL &
      nucleosome_signal < nucleosome_signal_UL & TSS.enrichment > TSSenrichment_LL &
      percent.mt < mt_UL & pct_reads_in_peaks > pct_reads_in_peaks_LL & blacklist_fraction < blacklist_fraction_UL &
      scrub.if_doublet == "Singlet")
  object_filt_doublet2 = object %>% subset(
    subset=
      nCount_RNA > nCount_RNA_LL & nCount_RNA < nCount_RNA_UL &
      nCount_ATAC > nCount_ATAC_LL &  nCount_ATAC < nCount_ATAC_UL &
      nFeature_RNA > nFeature_RNA_LL & nFeature_ATAC > nFeature_ATAC_LL &
      nucleosome_signal < nucleosome_signal_UL & TSS.enrichment > TSSenrichment_LL &
      percent.mt < mt_UL & pct_reads_in_peaks > pct_reads_in_peaks_LL & blacklist_fraction < blacklist_fraction_UL &
      scrub.if_doublet_percentile == "Singlet" )
  total_cell_left = ncol(object_filt)
  total_cell_left_doublet = ncol(object_filt_doublet)
  total_cell_left_doublet2 = ncol(object_filt_doublet2)
  total_filt = total_cell - total_cell_left 
  total_filt_doublet = total_cell - total_cell_left_doublet
  total_filt_doublet2 = total_cell - total_cell_left_doublet2
  filtStat = c(sample,nCountRNA_LL_filt,nCountRNA_UL_filt,nFeatureRNA_LL_filt,mt_UL_filt,
               nCountATAC_LL_filt,nCountATAC_UL_filt,nFeatureATAC_LL_filt,
               nucleasome_UL_filt,TSS_UL_filt,Frip_LL_filt,Blacklist_filt,
               total_filt, Doublet_filt,total_filt_doublet,Doublet_filt2,total_filt_doublet2,total_cell)
  df_filt_now = rbind(filtStat = filtStat )
  
  colnames(df_filt_now) = data.frame("Sample", "nCount_RNA_LL", "nCount_RNA_UL","nFeature_RNA_LL","percent_mt",
                         "nCount_ATAC_LL","nCount_ATAC_UL","nFeature_ATAC_LL",
                         "nucleosome_signal", "TSSenrichment",
                         "pct_fragments_in_peaks", "blacklist_fraction",
                         "total_filt","scrublet_doublet","total_with_doublet_filter","scrublet_doublet_percentile","total_with_doublet_filterpercentile","total")
  df_filter = rbind(df_filter,df_filt_now)
}

df_filter
df_threshold
write.table(df_filter, file = paste(plot_path,"/QC_stats_FilteredCell.tsv",sep=""), quote=FALSE,sep="\t", row.names=FALSE)
write.table(df_threshold, file = paste(plot_path,"/QC_stats_FilteredCell_threshold.tsv",sep=""), quote=FALSE,sep="\t", row.names=FALSE)


filt_object.list = list()
for (seob in samplelist){
  obj_name = paste0(seob,"_GEXATACobj_filt.rds")
  #print(obj_name)
  file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
  #print(sample_outdir_list[[seob]])
  object = readRDS(file_name) 
  filt_object.list[[seob]] = object
}

saveRDS(filt_object.list,paste0(output_path,"/seurat/3.object_list_GEXATAC_filt.rds"))

PlotQC(filt_object.list, prefix="post", outpath=plot_path,filter.doublets=FALSE)
