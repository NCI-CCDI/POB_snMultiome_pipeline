args=commandArgs(trailingOnly=T)
input_path = as.character(args[1])
working_dir = as.character(args[2])
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/Snake"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
cat(input_path,working_dir,"\n")
setwd(working_dir)
source('scripts/R/main/load_packages.R')

output_path = paste0(working_dir,"/result_wSoupX")
plot_path = paste0(output_path,"/qc_result_plot")
sample_plot_path =  paste0(output_path,"/sample_result_plot")

sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet))
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist

file_name = paste0(output_path,"/seurat/3.object_list_GEXATAC_filt.rds")
filt_object.list = readRDS(file_name )

macs_path = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/seuratv5_env/bin/macs2"
frag.files <- sapply(samplelist, function(project) paste0(input_path, project, '/outs/atac_fragments.tsv.gz'), simplify=FALSE) 

# call peaks using MACS2
macs_peaks.list = list()
for (s in samplelist) {
  obj = filt_object.list[[s]]
  peaks <- CallPeaks(obj, macs2.path = macs_path)
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  macs_peaks.list[[s]] = peaks
}

## now create a unified set of peaks to quantify in each dataset
# Create a GRangesList from the list of GRanges objects
new_granges <- GRangesList(macs_peaks.list)
# Used reduce function here
combined.peaks <- GenomicRanges::reduce(x =unlist(new_granges))
cat("the total number of peaks after Reduce is:",length(combined.peaks@ranges),"\n")

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
cat("the total number of peaks after Reduce and filtering is:",length(combined.peaks@ranges),"\n")
peak_file_name = paste0(output_path,"/seurat/3.ATAC_combined_peaks.rds")
saveRDS(combined.peaks,peak_file_name)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

#project = c("CCDI0268_TCU00389-0101-Rep1")
con_object.list <- sapply(samplelist, function(project) {
  object <- filt_object.list[[project]]
  fragments <- Fragments(object[['ATAC']])[[1]]
  cells.use <- GetFragmentData(fragments, slot='cells')               
  print(length(cells.use))
  print(ncol(object))
  
  object <- object %>% subset(cells=cells.use)
  fragpath <- paste0(input_path, "/", project, '/outs/atac_fragments.tsv.gz')
  
  atac_counts <- FeatureMatrix( 
    cells=cells.use,
    features=combined.peaks,
    fragments=fragments
  )
  object[['MACS']] <- CreateChromatinAssay(
    counts=atac_counts,
    fragments=fragpath,
    annotation=annotations
  )
  
  object
}, simplify=FALSE)

for (s in samplelist){
  obj = con_object.list[[s]]
  obj_file = paste0(s,"_GEXATACobj_filt_consensus.rds")
  saveRDS(obj,paste0(output_path, "/seurat/preprocessed/",obj_file))
  
}

saveRDS(con_object.list, paste0(output_path,"/seurat/3.object_list_filtGEXATAC_consnsus.rds"))

