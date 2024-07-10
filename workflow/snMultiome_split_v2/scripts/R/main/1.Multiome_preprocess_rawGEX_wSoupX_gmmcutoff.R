Sys.setenv(RETICULATE_PYTHON = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/reticulate-env/bin/python")

args=commandArgs(trailingOnly=T)
#args <- commandArgs()
print(paste0("args: ",args))
working_dir = as.character(args[1])
project_name = as.character(args[2])
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/Analysis_filtGEXATAC_split"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"

print(paste0("working_dir: ",working_dir))
print(paste0("project_name: ",project_name))
#cat(input_path,working_dir,"\n")
setwd(working_dir)
set.seed(123)
source('scripts/R/main/load_packages.R')
script_path = paste0(working_dir,"/scripts/R/utility")
output_path = paste0(working_dir,"/result_wSoupX")
subfolder_names = c("seurat","integration","qc_result_plot","integration_plot","seurat/preprocessed","seurat/merged","sample_result_plot")
#for (subfolder in subfolder_names) {
#  subfolder_path <- file.path(output_path, subfolder)
#  dir.create(subfolder_path, recursive = TRUE, showWarnings = FALSE)
#}

# Check if subfolders were created successfully
for (subfolder in subfolder_names) {
  subfolder_path <- file.path(output_path, subfolder)
  if (file.exists(subfolder_path)) {
    cat("Subfolder", subfolder, "created successfully.\n")
  } else {
    cat("Failed to create subfolder", subfolder, "\n")
  }
}

sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet)) 
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist
sample_outdir_list = list()

object.list = list()
cutoff_list = list()
for (sample in samplelist){
  obj_name = paste0(sample,"_GEXobj_raw.rds")
  print(obj_name)
  file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
  object = readRDS(file_name)
  object.list[[sample]]=object
  cutoff_list[["total"]][[sample]] = ncol(object)
  cutoff_list[["low_quality"]][[sample]] = sum(object@meta.data$scrub.if_doublet == "Low_quality")
  cutoff_list[["cutoff"]][[sample]] = unique(object@meta.data$scrub.cutoff_final)
  cutoff_list[["criteria"]][[sample]] = unique(object@meta.data$scrub.criteria)
  cutoff_list[["doublet"]][[sample]] = sum(object@meta.data$scrub.if_doublet == "Doublet")
  cutoff_list[["percentile"]][[sample]] = unique(object@meta.data$scrub.percentile)
  cutoff_list[["percentile_cutoff"]][[sample]] = unique(object@meta.data$scrub.cutoff_percentile)
  cutoff_list[["percentile_doublet"]][[sample]] = sum(object@meta.data$scrub.if_doublet_percentile =="Doublet")
  cutoff_list[["gmm_cutoff"]][[sample]] = unique(object@meta.data$scrub.cutoff_gmm)
  cutoff_list[["gmm_doublet"]][[sample]] = sum(object@meta.data$scrub.if_doublet_gmm =="Doublet")
  cutoff_list[["percentile_adj"]][[sample]] = unique(object@meta.data$scrub.percentile_adj)
  cutoff_list[["percentile_cutoff_adj"]][[sample]] = unique(object@meta.data$scrub.cutoff_percentile_adj)
  cutoff_list[["percentile_doublet_adj"]][[sample]] = sum(object@meta.data$scrub.if_doublet_percentile_adj =="Doublet")
  
}
saveRDS(object.list, paste0(output_path,"/seurat/1.object_list_rawGEX.rds"))
#saveRDS(object.list, paste0(output_path,"/seurat/1.object_list_rawraw.rds"))

#object.list <- Doublet_GMM_cutoff(object.list, project.name = project_name, outpath = output_path, filter.doublets = FALSE)


if (length(samplelist) != length(names(cutoff_list[["criteria"]]))) {
    stop("Length of samplelist does not match the number of criteria.")
} else {
  df <- data.frame(
    sample = unlist(samplelist),
    total_cell_number = unlist(cutoff_list[["total"]]),
    low_quality_cell_number = unlist(cutoff_list[["low_quality"]]),
    cutoff_criteria = unlist(cutoff_list[["criteria"]]),
    doublet_score_cutoff = unlist(cutoff_list[["cutoff"]]),  
    doublet_cell_number = unlist(cutoff_list[["doublet"]]), 
    percentile = unlist(cutoff_list[["percentile"]]),
    percentile_cutoff =  unlist(cutoff_list[["percentile_cutoff"]]),
    doublet_percentile = unlist(cutoff_list[["percentile_doublet"]]),
    percentile_adj = unlist(cutoff_list[["percentile_adj"]]),
    percentile_cutoff_adj =  unlist(cutoff_list[["percentile_cutoff_adj"]]),
    doublet_percentile_adj = unlist(cutoff_list[["percentile_doublet_adj"]]),
    gmm_cutoff = unlist(cutoff_list[["gmm_cutoff"]]),
    doublet_gmm = unlist(cutoff_list[["gmm_doublet"]]))
  print(df)
  write.csv(df,paste0(output_path,"/qc_result_plot/scrublet_doublet_filter_summary.csv"),row.names=F,quote =F)
}
# Extract sample names
#sample_names <- names(cutoff_list[["cutoff"]])

# Combine sample names with other elements of the list
#combined_list <- cbind(sample = sample_names, cutoff_list)

# Create dataframe using do.call
#df <- do.call(data.frame, combined_list)


#for (seob in samplelist){
#  obj_name = paste0(seob,"_GEXobj_raw.rds")
#  #print(obj_name)
#  file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
#  #print(sample_outdir_list[[seob]])
#  saveRDS(object.list[[seob]],file_name) 
#}



