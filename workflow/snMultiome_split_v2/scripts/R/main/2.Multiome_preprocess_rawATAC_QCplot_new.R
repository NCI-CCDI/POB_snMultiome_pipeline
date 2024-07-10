args=commandArgs(trailingOnly=T)
working_dir = as.character(args[1])
project_name = as.character(args[2])


#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/Snake"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
cat(working_dir,"\n")
setwd(working_dir)
library(RColorBrewer)
source('scripts/R/main/load_packages.R')

output_path = paste0(working_dir,"/result_wSoupX")
plot_path = paste0(output_path,"/qc_result_plot")
sample_plot_path =  paste0(output_path,"/sample_result_plot")
sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet))
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist
#for (s in samplelist) {
#  subfolder_path <- file.path(sample_plot_path, s)
#  dir.create(subfolder_path, recursive = TRUE, showWarnings = FALSE)
#}

newobject.list = list()
for (seob in samplelist){
  obj_name = paste0(seob,"_GEXATACobj_raw.rds")
  #print(obj_name)
  file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
  #print(sample_outdir_list[[seob]])
  object = readRDS(file_name) 
  newobject.list[[seob]] = object
}

saveRDS(newobject.list, paste0(output_path,"/seurat/2.object_list_rawGEXATAC.rds"))

PlotQC(newobject.list, prefix="before",outpath=plot_path,filter.doublets=FALSE)
