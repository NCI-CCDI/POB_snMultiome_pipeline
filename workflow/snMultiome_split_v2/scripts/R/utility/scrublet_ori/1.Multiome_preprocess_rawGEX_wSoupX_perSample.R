Sys.setenv(RETICULATE_PYTHON = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/reticulate-env/bin/python")
library(optparse)
# Create an option parser
option_parser <- OptionParser(usage = "Usage: %prog --arg1 <value1> --arg2 <value2>")

# Define the options
option_parser <- add_option(option_parser, "--input_path", type = "character", default = NULL,
                            help = "This is the path for the untared cellranger output")
option_parser <- add_option(option_parser, "--analysis_dir", type = "character", default = NULL,
                            help = "Current workding dir")
option_parser <- add_option(option_parser, "--project_name", type = "character", default = NULL,
                            help = "Project name")
option_parser <- add_option(option_parser, "--sample", type = "character", default = NULL,
                            help = "Sample")
option_parser <- add_option(option_parser, "--nCount_RNA_LL", type = "integer", default = 500,
                            help = "Lower Limit for nCount_RNA")
option_parser <- add_option(option_parser, "--nFeature_RNA_LL", type = "integer", default = 300,
                            help = "Lower Limit for nFeature_RNA")
option_parser <- add_option(option_parser, "--mt_UL", type = "integer", default = 10,
                            help = "Upper Limit for mt percentage")
option_parser <- add_option(option_parser, "--res", type = "numeric", default = 0.8,
                            help = "Selection of resolution for UMAP")
opt <- parse_args(option_parser)

working_dir = opt$analysis_dir
input_path = opt$input_path
project_name = opt$project_name
sample = opt$sample
nCount_RNA_LL = opt$nCount_RNA_LL
nFeature_RNA_LL = opt$nFeature_RNA_LL
mt_UL = opt$mt_UL
resolution = opt$res

#working_dir = "/data/CCRCCDI/analysis/ccrccdi5/analysis_YW"
#input_path = "/data/CCRCCDI/analysis/ccrccdi5/analysis_YW/cellranger_output"
#project_name="ccrccdi5"
#sample = "CCDI0054_1535"
print(paste0("inputpath: ",input_path))
print(paste0("working_dir: ",working_dir))
print(paste0("project_name: ",project_name))
print(paste0("nCount_RNA_LL:",nCount_RNA_LL))
print(paste0("nFeature_RNA_LL:",nFeature_RNA_LL))
print(paste0("mt_UL:",mt_UL))
print(paste0("sample:",sample))
#cat(input_path,working_dir,"\n")
setwd(working_dir)
set.seed(123)
source('scripts/R/main/load_packages.R')
script_path = paste0(working_dir,"/scripts/R/utility")
output_path = paste0(working_dir,"/result_wSoupX")
subfolder_names = c("seurat","integration","qc_result_plot","integration_plot","seurat/preprocessed","seurat/merged","sample_result_plot")
for (subfolder in subfolder_names) {
  subfolder_path <- file.path(output_path, subfolder)
  if (file.exists(subfolder_path)) {
    cat("Subfolder", subfolder, "created successfully before.\n")
  } else {
    cat("The subfolder does not exitst and create now:,", subfolder, "\n")
    dir.create(subfolder_path, recursive = TRUE, showWarnings = FALSE)

  }
}
sample_result_path = paste0(output_path,"/sample_result_plot/",sample)
dir.create(sample_result_path, recursive = TRUE, showWarnings = FALSE)

dir_list <- list.dirs(input_path, recursive = FALSE,full.names = FALSE)
sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet)) 
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist

# 0. the 10x hdf5 file contains both data types.
raw.input.list=list()
filtered.input.list = list()
raw.path <- paste0(input_path, "/", sample, '/outs/raw_feature_bc_matrix.h5')
filtered.path <- paste0(input_path, "/", sample, '/outs/filtered_feature_bc_matrix.h5')
raw.input.list[[sample]] <- Read10X_h5(raw.path)
filtered.input.list[[sample]] <- Read10X_h5(filtered.path)

soup_rate <- 0.20
env_folder = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/reticulate-env"
use_virtualenv(env_folder)
reticulate::source_python(paste0(script_path,'/scrublet_py_YW.py'))

dataset = sample
#samplelist
# 1. create object list
adj.matrix <- suppressWarnings(SoupCorrect(raw.input.list[[dataset]], filtered.input.list[[dataset]], contamination_rate=soup_rate))
object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=dataset)
adj.matrix_sparse  = as(adj.matrix, "CsparseMatrix")
write10xCounts(paste0(
  sample_result_path,"/10x_wSoupX"),adj.matrix_sparse,barcodes=colnames(adj.matrix),gene.id=rownames(adj.matrix))
object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=dataset)
rna_counts <- filtered.input.list[[dataset]]$`Gene Expression`
cat("The original total UMI for", dataset, sum(rna_counts),"\n")
cat("The after-SoupX total UMI is", sum(adj.matrix),"\n")
 

# Create Seurat object
object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')
object_filter_min =  object %>% subset(
  subset=
    nCount_RNA > nCount_RNA_LL  &
    nFeature_RNA > nFeature_RNA_LL & 
    percent.mt < mt_UL)
doublet_rate <- (ncol(object_filter_min) / 1000) * 0.008
cat("Doublet rate was set to: ",doublet_rate)
sampleID_labels <- sample(sample, size=ncol(object),replace=TRUE)
object$SampleID = sampleID_labels
for (i in 2:(ncol(sampleInfo)-1)){
  label_name = colnames(sampleInfo)[i]
  label_now = as.character(sampleInfo[sampleInfo$uniqueID==sample,label_name])
  meta_label = rep(label_now, each=ncol(object))
  object[[label_name]] = meta_label
}
object[['doublet_rate']] <- rep(doublet_rate,each=ncol(object))

#source("/data/CCRCCDI/software/Scripts/current_snakemake/snMultiome_split_v2/scripts/R/utility/scrublet_YW.r")
#object <- object %>% 
#  scrublet_YW(n_prin_comps=30, expected_doublet_rate=doublet_rate,sample=sample,path=working_dir) 
object <- object %>%
  doublet_Finder(n_prin_comps=30,expected_doublet_rate=doublet_rate,sample = sample, path=sample_result_path,umi_min =nCount_RNA_LL,ngene_min=nFeature_RNA_LL,mt_max=mt_UL,resolution=resolution)


object <- object %>% 
  scrublet_YW(n_prin_comps=30, expected_doublet_rate=doublet_rate,sample=sample,path=sample_result_path,umi_min =nCount_RNA_LL,ngene_min=nFeature_RNA_LL,mt_max=mt_UL) 

obj_name = paste0(sample,"_GEXobj_raw.rds")
print(obj_name)
file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
print(file_name)
saveRDS(object,file_name) 

