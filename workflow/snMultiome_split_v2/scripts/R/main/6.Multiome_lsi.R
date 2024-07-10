library(optparse)
# Create an option parser
option_parser <- OptionParser(usage = "Usage: %prog --arg1 <value1> --arg2 <value2>")

# Define the options
option_parser <- add_option(option_parser, "--analysis_dir", type = "character", default = NULL,
                            help = "This is the current working dir")
option_parser <- add_option(option_parser, "--species", type = "character", default = NULL,
                            help = "species")
option_parser <- add_option(option_parser, "--npcs", type = "integer", default = 30,
                            help = "number of pcs")
option_parser <- add_option(option_parser, "--threads", type = "integer", default = 8,
                            help = "number of threads")
option_parser <- add_option(option_parser, "--kanchor", type = "integer", default = 5,
                            help = "ATAC kanchor used for FindIntegrationAnchors")
option_parser <- add_option(option_parser, "--kfilter", type = "integer", default = 200,
                            help = "ATAC kfilter used for FindIntegrationAnchors")
option_parser <- add_option(option_parser, "--kweight", type = "integer", default = 100,
                            help = "ATAC kweight used for IntegrateEmbeddings")
#option_parser <- add_option(option_parser, "--res", type = "numeric", default = 0.8,
#                            help = "Description of argument 1")
opt <- parse_args(option_parser)

working_dir = opt$analysis_dir
species= opt$species
npcs_val = opt$npcs
threads = opt$threads
kanchor = opt$kanchor
kfilter = opt$kfilter
kweight = opt$kweight
#resolution = opt$res
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/Analysis_filtGEX"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
#species= "hg38"
#npcs_val = 50

cat(species,npcs_val,kanchor,kfilter,"----\n")
cat(working_dir,"----\n")
setwd(working_dir)
source('scripts/R/main/load_packages.R')
plan('multicore', workers=threads)
options(future.globals.maxSize=100 * 1000 * 1024^2)
options(future.rng.onMisuse="ignore")

output_path = paste0(working_dir,"/result_wSoupX")
plot_path = paste0(output_path,"/qc_result_plot")
integration_path = paste0(output_path,"/integration")
integration_plot_path = paste0(output_path,"/integration_plot")

sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet))
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist

object <- readRDS(paste0(integration_path,"/5.integrated_SeuratObj_RPCA.rds")) 

DefaultAssay(object) <- 'MACS'
object.list <- object %>% RunTFIDF() %>% SplitObject(split.by='SampleID')

object.list <- sapply(object.list, FindTopFeatures, min.cutoff='q5', simplify=FALSE)
features <- Reduce(union, lapply(object.list, VariableFeatures))

object <- object %>% FindTopFeatures(min.cutoff='q0') %>% RunSVD(features=features)

so_merge = object

object.list <- sapply(object.list, RunSVD, features=features, simplify=FALSE)

anchors <- object.list %>% FindIntegrationAnchors(anchor.features=features, reduction='rlsi', k.anchor=kanchor, k.filter = kfilter, dims=2:npcs_val)
integrated_object <- IntegrateEmbeddings(anchorset=anchors, reductions=object[['lsi']], new.reduction.name='integrated_lsi', dims.to.integrate=1:npcs_val, k.weight=kweight)
integrated_object<- RunUMAP(integrated_object, reduction = "integrated_lsi", dims = 2:npcs_val)
saveRDS(integrated_object, paste0(integration_path,"/6.integrated_SeuratObj_RLSI.rds"))

