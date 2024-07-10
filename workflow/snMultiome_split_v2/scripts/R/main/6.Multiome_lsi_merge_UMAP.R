library(optparse)
# Create an option parser
option_parser <- OptionParser(usage = "Usage: %prog --arg1 <value1> --arg2 <value2>")

# Define the options
option_parser <- add_option(option_parser, "--analysis_dir", type = "character", default = NULL,
                            help = "This is the current working dir")
#option_parser <- add_option(option_parser, "--species", type = "character", default = NULL,
#                            help = "species")
option_parser <- add_option(option_parser, "--npcs", type = "integer", default = 50,
                            help = "number of pcs")
option_parser <- add_option(option_parser, "--threads", type = "integer", default = 8,
                            help = "number of threads")
#option_parser <- add_option(option_parser, "--kanchor", type = "integer", default = 5,
#                            help = "ATAC kanchor used for FindIntegrationAnchors")
#option_parser <- add_option(option_parser, "--kfilter", type = "integer", default = 200,
#                            help = "ATAC kfilter used for FindIntegrationAnchors")
#option_parser <- add_option(option_parser, "--kweight", type = "integer", default = 100,
#                            help = "ATAC kweight used for IntegrateEmbeddings")
#option_parser <- add_option(option_parser, "--res", type = "numeric", default = 0.8,
#                            help = "Description of argument 1")
opt <- parse_args(option_parser)

working_dir = opt$analysis_dir
#species= opt$species
npcs_val = opt$npcs
threads = opt$threads
#resolution = opt$res
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/Analysis_filtGEX"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
#species= "hg38"
#npcs_val = 50

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
integrated_object<- readRDS(paste0(integration_path,"/6.integrated_SeuratObj_RLSI.rds"))
DefaultAssay(integrated_object) = "MACS"

####plot merged umap for ATAC
so_merge <- RunUMAP(so_merge, dims = 2:npcs_val, reduction = 'lsi')
p2<- DimPlot(so_merge, group.by = 'SampleID', split = "SampleID",pt.size = 0.1)
p1<- DimPlot(so_merge, group.by = 'SampleID',pt.size = 0.1)
p3 <-CoveragePlot(
  object = so_merge,
  group.by = 'SampleID',
  region = "chr14-100725000-100739000"
)

p5 <- DimPlot(integrated_object, group.by = 'SampleID', split = "SampleID",pt.size = 0.1)
p4 <- DimPlot(integrated_object, group.by = 'SampleID',pt.size = 0.1)


if (!is.null(p1)) {
  # Open a PDF device
  pdf(paste0(integration_plot_path,"/3.merged_umap_ATAC.pdf"),width=12, height=8)
  print(p1+ggtitle("Merged"))
  dev.off()
  pdf(paste0(integration_plot_path,"/4.integrated_umap_ATAC.pdf"),width=10, height=8)
  print(p4+ggtitle("Integrated"))
  dev.off()
  pdf(paste0(integration_plot_path,"/3.merged_umap_split_ATAC.pdf"),width=32, height=8)
  print(p2+ggtitle("Merged & Splited by sample") & theme(plot.title = element_text(size = 15)))
  dev.off()
  pdf(paste0(integration_plot_path,"/4.integrated_umap_split_ATAC.pdf"),width=28, height=8)
  print(p5+ggtitle("Integrated & Splited by sample") & theme(plot.title = element_text(size = 15)))
  dev.off()
#  pdf(paste0(integration_plot_path,"/2.peak_example_ATAC.pdf"),width=12, height=8)
#  print(p3+ggtitle(sample))
#  dev.off()
} else {
  cat("Error: Plot not created successfully.\n")
}



