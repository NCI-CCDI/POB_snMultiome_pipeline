library(optparse)
# Create an option parser
option_parser <- OptionParser(usage = "Usage: %prog --arg1 <value1> --arg2 <value2>")

# Define the options
option_parser <- add_option(option_parser, "--analysis_dir", type = "character", default = NULL,
                            help = "Current workding dir")
#option_parser <- add_option(option_parser, "--species", type = "character", default = NULL,
#                            help = "species")
option_parser <- add_option(option_parser, "--npcs", type = "integer", default = 30,
                            help = "Number of pcs")
option_parser <- add_option(option_parser, "--threads", type = "integer", default = 8,
                            help = "Number of threads")
option_parser <- add_option(option_parser, "--res", type = "numeric", default = 0.8,
                            help = "Selection of resolution for UMAP")
#option_parser <- add_option(option_parser, "--kanchor", type = "integer", default = 5,
#                            help = "GEX kanchor used for FindIntegrationAnchors")
#option_parser <- add_option(option_parser, "--kfilter", type = "integer", default = 200,
#                            help = "GEX kfilter used for FindIntegrationAnchors")
#option_parser <- add_option(option_parser, "--kweight", type = "integer", default = 100,
#                            help = "GEX kweight used for IntegrateEmbeddings")
opt <- parse_args(option_parser)

working_dir = opt$analysis_dir
species= opt$species
npcs_val = opt$npcs
threads = opt$threads
resolution = opt$res
#kanchor = opt$kanchor
#kfilter = opt$kfilter
#kweight = opt$kweight
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/Analysis_filtGEXATAC"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
#species= "hg38"
#npcs_val = 50

cat(working_dir,"\n",npcs_val,"\n",threads,"\n")

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


object.list = readRDS(paste0(output_path,"/seurat/4.object_list_filtGEXATAC_preprocess_addcellcycle.rds"))
prefix <- sampleInfo$CCDI
prefix
objGEX.merge <- merge(x = object.list[[1]], # first
                      y = object.list[-1], # other
                      add.cell.ids = prefix) # cell id as prefix


###get UMAP before and after inegration
noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA', 'doublet_scores', 'G2M.Score', 'S.Score')
DefaultAssay(objGEX.merge) = "RNA"
so_merged <- SCTransform(objGEX.merge, vars.to.regress = noise, assay = "RNA")
so_merged <- RunPCA(so_merged, npcs = npcs_val, verbose = F)
#so_merged <- RunTSNE(so_merged, reduction = "pca", dims = 1:30, perplexity = 30, max_iter = 1000,
#                     theta = 0.5, eta = 200, num_threads = 0)
#so_merged <- RunUMAP(so_merged, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30,
#                     n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
so_merged <- RunUMAP(so_merged, dims = 1:npcs_val, reduction = "pca", reduction.name = "umap.unintegrated")
so_merged <- FindNeighbors(so_merged, dims = 1:npcs_val, reduction = "pca")
so_merged <- FindClusters(so_merged, resolution = resolution, cluster.name = "unintegrated_clusters")

pdf(paste0 (integration_plot_path,"/1.merged_umap_GEX.pdf"),width=8, height= 6)
DimPlot(so_merged, reduction="umap.unintegrated",group.by = c("SampleID"))
dev.off()
pdf(paste0 (integration_plot_path,"/1.merged_umap_GEX_cluster_",resolution,".pdf"),width=8, height= 6)
DimPlot(so_merged, reduction="umap.unintegrated",group.by = c("unintegrated_clusters"))
dev.off()
pdf(paste0(integration_plot_path,"/1.merged_umap_GEX_split_by_sample.pdf"),width=24, height= 6)
DimPlot(so_merged, reduction="umap.unintegrated",group.by = c("unintegrated_clusters"), split.by=c("SampleID")) &
  theme(plot.title = element_text(size = 15),axis.text = element_text(size = 12))
dev.off()

so_integration = readRDS(paste0(integration_path,"/5.integrated_SeuratObj_RPCA.rds"))
DefaultAssay(so_integration) = "integratedrna"
so_integration <- ScaleData(so_integration,verbose=FALSE)
so_integration <- RunPCA(so_integration, npcs = npcs_val, verbose = F)
#so_merged <- RunTSNE(so_merged, reduction = "pca", dims = 1:30, perplexity = 30, max_iter = 1000,
#                     theta = 0.5, eta = 200, num_threads = 0)
#so_merged <- RunUMAP(so_merged, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30,
#                     n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
so_integration <- RunUMAP(so_integration, dims = 1:npcs_val, reduction = "pca", reduction.name = "umap.integrated")
so_integration <- FindNeighbors(so_integration, dims = 1:npcs_val, reduction = "pca")
so_integration <- FindClusters(so_integration, resolution = resolution, cluster.name = "integrated_clusters")

pdf(paste0 (integration_plot_path,"/2.integrated_umap_GEX.pdf"),width=8, height= 6)
DimPlot(so_integration, reduction="umap.integrated",group.by = c("SampleID"))
dev.off()
pdf(paste0 (integration_plot_path,"/2.integrated_umap_GEX_cluster_",resolution,".pdf"),width=8, height= 6)
DimPlot(so_integration, reduction="umap.integrated",group.by = c("integrated_clusters"))
dev.off()
pdf(paste0(integration_plot_path,"/2.integrated_umap_GEX_split_by_sample.pdf"),width=24, height= 6)
DimPlot(so_integration, reduction="umap.integrated",group.by = c("integrated_clusters"), split.by=c("SampleID")) &
  theme(plot.title = element_text(size = 15),axis.text = element_text(size = 12))
dev.off()


