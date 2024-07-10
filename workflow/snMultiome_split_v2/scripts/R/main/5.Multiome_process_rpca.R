library(optparse)
# Create an option parser
option_parser <- OptionParser(usage = "Usage: %prog --arg1 <value1> --arg2 <value2>")

# Define the options
option_parser <- add_option(option_parser, "--analysis_dir", type = "character", default = NULL,
                            help = "Current workding dir")
option_parser <- add_option(option_parser, "--species", type = "character", default = NULL,
                            help = "species")
option_parser <- add_option(option_parser, "--npcs", type = "integer", default = 30,
                            help = "Number of pcs")
option_parser <- add_option(option_parser, "--threads", type = "integer", default = 8,
                            help = "Number of threads")
option_parser <- add_option(option_parser, "--res", type = "numeric", default = 0.8,
                            help = "Selection of resolution for UMAP")
option_parser <- add_option(option_parser, "--kanchor", type = "integer", default = 5,
                            help = "GEX kanchor used for FindIntegrationAnchors")
option_parser <- add_option(option_parser, "--kfilter", type = "integer", default = 200,
                            help = "GEX kfilter used for FindIntegrationAnchors")
option_parser <- add_option(option_parser, "--kweight", type = "integer", default = 100,
                            help = "GEX kweight used for IntegrateEmbeddings")
opt <- parse_args(option_parser)

working_dir = opt$analysis_dir
species= opt$species
npcs_val = opt$npcs
threads = opt$threads
resolution = opt$res
kanchor = opt$kanchor
kfilter = opt$kfilter
kweight = opt$kweight
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/Analysis_filtGEXATAC"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
#species= "hg38"
#npcs_val = 50

cat(working_dir,"\n",species,"\n",npcs_val,"\n",threads,"\n")
cat(species,npcs_val,"\n"," kanchor:",kanchor," kfilter:",kfilter," kweight:",kweight,"----\n")
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

#object.list <- readRDS("/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/New_code2/result_sodo/3.object_list_filtGEXATAC_consnsus.rds")
object.list <- readRDS(paste0(output_path,"/seurat/3.object_list_filtGEXATAC_consnsus.rds"))

noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA', 'doublet_scores', 'G2M.Score', 'S.Score')

if(species=="hg38" || species == "hg19"){
  print("--proccesing human data")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
} else if(species=="mm10"){
  print("--proccesing mouse data")
  s.genes <- CONVERT_TO_HUMAN_GENELIST(cc.genes$s.genes)
  g2m.genes <- CONVERT_TO_HUMAN_GENELIST(cc.genes$g2m.genes)
}

object.list <- sapply(object.list, function(object) {
  DefaultAssay(object) <- 'RNA'
  genes <- rownames(object)
  
  object %>% 
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA") %>%
    CleanVarGenes(nHVG=5000)  %>% 
    ScaleData(features=genes,verbose=FALSE) %>%
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes,  set.ident = TRUE) %>%
    SCTransform(vars.to.regress=noise)
  
}, simplify=FALSE)

for (s in samplelist){
  object_now = object.list[[s]]
  object.list[[s]]$Sample = sampleInfo$CCDI[sampleInfo$uniqueID==s]
  obj_name = paste0(s,"_filtGEXATAC_preprocess_addcellcycle.rds")
  file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
  saveRDS(object_now, file_name)
}

saveRDS(object.list, paste0(output_path,"/seurat/4.object_list_filtGEXATAC_preprocess_addcellcycle.rds"))

#object.list = readRDS(paste0(output_path,"/seurat/4.object_list_filtGEXATAC_preprocess_addcellcycle.rds"))
features <- object.list %>% SelectIntegrationFeatures(nfeatures=2000)
writeLines(features, paste0(plot_path,"/selected_variable_features_2000_for_integration.txt"))

object.list <- sapply(object.list, function(object) {
  DefaultAssay(object) <- 'RNA'
  object %>% 
    ScaleData(features=features, vars.to.regress=noise, verbose=FALSE) %>%
    RunPCA(features=features, verbose=FALSE)
})

for (s in samplelist){
  object_now = object.list[[s]]
  obj_name = paste0(s,"_filtGEXATAC_preprocess_PCA.rds")
  file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
  saveRDS(object_now, file_name)
}

saveRDS(object.list, paste0(output_path,"/seurat/4.object_list_filtGEXATAC_preprocess_PCA.rds"))
#anchors <- object.list %>% FindIntegrationAnchors(anchor.features=features, reduction='rpca', k.anchor=7, scale=FALSE, dims=1:50)

anchors <- object.list %>% FindIntegrationAnchors(anchor.features=features, reduction='rpca', k.anchor=kanchor, scale=FALSE, dims=1:npcs_val)
object <- IntegrateData(anchorset=anchors, new.assay.name='integratedrna', k.weight=kweight, dims=1:npcs_val)
saveRDS(object, paste0(integration_path,"/5.integrated_SeuratObj_RPCA.rds"))

