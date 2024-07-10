args=commandArgs(trailingOnly=T)
input_path = as.character(args[1])
working_dir = as.character(args[2])
species= as.character(args[3])
npcs_val = as.character(args[4])
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/Snake"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
#species= "hg38"
#npcs_val = 50
resolution_list = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

cat(species,npcs_val,"----\n")
cat(input_path,working_dir,"----\n")
setwd(working_dir)
source('scripts/R/main/load_packages.R')
plan('multicore', workers=16)
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

#object2 <- readRDS(paste0(integration_path,"/6.integrated_SeuratObj_lsi_new.rds"))
object <- readRDS(paste0(integration_path,"/6.integrated_SeuratObj_RLSI.rds"))


##########integrate GEX and ATAC#############################
DefaultAssay(object) <- 'integratedrna'
VariableFeatures(object) <- rownames(object)

object <- object %>% 
  ScaleData(verbose=FALSE) %>% RunPCA(verbose=FALSE) %>% 
  FindNeighbors(reduction='pca', dims=1:npcs_val) 

object <- object %>% FindClusters(algorithm=3, resolution=0.2,cluster.name = "rna_clusters") 

object <- object  %>%  RunUMAP(reduction='pca', dims=1:npcs_val, reduction.name='umap.rna', reduction.key='rnaUMAP_')

object <- object %>% RunUMAP(reduction='integrated_lsi', dims=2:npcs_val, reduction.name='umap.atac', reduction.key='atacUMAP_')


#object_wnn <- object %>% 
#  FindMultiModalNeighbors(reduction.list=list('pca', 'integrated_lsi'), dims.list=list(1:50, 2:50)) %>% 
#  FindClusters(graph.name='wsnn', algorithm=3, resolution=resolutions,cluster.name="wnn_clusters") %>% 
#  
object_wnn <- object %>% FindMultiModalNeighbors(reduction.list=list('pca', 'integrated_lsi'), dims.list=list(1:50, 2:50))
for (res in resolution_list){
  object_wnn = FindClusters(object_wnn,graph.name='wsnn', algorithm=3, resolution=res)
}
object_wnn = FindClusters(object_wnn,graph.name='wsnn', algorithm=3, resolution=0.2,cluster.name="wnn_clusters") 
object_wnn = RunUMAP(object_wnn,nn.name='weighted.nn', reduction.name='wnn.umap', reduction.key='wnnUMAP_')

#object[['seurat_clusters']] <- NULL
object_wnn@meta.data[1:10,]
Idents(object_wnn) <- 'wnn_clusters'

RUN_SINGLEr_AVERAGE = function(obj,refFile,fineORmain){
  DefaultAssay(obj)="SCT"
  avg = AverageExpression(obj,assays = "SCT")
  avg = as.data.frame(avg)
  ref = refFile
  s = SingleR(test = as.matrix(avg),ref = ref,labels = ref[[fineORmain]])
  
  clustAnnot = s$labels
  names(clustAnnot) = colnames(avg)
  names(clustAnnot) = gsub("SCT.g","",names(clustAnnot))
  
  annotVect = clustAnnot[match(obj$wnn_clusters,names(clustAnnot))]
  names(annotVect) = colnames(obj)
  return(annotVect)
}

if(species == "hg38" || species == "hg19"){
  object_wnn$clustAnnot_HPCA_main <- RUN_SINGLEr_AVERAGE(object_wnn,celldex::HumanPrimaryCellAtlasData(),"label.main")
  object_wnn$clustAnnot_HPCA <-  RUN_SINGLEr_AVERAGE(object_wnn,celldex::HumanPrimaryCellAtlasData(),"label.fine")
  object_wnn$clustAnnot_BP_encode_main <-  RUN_SINGLEr_AVERAGE(object_wnn,celldex::BlueprintEncodeData(),"label.main")
  object_wnn$clustAnnot_BP_encode <-  RUN_SINGLEr_AVERAGE(object_wnn,celldex::BlueprintEncodeData(),"label.fine")
  object_wnn$clustAnnot_monaco_main <-  RUN_SINGLEr_AVERAGE(object_wnn,celldex::MonacoImmuneData(),"label.main")
  object_wnn$clustAnnot_monaco <- RUN_SINGLEr_AVERAGE(object_wnn,celldex::MonacoImmuneData(),"label.fine")
  object_wnn$clustAnnot_immu_cell_exp_main <-  RUN_SINGLEr_AVERAGE(object_wnn,celldex::DatabaseImmuneCellExpressionData(),"label.main")
  object_wnn$clustAnnot_immu_cell_exp <- RUN_SINGLEr_AVERAGE(object_wnn,celldex::DatabaseImmuneCellExpressionData(),"label.fine")
} else if(species == "mm10"){
  object_wnn$clustAnnot_immgen_main <-  RUN_SINGLEr_AVERAGE(object_wnn,celldex::ImmGenData(),"label.main")
  object_wnn$clustAnnot_immgen <- RUN_SINGLEr_AVERAGE(object_wnn,celldex::ImmGenData(),"label.fine")
  object_wnn$clustAnnot_mouseRNAseq_main <-  RUN_SINGLEr_AVERAGE(object_wnn,celldex::MouseRNAseqData(),"label.main")
  object_wnn$clustAnnot_mouseRNAseq <- RUN_SINGLEr_AVERAGE(object_wnn,celldex::MouseRNAseqData(),"label.fine")
}

saveRDS(object_wnn , paste0(integration_path,"/7.integrated_SeuratObj_WNN_anno0.2.rds"))
