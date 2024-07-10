args=commandArgs(trailingOnly=T)
input_path = as.character(args[1])
working_dir = as.character(args[2])
species= as.character(args[3])
npcs_val = as.character(args[4])
working_dir = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/Snake"
input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
species= "hg38"
npcs_val = 50
resolutions = 0.2

cat(species,npcs_val,"----\n")
cat(input_path,working_dir,"----\n")
setwd(working_dir)
source('scripts/R/main/load_packages.R')
plan('multicore', workers=10)
options(future.globals.maxSize=50 * 1000 * 1024^2)

output_path = paste0(working_dir,"/result_wSoupX")
plot_path = paste0(output_path,"/integration_plot")
integration_path = paste0(output_path,"/integration")

sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet))
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist
subfolder=samplelist
subfolder=samplelist
for (sub in subfolder){
  dir.create(paste0(scevanOut_path,"/",sub))
}

setwd(output_path)
#object2 <- readRDS(paste0(integration_path,"/6.integrated_SeuratObj_lsi_new.rds"))
object_filt.list = readRDS(paste0(output_path,"/3.object_list_filtGEXATAC_consnsus.rds"))
object.wnn <- readRDS(paste0(integration_path,"/7.integrated_SeuratObj_wnn_anno0.2.rds"))


####DimPlot with differnt resolution
cell.size = 0.2
resolution = 0.2
object = object.wnn


#ref_data =c("clustAnnot_HPCA_main","clustAnnot_HPCA","clustAnnot_BP_encode_main","clustAnnot_BP_encode","clustAnnot_monaco_main","clustAnnot_monaco", 
            "clustAnnot_immu_cell_exp_main","clustAnnot_immu_cell_exp") 
ref_data =c("clustAnnot_BP_encode_main","clustAnnot_BP_encode")
group = "wsnn_res.0.2"
#plot umap of RNA, ATAC, WNN
library(cowplot)
for (group in ref_data) {
  p1 <- DimPlot(object.wnn, reduction = "umap.rna", group.by = group, label = TRUE, label.size = 6,pt.size = 0.2, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(object.wnn, reduction = "umap.atac", group.by = group, label = TRUE, label.size = 6,pt.size = 0.2, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = group, label = TRUE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
  # Create a combined legend for p1, p2, and p3
  combined_legend <- cowplot::get_legend(p1 + p2 + p3)
  # Remove individual legends from p1, p2, and p3
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- p3 + theme(legend.position = "none")
  # Combine plots along with the combined legend
  combined_plot <- cowplot::plot_grid(
    p1,
    p2,
    p3,
    combined_legend, # Placeholder for the combined legend
    ncol = 4,
    rel_widths = c(8, 8, 8, 8) # Set relative widths so the legend occupies less space
  )
  combined_plot = combined_plot + theme(legend.text = element_text(size = 14)) +
     theme(axis.title = element_text(size = 16))
  combined_plot
  # Save the combined plot with the combined legend
  ggsave(
    plot = combined_plot,
    height = 8,
    width = 32,
    filename = paste0(plot_path, "/3.integrated_Multiome_3UMAP_", group, "_res0.2.pdf")
  )
}    

#plot cluster and annotation together resolution 0.2,0.5,0.8
group = "clustAnnot_BP_encode"
DimPlot(object.wnn, reduction = "wnn.umap", group.by = "wsnn_res.0.1", label = TRUE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = group, label = TRUE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
p5 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
p6 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
p7 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = "wsnn_res.0.2", label = FALSE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
combined_plot2 = p7 + p5 + p6 & theme(plot.title = element_text(hjust = 0.5))
ggsave(
  plot = combined_plot2,
  height = 8,
  width = 24,
  filename = paste0(plot_path, "/3.integrated_Multiome_UMAP_res0.20.50.8.pdf")
)
combined_plot3 = p4 + p7 & theme(plot.title = element_text(hjust = 0.5))
ggsave(
  plot = combined_plot3,
  height = 8,
  width = 16,
  filename = paste0(plot_path, "/3.integrated_Multiome_UMAP_res0.2.pdf")
)

p8 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = group,split="SampleID", label = TRUE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
p9 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = "wsnn_res.0.2",split="SampleID",label = TRUE, label.size = 6, pt.size = 0.2,repel = TRUE) + ggtitle("WNN")
combined_plot3 <- plot_grid(
  p8 + theme(plot.title = element_text(hjust = 0.5)), 
  p9 + theme(plot.title = element_text(hjust = 0.5)),
  ncol = 1
)
ggsave(
  plot = combined_plot3,
  height =16,  # Adjust height as needed
  width = 24,    # Adjust width as needed
  filename = paste0(plot_path, "/3.integrated_Multiome_UMAP_res0.2_splitbySample.pdf")
)

### plot cell cycle
options(ggrepel.max.overlaps = Inf) 
p10 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = "Phase", split.by = "SampleID",label = TRUE, label.size = 6,pt.size=0.2, repel = TRUE) + ggtitle("WNN")
p10  & theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p10,height =16,width =20, filename = paste0(plot_path, "/3.integrated_Multiome_wnnUMAP_cellcycle_res0.2_ncol2.pdf")
)
p11 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = "Phase", label = TRUE, label.size = 6,pt.size=0.2, repel = TRUE) + ggtitle("WNN")
p11  & theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = p11,height =8,width =8, filename = paste0(plot_path, "/3.integrated_Multiome_wnnUMAP_cellcycle_res0.2_one.pdf")
)
p12 <- DimPlot(object.wnn, 
               reduction = "wnn.umap", 
               group.by = "Phase", 
               split="SampleID",
               label = TRUE, 
               label.size = 6, 
               pt.size = 0.2, 
               repel = TRUE) +
  ggtitle("WNN") +
  theme(plot.title = element_text(hjust = 0.5)) #+
  facet_wrap(~ sampleID,ncol=2)

p12 <- DimPlot(object.wnn, reduction = "wnn.umap", group.by = "Phase", 
                        split.by = "SampleID",  # Split the plot by SampleID
                        label = TRUE, label.size = 3,  # Adjust label size
                        pt.size = 0.1, repel = TRUE) +
    ggtitle("Cell Cycle") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ SampleID, ncol = 2, scales = "free_x")  # Specify facet wrapping
p12
ggsave(plot = p12,height =8,width =8, filename = paste0(plot_path, "/3.integrated_Multiome_wnnUMAP_cellcycle_res0.2_ncol2.pdf")
)

####add module score in FeaturePlot
DefaultAssay(object.wnn)="RNA"
genes <- list (
  Fibroblast = c("COL1A1","COL1A2","COL3A1","ACTA2"),
  TcellNKcell = c("CD3D","CD3E","CD3G","CD247"),
  Bcell = c("CD79A","CD79B","IGHM","IGHD"),
  Thyrocytes = c("TG","EPCAM","KRT18","KRT19"),
  Myeloid = c("LYZ","S100A8","S100A9","CD14"),
  Endothelial = c("PECAM1","CD34","CDH5","VWF")
  )
  
colors <- c('lightgrey', 'navy')
object = object.wnn
DefaultAssay(object) = "RNA"
object_join = JoinLayers(object)
DefaultAssay(object_join) = "RNA"
DefaultLayer(object_join[["RNA"]]) <- 'data'
object_join <- object_join %>% AddModuleScore(features=genes,name=names(genes))
cell.size = 0.1
p <- object_join %>% 
    FeaturePlot(features=paste0(names(genes), seq(length(genes))), pt.size=cell.size, 
                min.cutoff='q10', max.cutoff='q90', ncol=3, cols=colors, reduction= "wnn.umap")
p
ggsave(plot=p, height=9, width=12, filename=paste0(plot_path, '/FeaturePlot_activity_module_umap_res0.2.pdf'))
  
#####more marker genes for Malignat cells
PTC_genes <- list (
  Malignant = c("S100A4","FN1","TMSB4X","APOE","CXCL14","TIMP1","LGALS3","IGFBP6",
                "APOC1","S100A10","TM4SF1","ZCCHC12","KRT19","NMU","SLC34A2","SERPINA1",
                "CST6","TACSTD2","LGALS1","S100A1","PRSS23","SLPI","S100A6","CXCL2","TMSB10"),
  Non_Malignant = c("TFF3","MT1G","TPO","MT1F","SLC26A7","TG","DIO2","SORBS2","MT1E","ID3","MT1X",
                    "SLC26A4","MT1H","ID4","IYD","CRABP1","GPX3","SLC26A4−AS1","FCGBP","PRDX1","GCSH",
                    "FHL1","SORD","MATN2","TXNL1")
)

DefaultLayer(object_join[["RNA"]]) <- 'data'
object_join <- object_join %>% AddModuleScore(features=PTC_genes,name=names(PTC_genes))
p2 <- object_join %>% 
  FeaturePlot(features=paste0(names(PTC_genes), seq(length(PTC_genes))),split="SampleID", pt.size=cell.size, 
              min.cutoff='q10', max.cutoff='q90', ncol=4, cols=colors, reduction= "wnn.umap")
p2
ggsave(plot=p2, height=9, width=16, filename=paste0(plot_path, '/FeaturePlot_malignant_module_umap_splitSample_res0.2.pdf'))

p3 <- object_join %>% 
  FeaturePlot(features=paste0(names(PTC_genes), seq(length(PTC_genes))), pt.size=cell.size, 
              min.cutoff='q10', max.cutoff='q90', ncol=3, cols=colors, reduction= "wnn.umap")
p3
ggsave(plot=p3, height=4.5, width=12, filename=paste0(plot_path, '/FeaturePlot_malignant_module_umap_res0.2.pdf'))





###ATAC plot##################################################################
DefaultAssay(object_wnn) ="MACS"
p1 <- CoveragePlot(
  object = object_wnn,
  region = "CDK8",
  features = "CDK8",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream =0,
  extend.downstream = 0
)
CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)

p3 <-CoveragePlot(
  object = object_wnn,
  group.by = 'SampleID',
  region = "chr14-100725000-100739000",
  expression.assay = "SCT",
  idents = idents.plot
)

####RNA plot
DefaultAssay(object_wnn) = "RNA"
FeaturePlot(object_wnn, features = markers, pt.size = 0.2,
            ncol = 4)
DefaultAssay(object_wnn) = "MACS"
pdf(paste0(integration_path,"/3.integrated_Multiome_3UMAP_",group,".pdf"),width=24, height=8)
p2 <- CoveragePlot(
  object = object_wnn,
  region = "MEG3",
  features = "MEG3",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)
dev.off()

patchwork::wrap_plots(p1, p2, ncol = 1)

p2 <- CoveragePlot(
  object = object_wnn,
  region = "CDK8",
  features = "CDK8",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)


