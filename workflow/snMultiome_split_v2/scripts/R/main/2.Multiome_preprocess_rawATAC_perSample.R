args=commandArgs(trailingOnly=T)
input_path = as.character(args[1])
working_dir = as.character(args[2])
project_name = as.character(args[3])
sample = as.character(args[4])
#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/Analysis_filtGEXATAC_split"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/Analysis_filtGEXATAC_split/cellranger_output"
#project_name = "ccdi4"
#sample = "CCDI0269_TCU00389-0101-Rep2"
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
for (s in samplelist) {
  subfolder_path <- file.path(sample_plot_path, s)
  dir.create(subfolder_path, recursive = TRUE, showWarnings = FALSE)
}

filtered.input.list=list()
filtered.path <- paste0(input_path, "/",sample, '/outs/filtered_feature_bc_matrix.h5')
filtered.input.list[[sample]] <- Read10X_h5(filtered.path)
metadata.path <- paste0(input_path, "/",sample, '/outs/per_barcode_metrics.csv')
obj_name = paste0(sample,"_GEXobj_raw.rds")
#print(obj_name)
file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
#print(sample_outdir_list[[seob]])
object = readRDS(file_name) 

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

frag.files <- paste0(input_path, "/",sample, '/outs/atac_fragments.tsv.gz')

dataset = sample
#ccdi = sampleInfo$CCDI[sampleInfo$uniqueID == sample]
#object <- RenameCells(object, add.cell.id = ccdi)
print(object)
atac_counts <- filtered.input.list[[dataset]]$Peaks
  # only keep the count in standard genomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)  # we'll only use peaks in standard chromosomes
atac_counts <- atac_counts[as.vector(grange.use), ]
dim(atac_counts)
metafile = metadata.path
metadata=read.csv(metafile,header = TRUE,row.names = 1)
print("now creat chromatin assay")
frag.file = frag.files
chrom_assay =CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = 'hg38',
        fragments = frag.file,
        min.cells=0,
        annotation = annotations
      )
print(chrom_assay)
object[["ATAC"]] <- chrom_assay # %>% 
DefaultAssay(object) ="ATAC"
  #Nucleosome signal score per cell
object = NucleosomeSignal(object)
  #TSS Enrichment score per cell
object = TSSEnrichment(object,fast=FALSE)
  #Blacklist ratio and fraction reads in peaks
object$blacklist_fraction = FractionCountsInRegion(object = object,assay = 'ATAC',regions = blacklist_hg38)
p1=DensityScatter(object, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
  #TSS enrihcment plot. Only works if TSS calculation above has fast=FALSE
object$tss_group =ifelse(object$TSS.enrichment>2,">2","<2")
p2=TSSPlot(object,group.by="tss_group")
  
  #Fragment length periodicity
object$nucleosome_group = ifelse(object$nucleosome_signal > 2, "NS>2","NS<=2")
p3=FragmentHistogram(object,group.by="nucleosome_group")
  # FragmentHistogram(object,group.by="tss_group")
  
  
  #calculate fraction of reads in peaks
total_fragments <- CountFragments(frag.file)
rownames(total_fragments) <- total_fragments$CB
object$fragments <- total_fragments[colnames(object), "frequency_count"]*2
  
object <- FRiP(
    object = object,
    assay = 'ATAC',
    total.fragments = 'fragments',
    col.name = 'pct_reads_in_peaks'
  )
object$FRiP_group =ifelse(object$pct_reads_in_peaks>0.15,">0.15","<=0.15")
  #VlnPlot(object = object,features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks'),pt.size = 0.0,ncol = 5)
  
p4= VlnPlot(
    object = object,
    features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks',
                 'nCount_RNA','nFeature_RNA','percent.mt','percent.rb'),
    pt.size = 0,
    ncol = 6,
    same.y.lims = FALSE
  ) &
  theme(
      plot.title = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )

library(cowplot)

title_plot <- ggdraw() +
  draw_label(sample, size = 14, x = 0.5, y = 0.5,fontface = 'bold')

#title_plot <- ggplot() +
#  ggtitle(sample) +  # Add your desired title here
#  theme_void()
combined_p4 <- title_plot + p4 +plot_layout(ncol = 1, heights = c(0.05, 0.95))

p5 <- FeatureScatter( object = object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme(legend.position = "none")
p6 <- FeatureScatter( object = object, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC")+
  theme(legend.position = "none")
p7 <- FeatureScatter( object = object, feature1 = "nCount_RNA", feature2 = "percent.mt")+
  theme(legend.position = "none")
p8 <- FeatureScatter( object = object, feature1 = "nCount_ATAC", feature2 = "percent.mt")+
  theme(legend.position = "none")

combined1 = p5 + p7 +p6 + p8 + plot_layout(ncol = 2)
combined_scatter = title_plot + combined1 + plot_layout(ncol = 1, heights = c(0.05, 0.95))

  #Set up outliers upper filter on peaks (3 median absolute deviations)
#peaksUpper=outliers_mad(object$nCount_ATAC,threshold=3)$limits[2]
#peaksLower = 500
#cat(dataset,"ncount_ATAC threshold:",peaksLower,"-", peaksUpper,"\n" )
#object$nCountATAC_group = "Normal"
#object$nCountATAC_group[object$nCount_ATAC <= peaksLower] = "Low"
#object$nCountATAC_group[object$nCount_ATAC >= peaksUpper] = "Outlier"
  
  # save the raw seurat object
  #obj_name = paste0(dataset,"_GEXATACobj_raw.rds")
  #saveRDS(object,paste0(sample_outdir_list[[dataset]],"/SeuratObj/",obj_name))
# newobject.list[[dataset]] = object
  
if (!is.null(p1)) {
    # Open a PDF device
    path_now = paste0(sample_plot_path,"/",dataset)
    pdf(paste0(path_now,"/1.qc_before_filtering.pdf"),width=8,height = 6)
    print(p1+ggtitle(dataset))
    print(p2+ggtitle(dataset))
    print(p3+ggtitle(dataset))
    dev.off()
    pdf(paste0(path_now,"/1.qc_before_filtering_violinPlot.pdf"),width=12,height = 8)
    print(combined_p4)
    dev.off()
    pdf(paste0(path_now,"/1.qc_before_filtering_scatterPlot.pdf"),width=8,height = 8)
    print(combined_scatter)
    dev.off()
} else {
    cat("Error: Plot not created successfully.\n")
}
  

obj_name = paste0(sample,"_GEXATACobj_raw.rds")
  #print(obj_name)
file_name = paste0(output_path,"/seurat/preprocessed/",obj_name)
  #print(sample_outdir_list[[seob]])
saveRDS(object,file_name) 
