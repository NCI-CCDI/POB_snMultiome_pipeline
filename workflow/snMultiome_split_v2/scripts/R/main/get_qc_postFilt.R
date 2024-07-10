args=commandArgs(trailingOnly=T)
working_dir = as.character(args[1])

#working_dir = "/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/Snake"
#input_path = "/data/CCRCCDI/analysis/ccrccdi4/02_PrimaryAnalysisOutput"
cat(working_dir,"\n")
setwd(working_dir)
source('scripts/R/main/load_packages.R')

output_path = paste0(working_dir,"/result_wSoupX")
plot_path = paste0(output_path,"/qc_result_plot")
sample_plot_path =  paste0(output_path,"/sample_result_plot")

sampleSheet = "assets/input_manifest_cellranger.csv"
sampleInfo = read.csv(paste0(working_dir,"/",sampleSheet))
samplelist = as.list(sampleInfo$uniqueID)
names(samplelist) = samplelist

#print(sample_outdir_list)
object.list <- readRDS(paste0(output_path,"/seurat/3.object_list_GEXATAC_filt.rds"))
df <- data.frame("Sample" = character(), "Cell Number" = numeric(), "Total fragment" = numeric(), 
                      "Mean fragment per cell" = numeric(), "Median fragment per cell"= numeric(), 
                      "Peaks per cell"= numeric(), "Total Peaks"= numeric(),stringsAsFactors = FALSE)

for (s in samplelist) {
  sample = s
  object = object.list[[s]]
  DefaultAssay(object) = "RNA"
  total_cells <- length(object$nCount_RNA)
  total_reads <- sum(object$nCount_RNA)
  mean_reads_per_cell <- sum(object$nCount_RNA)/length(object$nCount_RNA)
  mean_reads_per_cell <- format(round(mean_reads_per_cell, 2), nsmall = 2)
  median_reads_per_cell <- median(object$nCount_RNA)
  median_genes_detected <- median(object$nFeature_RNA)
  dd = object@assays[["RNA"]]@layers[["counts"]]
  cell_per_gene=rowSums(dd>0)
  total_gene_detected = length(cell_per_gene[cell_per_gene>0])
  qcstat = c(sample, total_cells, total_reads, mean_reads_per_cell, median_reads_per_cell, median_genes_detected,total_gene_detected)
  cat(qcstat, "\n")
  df_now = rbind(qcstat=qcstat)
  colnames(df_now) = c("Sample", "Cell  Number", "Total UMI", "Mean UMI per cell", "Median UMI per cell", "Median Gene Detected", "Total Gene Detected")
  df = rbind(df,df_now)
}

#df
#write.table(df, file = paste(plot_path,"/GEX_QC_stats_preFilt.tsv",sep=""), quote=FALSE,sep="\t", row.names=FALSE)


df_ATAC <- data.frame("Sample" = character(), "Cell Number" = numeric(), "Total Fragments" = numeric(), 
                 "Mean Fragments per cell" = numeric(), "Median Fragments per cell"= numeric(), 
                 "Median Peak number"= numeric(),  "Total Peak number"= numeric(),stringsAsFactors = FALSE)

for (s in samplelist) {
  sample = s
  object = object.list[[s]]
  DefaultAssay(object) = "ATAC"
  total_cells <- length(object$nCount_ATAC)
  total_reads <- sum(object$nCount_ATAC)
  mean_reads_per_cell <- sum(object$nCount_ATAC)/length(object$nCount_ATAC)
  mean_reads_per_cell <- format(round(mean_reads_per_cell, 2), nsmall = 2)
  median_reads_per_cell <- median(object$nCount_ATAC)
  median_genes_detected <- median(object$nFeature_ATAC)
  dd = object@assays[["ATAC"]]@counts
  cell_per_gene=rowSums(dd>0)
  total_gene_detected = length(cell_per_gene[cell_per_gene>0])
  qcstat = c(sample, total_cells, total_reads, mean_reads_per_cell, median_reads_per_cell, median_genes_detected,total_gene_detected)
  cat(qcstat, "\n")
  df_now = rbind(qcstat=qcstat)
  colnames(df_now) = c("Sample", "Cell  Number", "Total Fragments", "Mean Fragments per cell", "Median Fragments per cell", "Median Peak number", "Total Peak number")
  df_ATAC = rbind(df_ATAC,df_now)
}
#df_ATAC
##rite.table(df, file = paste(plot_path,"/GEX_QC_stats_preFilt.tsv",sep=""), quote=FALSE,sep="\t", row.names=FALSE)

df_all = cbind(df,df_ATAC)
df_all
write.table(df_all, file = paste(plot_path,"/QC_stats_postFilt.tsv",sep=""), quote=FALSE,sep="\t", row.names=FALSE)

