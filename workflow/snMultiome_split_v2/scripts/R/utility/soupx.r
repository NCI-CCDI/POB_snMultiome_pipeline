

SoupCorrect <- function(raw.matrix, filt.matrix, contamination_rate=NULL) {
  
  srat  <- CreateSeuratObject(counts=filt.matrix[['Gene Expression']])
  soup.channel  <- SoupChannel(raw.matrix[['Gene Expression']], filt.matrix[['Gene Expression']])
  
  srat <- srat %>% 
    NormalizeData() %>%
    FindVariableFeatures() %>% 
    ScaleData(verbose=FALSE) %>%
    RunPCA(verbose=FALSE) %>% 
    RunUMAP(dims=1:30, verbose=FALSE) %>% 
    FindNeighbors(dims=1:30, verbose=FALSE) %>% 
    FindClusters(verbose=FALSE, algorithm=3, resolution=0.8)
  
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  
  if (is.null(contamination_rate)) {
    soup.channel  <- autoEstCont(soup.channel, forceAccept=TRUE, doPlot=FALSE)
  } else {
    soup.channel <- setContaminationFraction(soup.channel, contamination_rate, forceAccept=TRUE)
  }

  adj.matrix  <- adjustCounts(soup.channel)

  return(adj.matrix)
  
}

  