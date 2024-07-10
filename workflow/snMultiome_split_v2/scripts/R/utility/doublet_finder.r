doublet_Finder <-function(dfso, n_prin_comps=30,expected_doublet_rate=0.06, sample="CCDI", path = "/data/CCRCCDI", umi_min =0,ngene_min=0,mt_max=100,resolution=0.8){
  library(DoubletFinder)
  #dfso=object
  #umi_min =500
  #ngene_min=300
  #mt_max=10
  #resolution = 0.2
  #npcs = 30
  #expected_doublet_rate=doublet_rate
  dfso_ori = dfso
  npcs = n_prin_comps
  cat("doublet_rate:",expected_doublet_rate,", min umi:",umi_min,", min gene:",ngene_min,", mt:", mt_max,", res:", resolution,",npcs:",npcs,"\n")
  #doublet_rate = 0.008 * ncol(dfso_ori)/1000


  dfso <- subset(dfso, subset = nFeature_RNA >ngene_min  & percent.mt < mt_max & nCount_RNA > umi_min) 
  noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA')
  dfso <- SCTransform(dfso,vars.to.regress=noise,verbose = FALSE)
  dfso <- RunPCA(dfso, pc.genes = dfso@var.genes,npcs =npcs )
  dfso  <- FindNeighbors(object = dfso , dims = 1:npcs,  reduction = "pca", verbose = FALSE)
  dfso  <- FindClusters(object = dfso, resolution = resolution, verbose = FALSE) 
  dfso <- RunUMAP(dfso,dims=1:npcs,reduction="pca")
  
  sweep.res.list_kidney <- paramSweep(dfso,PCs = 1:npcs, sct = T)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK_value <- as.numeric(as.character(bcmvn_kidney$pK[bcmvn_kidney$BCmetric == max(bcmvn_kidney$BCmetric)])) #get optimal pK
  
  
  pdf(paste0(path,"/",sample,"_DoubletFinder_pKvalue.pdf"))
  par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  plot(x = as.numeric(as.character(bcmvn_kidney$pK)), y = bcmvn_kidney$BCmetric, pch = 16,type="b", col = "blue",lty=1)
  abline(v=pK_value,lwd=2,col='red',lty=2)
  title(paste0(sample,"\nThe BCmvn distributions"))
  text(pK_value,max(bcmvn_kidney$BCmetric),as.character(pK_value),pos = 4,col = "red")
  dev.off()
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(dfso$seurat_clusters)
  perc = 0.008 * (length(colnames(dfso))/1000)
  nExp_poi <- round(expected_doublet_rate*ncol(dfso))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  pN_value = 0.25
  dfso <- doubletFinder(dfso, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE,PCs = 1:npcs,sct = T)
  p1 = DimPlot(object = dfso, pt.size=0.05, reduction = "umap",group.by="seurat_clusters", label = F)
  pAAN_value=tail(names(dfso@meta.data),2)[1]
  #pAAN_value = paste0("pANN_",pN_value,"_",pK_value,"_",nExp_poi)
  dfso <- doubletFinder(dfso, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pAAN_value,PCs = 1:npcs,sct = T)
  label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi.adj)
  dfso@meta.data$DoubletFinder_if.doublet <- dfso@meta.data[, label]
  db_cell = rownames(dfso@meta.data[dfso@meta.data$DoubletFinder_if.doublet=="Doublet",])
  p2 = DimPlot(object = dfso, pt.size=0.05, reduction = "umap",cells.highlight = db_cell, cols.highlight = "red", cols = "gray",label = F)+
    NoLegend() +
    ggtitle(paste0(sample,"_doublet"))
  
  if (length(dev.list()) > 0) {
    dev.off()
  }
  
  print("now plot UMAP!!!!!!!!!!!!")
  pdf(paste0(path,"/",sample,"_DoubletFinder_UMAP.pdf"),width=12, height= 6)
  plot(p2+p1)
  dev.off()
  
  d <- data.table(cells=colnames(dfso),DoubletFinder_doublet_scores = dfso@meta.data[,pAAN_value],DoubletFinder_if.doublet=dfso@meta.data$DoubletFinder_if.doublet)
  #head(d)
  d2 <- data.table(cells=colnames(dfso_ori))
  head(d2)
  d3 <- left_join(d2, d, by = "cells")
  d3 <- d3 %>%
    mutate(DoubletFinder_if.doublet = ifelse(is.na(DoubletFinder_if.doublet), "Low_quality",DoubletFinder_if.doublet))
  doublet_scores <- d3$DoubletFinder_doublet_scores
  names(doublet_scores) <- d3$cells
  if.doublet <- d3$DoubletFinder_if.doublet
  names(if.doublet) <- d3$cells
  dfso_ori <- dfso_ori %>% AddMetaData(metadata=doublet_scores, col.name='DoubletFinder_doublet_scores') %>%
    AddMetaData(metadata=if.doublet, col.name='DoubletFinder_if.doublet')
  doublet_rate_adj = expected_doublet_rate* (1-homotypic.prop)
  doublet_rate_adj = rep(doublet_rate_adj,each=ncol(dfso_ori))
  dfso_ori[["doublet_rate_adj"]] = doublet_rate_adj
  print(head(dfso_ori@meta.data[1:2,]))
  
  return(dfso_ori)
}
