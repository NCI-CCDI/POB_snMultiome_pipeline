
Doublet_GMM_cutoff <- function(object.list, project.name='', outpath='', filter.doublets=TRUE) {
  newobj_list = list()
  samplelist = names(object.list)
  cutoff_list=list()
  for (s in samplelist) {
    print(s)
    m <- copy(object.list[[s]]@meta.data)
   cutoff <- NormalMixCutoff(mixtools::normalmixEM(m[, "doublet_scores"], k=2))
   if(length(cutoff) >1) {
     cutoff = max(cutoff)
   }
   cutoff_list[["cutoff"]][[s]]=cutoff
   print(cutoff)
   setDT(m, keep.rownames='cells')
   #print(dim(m))
   m[, `:=` (
     project=rep(project.name, nrow(m)),
     predicted_gmm_doublets=fifelse(doublet_scores < cutoff, 'singlet', 'doublet')
   )]
   project <- m[, project]; names(project) <- m[, cells]
   if.doublet.gmm <- m[,predicted_gmm_doublets]; names(if.doublet.gmm) = m[,cells]
   object = object.list[[s]]
   object <- object %>% AddMetaData(metadata=project, col.name='project') %>%
     AddMetaData(metadata=if.doublet.gmm, col.name='if.doublet.gmm')
   newobj_list[[s]] = object
   print(object@meta.data[1:2,])
   cutoff_list[["doublet"]][[s]]=sum(object@meta.data$if.doublet.gmm=="doublet")
   cutoff_list[["total"]][[s]]=nrow(object@meta.data)
  }   
  print(cutoff_list)
  #cat("doublet cutoff for each sample is ", unlist(cutoff_list))
  df <- data.frame(
    sample = unlist(samplelist),
    gmm_doublet_score_cutoff = unlist(cutoff_list[["cutoff"]]),  
    doublet_cell_number = unlist(cutoff_list[["doublet"]]),  
    total_cell_number = unlist(cutoff_list[["total"]]))
  print(df)
  write.csv(df,paste0(output_path,"/qc_result_plot/2.doublet_filter_threshold_gmm.csv"),row.names=F,quote =F)
  return(newobj_list)

}



