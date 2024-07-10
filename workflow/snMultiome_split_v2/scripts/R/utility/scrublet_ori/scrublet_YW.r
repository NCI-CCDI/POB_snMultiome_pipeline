
scrublet_YW <- function(object, min_counts=2, min_cells=3, expected_doublet_rate=0.06, min_gene_variability_pctl=85, n_prin_comps=30, sim_doublet_ratio=2,sample="ccdi",path="/data/CCRCCDI", umi_min = 0, ngene_min= 0, mt_max = 100,n_neighbors=NULL) {
    #umi_min = 500
    #ngene_min = 300
    #mt_max = 10
    if (class(object) != 'Seurat') stop('object is not of type "Seurat"') 
    object_ori = object
    cat("filter umi:", umi_min,"\n",
        "filter ngene:", ngene_min, "\n",
        "filter mt:", mt_max, "\n")
    object = object %>% subset(
      subset=
        nCount_RNA > umi_min &
        nFeature_RNA > ngene_min &
        percent.mt < mt_max)
    
    #reticulate::source_python('/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/New_code2/scripts/utility/scrublet_py.py')
    reticulate::source_python('/data/CCRCCDI/software/Scripts/current_snakemake/snMultiome_split_v2/scripts/R/utility/scrublet_py_YW.py')
    DefaultAssay(object) = "RNA"
    X <- as(t(as.matrix(GetAssayData(object=object, layer='counts'))), 'TsparseMatrix')
    val <- X@x
    i <- as.integer(X@i)
    j <- as.integer(X@j)
    dim <- as.integer(X@Dim)
    #print(X[1:10,1:5])
    if (is.null(n_neighbors)) n_neighbors <- round(0.5 * sqrt(nrow(X)))
    cat("n_neighbors:",n_neighbors)

    scrublet_py_args <- c(
      list(
        i=i, j=j, val=val, dim=dim,
        expected_doublet_rate=doublet_rate, min_counts=2, 
        min_cells=3, min_gene_variability_pctl=85, 
        n_prin_comps=30, sim_doublet_ratio=2, n_neighbors=n_neighbors,
        sample=sample,path=path
      )
    )
    
    scrublet_res <- do.call(scrublet_py_YW, scrublet_py_args)
    #(doublet_scores, predicted_doublets,cutoff_auto,cutoff_gmm1, cutoff_percentile1,cutoff_final,scr_cri,percentile)
    names(scrublet_res) <- c('doublet_scores', 'scrub.auto_doublet','scrub.cutoff_auto','scrub.cutoff_gmm','scrub.cutoff_percentile','scrub.cutoff_final','scrub.criteria','scrub.percentile')
    #print(scrublet_res)

    d <- data.table(cells=rownames(X), doublet_scores=as.numeric(scrublet_res$doublet_scores),scrub.auto_doublet=scrublet_res$scrub.auto_doublet)
    #head(d)
    d2 <- data.table(cells=colnames(object_ori))
    head(d2)
    d3 <- left_join(d2, d, by = "cells")
    d3 <- d3 %>%
      mutate(scrub.auto_doublet = ifelse(is.na(scrub.auto_doublet), "Low_quality", ifelse(scrub.auto_doublet, "Doublet", "Singlet")))
    doublet_scores <- d3$doublet_scores
    names(doublet_scores) <- d3$cells
    cat("the length of doublet_scores:",length(doublet_scores),"\n")
    
    scrub.auto_doublet <- d3$scrub.auto_doublet
    names(scrub.auto_doublet) <- d3$cells
    cat("the lenfth of scrub.auto_doublet",length(scrub.auto_doublet),"\n")
    
    object_ori <- object_ori %>% AddMetaData(metadata=doublet_scores, col.name='doublet_scores') %>%
      AddMetaData(metadata=scrub.auto_doublet, col.name='scrub.auto_doublet')

    for (i in 3:length(names(scrublet_res))){
      label_name = names(scrublet_res)[i]
      label_now = as.character(scrublet_res[[label_name]])
      #print(label_name)
      #print(label_now)
      meta_label = rep(label_now, each=ncol(object_ori))
      #print(head(meta_label))
      object_ori[[label_name]] = meta_label
    } 
    
    object_ori@meta.data <- object_ori@meta.data %>%
      mutate(scrub.if_doublet_gmm = ifelse(is.na(doublet_scores), "Low_quality", 
                                           ifelse(doublet_scores < scrublet_res$scrub.cutoff_gmm, "Singlet", "Doublet"))) %>%
      mutate(scrub.if_doublet_percentile = ifelse(is.na(doublet_scores), "Low_quality", 
                                                  ifelse(doublet_scores < scrublet_res$scrub.cutoff_percentile, "Singlet", "Doublet"))) %>%
      mutate(scrub.if_doublet = ifelse(is.na(doublet_scores), "Low_quality", 
                                       ifelse(doublet_scores < scrublet_res$scrub.cutoff_final, "Singlet", "Doublet")))
    
    table(object_ori@meta.data$scrub.if_doublet)
    print(head(object_ori@meta.data[1:2,]))
    
    return(object_ori)
  
}

