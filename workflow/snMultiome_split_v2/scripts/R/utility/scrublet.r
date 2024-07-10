
scrublet <- function(object, min_counts=2, min_cells=3, expected_doublet_rate=0.06, min_gene_variability_pctl=85, n_prin_comps=30, sim_doublet_ratio=2,n_neighbors=NULL) {
 
    if (class(object) != 'Seurat') stop('object is not of type "Seurat"') 
    #reticulate::source_python('/data/CCRCCDI/analysis/ccrccdi4/03_Secondary_analysis/Multiome_analysis_Rproject/New_code2/scripts/utility/scrublet_py.py')
    reticulate::source_python('/data/CCRCCDI/software/Scripts/current_snakemake/snMultiome_split_v2/scripts/R/utility/scrublet_py.py')
    X <- as(t(as.matrix(GetAssayData(object=object, layer='counts'))), 'TsparseMatrix')
    
    val <- X@x
    i <- as.integer(X@i)
    j <- as.integer(X@j)
    dim <- as.integer(X@Dim)

    if (is.null(n_neighbors)) n_neighbors <- round(0.5 * sqrt(nrow(X)))

    scrublet_py_args <- c(
        list(
                i=i, j=j, val=val, dim=dim,
                expected_doublet_rate=expected_doublet_rate, min_counts=min_counts, 
                min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, 
                n_prin_comps=n_prin_comps, sim_doublet_ratio=sim_doublet_ratio, n_neighbors=n_neighbors
            )
        )

    scrublet_res <- do.call(scrublet_py_YW, scrublet_py_args)
    names(scrublet_res) <- c('doublet_scores', 'predicted_doublets')

    #d <- data.table(cells=rownames(X), doublet_scores=as.numeric(scrublet_res$doublet_scores))

    #doublet_scores <- d[, doublet_scores]
    #names(doublet_scores) <- d[, cells]


    #object <- object %>% AddMetaData(metadata=doublet_scores, col.name='doublet_scores')
    
    d <- data.table(cells=rownames(X), doublet_scores=as.numeric(scrublet_res$doublet_scores), if.doublet=scrublet_res$predicted_doublets)
    doublet_scores <- d[, doublet_scores]
    names(doublet_scores) <- d[, cells]
    #print(head(doublet_scores))
    #cat("now is in result",head(doublet_scores),"\n")
    
    if.doublet <- d[, if.doublet]
    names(if.doublet) <- d[, cells]
    #cat("now is in result",head(if.doublet),"\n")
    
    object <- object %>% AddMetaData(metadata=doublet_scores, col.name='doublet_scores') %>%
      AddMetaData(metadata=if.doublet, col.name='if.doublet')
    #print(head(object@meta.data[1:2,]))
                    
    return(object)
  
}

