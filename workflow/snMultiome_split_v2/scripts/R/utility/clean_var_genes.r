
require(Seurat)
require(data.table)

CleanVarGenes <- function(object, nHVG=2000) {
    
    if (class(object) != 'Seurat') stop('please provide a Seurat object...')
    
    g <- data.table(genes=rownames(object))
    
    genes <- object %>% 
        subset(features=g[!(genes %like% '^MT' | genes %like% '^RP[LS]'), genes]) %>% 
        NormalizeData() %>% FindVariableFeatures(nfeatures=nHVG) %>% VariableFeatures()

    VariableFeatures(object) <- genes

    return(object)
}