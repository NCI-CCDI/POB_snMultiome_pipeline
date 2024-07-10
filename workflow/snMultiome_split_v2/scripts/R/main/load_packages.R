

PKGS <- c('Seurat', 'Signac', 'EnsDb.Hsapiens.v86','dplyr', 'tidyverse','ggplot2', 'SoupX', 'GenomeInfoDb', 
          'AnnotationHub', 'ensembldb', 'data.table', 'reticulate', 'biovizBase', 'rootSolve','patchwork',
          'Routliers','biovizBase','mixtools','ggpubr','SingleR','DropletUtils')


invisible(sapply(PKGS, require, character.only=TRUE))
invisible(lapply(list.files('scripts/R/utility', full.names=TRUE, pattern='\\.r$'), source))
invisible(lapply(list.files('scripts/R/utility', full.names=TRUE, pattern='\\.R$'), source))

#data_path <- 'data/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/Multiome/'


setDTthreads(threads=1)
