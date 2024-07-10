
PlotQC <- function(object.list, prefix='', outpath='', filter.doublets=TRUE) {
    m <- rbindlist(lapply(object.list, function(object) {
      m <- copy(object@meta.data); setDT(m, keep.rownames='cells')
    }))
    
    #noise <- c('nCount_ATAC', 'nCount_RNA', 'nFeature_ATAC', 'nFeature_RNA', 'percent.mt', 'percent.rb', 'doublet_scores', 'nucleosome_signal', 'TSS.enrichment')
    features1 = c('nCount_ATAC','nFeature_ATAC', 'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal', 'pct_reads_in_peaks','nCount_RNA','nFeature_RNA','percent.mt','doublet_scores')
    #head(object.list[[5]]@meta.data)
    print(unique(m$SampleID))
    print(unique(m$CCDI))
    num_samples <- length(object.list)
    print(num_samples)
    library(RColorBrewer)
    palette <- rainbow(num_samples)
    expected_order <- names(object.list)
    m$SampleID <- factor(m$SampleID, levels = expected_order)
    plot.list1 <- lapply(features1, function(feature) {
      m %>% 
        ggplot(aes(x = SampleID, y = !!as.name(feature), fill = SampleID)) + 
        geom_violin(color = 'gray30') + 
        scale_fill_manual(values = palette) +  # Set custom color palette
        labs(title = feature) +  # Set the feature as the title
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust=0,size = 12, face = "bold",color = "black"),  # Rotate x-axis labels by 90 degrees
              axis.title.x = element_blank(),  # Remove x-axis title
              axis.text.y = element_text( size = 14, face = "bold",color = "black"),
              axis.title.y = element_blank(),  # Remove y-axis title
              axis.line = element_line(color = "black"),  # Change axis line color to black
              axis.ticks = element_line(color = "black"),  # Change axis ticks color to black
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"))
    })
    
    p1 <- ggpubr::ggarrange(plotlist=plot.list1, legend='none', align='hv', ncol=3, nrow=4)
    names(object.list)
    name_lengths <- nchar(names(object.list))
    max_len = max(name_lengths)
    width =10
    height = 15
    if (num_samples >4 ) {
      width = 10 + (num_samples - 4) * 0.5
    }
    if (max_len>10) {
      height = 15 + (max_len - 10) * 0.5
    }
    width
    height
    file_name1 = paste0("qc_",prefix,"_filtering.pdf")
    ggsave(plot=p1, width=width, height=height, filename=paste0(outpath, '/',file_name1))
    
    plot.list <- lapply(features1, function(feature) {
      m %>% 
        ggplot(aes(x = SampleID, y = !!as.name(feature), fill = SampleID)) + 
        geom_violin(color = 'gray30') + 
        geom_point(color = 'gray70', size = 0.001,alpha = 0.5, position = position_jitter(width = 0.15)) +
        scale_fill_manual(values = palette) +  # Set custom color palette
        labs(title = feature) +  # Set the feature as the title
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 12, face = "bold",color = "black"),  # Rotate x-axis labels by 90 degrees
              axis.title.x = element_blank(),  # Remove x-axis title
              axis.text.y = element_text( size = 14, face = "bold",color = "black"),
              axis.title.y = element_blank(),  # Remove y-axis title
              axis.line = element_line(color = "black"),  # Change axis line color to black
              axis.ticks = element_line(color = "black"),  # Change axis ticks color to black
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
              plot.margin = margin(0.2, 0.2, 1.5, 0.2, "cm"))  # Adjust plot margin
    })

    p <- ggpubr::ggarrange(plotlist=plot.list, legend='none', align='hv', ncol=3, nrow=4)
    file_name2 = paste0("qc_",prefix,"_filtering_point.png")
    ggsave(plot=p, width=width, height=height,  units = "in", dpi=300, filename=paste0(outpath, '/', file_name2))
}



