library(dplyr)
library(ggplot2)

filter_seurat <- function(seurat_obj, project) {
        # seurat_obj@meta.data$orig.ident <- sapply(seurat_obj@meta.data$orig.ident, function(i) substr(i,1,3))
        seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
        
        # Compute percent mito ratio
        seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
        seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
        
        metadata <- seurat_obj@meta.data
        metadata$cells <- rownames(metadata)
        metadata <- metadata %>%
                dplyr::rename(nUMI = nCount_RNA,
                              nGene = nFeature_RNA)
        
        metadata$sample <- NA
        metadata$sample <- metadata$orig.ident
        # metadata$sample[which(str_detect(metadata$cells, "^N_"))] <- "seurat_obj"
        # metadata$sample[which(str_detect(metadata$cells, "^S_"))] <- "seurat_obj"
        
        seurat_obj@meta.data <- metadata
        
        # Visualize the number of cell counts per sample
        metadata %>% 
                ggplot(aes(x=sample, fill=sample)) + 
                geom_bar() +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                theme(plot.title = element_text(hjust=0.5, face="bold")) +
                ggtitle("NCells")
        ggsave("outputs/test_cell_count.png")
        # Visualize the number UMIs/transcripts per cell
        metadata %>% 
                ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
                geom_density(alpha = 0.2) + 
                scale_x_log10() + 
                theme_classic() +
                ylab("Cell density") +
                geom_vline(xintercept = 500)
        ggsave("outputs/test_UMI_density.png")
        
        # Visualize the distribution of genes detected per cell via boxplot
        metadata %>% 
                ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
                geom_boxplot() + 
                theme_classic() +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                theme(plot.title = element_text(hjust=0.5, face="bold")) +
                ggtitle("NCells vs NGenes")
        ggsave("outputs/test_gene_count_boxplot.png")
        
        # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
        metadata %>% 
                ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
                geom_point() + 
                scale_colour_gradient(low = "gray90", high = "black") +
                stat_smooth(method=lm) +
                scale_x_log10() + 
                scale_y_log10() + 
                theme_classic() +
                geom_vline(xintercept = 500) +
                geom_hline(yintercept = 250) +
                facet_wrap(~sample)
        ggsave("outputs/test_nGene~nUMI.png")
        # Visualize the distribution of mitochondrial gene expression detected per cell
        metadata %>% 
                ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
                geom_density(alpha = 0.2) + 
                scale_x_log10() + 
                theme_classic() +
                geom_vline(xintercept = 0.2)
        ggsave("outputs/test_mito_distribution.png")
        # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
        metadata %>%
                ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
                geom_density(alpha = 0.2) +
                theme_classic() +
                geom_vline(xintercept = 0.8)
        ggsave("outputs/test_gene_complexity.png")
        
        # Filter out low-quality cells
        filtered_seurat <- subset(x = seurat_obj, 
                                  subset=
                                          (nGene >= 250) & 
                                          (log10GenesPerUMI > 0.80) & 
                                          (mitoRatio < 0.20))
        
        # Gene level filtering 
        counts <- GetAssayData(object = filtered_seurat, slot = "counts")
        nonzero <- counts > 0
        keep_genes <- Matrix::rowSums(nonzero) >= 10
        filtered_counts <- counts[keep_genes, ]
        filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data, project = project)
        
        
        
        return(filtered_seurat)
}