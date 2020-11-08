library(Seurat)
library(dplyr)
library(ggplot2)
source("filter_seurat.R")

read_filter <- function(sample_list, project) {
        dir_list <- paste0(dir_list, "_filtered_feature_bc_matrix")
        seurat_list = list()
        for (file in dir_list){
                seurat_data <- Read10X(data.dir = file)
                seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                                 min.features = 100, 
                                                 project = file)
                seurat_list[[file]] <- seurat_obj
        }
        merged_seurat <- merge(object_list[[1]], y = tail(seurat_list,-1), add.cell.ids = sample_lists, project = project)
        filtered_seurat <- filter_seurat(merged_seurat)

        return(merged_seurat)
}