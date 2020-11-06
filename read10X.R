library(Seurat)

read_input <- function(sample_list, project) {
        dir_list <- paste0(dir_list, "_filtered_feature_bc_matrix")
        seurat_list = list()
        for (file in dir_list){
                seurat_data <- Read10X(data.dir = dir_list)
                seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                                 min.features = 100, 
                                                 project = file)
                seurat_list[[file]] <- seurat_obj
        }
        merged_seurat <- merge(object_list[[1]], y = tail(seurat_list,-1), add.cell.ids = c("VS5", "VS8"), project = project)
        return(merged_seurat)
}