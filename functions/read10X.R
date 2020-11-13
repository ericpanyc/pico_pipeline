read_filter <- function(sample_list, project) {
        
        seurat_list = list()
        for (sample in sample_list){
                seurat_data <- Read10X(data.dir = paste0("../data/",sample, "_filtered_feature_bc_matrix"))
                seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                                 min.features = 100,
                                                 project = sample)
                seurat_list[[sample]] <- seurat_obj
        }
        if (length(sample_list) > 1) {
                print("start to merge samples")
                merged_seurat <- merge(seurat_list[[1]], y = tail(seurat_list,-1), add.cell.ids = sample_list, project = project)
        } else {
                print("skip merge")
                merged_seurat <- seurat_list[[1]]
        }
        return(merged_seurat)

        
}