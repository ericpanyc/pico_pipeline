library(Seurat)

read_input <- function(dir_list) {
        object_list <- paste0(dir_list, "_filtered_feature_bc_matrix")
        for (file in object_list){
                seurat_data <- Read10X(data.dir = paste0("data/", file))
                seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                                 min.features = 100, 
                                                 project = file)
                assign(file, seurat_obj)
        }
        merge() <- merge(object_list[1], y = ob, add.cell.ids = c("VS5", "VS8"), project = "naive")
        
}