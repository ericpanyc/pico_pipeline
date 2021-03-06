# pico_pipeline
This command line tool will take a list of 10x single-cell RNA sequencing results as input (up to 5 samples for now), and will run preliminary filtering and processing on the input data. It will generate a set fo QC plots for you to check the 10x results and a filtered seurat object in the form of RDS for the downstream analysis at the outputs/ directory<br/>

Usage:<br/>
--sample1~5 {XYZ1} to pass in up to 5 samples of 10x results, should be stored in the data/ folder under the parent directory of the main.R, and the suffix of the input data should be _filtered_feature_bc_matrix<br/>
--project {control} to specify the group of the input data come from<br/>
--nGene {200} to give the minimum threshold of number of genes per cell for filtering, default value is 250<br/>
--geneUMI {0.7} to give the minimum threshold of log10GenesPerUMI for filtering, default value is 0.8<br/>
--mitoRatio {0.2} to give the maximum threshold of mito genes ratio for filtering, default value is 0.2<br/>

Users should install the following packages on R 3.6.3 up front:<br/>
* optparse_1.6.6 
* mice_3.11.0
* miceadds_3.10-28
* Seurat_3.1.5
* dplyr_1.0.0
* ggplot2_3.3.1
