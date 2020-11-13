# pico_pipeline
This command line tool will take a list of 10x single-cell RNA sequencing results as input (up to 5 samples for now), and will run preliminary filtering and processing on the input data. It will generate a set fo QC plots for you to check the 10x results and a filtered seurat object in the form of RDS for the downstream analysis at the outputs/ directory. \n

Usage:
--sample1~5 <XYZ1> to pass in up to 5 samples of 10x results, should be stored in the data/ folder under the same directory as the main.R, and the suffix of the input data should be _filtered_feature_bc_matrix. \n
--project <control> to specify the group of the input data come from \n
--nGene <200> to give the minimum threshold of number of genes per cell for filtering, default value is 250 \n
--geneUMI <0.7> to give the minimum threshold of log10GenesPerUMI for filtering, default value is 0.8 \n
--mitoRatio <0.2> to give the maximum threshold of mito genes ratio for filtering, default value is 0.2 \n
