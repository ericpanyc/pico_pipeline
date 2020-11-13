suppressMessages(library(optparse))
suppressMessages(library(miceadds))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

path <- "~/Downloads/usc_folder/kgp/seurat/pico_pipeline/functions"
source.all(path, grepstring="\\.R",  print.source=F, file_sep="__")

option_list = list(
        make_option("--sample1", type="character", default=NULL, 
                    help="sample 1 name", metavar="character"),
        make_option("--sample2", type="character", default=NULL, 
                    help="sample 2 name", metavar="character"),
        make_option("--sample3", type="character", default=NULL, 
                    help="sample 3 name", metavar="character"),
        make_option("--sample4", type="character", default=NULL, 
                    help="sample 4 name", metavar="character"),
        make_option("--sample5", type="character", default=NULL, 
                    help="sample 5 name", metavar="character"),
        make_option("--project", type="character", default="project", 
                    help="project name", metavar="character"),
        make_option("--nGene", type="integer", default=250, 
                    help="threshold for number of genes", metavar="number"),
        make_option("--geneUMI", type="double", default=0.8, 
                    help="threshold for log10GenesPerUMI", metavar="number"),
        make_option("--mitoRatio", type="double", default=0.2, 
                    help="threshold for mito gene ratio", metavar="number")
        
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sample1)){
        print_help(opt_parser)
        stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


sample_list = c(opt$sample1,opt$sample2,opt$sample3,opt$sample4,opt$sample5)
sample_names <- paste(sample_list,collapse = "_")
group <- opt$project
nGene <- opt$nGene
geneUMI <- opt$geneUMI
mitoR <- opt$mitoRatio




print("arguments parsed successfully, start to load data")
merged_seurat <- read_filter(sample_list, group)

print("seurat object created successfully, start to filter data")
filtered_seurat <- filter_seurat(merged_seurat, group, sample_names, nGene, geneUMI, mitoR)

print("seurat object filtered successfully, start to save data")
saveRDS(filtered_seurat, paste0("outputs/",sample_names,"_filtered.RDS"))

print("complete!")

