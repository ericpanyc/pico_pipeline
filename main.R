suppressMessages(library(optparse))
suppressMessages(library(miceadds))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

path <- "~/Downloads/usc_folder/kgp/seurat/pico_pipeline/functions"
source.all(path, grepstring="\\.R",  print.source=F, file_sep="__")

option_list = list(
        make_option(c("-a", "--sample1"), type="character", default=NULL, 
                    help="sample 1 name", metavar="character"),
        make_option(c("-b", "--sample2"), type="character", default=NULL, 
                    help="sample 2 name", metavar="character"),
        make_option(c("-c", "--sample3"), type="character", default=NULL, 
                    help="sample 3 name", metavar="character"),
        make_option(c("-d", "--sample4"), type="character", default=NULL, 
                    help="sample 4 name", metavar="character"),
        make_option(c("-e", "--sample5"), type="character", default=NULL, 
                    help="sample 5 name", metavar="character"),
        make_option(c("-p", "--project"), type="character", default="project", 
                    help="project name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sample1)){
        print_help(opt_parser)
        stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


sample_list = c(opt$sample1,opt$sample2,opt$sample3,opt$sample4,opt$sample5)
group <- opt$project





merged_seurat <- read_filter(sample_list, group)
filtered_seurat <- filter_seurat(merged_seurat, group)

print("complete!")

