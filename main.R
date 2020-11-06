#!/usr/local/bin Rscript
library(optparse)
source("read10X.R")

option_list = list(
        make_option(c("-s1", "--sample1"), type="character", default=NULL, 
                    help="sample 1 name", metavar="character"),
        make_option(c("-s2", "--sample2"), type="character", default=NULL, 
                    help="sample 2 name", metavar="character"),
        make_option(c("-p", "--project"), type="character", default=NULL, 
                    help="project name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
        print_help(opt_parser)
        stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

sample_list = c(opt$sample1, opt$sample2)

merged_seurat <- read_input(sample_list, opt$project)

