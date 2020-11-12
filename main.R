#!/usr/local/bin/env Rscript --vanilla
library(optparse)
source("read10X.R")

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
# sample_list = 1:3
if (length(sample_list) > 1) {
        is_merge = T
} else {
        is_merge = F
}

print("arguments parsed successfully, start to load data")

merged_seurat <- read_filter(sample_list, opt$project, is_merge)

print("complete!")

