#!/usr/local/bin Rscript
library(optparse)
source("read10X.R")

option_list = list(
        make_option(c("-d1", "--dir1"), type="character", default=NULL, 
                    help="data directory 1", metavar="character"),
        make_option(c("-d2", "--dir2"), type="character", default=NULL, 
                    help="ata directory 2", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
        print_help(opt_parser)
        stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

input_data_list = c(opt$dir1, opt$dir2)

read_input(input_data_list)

