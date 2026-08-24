#! /usr/bin/env Rscript

# Load necessary library
library(phytools)

# Load necessary library
if (!require("optparse", quietly = TRUE)) {
    stop("The 'optparse' package is required but not installed. Please install it using install.packages('optparse').")
}

# 1. Configure command-line arguments
option_list = list(
 make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input .tree", metavar="archivo"),
 make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output .png", metavar="archivo"),
 make_option(c("-r", "--rooted"), action="store_true", default=FALSE, 
             help="Plot the tree as rooted"),
 make_option(c("-n", "--node_labels"), action="store_true", default=FALSE, 
             help="Plot node labels"),
 make_option(c("-s", "--size"), type="integer", default=800, 
              help="Size of the output image [default %default]", metavar="size")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input) || is.null(opt$output)){
  print_help(opt_parser)
  stop("Input and output files are required.", call.=FALSE)
}


## Script logic
size = opt$size

tree = read.tree(file=opt$input)

if (!opt$rooted){
    
    tree = unroot(tree)

}

png(opt$output,width=size,height=size)

plotTree(tree, fsize=0.6, lwd=1, ftype="i")

if (opt$node_labels){

    nodelabels()
    
}

dev.off()
