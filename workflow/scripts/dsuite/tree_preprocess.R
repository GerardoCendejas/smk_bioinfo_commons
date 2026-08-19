#!/usr/bin/env Rscript

# Load necessary library
if (!require("optparse", quietly = TRUE)) {
    stop("The 'optparse' package is required but not installed. Please install it using install.packages('optparse').")
}

# 1. Configure command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input .tre", metavar="archivo"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output .tre file", metavar="archivo"),
  make_option(c("-g","--outgroup"), type="character", default=NULL,
              help="Outgroup names comma separated", metavar="string")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Verify required argument
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$outgroup)){
  print_help(opt_parser)
  stop("❌ Error: Please provide arguments: --input, --output and --outgroup", call.=FALSE)
}

# 2. Load data
message(paste("Reading:", opt$input))
input_tree = opt$input
output_tree = opt$output
outgroup_pops = unlist(strsplit(opt$outgroup, split=","))

suppressPackageStartupMessages(library(ape))

# Read tree
if (!file.exists(input_tree)) stop(paste("File not found:", input_tree))
tree <- read.tree(input_tree)

# Check if outgroup populations are in the tree
missing <- outgroup_pops[!outgroup_pops %in% tree$tip.label]
if (length(missing) > 0) {
  warning(paste("Warning: These populations not in tree:", paste(missing, collapse=", ")))
}

valid_outgroups <- outgroup_pops[outgroup_pops %in% tree$tip.label]

if (length(valid_outgroups) == 0) {
  stop("Error: No valid outgroup populations found in tree. Please check your input.")
}

keeper <- valid_outgroups[1]
droppers <- valid_outgroups[2:length(valid_outgroups)]


if (length(droppers) > 0) {
  tree <- drop.tip(tree, droppers)
  message(paste("Removed:", paste(droppers, collapse=", ")))
}

idx <- which(tree$tip.label == keeper)
tree$tip.label[idx] <- "Outgroup"
message(paste("Tip renamed '", keeper, "' a 'Outgroup'", sep=""))

tree$node.label <- NULL

write.tree(tree, file = output_tree)
message(paste("Tree saved in:", output_tree))
