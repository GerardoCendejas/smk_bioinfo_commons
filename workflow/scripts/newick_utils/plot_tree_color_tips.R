#! /usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
    library(optparse)
    library(ggtree)
    library(ggplot2)
    library(ape)
    library(phytools)
})

# Config command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input .tree (Newick)", metavar="archivo"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output image (.png, .pdf)", metavar="archivo"),
  make_option(c("-m", "--map"), type="character", default=NULL, 
              help="CSV file without header (col1: ID, col2: Group). Optional.", metavar="archivo"),
  make_option(c("-r", "--rooted"), action="store_true", default=TRUE, 
              help="Keep tree rooted (default is to unroot)"),
  make_option(c("-n", "--node_labels"), action="store_true", default=FALSE, 
              help="Plot node labels (bootstrap/support)"),
  make_option(c("-s", "--size"), type="integer", default=10, 
              help="Size of the output image [default %default]", metavar="size"),
  make_option(c("-d", "--dpi"), type="integer", default=300, 
              help="DPI of the output image [default %default]", metavar="dpi")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input) || is.null(opt$output)){
  print_help(opt_parser)
  stop("Input and output files are required.", call.=FALSE)
}

# Save size of image

size = opt$size
dpi = opt$dpi

# Read the tree
tree <- read.tree(file=opt$input)

# Rooting logic
if (!opt$rooted){
    tree <- unroot(tree)
}

# Initialize ggtree plot
p <- ggtree(tree, size=0.8) 

# Coloring tips based on map (if provided)
if (!is.null(opt$map)) {
    
    # With map (color tips)
    # Read CSV mapping file
    # Species_name,Group
    info <- read.csv(opt$map, header=FALSE, stringsAsFactors=FALSE)
    colnames(info) <- c("label", "group") # Naming columns for ggtree join
    
    # Join the mapping info to the tree data
    p <- p %<+% info + 
         # Color the tip labels based on group
         geom_tiplab(size=3, align=TRUE, linesize=0.5) +
         # Add colored points at tips
         geom_tippoint(aes(color=group), size=3) +
         # Legend title
         scale_color_discrete(name="") +
         theme(legend.position="top")

} else {

    # Without map (default tip labels)
    p <- p + geom_tiplab(size=3)
}

# Node labels (if requested)
if (opt$node_labels){
    p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-0.3, size=2)
}

p <- p +
    theme_tree2() 

tryCatch({
    # Obtenemos la altura máxima de cualquier punta o nodo
    max_x <- max(phytools::nodeHeights(tree))
}, error = function(e) {
    # Respaldo de emergencia si phytools falla totalmente
    message("Resort to ape depth calculation")
    max_x <- max(ape::nodeDepth.edgelen(tree))
})
p <- p + xlim(0, max_x * 1.3)

# Save the plot
ggsave(filename=opt$output, plot=p, width=size, height=size, dpi=dpi)
