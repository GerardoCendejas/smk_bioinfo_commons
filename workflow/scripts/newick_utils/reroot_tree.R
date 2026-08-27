#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input tree file (.tree, .newick)", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Output rerooted tree", metavar = "file"),
  make_option(c("-g", "--outgroup"), type = "character", default = NULL, 
              help = "Path to a text file with tip labels (one per line) OR a comma-separated string", 
              metavar = "file_or_string")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$outgroup)) {
  print_help(opt_parser)
  stop("Input, Output and Outgroup are mandatory.")
}

if (file.exists(opt$outgroup)) {
  message(paste("Reading file:", opt$outgroup))
  outgroup_vec <- readLines(opt$outgroup)
} else {
  message("Reading names.")
  outgroup_vec <- unlist(strsplit(opt$outgroup, ","))
}

# Clean up outgroup vector
outgroup_vec <- trimws(outgroup_vec)
outgroup_vec <- outgroup_vec[outgroup_vec != ""]

# --- Script logic ---

tree <- read.tree(file = opt$input)

if (is.null(tree$edge.length)) {
  message("No branch lengths detected. Assigning unit lengths for rerooting.")
  tree$edge.length <- rep(1.0, nrow(tree$edge))
  has_orig_lengths <- FALSE
} else {
  has_orig_lengths <- TRUE
  if (any(is.na(tree$edge.length))) {
    message("Sanitizing branch lengths (converting NA/NaN to 0).")
    tree$edge.length[is.na(tree$edge.length)] <- 0
  }
}

missing_tips <- outgroup_vec[!outgroup_vec %in% tree$tip.label]

if (length(missing_tips) > 0) {
  if (length(missing_tips) == length(outgroup_vec)) {
    stop("ERROR: No tip labels from the outgroup are present in the tree.")
  } else {
    warning(paste("Omitting", length(missing_tips), "species from outgroup list since they are not in the tree."))
    outgroup_vec <- outgroup_vec[outgroup_vec %in% tree$tip.label]
  }
}

tryCatch({
  # Reroot using ape::root to handle paraphyletic/grade outgroups robustly
  rerooted_tree <- ape::root(
    tree, 
    outgroup = outgroup_vec, 
    resolve.root = TRUE, 
    edgelabel = TRUE
  )

  # Cleanup branch lengths if original lacked them or introduced NAs
  if (!has_orig_lengths) {
    rerooted_tree$edge.length <- NULL
  } else if (!is.null(rerooted_tree$edge.length)) {
    rerooted_tree$edge.length[is.na(rerooted_tree$edge.length)] <- 0
  }

  write.tree(rerooted_tree, file = opt$output)
  message(paste("Successfully rerooted to:", opt$output))

}, error = function(e) {
  stop(paste("ERROR during rerooting:", e$message))
})
