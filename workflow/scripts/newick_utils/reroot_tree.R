#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
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
  # 1. Start from an unrooted tree to avoid root placement artifacts
  unrooted_tr <- unroot(tree)
  
  ingroup_vec <- setdiff(unrooted_tr$tip.label, outgroup_vec)
  n_tips <- length(unrooted_tr$tip.label)

  if (length(outgroup_vec) == 1) {
    target_node <- which(unrooted_tr$tip.label == outgroup_vec)
    edge_idx <- which(unrooted_tr$edge[, 2] == target_node)
    branch_len <- unrooted_tr$edge.length[edge_idx]
    if (is.na(branch_len)) branch_len <- 0
    rerooted_tree <- reroot(unrooted_tr, node = target_node, position = branch_len / 2)

  } else {
    # 2. Check bipartitions to find the branch that isolates the outgroup
    # Propagate clades from downstream nodes
    target_edge <- NULL
    target_node <- NULL
    
    # Check all internal nodes and edges
    for (i in 1:nrow(unrooted_tr$edge)) {
      desc_node <- unrooted_tr$edge[i, 2]
      
      # Get all tips descending from desc_node
      desc_tips <- unrooted_tr$tip.label[propagate_tips <- if (desc_node <= n_tips) desc_node else {
        descendants <- phytools::getDescendants(unrooted_tr, desc_node)
        descendants[descendants <= n_tips]
      }]
      
      # Check if this edge exactly separates outgroup or ingroup
      if (setequal(desc_tips, outgroup_vec)) {
        target_node <- desc_node
        target_edge <- i
        break
      } else if (setequal(desc_tips, ingroup_vec)) {
        # The parent node connects to the outgroup stem
        target_node <- desc_node
        target_edge <- i
        break
      }
    }

    if (!is.null(target_edge)) {
      branch_len <- unrooted_tr$edge.length[target_edge]
      if (is.na(branch_len)) branch_len <- 0
      rerooted_tree <- reroot(unrooted_tr, node = target_node, position = branch_len / 2)
    } else {
      # Fallback: Reroot on ingroup MRCA first, then resolve outgroup MRCA
      ingroup_mrca <- getMRCA(unrooted_tr, ingroup_vec)
      temp_tr <- reroot(unrooted_tr, node = ingroup_mrca, position = 0)
      
      out_mrca <- getMRCA(temp_tr, outgroup_vec)
      edge_idx <- which(temp_tr$edge[, 2] == out_mrca)
      b_len <- temp_tr$edge.length[edge_idx]
      if (is.na(b_len)) b_len <- 0
      
      rerooted_tree <- reroot(temp_tr, node = out_mrca, position = b_len / 2)
    }
  }

  # Clean branch lengths if original was purely topological
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
