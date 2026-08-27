#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(optparse)
})

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input tree file (.tree, .newick)", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output rerooted tree", metavar="file"),
  make_option(c("-g", "--outgroup"), type="character", default=NULL, 
              help="Path to a text file with tip labels (one per line) OR a comma-separated string", 
              metavar="file_or_string")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output) || is.null(opt$outgroup)){
  print_help(opt_parser)
  stop("Input, Output and Outgroup are mandatory.")
}

if (file.exists(opt$outgroup)) {
  # First case: it's a file
  message(paste("Reading file:", opt$outgroup))
  outgroup_vec <- readLines(opt$outgroup)
} else {
  # Second case: it's a comma-separated string
  message("Reading names.")
  outgroup_vec <- unlist(strsplit(opt$outgroup, ","))
}

# Clean up outgroup vector
outgroup_vec <- trimws(outgroup_vec)
outgroup_vec <- outgroup_vec[outgroup_vec != ""] # Remove empty strings

# --- Script logic ---

tree <- read.tree(file = opt$input)

if (is.null(tree$edge.length)) {
    message("No branch lengths detected. Assigning unit lengths for rerooting.")
    tree$edge.length <- rep(1.0, nrow(tree$edge))
    has_orig_lengths <- FALSE
} else {
    has_orig_lengths <- TRUE
    n_tips <- length(tree$tip.label)
    edge_to_tips <- tree$edge[, 2] <= n_tips
    
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
                                        # 2. Identify the target node (MRCA of the outgroup)
    if (length(outgroup_vec) == 1) {
                                        # If outgroup is a single species, find its index
        target_node <- which(tree$tip.label == outgroup_vec)
    } else {
                                        # If outgroup is a clade, find the MRCA node
        target_node <- getMRCA(tree, outgroup_vec)
    }

    root_node <- length(tree$tip.label) + 1

    if (is.null(target_node)) {
        stop("Could not resolve MRCA for the provided outgroup taxa.")
    }

    if (target_node == root_node) {
        message("The selected outgroup is already at the root. No rerooting needed.")
        rerooted_tree <- tree
    } else {
        
                                        # 3. Find the edge connecting this node to the rest of the tree
        edge_index <- which(tree$edge[, 2] == target_node)

        if (length(edge_index) == 0) {
            stop("The selected outgroup node appears to be the root already or has no upstream branch.")
        }

                                        # 4. Get branch length and calculate midpoint
        branch_len <- tree$edge.length[edge_index]
        
                                        # Asegurar que branch_len sea numérico para la división
        if (is.na(branch_len)) branch_len <- 0
        
        split_position <- branch_len / 2
        
                                        # 5. Reroot using phytools
        rerooted_tree <- reroot(tree, node = target_node, position = split_position)
        
                                        # 6. Limpieza final antes de guardar
                                        # Si el árbol original era topológico, borramos las ramas. 
                                        # De lo contrario, nos aseguramos de que no queden NaNs residuales.
        if (!has_orig_lengths) {
            rerooted_tree$edge.length <- NULL
        } else {
            if (any(is.na(rerooted_tree$edge.length))) {
                rerooted_tree$edge.length[is.na(rerooted_tree$edge.length)] <- 0
            }
        }
        
                                        # 7. Save output
        write.tree(rerooted_tree, file = opt$output)
        message(paste("Successfully rerooted to:", opt$output))
    }    
}, error = function(e) {
    stop(paste("ERROR during rerooting:", e$message))
})
