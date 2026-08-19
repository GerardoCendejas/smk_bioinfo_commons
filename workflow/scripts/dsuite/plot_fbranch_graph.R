#!/usr/bin/env Rscript

# Load necessary library
if (!require("optparse", quietly = TRUE)) {
    stop("The 'optparse' package is required but not installed. Please install it using install.packages('optparse').")
}

# 1. Configure command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input .txt", metavar="archivo"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output .png file", metavar="archivo"),
  make_option(c("-p", "--p_value"), type="double", default=0.01,
              help="P-value threshold [default= %defeautl]", metavar="double"),
  make_option(c("-g","--graph_script"), type="character", default="igraphplot2.R",
              help="Path to modified igraph plotting script (default: igraphplot2.R)", metavar="archivo")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(igraph)
source(opt$graph_script)

environment(plot.igraph2) <- asNamespace('igraph')
environment(igraph.Arrows2) <- asNamespace('igraph')

input_file <- opt$input
p_threshold <- opt$p_value
output_file <- opt$output

# Read the fbranch file
read_dsuite_section <- function(file, section_start_pattern, section_end_pattern=NULL) {
  lines <- readLines(file)
  
  # Find start
  if (is.null(section_start_pattern)) {
    start_idx <- 1 # First matrix
  } else {
    start_idx <- grep(section_start_pattern, lines) + 1
  }
  
  # Locate end
  if (is.null(section_end_pattern)) {
    next_headers <- grep("^#", lines)
    next_headers <- next_headers[next_headers > start_idx]
    if(length(next_headers) > 0) {
      end_idx <- next_headers[1] - 2 
    } else {
      end_idx <- length(lines)
    }
  } else {
    end_idx <- length(lines)
  }
  
  # na.strings = "nan" Important to avoid parsing issues with non-significant values
  text_block <- paste(lines[start_idx:end_idx], collapse="\n")
  data <- read.table(text = text_block, header = TRUE, sep = "\t", check.names = FALSE, na.strings = c("nan", "NaN", "-nan"))
  return(data)
}

# Read f-branch matrix
data_f <- read_dsuite_section(input_file, NULL)

# Read p-values matrix
data_p <- read_dsuite_section(input_file, "# p-values:")

# Clean and align matrices
clean_matrix <- function(df) {
  rownames(df) <- df$branch_descendants
  
  df <- df[,-c(1,2,3)]
  
  df <- df[match(colnames(df), rownames(df)), ]
  
  return(as.matrix(df))
}

mat_f <- clean_matrix(data_f)
mat_p <- clean_matrix(data_p)

# Parse Nas
mat_f[is.na(mat_f)] <- 0
mat_p[is.na(mat_p)] <- 1 # Si no hay p-value, asumimos no significativo (1)

# Remove edges that are not significant
mat_f[mat_p > p_threshold] <- 0

g <- graph_from_adjacency_matrix(mat_f, mode="directed", weighted=TRUE)

# Remove isolated nodes
g <- delete.vertices(g, which(degree(g) == 0))

vert <- names(V(g))
w <- E(g)$weight

island <- c()
species <- c()

for(i in 1:length(vert)){
  name <- strsplit(vert[i], split="_")[[1]]
  if(length(name) >= 3) {
    tmp_isl <- paste(name[3:length(name)], collapse = "_")
    tmp_sp <- paste(name[1:2], collapse = "_")
  } else {
    tmp_isl <- "Unknown"
    tmp_sp <- vert[i]
  }
  island <- c(island, tmp_isl)
  species <- c(species, tmp_sp)
}

V(g)$island <- island
V(g)$species <- species
V(g)$color <- as.factor(species)

# Calculate size by degree (number of connections)
degrees <- igraph::degree(g)
# Normalize sizes to a reasonable range (e.g., 5 to 20)
max_deg <- max(degrees)
if(max_deg == 0) max_deg <- 1
sizes <- degrees/max_deg * 20

# --- PLOT ---

png(output_file, 1200, 1200)

if(length(V(g)) > 0) {
  layout <- layout_with_lgl(g)
  
  max_w <- max(w)
  if(max_w == 0) max_w <- 1
  edge_alphas <- w/max_w
  
  plot.igraph2(g,
               edge.arrow.size = (E(g)$weight * 2) * 5,
               edge.curved = 0.2,
               edge.color = rgb(0, 0, 0, edge_alphas),
               edge.width = (E(g)$weight * 2) * 20,
               vertex.label.color = "black",
               vertex.size = sizes,
               layout = layout,
               main = paste("F-branch (p <", p_threshold, ")"))
} else {
  plot(0, type='n', axes=FALSE, ann=FALSE)
  text(1, 1, "No significant edges found")
}

dev.off()
