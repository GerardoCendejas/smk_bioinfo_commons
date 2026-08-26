#!/usr/bin/env Rscript

# Load libraries silently to keep logs clean
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra) 
  library(RColorBrewer)
})

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = ".",
              help = "Directory containing the .Q files (e.g., 2.Q, 3.Q) [default: %default]", metavar = "DIR"),
  make_option(c("-f","--fam"), type = "character", default = NULL,
              help = "Full path to the .fam file (used for individual IDs)", metavar = "FILE"),
  make_option(c("-p","--pop_map"), type = "character", default = NULL,
              help = "Population map file (Col1: ID, Col2: Species/Population)", metavar = "FILE"),
  make_option(c("-s","--sort"), type = "character", default = NULL,
              help = "File containing the desired order of populations (one per line)", metavar = "FILE"),
  make_option(c("--min"), type = "integer", default = 2,
              help = "Minimum K to plot [default: %default]", metavar = "INT"),
  make_option(c("--max"), type = "integer", default = 5,
              help = "Maximum K to plot [default: %default]", metavar = "INT"),
  make_option(c("-o", "--out"), type = "character", default = "admixture_plot.png",
              help = "Output filename (.png) [default: %default]", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$fam) || is.null(opt$pop_map) || is.null(opt$sort)) {
  print_help(opt_parser)
  stop("Error: Missing mandatory arguments (--fam, --meta, --pops).", call. = FALSE)
}

if (!file.exists(opt$fam)) stop(paste("File not found:", opt$fam))
if (!file.exists(opt$meta)) stop(paste("File not found:", opt$meta))
if (!file.exists(opt$pops)) stop(paste("File not found:", opt$pops))


# Read metadata (ID Species)
records <- read.table(opt$meta, header = FALSE, stringsAsFactors = FALSE)
colnames(records)[1:2] <- c("Ind", "Pop")

# Read population order
populations_order <- read.table(opt$pops, header = FALSE, stringsAsFactors = FALSE)$V1

# Read .fam file to match IDs with Admixture rows
fam_data <- read.table(opt$fam, stringsAsFactors = FALSE)
ind_ids <- fam_data$V2

# Plot function for each K
make_admixture_plot <- function(k, input_dir, show_x = FALSE, spacing = 1) {

  # Construct filename: just "K.Q" (e.g., "3.Q") inside the input directory
  q_filename <- file.path(input_dir, paste0(k, ".Q"))
  
  if (!file.exists(q_filename)) {
    warning(paste("File not found:", q_filename, "- Skipping K =", k))
    return(NULL)
  }
  
  # Read Q matrix
  data <- read.table(q_filename)
  
  # Assign IDs
  data$ind <- ind_ids
  # Match with metadata
  data$Species <- records$Pop[match(data$ind, records$Ind)]
  
  # Filter individuals missing from metadata
  data <- data[!is.na(data$Species), ]
  
  # Filter for target populations and set factor levels for ordering
  data <- data %>% 
    filter(Species %in% populations_order) %>%
    mutate(Species = factor(Species, levels = populations_order))
  
  # Identify ancestry columns (V1, V2...)
  ancestry_cols <- grep("^V\\d+", names(data), value = TRUE)
  
  # Calculate means per species
  species_means <- data %>% 
    group_by(Species) %>% 
    summarise(across(all_of(ancestry_cols), mean, na.rm = TRUE), .groups = "drop")
  
  # Identify dominant ancestry for each species
  species_dominant_ancestry <- species_means %>% 
    pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "MeanProp") %>% 
    group_by(Species) %>% 
    slice_max(order_by = MeanProp, n = 1, with_ties = FALSE) %>% 
    select(Species, DominantAncestry = Ancestry)
  
  # Log to console
  cat(sprintf("\n--- Average ancestry proportions for K = %d ---\n", k))
  print(species_means %>% column_to_rownames("Species") %>% round(3))
  
  # Sort individuals: Primary by Species, Secondary by their dominant ancestry proportion
  data <- left_join(data, species_dominant_ancestry, by = "Species")
  data <- data %>% 
    rowwise() %>% 
    mutate(DomValue = get(DominantAncestry)) %>% 
    ungroup() %>% 
    arrange(Species, desc(DomValue))
  
  # --- Spacing (Gaps between populations) ---
  # Create a blank row
  blank_row <- as_tibble(matrix(NA_real_, nrow = spacing, ncol = length(ancestry_cols)), .name_repair = "minimal")
  colnames(blank_row) <- ancestry_cols
  
  # Insert blank rows between groups
  df_with_gaps <- data %>% 
    group_split(Species) %>% 
    purrr::map_dfr(~bind_rows(.x, blank_row)) %>% 
    mutate(RowID = row_number()) %>% 
    tidyr::fill(Species, .direction = "downup") # Fill Species in blank rows to avoid grouping errors
  
  # Pivot to long format for ggplot
  df_long <- df_with_gaps %>% 
    pivot_longer(cols = all_of(ancestry_cols), names_to = "Ancestry", values_to = "Proportion") %>% 
    filter(!is.na(Proportion)) # Remove NAs (gaps will be rendered as white space)
  
  # --- Color Palette ---
  # Use Set1, interpolate if K > 9
  if (k <= 9) {
    colors <- brewer.pal(max(3, k), "Set1")[1:k] # max(3) prevents errors if k < 3
  } else {
    colors <- colorRampPalette(brewer.pal(9, "Set1"))(k)
  }
  
  # --- Generate Plot ---
  p <- ggplot(df_long, aes(x = RowID, y = Proportion, fill = Ancestry)) + 
    geom_bar(stat = "identity", width = 1) + 
    scale_fill_manual(values = colors) +
    theme_minimal(base_size = 12) + 
    labs(y = paste0("K = ", k)) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.ticks = element_blank(), 
      panel.grid = element_blank(),
      axis.title.x = element_blank(), 
      axis.text.y = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold", angle = 0, vjust = 0.5), # K label orientation
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0, "cm")
    )
  
  # Handle X-axis labels (only for the last plot)
  if (show_x) {
    species_positions <- df_long %>% 
      group_by(Species) %>% 
      summarise(mid = mean(RowID, na.rm = TRUE), .groups = "drop") %>% 
      filter(Species %in% populations_order)
    
    p <- p + 
      theme(
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        plot.margin = margin(0, 0, 0.5, 0, "cm") # Extra space at bottom
      ) + 
      scale_x_continuous(
        breaks = species_positions$mid, 
        labels = species_positions$Species, 
        expand = expansion(add = 0)
      )
  } else {
    p <- p + 
      theme(axis.text.x = element_blank()) + 
      scale_x_continuous(expand = expansion(add = 0))
  }
  
  return(p)
}

# --- 5. Main Loop ---

plot_list <- list()

# Ensure K values are integers
start_k <- as.integer(opt$min)
end_k <- as.integer(opt$max)

cat(sprintf("Generating plots from K=%d to K=%d...\n", start_k, end_k))

for (k in start_k:end_k) {
  # Show X-axis labels only for the last K
  is_last <- (k == end_k)
  
  g_plot <- make_admixture_plot(k, opt$dir, show_x = is_last)
  
  if (!is.null(g_plot)) {
    plot_list[[length(plot_list) + 1]] <- ggplotGrob(g_plot)
  }
}

# --- 6. Save Final Image ---

if (length(plot_list) > 0) {
  # Calculate dynamic height: base height + extra per plot
  h <- 300 * length(plot_list) + 300
  
  cat(paste("Saving combined plot to:", opt$out, "\n"))
  
  png(file = opt$out, width = 2400, height = h, res = 200)
  grid.newpage()
  # Use grid.draw with rbind (from gridExtra) to stack plots perfectly aligned
  grid.draw(do.call(rbind, plot_list)) 
  dev.off()
  
  message("Done successfully!")
} else {
  stop("No plots were generated. Please check if .Q files exist in the specified directory.")
}
