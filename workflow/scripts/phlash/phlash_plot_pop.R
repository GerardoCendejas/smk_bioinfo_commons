#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

option_list <- list(
  make_option(c("-i", "--inputs"), type = "character", default = NULL, 
              help = "Comma-separated list of input TSV files", metavar = "files"),
  make_option(c("-p", "--target"), type = "character", default = NULL, 
              help = "Target Population(s) OR Species name(s)", metavar = "name"),
  make_option(c("-m", "--mapping"), type = "character", default = NULL, 
              help = "Mapping TSV file (Col1: Population, Col2: Species)", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = "phlash_dynamic_plot.png", 
              help = "Output PNG file name [default = %default]", metavar = "file"),
  make_option(c("-g", "--generation_time"), type = "numeric", default = 1.0, 
              help = "Generation time in years [default = %default]", metavar = "number"),
  make_option(c("-s", "--start_generations"), type = "numeric", default = 0, 
              help = "Filter out generations before this value [default = %default]", metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

if (is.null(args$inputs) || is.null(args$mapping) || is.null(args$target)) {
  print_help(opt_parser)
  stop("Error: --inputs, --target, and --mapping are required.", call. = FALSE)
}

# 1. Load mapping
mapping_df <- read.delim(args$mapping, header = FALSE, col.names = c("Population", "Species"), stringsAsFactors = FALSE)

# 2. Advanced Logic: Handle mixtures of Species and Populations
input_targets <- unlist(strsplit(args$target, ","))
target_pops <- c()

# Determine if we should color by Population or Species
# If the user provides only ONE species, we usually want to see the different populations (Islands)
if (length(input_targets) == 1 && input_targets %in% mapping_df$Species) {
    target_pops <- mapping_df$Population[mapping_df$Species == input_targets]
    color_by <- "Population"
    legend_title <- paste("Islands (", input_targets, ")")
} else {
    # If multiple items, resolve each one
    for (t in input_targets) {
        if (t %in% mapping_df$Species) {
            # It's a species name: add all its populations
            target_pops <- c(target_pops, mapping_df$Population[mapping_df$Species == t])
        } else {
            # It's a population name: add it directly
            target_pops <- c(target_pops, t)
        }
    }
    color_by <- "Species"
    legend_title <- "Species"
}

target_pops <- unique(target_pops)

# 3. Process Files
input_files <- unlist(strsplit(args$inputs, ","))
all_data <- data.frame()

for (fpath in input_files) {
  # Get the filename without extension (e.g., Geospiza_scandens_Daphne)
  pop_id <- tools::file_path_sans_ext(basename(fpath))

  # Check if this file is one of our resolved targets
  if (!(pop_id %in% target_pops)) next
  
  species_id <- mapping_df$Species[mapping_df$Population == pop_id]
  if (length(species_id) == 0) next

  # Read TSV (phlash output)
  df <- read.delim(fpath, header = TRUE, check.names = FALSE)
  filtered_df <- df[df[[1]] >= args$start_generations, , drop = FALSE]
  
  if (nrow(filtered_df) > 0) {
    # Calculate stats
    tmp_data <- data.frame(
      Time = filtered_df[[1]] * args$generation_time,
      Ne_median = apply(filtered_df[, -1, drop = FALSE], 1, median),
      Ne_low = apply(filtered_df[, -1, drop = FALSE], 1, quantile, probs = 0.05),
      Ne_high = apply(filtered_df[, -1, drop = FALSE], 1, quantile, probs = 0.95),
      Population = pop_id,
      Species = species_id
    )
    all_data <- rbind(all_data, tmp_data)
  }
}

if (nrow(all_data) == 0) {
    stop(paste("Error: No data found. \nTargets looked for:", paste(target_pops, collapse=", ")))
}

# 4. Plotting
p <- ggplot(all_data, aes(x = Time, group = Population, color = !!sym(color_by), fill = !!sym(color_by))) +
  geom_ribbon(aes(ymin = Ne_low, ymax = Ne_high), alpha = 0.15, color = NA) +
  geom_line(aes(y = Ne_median), linewidth = 0.8) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(
    x = ifelse(args$generation_time == 1, "Time (Generations)", "Time (Years)"),
    y = expression(Effective ~ Population ~ Size ~ (N[e])),
    title = "Demographic History (phlash)",
    color = legend_title,
    fill = legend_title
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(args$output, plot = p, width = 10, height = 7, dpi = 300)
cat(paste("Success! Found", length(unique(all_data$Population)), "populations. Plot saved to:", args$output, "\n"))
