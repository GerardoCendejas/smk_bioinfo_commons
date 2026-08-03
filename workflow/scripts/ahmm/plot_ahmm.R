#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggtern)
  library(RColorBrewer)
})

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, help="Folder with .posterior files"),
  make_option(c("-o", "--outdir"), type="character", default="plots_ancestry", help="Output folder"),
  make_option(c("-l", "--label"), type="character", default="genomewide", help="Label for titles")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("Error: --input_dir is mandatory.", call.=FALSE)
}

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# --- 1. DATA LOADING ---
files <- list.files(opt$input_dir, pattern = "\\.posterior$", full.names = TRUE)
if (length(files) == 0) stop("No .posterior files found.")

message(paste("Reading", length(files), "files..."))

read_posterior <- function(f) {
  df <- read.table(f, header = TRUE, check.names = FALSE)
  colnames(df) <- c("chrom", "position", "p20", "p11", "p02")
  df$sample <- gsub("\\.posterior$", "", basename(f))
  return(df)
}

all_data <- map_df(files, read_posterior)

# Mean Ancestry Calculation (0 to 1)
# 1.0 = pop_a (2,0) | 0.5 = het (1,1) | 0.0 = pop_b (0,2)
all_data <- all_data %>%
  mutate(mean_ancestry = (2*p20 + 1*p11) / 2)

# --- 2. TERNARY PLOTS (Points and Density) ---
# Mapping: x=Left(p20), y=Top(p11), z=Right(p02)
p_tern_base <- ggtern(data = all_data, aes(x = p20, y = p11, z = p02)) +
  geom_point(alpha = 0.1, size = 0.4, color = "darkslategrey") +
  labs(T = "het (1,1)", L = "pop_a (2,0)", R = "pop_b (0,2)",
       title = paste("Ternary Distribution:", opt$label)) +
  theme_rgbw() +
  theme_showarrows() +
  theme(
    tern.axis.arrow.text.T = element_text(color = "#FB8C00", face = "bold"),
    tern.axis.arrow.text.L = element_text(color = "#D32F2F", face = "bold"),
    tern.axis.arrow.text.R = element_text(color = "#1976D2", face = "bold")
  )

# A. Save Ternary with only Points
ggsave(file.path(opt$outdir, "ternary_points.png"), p_tern_base, width = 8, height = 8, dpi = 300)

# B. Save Ternary with Density (exactly as plot_qs.R style)
p_tern_density <- p_tern_base +
  stat_density_tern(
    geom = 'polygon',
    aes(fill = ..level..),
    base = "identity",
    colour = "white",
    alpha = 0.5
  ) +
  scale_fill_viridis_c(option = "magma", name = "Density") +
  guides(fill = guide_colourbar(order = 1))

ggsave(file.path(opt$outdir, "ternary_density.png"), p_tern_density, width = 8, height = 8, dpi = 300)

# --- 3. DENSITY PLOTS (Overall and Per Chromosome) ---
density_data <- all_data %>%
  pivot_longer(cols = c(p20, p11, p02), names_to = "state", values_to = "prob")

# Helper for density style
plot_density_style <- function(df, title_text) {
  ggplot(df, aes(x = prob, fill = state)) +
    geom_density(alpha = 0.4, color = "white") +
    scale_fill_manual(values = c("p20" = "#D32F2F", "p11" = "#FB8C00", "p02" = "#1976D2"),
                      labels = c("p02"="pop_b (0,2)", "p11"="het (1,1)", "p20"="pop_a (2,0)")) +
    labs(title = title_text, x = "Posterior Probability", y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Overall
ggsave(file.path(opt$outdir, "overall_density.png"), 
       plot_density_style(density_data, paste("Overall Probability Density:", opt$label)), 
       width = 10, height = 6)

# Per Chromosome
p_dens_chr <- plot_density_style(density_data, paste("Density by Chromosome:", opt$label)) + 
  facet_wrap(~chrom, scales = "free_y")
ggsave(file.path(opt$outdir, "chromosome_densities.png"), p_dens_chr, width = 16, height = 12)

# --- 4. GENOMIC MANHATTAN PLOT (Mean Ancestry) ---
p_genomic <- ggplot(all_data, aes(x = position/1e6, y = mean_ancestry, color = sample)) +
  geom_point(alpha = 0.3, size = 0.3) + 
  facet_wrap(~chrom, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(title = paste("Genomic Ancestry Track:", opt$label),
       subtitle = "Y=1.0: pop_a (2,0) | Y=0.5: het (1,1) | Y=0.0: pop_b (0,2)",
       x = "Position (Mb)", y = "Mean Ancestry") +
  theme_minimal() + 
  theme(legend.position = "none")

ggsave(file.path(opt$outdir, "genomic_ancestry_manhattan.png"), p_genomic, width = 20, height = 12)
