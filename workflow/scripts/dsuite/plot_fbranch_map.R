#!/usr/bin/env Rscript

# Load necessary libraries
if (!require("optparse", quietly = TRUE)) install.packages("optparse")
suppressPackageStartupMessages({
  library(optparse)
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(dplyr)
  library(tidyr)
})

# 1. Configure command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input fbranch .txt matrix", metavar="file"),
  make_option(c("-c", "--coords"), type="character", default=NULL, 
              help="Input CSV file with coordinates (location,lon,lat)", metavar="file"),
  make_option(c("-o", "--out_prefix"), type="character", default="fbranch_map",
              help="Output prefix for the two generated .png files", metavar="prefix"),
  make_option(c("-p", "--p_value"), type="double", default=0.01,
              help="P-value threshold [default= %default]", metavar="double")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_file <- opt$input
coords_file <- opt$coords
p_threshold <- opt$p_value
out_prefix <- opt$out_prefix

if (is.null(input_file) || is.null(coords_file)) {
  stop("Error: Both input matrix (-i) and coordinates CSV (-c) must be provided.")
}

# 2. Read Island coordinates
island_coords <- read.table(coords_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(island_coords)[colnames(island_coords) == "location"] <- "Island"

# 3. Read the fbranch file matrices
read_dsuite_section <- function(file, section_start_pattern, section_end_pattern=NULL) {
  lines <- readLines(file)
  start_idx <- if (is.null(section_start_pattern)) 1 else grep(section_start_pattern, lines) + 1
  
  if (is.null(section_end_pattern)) {
    next_headers <- grep("^#", lines)
    next_headers <- next_headers[next_headers > start_idx]
    end_idx <- if(length(next_headers) > 0) next_headers[1] - 2 else length(lines)
  } else {
    end_idx <- length(lines)
  }
  
  text_block <- paste(lines[start_idx:end_idx], collapse="\n")
  data <- read.table(text = text_block, header = TRUE, sep = "\t", check.names = FALSE, na.strings = c("nan", "NaN", "-nan"))
  return(data)
}

data_f <- read_dsuite_section(input_file, NULL)
data_p <- read_dsuite_section(input_file, "# p-values:")

# 4. Clean matrices and convert to long format
clean_to_long <- function(df, value_name) {
  rownames(df) <- df$branch_descendants
  df <- df[, -c(1,2,3)] # Remove branch, branch_descendants, Outgroup
  df$source <- rownames(df)
  df_long <- pivot_longer(df, cols = -source, names_to = "target", values_to = value_name)
  return(df_long)
}

edges_f <- clean_to_long(data_f, "f_val")
edges_p <- clean_to_long(data_p, "p_val")
edges <- merge(edges_f, edges_p, by = c("source", "target"))

# 5. Filter significant edges and keep ONLY terminal branches
edges <- edges %>%
  filter(!is.na(f_val), !is.na(p_val)) %>%
  filter(f_val > 0) %>%
  filter(p_val <= p_threshold) %>%
  filter(!grepl(",", source)) # Terminal branches lack commas

# 6. Extract nodes and parse taxonomy
all_nodes <- unique(c(edges$source, edges$target))
nodes_df <- data.frame(node = all_nodes, stringsAsFactors = FALSE)

parse_node <- function(n) {
  parts <- strsplit(n, "_")[[1]]
  if (length(parts) >= 3) {
    sp <- paste(parts[1:2], collapse = "_")
    isl <- paste(parts[3:length(parts)], collapse = "_")
  } else {
    sp <- n
    isl <- "Unknown"
  }
  return(c(species = sp, island = isl))
}

parsed <- t(sapply(nodes_df$node, parse_node))
nodes_df$species <- parsed[, "species"]
nodes_df$island <- parsed[, "island"]

# Distribute nodes sharing the same island in a small radius
nodes_df <- merge(nodes_df, island_coords, by.x = "island", by.y = "Island", all.x = TRUE)
nodes_df <- nodes_df[!is.na(nodes_df$lon), ]

nodes_df <- nodes_df %>%
  group_by(island) %>%
  mutate(
    n_nodes = n(),
    radius = ifelse(n_nodes > 1, 0.08, 0),
    angle = ifelse(n_nodes > 1, seq(0, 2 * pi, length.out = n() + 1)[-(n() + 1)], 0),
    plot_lon = lon + radius * cos(angle),
    plot_lat = lat + radius * sin(angle)
  ) %>%
  ungroup()

# 7. Merge coordinates, calculate non-linear scaling, and classify event type
edges <- edges %>%
  inner_join(nodes_df %>% select(node, plot_lon, plot_lat, species), by = c("source" = "node")) %>%
  rename(source_lon = plot_lon, source_lat = plot_lat, source_sp = species) %>%
  inner_join(nodes_df %>% select(node, plot_lon, plot_lat, species), by = c("target" = "node")) %>%
  rename(target_lon = plot_lon, target_lat = plot_lat, target_sp = species) %>%
  mutate(
    curve_val = ifelse(source_lon < target_lon, 0.2, -0.2),
    # Non-linear transformation to exaggerate large values and suppress noise
    f_nl = (f_val / max(f_val))^2,
    # Classify event
    event_type = ifelse(source_sp == target_sp, "Intraspecific", "Interspecific")
  )

# 8. Prepare spatial data
sf::sf_use_s2(FALSE)
world <- ne_countries(scale = "large", returnclass = "sf")
bbox <- st_bbox(c(xmin = -92.5, xmax = -88, ymin = -2.2, ymax = 1.8), crs = st_crs(world))
galapagos <- st_crop(world, bbox)
galapagos <- st_make_valid(galapagos)

island_labels <- island_coords %>% filter(Island %in% nodes_df$island)

# ==============================================================================
# 9. PLOT 1: Magnitude Focus (Dark lines, exaggerated alpha/width)
# ==============================================================================
out_file_mag <- paste0(out_prefix, ".png")
png(out_file_mag, width = 1400, height = 1200, res = 120)

ggplot() +
  geom_sf(data = galapagos, fill = "#f0f0f0", color = "black", size=1) +
  # Lines scaled by f_nl (the squared, non-linear metric)
  geom_curve(data = edges, 
             aes(x = source_lon, y = source_lat, 
                 xend = target_lon, yend = target_lat, 
                 linewidth = f_nl, alpha = f_nl),
             color = "#333333", curvature = 0.2, 
             arrow = arrow(length = unit(0.015, "npc"), type = "closed")) +
  geom_point(data = nodes_df, aes(x = plot_lon, y = plot_lat, fill = species), 
             shape = 21, size = 5, color = "white", stroke = 1.2) +
  geom_text(data = island_labels, aes(x = lon, y = lat, label = Island),
            nudge_y = -0.15, size = 3, fontface = "bold.italic") +
  scale_linewidth_continuous(range = c(0.1, 4.5), guide = "none") +
  scale_alpha_continuous(range = c(0.1, 0.95), guide = "none") +
  scale_fill_brewer(palette = "Paired", name = "Species") + 
  coord_sf(crs = st_crs(4326), xlim = c(-92.5, -88), ylim = c(-2.2, 1.8), expand = FALSE) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#ffffff"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 18)) +
  labs(title = paste("Introgression Magnitude (p <", p_threshold, ")"),
       subtitle = "Non-linear scaling emphasizes major gene flow events.")

invisible(dev.off())

# ==============================================================================
# 10. PLOT 2: Event Type Focus (Intra vs Inter-specific)
# ==============================================================================
out_file_type <- paste0(out_prefix, "_2.png")
png(out_file_type, width = 1400, height = 1200, res = 120)

ggplot() +
  geom_sf(data = galapagos, fill = "#f0f0f0", color = "black", size=1) +
  geom_curve(data = edges, 
             aes(x = source_lon, y = source_lat, 
                 xend = target_lon, yend = target_lat, 
                 linewidth = f_nl, alpha = f_nl, color = event_type),
             curvature = 0.2, 
             arrow = arrow(length = unit(0.015, "npc"), type = "closed")) +
  geom_point(data = nodes_df, aes(x = plot_lon, y = plot_lat, fill = species), 
             shape = 21, size = 5, color = "white", stroke = 1.2) +
  geom_text(data = island_labels, aes(x = lon, y = lat, label = Island),
            nudge_y = -0.15, size = 3, fontface = "bold.italic") +
  scale_linewidth_continuous(range = c(0.1, 4.5), guide = "none") +
  scale_alpha_continuous(range = c(0.1, 0.95), guide = "none") +
  # Custom colors for event types
  scale_color_manual(values = c("Intraspecific" = "#ff7f00", "Interspecific" = "#984ea3"), 
                     name = "Event Type") +
  scale_fill_brewer(palette = "Paired", name = "Species") + 
  coord_sf(crs = st_crs(4326), xlim = c(-92.5, -88), ylim = c(-2.2, 1.8), expand = FALSE) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#ffffff"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 18)) +
  labs(title = paste("Intra vs Inter-specific Introgression (p <", p_threshold, ")"),
       subtitle = "Colors denote if introgression occurred within the same species or across species.")

invisible(dev.off())
