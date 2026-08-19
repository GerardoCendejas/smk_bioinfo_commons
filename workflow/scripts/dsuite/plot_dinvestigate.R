#!/usr/bin/env Rscript

library("ggplot2")
library("tidyverse")
library("zoo")

# Load necessary library
if (!require("optparse", quietly = TRUE)) {
    stop("The 'optparse' package is required but not installed. Please install it using install.packages('optparse').")
}

# 1. Configure command-line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input file from dinvestigate", metavar="archivo"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output plot dir", metavar="archivo"),
  make_option(c("-c","--chroms"), type="character", default=NULL,
              help="Chromosome order file", metavar="string")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Verify required argument
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$chroms)){
  print_help(opt_parser)
  stop("❌ Error: Please provide arguments: --input, --output and --chroms", call.=FALSE)
}

data_file = opt$input
out_name = opt$output
chr_file = opt$chroms

out_name <- sub("/$", "", out_name)

## Read data
data = read.table(data_file, sep = "\t", header = TRUE)
chrs = as.vector(unique(data$chr))
data$mid = (data$windowStart + data$windowEnd)/2
chr_order = read.table(chr_file, header = FALSE)
stats = c("D", "f_d", "f_dM", "d_f")
data$chr = factor(data$chr, levels = chr_order$V1)
data = data[order(data$chr), ]

## Cumulative positions
data$cumulative_mid = NA
cumulative_length = 0
chrs_ordered = levels(data$chr) 

for (chr_name in chrs_ordered) {
    chr_mids = data$mid[which(data$chr == chr_name)]
    data$cumulative_mid[which(data$chr == chr_name)] = chr_mids + cumulative_length
    cumulative_length = cumulative_length + max(data$windowEnd[which(data$chr == chr_name)])
}

axis_labels = data %>%
    group_by(chr) %>%
    summarize(chr_center = (min(cumulative_mid) + max(cumulative_mid)) / 2)

n_chrs = length(chrs_ordered)
color_list = rep(c("grey30", "grey"), length.out = n_chrs)

# --- Main loop for plots (remains as is) ---
all_outliers_list <- list()

for(j in 1:4){
    m = mean(data[,stats[j]], na.rm = TRUE)
    s = sd(data[,stats[j]], na.rm = TRUE)
    is_outlier = (data[,stats[j]] > (m + 3*s)) | (data[,stats[j]] < (m - 3*s))
    
    # Store outliers for each stat
    all_outliers_list[[stats[j]]] <- data[is_outlier, c("chr", "windowStart", "windowEnd", stats[j])]

    # Chromosome plots
    for(i in 1:length(chrs)){
        chr = chrs[i]
        p_data = data[which(data$chr==chr),]
        p_outlier = is_outlier[which(data$chr==chr)]

        p = ggplot(p_data, aes_string(x = "mid/1e6", y = stats[j])) +
            geom_point() +
            geom_hline(yintercept = 0, linetype = "solid", color = "black") +
            geom_hline(yintercept = m, linetype = "dashed", color = "red") +
            geom_hline(yintercept = m + 3*s, linetype = "dotted", color = "red") +
            geom_hline(yintercept = m - 3*s, linetype = "dotted", color = "red") +
            geom_point(data = p_data[p_outlier, ], aes_string(x = "mid/1e6", y = stats[j]), color = "blue") +
            labs(title = paste("Chromosome", chr, "-", stats[j]), x = "Position (Mb)", y = stats[j]) +
            theme_classic() + theme(legend.position = "none")

        ggsave(file.path(out_name,paste0(chr,"_",stats[j],".png")), plot = p, width = 16, height = 8, dpi = 300, bg = "white")
    }

    # Combined plot
    p = ggplot(data, aes_string(x = "cumulative_mid/1e6", y = stats[j])) +
        geom_point(aes_string(color = "chr"), alpha = 0.8) +
        scale_color_manual(values = color_list) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        geom_hline(yintercept = m, linetype = "dashed", color = "red") +
        geom_hline(yintercept = m + 3 * s, linetype = "dotted", color = "red") +
        geom_hline(yintercept = m - 3 * s, linetype = "dotted", color = "red") +
        geom_point(data = data[is_outlier, ], aes_string(x = "cumulative_mid/1e6", y = stats[j]), color = "blue", size = 1.5) +
        scale_x_continuous(breaks = axis_labels$chr_center / 1e6, labels = axis_labels$chr) +
        labs(title = paste("Genome-wide -", stats[j]), x = "Chromosome", y = stats[j]) +
        theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    ggsave(file.path(out_name, paste0("genome_", stats[j], ".png")), plot = p, width = 16, height = 8, dpi = 300, bg = "white")
}

# --- NEW LOGIC FOR LOCI EXTRACTION (d_f based) ---

# 1. Take only the d_f outliers
fd_sig_windows <- all_outliers_list[["d_f"]] %>%
    arrange(chr, windowStart) %>%
    # Prefijo p/n según el valor de f_d
    mutate(prefix = ifelse(d_f > 0, "p", "n")) %>%
    # Columna chrom sin "chr" (ej: 1, 1A, 2)
    mutate(chrom_num = gsub("chr", "", as.character(chr))) %>%
    # Locus_id manteniendo "chr" (ej: p-chr12-13476085)
    mutate(locus_id = paste0(prefix, "-", chr, "-", windowStart))

# Save ALL significant windows for d_f with their values
write_tsv(fd_sig_windows, file.path(out_name, "significant_windows_all.tsv"))

# 2. Function to merge regions
merge_regions <- function(df, gap = 500) {
    if(nrow(df) == 0) return(data.frame())
    df %>%
        group_by(chr) %>%
        # Unimos ventanas según el gap
        mutate(is_new = windowStart > (lag(windowEnd, default = 0) + gap)) %>%
        mutate(region_id = cumsum(is_new)) %>%
        group_by(chr, region_id) %>%
        summarize(
            start = min(windowStart),
            end = max(windowEnd),
            n_windows = n(),
            avg_df = mean(d_f),
            .groups = "drop"
        ) %>%
        # Creamos el ID secuencial: p-chr1-1-w1, p-chr1-2-w1, etc.
        group_by(chr) %>%
        mutate(idx = row_number()) %>%
        mutate(prefix = ifelse(avg_df > 0, "p", "n")) %>%
        mutate(chrom_num = gsub("chr", "", as.character(chr))) %>%
        mutate(locus_id = paste0(prefix, "-", chr, "-", idx, "-w", n_windows)) %>%
        ungroup() %>%
        select(locus_id, chrom = chrom_num, start, end, n_windows)
}

significant_regions_collapsed <- merge_regions(fd_sig_windows)

# Guardamos el archivo final para la extracción con Snakemake
write_tsv(significant_regions_collapsed, file.path(out_name, "significant_regions_collapsed.tsv"))

# Success marker
success_file <- file.path(out_name, "success.log")
write(paste("Dinvestigate plotting and d_f loci extraction completed at:", Sys.time()), file = success_file)
