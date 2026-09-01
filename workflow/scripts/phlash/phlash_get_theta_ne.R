#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))

option_list <- list(
    make_option(c("-m", "--mutation"), type = "double", default = NULL,
                help = "Mutation rate (m) per site per generation", metavar = "number"),
    make_option(c("-o", "--output"), type = "character", default = "ne_results.csv",
                help = "Output CSV file name [default %default]", metavar = "character")
)

parser <- OptionParser(
    usage = "usage: %prog [options] file1.log file2.log ...",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = TRUE)

opt <- arguments$options
log_files <- arguments$args

if (is.null(opt$mutation)) {
    print_help(parser)
    stop("Error: You must provide the mutation rate using the -m flag.", call. = FALSE)
}

if (length(log_files) == 0) {
    stop("Error: No .log files provided. Please specify the input files.", call. = FALSE)
}

results_list <- list()

for (f in log_files) {
    if (!file.exists(f)) {
        warning(paste("File not found, skipping:", f))
        next
    }
    
    lines <- readLines(f, warn = FALSE)
    
    match_line <- grep("Scaled mutation rate Θ=", lines, value = TRUE)
    
    if (length(match_line) > 0) {
        theta_str <- str_extract(match_line[1], "(?<=Θ=)[0-9.]+(e-?[0-9.]+)?")
        theta_val <- as.numeric(theta_str)
        
        if (!is.na(theta_val)) {
            ne_val <- theta_val / (4 * opt$mutation)
            
            pop_name <- tools::file_path_sans_ext(basename(f))
            
            results_list[[f]] <- data.frame(
                population = pop_name,
                theta = theta_val,
                Ne = ne_val,
                stringsAsFactors = FALSE
            )
        }
    } else {
        warning(paste("Pattern 'Scaled mutation rate Θ=' not found in:", f))
    }
}

if (length(results_list) > 0) {
    final_df <- do.call(rbind, results_list)
    write.table(final_df, opt$output, sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat(sprintf("\nSuccessfully processed %d files.\n", nrow(final_df)))
    cat(sprintf("Results saved to: %s\n", opt$output))
} else {
    cat("\nNo valid data found in the provided files.\n")
}
