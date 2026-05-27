#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(tidyverse)
    library(fs)
})

option_list <- list(
    make_option(c("-i", "--input-dir"), dest = "input_dir", type = "character", help = "Directory containing KneadData log files"),
    make_option(c("-o", "--output"), dest = "output", type = "character", help = "Output summary TSV file")
)

parser <- OptionParser(
    option_list = option_list,
    description = paste(
        "Summarise KneadData logs into one TSV table.",
        "The table includes raw reads, trimmed reads, clean reads, singletons,",
        "and percentages for trimming, host removal and final retention.",
        sep = " "
    )
)

opt <- parse_args(parser)

if (is.null(opt$input_dir)) stop("Missing required option: --input-dir")
if (is.null(opt$output)) stop("Missing required option: --output")
if (!dir_exists(opt$input_dir)) stop("Input directory does not exist: ", opt$input_dir)

log_files <- dir_ls(opt$input_dir, regexp = "\\.log$")
if (length(log_files) == 0) stop("No KneadData log files found in: ", opt$input_dir)

extract_count <- function(lines, file_pattern) {
    hit <- lines[str_detect(lines, fixed(file_pattern))]
    if (length(hit) == 0) return(NA_real_)
    str_match(hit[1], ":\\s*([0-9]+\\.?[0-9]*)")[, 2] %>% as.numeric()
}

parse_log <- function(path) {
    run <- path_file(path) %>% str_remove("\\.log$")
    lines <- read_lines(path)

    tibble(
        run_accession = run,
        initial_r1 = extract_count(lines, "_R1.fastq"),
        initial_r2 = extract_count(lines, "_R2.fastq"),
        trimmed_paired_r1 = extract_count(lines, paste0(run, ".trimmed.1.fastq")),
        trimmed_paired_r2 = extract_count(lines, paste0(run, ".trimmed.2.fastq")),
        trimmed_single_r1 = extract_count(lines, paste0(run, ".trimmed.single.1.fastq")),
        trimmed_single_r2 = extract_count(lines, paste0(run, ".trimmed.single.2.fastq")),
        clean_paired_r1 = extract_count(lines, paste0(run, "_paired_1.fastq")),
        clean_paired_r2 = extract_count(lines, paste0(run, "_paired_2.fastq")),
        clean_unmatched_r1 = extract_count(lines, paste0(run, "_unmatched_1.fastq")),
        clean_unmatched_r2 = extract_count(lines, paste0(run, "_unmatched_2.fastq"))
    )
}

kneaddata_summary <- map_dfr(log_files, parse_log) %>%
    mutate(
        raw_pairs = pmin(initial_r1, initial_r2, na.rm = TRUE),
        trimmed_pairs = pmin(trimmed_paired_r1, trimmed_paired_r2, na.rm = TRUE),
        clean_pairs = pmin(clean_paired_r1, clean_paired_r2, na.rm = TRUE),
        trimmed_singletons = trimmed_single_r1 + trimmed_single_r2,
        clean_singletons = clean_unmatched_r1 + clean_unmatched_r2,
        pairs_removed_by_trim = raw_pairs - trimmed_pairs,
        pairs_removed_by_host = trimmed_pairs - clean_pairs,
        pairs_removed_total = raw_pairs - clean_pairs,
        pct_pairs_removed_by_trim = 100 * pairs_removed_by_trim / raw_pairs,
        pct_pairs_removed_by_host = 100 * pairs_removed_by_host / trimmed_pairs,
        pct_pairs_retained_final = 100 * clean_pairs / raw_pairs
    ) %>%
    arrange(run_accession)

dir_create(path_dir(opt$output))
write_tsv(kneaddata_summary, opt$output)

message("Written KneadData summary table to: ", opt$output)
