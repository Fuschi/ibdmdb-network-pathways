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

extract_read_count <- function(lines, label) {
    pattern <- paste0("READ COUNT: ", label, " :")
    hit <- lines[str_detect(lines, fixed(pattern))]

    if (length(hit) == 0) return(NA_real_)

    str_match(hit[1], "\\):\\s*([0-9]+\\.?[0-9]*)\\s*$")[, 2] %>%
        as.numeric()
}

parse_log <- function(path) {
    run <- path_file(path) %>% str_remove("\\.log$")
    lines <- read_lines(path)

    tibble(
        run_accession = run,
        raw_pair1 = extract_read_count(lines, "raw pair1"),
        raw_pair2 = extract_read_count(lines, "raw pair2"),
        trimmed_pair1 = extract_read_count(lines, "trimmed pair1"),
        trimmed_pair2 = extract_read_count(lines, "trimmed pair2"),
        trimmed_orphan1 = extract_read_count(lines, "trimmed orphan1"),
        trimmed_orphan2 = extract_read_count(lines, "trimmed orphan2"),
        final_pair1 = extract_read_count(lines, "final pair1"),
        final_pair2 = extract_read_count(lines, "final pair2"),
        final_orphan1 = extract_read_count(lines, "final orphan1"),
        final_orphan2 = extract_read_count(lines, "final orphan2")
    )
}

kneaddata_summary <- map_dfr(log_files, parse_log) %>%
    mutate(
        raw_pairs = pmin(raw_pair1, raw_pair2, na.rm = TRUE),
        trimmed_pairs = pmin(trimmed_pair1, trimmed_pair2, na.rm = TRUE),
        clean_pairs = pmin(final_pair1, final_pair2, na.rm = TRUE),
        trimmed_singletons = coalesce(trimmed_orphan1, 0) + coalesce(trimmed_orphan2, 0),
        clean_singletons = coalesce(final_orphan1, 0) + coalesce(final_orphan2, 0),
        pairs_removed_by_trim = raw_pairs - trimmed_pairs,
        pairs_removed_after_trim = trimmed_pairs - clean_pairs,
        pairs_removed_total = raw_pairs - clean_pairs,
        pct_pairs_removed_by_trim = 100 * pairs_removed_by_trim / raw_pairs,
        pct_pairs_removed_after_trim = 100 * pairs_removed_after_trim / trimmed_pairs,
        pct_pairs_retained_final = 100 * clean_pairs / raw_pairs
    ) %>%
    arrange(run_accession)

dir_create(path_dir(opt$output))
write_tsv(kneaddata_summary, opt$output)

message("Written KneadData summary table to: ", opt$output)
