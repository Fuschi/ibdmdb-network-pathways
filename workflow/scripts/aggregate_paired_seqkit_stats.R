#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(janitor)
  library(fs)
})

option_list <- list(
  make_option(c("-i", "--input-dir"), dest = "input_dir", type = "character", help = "Directory containing per-FASTQ SeqKit TSV files"),
  make_option(c("-o", "--output"), dest = "output", type = "character", help = "Output aggregated TSV file")
)

parser <- OptionParser(
  option_list = option_list,
  description = paste(
    "Aggregate per-FASTQ SeqKit stats into one paired-end table.",
    "Input files must be named as <run>_R1.tsv and <run>_R2.tsv.",
    "R1/R2 columns are kept with suffixes _r1 and _r2.",
    "Aggregated paired-end columns use suffix _r12:",
    "counts and lengths are summed;",
    "min/max lengths use min/max;",
    "average length, GC%, Q20%, Q30% and average quality are weighted by sequence length or read count;",
    "Q1/Q2/Q3 and N50 are averaged;",
    "paired_reads_match checks whether R1 and R2 have the same number of reads.",
    sep = " "
  )
)

opt <- parse_args(parser)

if (is.null(opt$input_dir)) stop("Missing required option: --input-dir")
if (is.null(opt$output)) stop("Missing required option: --output")
if (!dir_exists(opt$input_dir)) stop("Input directory does not exist: ", opt$input_dir)

seqkit_files <- dir_ls(opt$input_dir, regexp = "_R[12]\\.tsv$")
if (length(seqkit_files) == 0) stop("No SeqKit TSV files found in: ", opt$input_dir)

message("Found ", length(seqkit_files), " SeqKit TSV files")

parse_seqkit_filename <- function(path) {
  file_name <- path_file(path)
  sample_name <- str_remove(file_name, "\\.tsv$")
  parsed <- str_match(sample_name, "^(.*)_R([12])$")
  if (any(is.na(parsed))) stop("Cannot parse run/mate from file name: ", file_name)
  tibble(run_accession = parsed[, 2], mate = paste0("R", parsed[, 3]))
}

read_seqkit_file <- function(path) {
  id <- parse_seqkit_filename(path)
  read_tsv(path, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(run_accession = id$run_accession, mate = id$mate, .before = 1)
}

has_cols <- function(data, cols) all(cols %in% names(data))

weighted_mean_pair <- function(x1, x2, w1, w2) {
  if_else(is.na(w1) | is.na(w2) | (w1 + w2) == 0, NA_real_, ((x1 * w1) + (x2 * w2)) / (w1 + w2))
}

seqkit_long <- map_dfr(seqkit_files, read_seqkit_file)

seqkit_wide <- seqkit_long %>%
  pivot_wider(id_cols = run_accession, names_from = mate, values_from = -c(run_accession, mate), names_glue = "{.value}_{str_to_lower(mate)}")

if (has_cols(seqkit_wide, c("num_seqs_r1", "num_seqs_r2"))) {
  seqkit_wide <- seqkit_wide %>%
    mutate(paired_reads_match = num_seqs_r1 == num_seqs_r2,
           num_read_pairs_r12 = if_else(paired_reads_match, num_seqs_r1, NA_real_),
           num_seqs_r12 = num_seqs_r1 + num_seqs_r2)
}

if (has_cols(seqkit_wide, c("sum_len_r1", "sum_len_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(sum_len_r12 = sum_len_r1 + sum_len_r2)
}

if (has_cols(seqkit_wide, c("sum_n_r1", "sum_n_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(sum_n_r12 = sum_n_r1 + sum_n_r2)
}

if (has_cols(seqkit_wide, c("sum_gap_r1", "sum_gap_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(sum_gap_r12 = sum_gap_r1 + sum_gap_r2)
}

if (has_cols(seqkit_wide, c("min_len_r1", "min_len_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(min_len_r12 = pmin(min_len_r1, min_len_r2, na.rm = TRUE))
}

if (has_cols(seqkit_wide, c("max_len_r1", "max_len_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(max_len_r12 = pmax(max_len_r1, max_len_r2, na.rm = TRUE))
}

if (has_cols(seqkit_wide, c("avg_len_r1", "avg_len_r2", "num_seqs_r1", "num_seqs_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(avg_len_r12 = weighted_mean_pair(avg_len_r1, avg_len_r2, num_seqs_r1, num_seqs_r2))
}

if (has_cols(seqkit_wide, c("avg_qual_r1", "avg_qual_r2", "sum_len_r1", "sum_len_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(avg_qual_r12 = weighted_mean_pair(avg_qual_r1, avg_qual_r2, sum_len_r1, sum_len_r2))
}

if (has_cols(seqkit_wide, c("gc_percent_r1", "gc_percent_r2", "sum_len_r1", "sum_len_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(gc_percent_r12 = weighted_mean_pair(gc_percent_r1, gc_percent_r2, sum_len_r1, sum_len_r2))
}

if (has_cols(seqkit_wide, c("q20_percent_r1", "q20_percent_r2", "sum_len_r1", "sum_len_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(q20_percent_r12 = weighted_mean_pair(q20_percent_r1, q20_percent_r2, sum_len_r1, sum_len_r2))
}

if (has_cols(seqkit_wide, c("q30_percent_r1", "q30_percent_r2", "sum_len_r1", "sum_len_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(q30_percent_r12 = weighted_mean_pair(q30_percent_r1, q30_percent_r2, sum_len_r1, sum_len_r2))
}

if (has_cols(seqkit_wide, c("q1_r1", "q1_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(q1_mean_r12 = rowMeans(pick(q1_r1, q1_r2), na.rm = TRUE))
}

if (has_cols(seqkit_wide, c("q2_r1", "q2_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(q2_mean_r12 = rowMeans(pick(q2_r1, q2_r2), na.rm = TRUE))
}

if (has_cols(seqkit_wide, c("q3_r1", "q3_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(q3_mean_r12 = rowMeans(pick(q3_r1, q3_r2), na.rm = TRUE))
}

if (has_cols(seqkit_wide, c("n50_r1", "n50_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(n50_mean_r12 = rowMeans(pick(n50_r1, n50_r2), na.rm = TRUE))
}

if (has_cols(seqkit_wide, c("n50_num_r1", "n50_num_r2"))) {
  seqkit_wide <- seqkit_wide %>% mutate(n50_num_r12 = n50_num_r1 + n50_num_r2)
}

seqkit_wide <- seqkit_wide %>%
  arrange(run_accession) %>%
  relocate(any_of(c("run_accession", "paired_reads_match", "num_read_pairs_r12", "num_seqs_r12", "sum_len_r12",
                    "avg_len_r12", "min_len_r12", "max_len_r12", "gc_percent_r12", "q20_percent_r12",
                    "q30_percent_r12", "avg_qual_r12", "sum_n_r12", "sum_gap_r12")))

dir_create(path_dir(opt$output))
write_tsv(seqkit_wide, opt$output)

message("Written aggregated table to: ", opt$output)
