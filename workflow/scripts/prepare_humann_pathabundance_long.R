library(tidyverse)

input_files <- snakemake@input

read_humann_pathabundance <- function(file) {
  run <- basename(file) |>
    stringr::str_remove("_4_pathabundance\\.tsv$")
  
  readr::read_tsv(file, show_col_types = FALSE) |>
    dplyr::rename(feature = 1, abundance = 2) |>
    dplyr::mutate(run = run) |>
    dplyr::select(run, feature, abundance)
}

humann_long <- purrr::map_dfr(input_files, read_humann_pathabundance) |>
  dplyr::mutate(
    abundance = as.numeric(abundance),
    
    is_stratified = stringr::str_detect(feature, "\\|"),
    feature_main = stringr::str_remove(feature, "\\|.*$"),
    
    taxon_raw = dplyr::if_else(
      is_stratified,
      stringr::str_replace(feature, "^.*\\|", ""),
      NA_character_
    ),
    
    is_unmapped = feature_main == "UNMAPPED",
    is_unintegrated = feature_main == "UNINTEGRATED",
    is_special = is_unmapped | is_unintegrated,
    is_pathway = !is_special,
    
    pathway_id = dplyr::case_when(
      is_pathway ~ stringr::str_extract(feature_main, "^[^:]+"),
      TRUE ~ NA_character_
    ),
    
    pathway_name = dplyr::case_when(
      is_pathway & stringr::str_detect(feature_main, ": ") ~
        stringr::str_replace(feature_main, "^[^:]+: ", ""),
      TRUE ~ NA_character_
    ),
    
    species_raw = dplyr::case_when(
      is_stratified & stringr::str_detect(taxon_raw, "^s__") ~ taxon_raw,
      is_stratified & taxon_raw == "unclassified" ~ "unclassified",
      TRUE ~ NA_character_
    ),
    
    species = dplyr::case_when(
      species_raw == "unclassified" ~ "unclassified",
      !is.na(species_raw) ~ species_raw |>
        stringr::str_remove("^s__") |>
        stringr::str_remove("\\.t__.*$") |>
        stringr::str_replace_all("_", " "),
      TRUE ~ NA_character_
    ),
    
    sgb = dplyr::case_when(
      !is.na(species_raw) & stringr::str_detect(species_raw, "\\.t__") ~
        stringr::str_extract(species_raw, "t__.*$") |>
        stringr::str_remove("^t__"),
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::select(
    run,
    feature,
    feature_main,
    pathway_id,
    pathway_name,
    taxon_raw,
    species,
    sgb,
    is_pathway,
    is_special,
    is_unmapped,
    is_unintegrated,
    is_stratified,
    abundance
  )

readr::write_tsv(humann_long, snakemake@output[["long"]])
