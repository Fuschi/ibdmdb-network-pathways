library(tidyverse)

# -----------------------------------------------------------------------------
# Input and output files from Snakemake
# -----------------------------------------------------------------------------

metaphlan_merged_file <- snakemake@input[["merged"]]

species_abundance_file <- snakemake@output[["abundance"]]
species_taxonomy_file <- snakemake@output[["taxonomy"]]

# -----------------------------------------------------------------------------
# Read merged MetaPhlAn table
# -----------------------------------------------------------------------------

metaphlan_merged <- read_tsv(
    metaphlan_merged_file,
    comment = "#",
    show_col_types = FALSE
)

# -----------------------------------------------------------------------------
# Keep only species-level profiles
# -----------------------------------------------------------------------------

metaphlan_species <- metaphlan_merged %>%
    filter(str_detect(clade_name, "s__")) %>%
    filter(!str_detect(clade_name, "t__")) %>%
    mutate(species = str_extract(clade_name, "s__[^|]+$")) %>%
    mutate(species = str_remove(species, "^s__"))

# -----------------------------------------------------------------------------
# Build abundance matrix: samples x species
# -----------------------------------------------------------------------------

species_abundance_matrix <- metaphlan_species %>%
    select(-clade_name) %>%
    column_to_rownames("species") %>%
    t() %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("sample_id") %>%
    as_tibble() %>%
    mutate(sample_id = str_remove(sample_id, "\\.metaphlan$")) %>%
    arrange(sample_id)

# -----------------------------------------------------------------------------
# Build taxonomy matrix: species x taxonomy ranks
# -----------------------------------------------------------------------------

species_taxonomy_matrix <- metaphlan_species %>%
    select(clade_name) %>%
    separate(
        col = clade_name,
        sep = "\\|",
        fill = "right",
        into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
        remove = FALSE
    ) %>%
    mutate(
        domain = str_remove(domain, "^k__"),
        phylum = str_remove(phylum, "^p__"),
        class  = str_remove(class,  "^c__"),
        order  = str_remove(order,  "^o__"),
        family = str_remove(family, "^f__"),
        genus  = str_remove(genus,  "^g__"),
        species = str_remove(species, "^s__")
    ) %>%
    relocate(species, genus, family, order, class, phylum, domain, clade_name) %>%
    distinct()

# -----------------------------------------------------------------------------
# Sanity check: abundance columns and taxonomy rows must be aligned
# -----------------------------------------------------------------------------

if (!all(colnames(species_abundance_matrix)[-1] == species_taxonomy_matrix$species)) {
    stop("Species names in the abundance matrix and taxonomy matrix are not aligned.")
}

# -----------------------------------------------------------------------------
# Write output files
# -----------------------------------------------------------------------------

write_tsv(species_abundance_matrix, species_abundance_file)
write_tsv(species_taxonomy_matrix, species_taxonomy_file)
