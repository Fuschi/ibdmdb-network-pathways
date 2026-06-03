# -----------------------------------------------------------------------------
# Taxonomic profiling with MetaPhlAn 4.2.4
# -----------------------------------------------------------------------------

rule profile_metaphlan:
    input:
        r1 = WORK + "/clean/{run}_R1.fastq.gz",
        r2 = WORK + "/clean/{run}_R2.fastq.gz"
    output:
        profile = WORK + "/metaphlan/{run}.metaphlan.tsv",
        mapout = WORK + "/metaphlan/{run}.metaphlan.mapout.bz2"
    log:
        WORK + "/logs/metaphlan/{run}.log"
    benchmark:
        WORK + "/benchmarks/metaphlan/{run}.tsv"
    threads: 12
    resources:
        cpus_per_task = 12,
        mem_mb = 32000,
        runtime = 120
    container:
        CONTAINERS["metaphlan"]
    params:
        db_dir = REFERENCES["metaphlan"]["db_dir"],
        index = REFERENCES["metaphlan"]["index"]
    shell:
        """
        metaphlan \
            {input.r1},{input.r2} \
            --input_type fastq \
            --db_dir {params.db_dir} \
            -x {params.index} \
            --offline \
            --sample_id {wildcards.run} \
            --mapout {output.mapout} \
            --nproc {threads} \
            -t rel_ab_w_read_stats \
            -o {output.profile} \
            > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# Merge MetaPhlAn taxonomic profiles
# -----------------------------------------------------------------------------

rule merge_metaphlan_profiles:
    input:
        expand(WORK + "/metaphlan/{run}.metaphlan.tsv", run=RUNS)
    output:
        "results/metaphlan/metaphlan_merged.tsv",
    log:
        WORK + "/logs/metaphlan_merge/merge_metaphlan_profiles.log"
    benchmark:
        WORK + "/benchmarks/metaphlan_merge/merge_metaphlan_profiles.tsv"
    threads: 1
    resources:
        cpus_per_task = 1,
        mem_mb = 4000,
        runtime = 30
    container:
        CONTAINERS["metaphlan"]
    shell:
        """
        merge_metaphlan_tables.py \
            {input} \
            -o {output} \
            --overwrite \
            > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# Format MetaPhlAn species profiles
# -----------------------------------------------------------------------------

rule format_metaphlan_species:
    input:
        merged = "results/metaphlan/metaphlan_merged.tsv"
    output:
        abundance = "results/metaphlan/metaphlan_species_abundance_matrix.tsv",
        taxonomy = "results/metaphlan/metaphlan_species_taxonomy_matrix.tsv"
    log:
        WORK + "/logs/metaphlan_merge/format_metaphlan_species.log"
    benchmark:
        WORK + "/benchmarks/metaphlan_merge/format_metaphlan_species.tsv"
    threads: 1
    resources:
        cpus_per_task = 1,
        mem_mb = 4000,
        runtime = 30
    container:
        CONTAINERS["r_postprocess"]
    script:
        "../scripts/format_metaphlan_species.R"

