# -----------------------------------------------------------------------------
# Run HUMAnN functional profiling on clean reads using precomputed MetaPhlAn
# taxonomic profiles.
# -----------------------------------------------------------------------------

rule profile_humann:
    input:
        r1 = WORK + "/clean/{run}_R1.fastq.gz",
        r2 = WORK + "/clean/{run}_R2.fastq.gz",
        metaphlan = WORK + "/metaphlan/{run}.metaphlan.tsv"
    output:
        genefamilies = WORK + "/humann/{run}_2_genefamilies.tsv",
        reactions = WORK + "/humann/{run}_3_reactions.tsv",
        pathabundance = WORK + "/humann/{run}_4_pathabundance.tsv"
    log:
        WORK + "/logs/humann/{run}.log"
    benchmark:
        WORK + "/benchmarks/humann/{run}.tsv"
    threads: 16
    resources:
        cpus_per_task = 16,
        mem_mb = 64000,
        runtime = 120
    params:
        container = CONTAINERS["humann"],
        db_sqfs = REFERENCES["humann_db_sqfs"],
        outdir = WORK + "/humann",
        tmp_input = WORK + "/humann/{run}.input.fastq.gz",
        db_mount = "/humann_db"
    shell:
        """
        cat {input.r1} {input.r2} > {params.tmp_input}

        apptainer exec \
            --cleanenv \
            --bind /scratch:/scratch \
            --bind /archive:/archive \
            --bind workflow/bin:/humann_wrappers \
            --bind {params.db_sqfs}:{params.db_mount}:image-src=/ \
            {params.container} \
            humann \
                --input {params.tmp_input} \
                --input-format fastq.gz \
                --output {params.outdir} \
                --output-basename {wildcards.run} \
                --threads {threads} \
                --taxonomic-profile {input.metaphlan} \
                --metaphlan /humann_wrappers \
                --memory-use maximum \
                --nucleotide-database {params.db_mount}/chocophlan/chocophlan \
                --protein-database {params.db_mount}/uniref/uniref \
                --utility-database {params.db_mount}/utility_mapping/utility_mapping \
                --pathways-database {params.db_mount}/utility_mapping/utility_mapping/metacyc_reactions_level4ec_only.uniref.bz2,{params.db_mount}/utility_mapping/utility_mapping/metacyc_pathways_structured_filtered_v24_subreactions \
                --remove-temp-output \
                --log-level DEBUG \
                --o-log {log}

        rm -f {params.tmp_input}
        """

# -----------------------------------------------------------------------------
# Prepare a long-format HUMAnN metabolic pathway abundance table.
# -----------------------------------------------------------------------------
rule prepare_humann_pathabundance_long:
    input:
        expand(WORK + "/humann/{run}_4_pathabundance.tsv", run=RUNS)
    output:
        long = "results/humann_metabolic_pathway_abundance_long.tsv"
    log:
        "workflow/work/logs/humann/prepare_pathabundance_long.log"
    benchmark:
        "workflow/work/benchmarks/humann/prepare_pathabundance_long.tsv"
    threads: 1
    resources:
        cpus_per_task = 1,
        mem_mb = 4000,
        runtime = 30
    container:
        CONTAINERS["r_postprocess"]
    script:
        "../scripts/prepare_humann_pathabundance_long.R"
