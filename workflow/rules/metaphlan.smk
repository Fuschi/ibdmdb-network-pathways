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
        mem_mb = 24000,
        runtime = 240
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
