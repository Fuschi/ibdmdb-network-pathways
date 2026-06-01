# -----------------------------------------------------------------------------
# Collect clean FASTQ statistics with SeqKit
# -----------------------------------------------------------------------------

rule collect_clean_stats:
    input:
        WORK + "/clean/{run}_R{mate}.fastq.gz"
    output:
        WORK + "/qc/clean/seqkit/{run}_R{mate}.tsv"
    log:
        WORK + "/logs/clean/seqkit/{run}_R{mate}.log"
    benchmark:
        WORK + "/benchmarks/clean/seqkit/{run}_R{mate}.tsv"
    threads: 2
    resources:
        cpus_per_task = 2,
        mem_mb = 2000,
        runtime = 30
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        seqkit stats \
            --all \
            --tabular \
            --threads {threads} \
            {input} \
            > {output} \
            2> {log}
        """

# -----------------------------------------------------------------------------
# Aggregate clean paired-end SeqKit statistics
# -----------------------------------------------------------------------------

rule aggregate_clean_seqkit_stats:
    input:
        expand(WORK + "/qc/clean/seqkit/{run}_R{mate}.tsv", run=RUNS, mate=["1", "2"])
    output:
        "results/clean_paired_seqkit_stats.tsv"
    log:
        WORK + "/logs/clean/seqkit/aggregate_paired_seqkit_stats.log"
    benchmark:
        WORK + "/benchmarks/clean/seqkit/aggregate_paired_seqkit_stats.tsv"
    resources:
        mem_mb = 2000,
        runtime = 30
    container:
        CONTAINERS["r_postprocess"]
    shell:
        """
        Rscript workflow/script/aggregate_paired_seqkit_stats.R \
            --input-dir workflow/work/qc/clean/seqkit \
            --output {output} \
            > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# Aggregate KneadData log statistics
# -----------------------------------------------------------------------------

rule aggregate_kneaddata_logs:
    input:
        expand(WORK + "/logs/clean/{run}.log", run=RUNS)
    output:
        "results/kneaddata_summary.tsv"
    log:
        WORK + "/logs/clean/kneaddata_summary.log"
    benchmark:
        WORK + "/benchmarks/clean/kneaddata_summary.tsv"
    resources:
        mem_mb = 2000,
        runtime = 30
    container:
        CONTAINERS["r_postprocess"]
    shell:
        """
        Rscript workflow/script/summarise_kneaddata_logs.R \
            --input-dir workflow/work/logs/clean \
            --output {output} \
            > {log} 2>&1
        """

