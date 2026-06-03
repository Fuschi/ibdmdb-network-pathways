# -----------------------------------------------------------------------------
# Collect raw FASTQ statistics with SeqKit
# -----------------------------------------------------------------------------

rule collect_raw_stats:
    input:
        WORK + "/raw/{run}_R{mate}.fastq.gz"
    output:
        WORK + "/qc/raw/seqkit/{run}_R{mate}.tsv"
    log:
        WORK + "/logs/raw/seqkit/{run}_R{mate}.log"
    benchmark:
        WORK + "/benchmarks/raw/seqkit/{run}_R{mate}.tsv"
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
# Aggregate raw paired-end SeqKit statistics
# -----------------------------------------------------------------------------

rule aggregate_raw_seqkit_stats:
    input:
        expand(WORK + "/qc/raw/seqkit/{run}_R{mate}.tsv", run=RUNS, mate=["1", "2"])
    output:
        "results/raw_paired_seqkit_stats.tsv"
    log:
        WORK + "/logs/raw/seqkit/aggregate_paired_seqkit_stats.log"
    benchmark:
        WORK + "/benchmarks/raw/seqkit/aggregate_paired_seqkit_stats.tsv"
    resources:
        mem_mb = 2000,
        runtime = 30
    container:
        CONTAINERS["r_postprocess"]
    shell:
        """
        Rscript workflow/scripts/aggregate_paired_seqkit_stats.R \
            --input-dir workflow/work/qc/raw/seqkit \
            --output {output} \
            > {log} 2>&1
        """
