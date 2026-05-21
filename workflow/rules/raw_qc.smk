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
    resources:
        mem_mb=2000,
        runtime=30,
        constraint="blade"
    container:
        CONTAINERS["seqkit"]
    shell:
        """
        seqkit stats \
            --all \
            --tabular \
            {input} \
            > {output} \
            2> {log}
        """
