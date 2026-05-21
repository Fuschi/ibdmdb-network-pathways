# -----------------------------------------------------------------------------
# Helper function
# -----------------------------------------------------------------------------

def fastq_url(wildcards):
    fastq_ftp = metadata.loc[metadata["run_accession"] == wildcards.run, "fastq_ftp"].iloc[0]
    fastq_url = fastq_ftp.split(";")[int(wildcards.mate) - 1]

    return "ftp://" + fastq_url


# -----------------------------------------------------------------------------
# Rule
# -----------------------------------------------------------------------------

rule download_fastq:
    output:
        WORK + "/raw/{run}_R{mate}.fastq.gz"
    params:
        url=fastq_url
    log:
        WORK + "/logs/download/{run}_R{mate}.log"
    benchmark:
        WORK + "/benchmarks/download/{run}_R{mate}.tsv"
    resources:
        mem_mb=10000,
        runtime=60,
        constraint="blade"
    shell:
        """
        wget \
            --continue \
            --tries=3 \
            --output-document {output} \
            {params.url} \
            > {log} 2>&1
        """
