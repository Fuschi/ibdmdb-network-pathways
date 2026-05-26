# -----------------------------------------------------------------------------
# Clean paired-end FASTQ reads with KneadData
# -----------------------------------------------------------------------------

rule clean_reads_kneaddata:
    input:
        r1 = WORK + "/raw/{run}_R1.fastq.gz",
        r2 = WORK + "/raw/{run}_R2.fastq.gz"
    output:
        r1 = WORK + "/clean/{run}_R1.fastq.gz",
        r2 = WORK + "/clean/{run}_R2.fastq.gz"
    params:
        clean = WORK + "/clean/kneaddata",
        scratch = WORK + "/clean/tmp_{run}",
        db = REFERENCES["human_kneaddata"],
        prefix = "{run}"
    log:
        WORK + "/logs/clean/kneaddata/{run}.log"
    benchmark:
        WORK + "/benchmarks/clean/kneaddata/{run}.tsv"
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 240
    container:
        CONTAINERS["kneaddata"]
    shell:
        """
        mkdir -p {params.clean} {params.scratch}

        kneaddata \
            --input1 {input.r1} \
            --input2 {input.r2} \
            --reference-db {params.db} \
            --output {params.clean} \
            --scratch {params.scratch} \
            --output-prefix {params.prefix} \
            --sequencer-source NexteraPE \
            --quality-scores phred33 \
            --threads {threads} \
            --remove-intermediate-output \
            --log-level INFO \
            --log {log}

        pigz -p {threads} {params.clean}/{wildcards.run}_paired_1.fastq
        pigz -p {threads} {params.clean}/{wildcards.run}_paired_2.fastq

        mv {params.clean}/{wildcards.run}_paired_1.fastq.gz {output.r1}
        mv {params.clean}/{wildcards.run}_paired_2.fastq.gz {output.r2}

        rm -rf {params.scratch}
        """

