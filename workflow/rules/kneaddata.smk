# -----------------------------------------------------------------------------
# Clean paired-end FASTQ reads with KneadData
# -----------------------------------------------------------------------------

rule clean_reads_kneaddata:
    input:
        r1 = WORK + "/raw/{run}_R1.fastq.gz",
        r2 = WORK + "/raw/{run}_R2.fastq.gz"
    output:
        r1 = WORK + "/clean/{run}_R1.fastq.gz",
        r2 = WORK + "/clean/{run}_R2.fastq.gz",
        unmatched = WORK + "/clean/{run}_unmatched.tar.gz",
        contam = WORK + "/clean/{run}_human_contam.tar.gz"
    params:
        clean = WORK + "/clean",
        scratch = WORK + "/clean/tmp_{run}",
        db = "/archive/extra/daniel.remondini2/bioinf_db/databases/kneaddata_human_hg39_T2T_bowtie2_v0.1/database"
    log:
        WORK + "/logs/clean/{run}.log"
    benchmark:
        WORK + "/benchmarks/clean/{run}.tsv"
    threads: 12
    resources:
        cpus_per_task = 12,
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
            --output-prefix {wildcards.run} \
            --sequencer-source NexteraPE \
            --quality-scores phred33 \
            --threads {threads} \
            --log-level INFO \
            --log {log}

        pigz -p {threads} {params.clean}/{wildcards.run}_paired_1.fastq
        pigz -p {threads} {params.clean}/{wildcards.run}_paired_2.fastq

        mv {params.clean}/{wildcards.run}_paired_1.fastq.gz {output.r1}
        mv {params.clean}/{wildcards.run}_paired_2.fastq.gz {output.r2}

        tar -czf {output.unmatched} \
            -C {params.clean} \
            {wildcards.run}_unmatched_1.fastq \
            {wildcards.run}_unmatched_2.fastq

        tar -czf {output.contam} \
            -C {params.scratch} \
            {wildcards.run}_hg_39_bowtie2_paired_contam_1.fastq \
            {wildcards.run}_hg_39_bowtie2_paired_contam_2.fastq \
            {wildcards.run}_hg_39_bowtie2_unmatched_1_contam.fastq \
            {wildcards.run}_hg_39_bowtie2_unmatched_2_contam.fastq

        rm -f {params.clean}/{wildcards.run}_unmatched_1.fastq
        rm -f {params.clean}/{wildcards.run}_unmatched_2.fastq
        rm -rf {params.scratch}
        """

