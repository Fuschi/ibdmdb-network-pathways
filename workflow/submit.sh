#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --account=applicata
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=workflow/work/logs/controller/%x_%j.out
#SBATCH --error=workflow/work/logs/controller/%x_%j.err
#SBATCH --constraint=blade

set -euo pipefail

cd /scratch/applicata/alessandro.fuschi2/ibdmdb-network-pathways

source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-slurm

mkdir -p workflow/work/logs/controller

snakemake --profile workflow/profiles/slurm &> workflow/work/logs/controller/snakemake_${SLURM_JOB_ID}.log

# test on a single file
#snakemake --profile workflow/profiles/slurm \
#    workflow/work/raw/SRR5935751_R1.fastq.gz \
#    &> workflow/work/logs/controller/snakemake_${SLURM_JOB_ID}.log

