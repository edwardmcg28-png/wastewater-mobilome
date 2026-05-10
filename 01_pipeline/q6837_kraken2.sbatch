#!/bin/bash
#SBATCH --job-name=q6837_kraken2
#SBATCH --array=1-43%8
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/kraken2_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/kraken2_%A_%a.err

set -euo pipefail
unset PYTHONHOME PYTHONPATH PYTHONUSERBASE PYTHONSTARTUP
export PYTHONNOUSERSITE=1

module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh
conda activate /scratch/user/uqnzhai/conda_envs/kraken2_bracken

TSV=/QRISdata/Q6837/samples.tsv
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TSV" | cut -f1)
R1=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_1.fastq.gz
R2=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_2.fastq.gz
DB=/QRISdata/Q6636/metawrap/kraken2_db
OUTDIR=/QRISdata/Q6837/kraken2/${SAMPLE}
DONE="${OUTDIR}/${SAMPLE}.bracken.species.txt"

[[ ! -s "$R1" ]] && echo "[SKIP] No reads: $SAMPLE" && exit 0
[[ -s "$DONE" ]] && echo "[SKIP] Done: $SAMPLE" && exit 0

mkdir -p "$OUTDIR"
echo "START: $SAMPLE  $(date)"

kraken2 \
    --db "$DB" \
    --threads ${SLURM_CPUS_PER_TASK} \
    --paired \
    --gzip-compressed \
    --report "${OUTDIR}/${SAMPLE}.kraken2.report" \
    --output "${OUTDIR}/${SAMPLE}.kraken2.output" \
    "$R1" "$R2"

bracken \
    -d "$DB" \
    -i "${OUTDIR}/${SAMPLE}.kraken2.report" \
    -o "${OUTDIR}/${SAMPLE}.bracken.species.txt" \
    -r 100 \
    -l S \
    -t 10

bracken \
    -d "$DB" \
    -i "${OUTDIR}/${SAMPLE}.kraken2.report" \
    -o "${OUTDIR}/${SAMPLE}.bracken.genus.txt" \
    -r 100 \
    -l G \
    -t 10

echo "DONE: $SAMPLE  $(date)"
