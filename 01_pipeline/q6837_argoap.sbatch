#!/bin/bash
#SBATCH --job-name=q6837_argoap
#SBATCH --array=1-43%8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/argoap_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/argoap_%A_%a.err

module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh
conda activate /scratch/user/uqnzhai/conda_envs/args_oap

TSV=/QRISdata/Q6837/samples.tsv
SAMPLE=$(awk "NR==$((SLURM_ARRAY_TASK_ID+1))" "$TSV" | cut -f1)
R1=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_1.fastq.gz
R2=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_2.fastq.gz
OUTDIR=/QRISdata/Q6837/argoap/${SAMPLE}
DONE="${OUTDIR}/stage_two/${SAMPLE}.gene.tab.txt"

[[ -s "$DONE" ]] && echo "[SKIP] $SAMPLE" && exit 0

mkdir -p "${OUTDIR}/input"
ln -sf "$R1" "${OUTDIR}/input/${SAMPLE}_1.fq.gz"
ln -sf "$R2" "${OUTDIR}/input/${SAMPLE}_2.fq.gz"

echo "START: $SAMPLE  $(date)"

args_oap stage_one \
    -i "${OUTDIR}/input" \
    -o "${OUTDIR}/stage_one" \
    -t ${SLURM_CPUS_PER_TASK} \
    -f fq

args_oap stage_two \
    -i "${OUTDIR}/stage_one" \
    -o "${OUTDIR}/stage_two" \
    -t ${SLURM_CPUS_PER_TASK}

echo "DONE: $SAMPLE  $(date)"
