#!/bin/bash
set -euo pipefail
unset PYTHONHOME PYTHONPATH PYTHONUSERBASE PYTHONSTARTUP
export PYTHONNOUSERSITE=1
#SBATCH --job-name=q6837_rgi
#SBATCH --array=1-43%10
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/rgi_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/rgi_%A_%a.err

module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh
conda activate /scratch/user/uqnzhai/conda_envs/rgi

TSV=/QRISdata/Q6837/samples.tsv
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TSV" | cut -f1)
CONTIGS=/QRISdata/Q6837/assemblies/${SAMPLE}_contigs.fa
OUTDIR=/QRISdata/Q6837/rgi

[[ ! -s "$CONTIGS" ]] && echo "[SKIP] No contigs: $SAMPLE" && exit 0
[[ -s "${OUTDIR}/${SAMPLE}.txt" ]] && echo "[SKIP] Done: $SAMPLE" && exit 0

echo "START: $SAMPLE  $(date)"
rgi main \
    --input_sequence "$CONTIGS" \
    --output_file "${OUTDIR}/${SAMPLE}" \
    --input_type contig \
    --alignment_tool DIAMOND \
    --num_threads ${SLURM_CPUS_PER_TASK} \
    --include_loose \
    --clean

echo "DONE: $SAMPLE  $(date)"
