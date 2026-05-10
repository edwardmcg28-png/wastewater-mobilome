#!/bin/bash
set -euo pipefail
unset PYTHONHOME PYTHONPATH PYTHONUSERBASE PYTHONSTARTUP
export PYTHONNOUSERSITE=1
#SBATCH --job-name=q6837_genomad
#SBATCH --array=1-43%8
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/genomad_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/genomad_%A_%a.err

module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh
conda activate /scratch/user/uqnzhai/conda_envs/genomad

export PATH=/scratch/user/uqnzhai/conda_envs/genomad/bin:$PATH
export LD_LIBRARY_PATH=/scratch/user/uqnzhai/conda_envs/genomad/lib:$LD_LIBRARY_PATH
DB=/scratch/user/uqnzhai/genomad_db/genomad_db

TSV=/QRISdata/Q6837/samples.tsv
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TSV" | cut -f1)
CONTIGS=/QRISdata/Q6837/assemblies/${SAMPLE}_contigs.fa
OUTDIR=/QRISdata/Q6837/genomad/${SAMPLE}
DONE="${OUTDIR}/${SAMPLE}_contigs_summary/${SAMPLE}_contigs_plasmid_summary.tsv"

[[ ! -s "$CONTIGS" ]] && echo "[SKIP] No contigs: $SAMPLE" && exit 0
[[ -s "$DONE" ]] && echo "[SKIP] Done: $SAMPLE" && exit 0

echo "START: $SAMPLE  $(date)"
genomad end-to-end \
    --cleanup \
    --threads ${SLURM_CPUS_PER_TASK} \
    "$CONTIGS" "$OUTDIR" "$DB"

echo "DONE: $SAMPLE  $(date)"
