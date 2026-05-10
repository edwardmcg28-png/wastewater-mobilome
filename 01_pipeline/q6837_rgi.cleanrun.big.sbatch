#!/bin/bash
#SBATCH --job-name=q6837_rgi
#SBATCH --array=1-43%8
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/rgi_big_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/rgi_big_%A_%a.err

set -euo pipefail
unset PYTHONHOME PYTHONPATH PYTHONUSERBASE PYTHONSTARTUP
export PYTHONNOUSERSITE=1

module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh

ENV=/scratch/user/uqnzhai/conda_envs/rgi
TSV=/QRISdata/Q6837/samples.tsv
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TSV" | cut -f1)

CONTIGS=/QRISdata/Q6837/assemblies/${SAMPLE}_contigs.fa
OUTDIR=/QRISdata/Q6837/rgi
OUT=${OUTDIR}/${SAMPLE}.txt

[[ ! -s "$CONTIGS" ]] && echo "[SKIP] No contigs: $SAMPLE" && exit 0
[[ -s "$OUT" ]] && echo "[SKIP] Done: $SAMPLE" && exit 0

mkdir -p "$OUTDIR"
echo "START: $SAMPLE  $(date)"

# 用 conda run 避免 activate 污染
conda run -p "$ENV" --no-capture-output rgi main \
  --input_sequence "$CONTIGS" \
  --output_file "${OUTDIR}/${SAMPLE}" \
  --input_type contig \
  --alignment_tool DIAMOND \
  --num_threads ${SLURM_CPUS_PER_TASK} \
  --include_loose \
  --clean

# 验收：必须有非空 txt
[[ -s "$OUT" ]] || { echo "[ERROR] missing output: $OUT"; exit 2; }

echo "DONE: $SAMPLE  $(date)"
