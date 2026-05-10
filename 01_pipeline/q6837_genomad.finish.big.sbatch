#!/bin/bash
#SBATCH --job-name=q6837_ge_finish
#SBATCH --array=1-43%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=96:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/genomad_finish_big_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/genomad_finish_big_%A_%a.err

set -euo pipefail
unset PYTHONHOME PYTHONPATH PYTHONUSERBASE PYTHONSTARTUP
export PYTHONNOUSERSITE=1

module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh

ENV=/scratch/user/uqnzhai/conda_envs/genomad
TSV=/QRISdata/Q6837/samples.tsv
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$TSV" | cut -f1)

CONTIGS=/QRISdata/Q6837/assemblies/${SAMPLE}_contigs.fa
OUT=/QRISdata/Q6837/genomad/${SAMPLE}

[[ ! -s "$CONTIGS" ]] && echo "[SKIP] No contigs: $SAMPLE" && exit 0

mkdir -p "$OUT"
echo "START: $SAMPLE  $(date)"

# 如果已经有任何 summary 输出（避免重复跑）
if find "$OUT" -path "*_contigs_summary/*" -type f -size +0c 2>/dev/null | grep -q . ; then
  echo "[SKIP] summary already exists: $SAMPLE"
  exit 0
fi

# 只补跑后半段（nn / aggregated / summary）
# 注意：这些子命令不需要 DB 参数
conda run -p "$ENV" --no-capture-output genomad nn-classification \
  -t ${SLURM_CPUS_PER_TASK} \
  --batch-size 64 \
  "$CONTIGS" "$OUT"

conda run -p "$ENV" --no-capture-output genomad aggregated-classification \
  "$CONTIGS" "$OUT"

conda run -p "$ENV" --no-capture-output genomad summary \
  "$CONTIGS" "$OUT"

# 验收：summary 目录必须出现非空文件
if ! find "$OUT" -path "*_contigs_summary/*" -type f -size +0c 2>/dev/null | grep -q . ; then
  echo "[ERROR] summary missing after run: $SAMPLE"
  exit 3
fi

echo "DONE: $SAMPLE  $(date)"
