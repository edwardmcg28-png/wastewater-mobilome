#!/bin/bash
#SBATCH --job-name=q6837_argoap
#SBATCH --array=1-43%8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/argoap_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/argoap_%A_%a.err

set -euo pipefail

module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh
conda activate /scratch/user/uqnzhai/conda_envs/args_oap

# 防污染
hash -r
unset PYTHONPATH

TSV=/QRISdata/Q6837/samples.tsv
# 默认TSV第1行是header
SAMPLE=$(awk -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1; exit}' "$TSV")

if [[ -z "${SAMPLE:-}" ]]; then
  echo "[ERROR] SAMPLE empty. Check samples.tsv header/row count. TASK=${SLURM_ARRAY_TASK_ID}"
  exit 2
fi

R1=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_1.fastq.gz
R2=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_2.fastq.gz
OUTDIR=/QRISdata/Q6837/argoap/${SAMPLE}
INDIR=${OUTDIR}/input
STAGE1=${OUTDIR}/stage_one
STAGE2=${OUTDIR}/stage_two
DONE="${STAGE2}/${SAMPLE}.gene.tab.txt"

echo "============================================================"
echo "START: ${SAMPLE}   $(date)"
echo "HOST:  $(hostname)"
echo "CPUS:  ${SLURM_CPUS_PER_TASK}   MEM: ${SLURM_MEM_PER_NODE:-NA}"
echo "ENV:   /scratch/user/uqnzhai/conda_envs/args_oap"
echo "PATH:  $PATH"
echo "python: $(which python)"
python -V
echo "pip:    $(which pip)"
pip -V
echo "args_oap: $(command -v args_oap || echo MISSING)"
echo "diamond:  $(command -v diamond || echo MISSING)"
echo "bwa:      $(command -v bwa || echo MISSING)"
echo "samtools: $(command -v samtools || echo MISSING)"
echo "blastp:   $(command -v blastp || echo MISSING)"
blastp -version | head -n 2 || true
diamond version || true
echo "============================================================"

# 必要命令检查
command -v args_oap >/dev/null 2>&1 || { echo "[ERROR] args_oap not found in env. Abort."; exit 10; }
command -v diamond  >/dev/null 2>&1 || { echo "[ERROR] diamond not found. Abort."; exit 11; }
command -v bwa      >/dev/null 2>&1 || { echo "[ERROR] bwa not found. Abort."; exit 12; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found. Abort."; exit 13; }
command -v blastp   >/dev/null 2>&1 || { echo "[ERROR] blastp not found. Abort."; exit 14; }

# 输入检查
[[ -s "$R1" ]] || { echo "[ERROR] Missing R1: $R1"; exit 20; }
[[ -s "$R2" ]] || { echo "[ERROR] Missing R2: $R2"; exit 21; }

# 已完成则跳过
if [[ -s "$DONE" ]]; then
  echo "[SKIP] Found done file: $DONE"
  echo "DONE: ${SAMPLE}   $(date)"
  exit 0
fi

mkdir -p "$INDIR" "$OUTDIR"

# 同时创建两套命名，最大化兼容性
ln -sf "$R1" "${INDIR}/${SAMPLE}_1.fq.gz"
ln -sf "$R2" "${INDIR}/${SAMPLE}_2.fq.gz"
ln -sf "$R1" "${INDIR}/${SAMPLE}_clean_1.fastq.gz"
ln -sf "$R2" "${INDIR}/${SAMPLE}_clean_2.fastq.gz"

echo "[INFO] Input links:"
ls -lh "$INDIR"

echo "[RUN] stage_one"
args_oap stage_one \
  -i "$INDIR" \
  -o "$STAGE1" \
  -t "${SLURM_CPUS_PER_TASK}" \
  -f fq

# 强制验收：metadata 必须存在且非空
META="${STAGE1}/metadata.txt"
if [[ ! -s "$META" ]]; then
  echo "[ERROR] stage_one finished but metadata.txt missing/empty: $META"
  echo "[DEBUG] stage_one listing (top 200):"
  find "$STAGE1" -maxdepth 2 -type f -ls | head -n 200 || true
  exit 30
fi

echo "[INFO] metadata.txt lines:"
wc -l "$META" || true
head -n 5 "$META" || true

echo "[RUN] stage_two"
args_oap stage_two \
  -i "$STAGE1" \
  -o "$STAGE2" \
  -t "${SLURM_CPUS_PER_TASK}"

# 输出验收
if [[ ! -s "$DONE" ]]; then
  echo "[ERROR] stage_two finished but done file missing/empty: $DONE"
  echo "[DEBUG] stage_two listing (top 200):"
  find "$STAGE2" -maxdepth 2 -type f -ls | head -n 200 || true
  exit 40
fi

echo "DONE: ${SAMPLE}   $(date)"
