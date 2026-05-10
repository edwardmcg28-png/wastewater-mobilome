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

hash -r
unset PYTHONPATH

TSV=/QRISdata/Q6837/samples.tsv
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
DONE="${STAGE2}/tpm.gene.txt"

LOCKDIR=/QRISdata/Q6837/argoap/.db_lock
READY=/QRISdata/Q6837/argoap/.db_ready

echo "============================================================"
echo "START: ${SAMPLE}   $(date)"
echo "HOST:  $(hostname)"
echo "CPUS:  ${SLURM_CPUS_PER_TASK}"
echo "ENV:   /scratch/user/uqnzhai/conda_envs/args_oap"
echo "python: $(which python)"; python -V
echo "args_oap: $(command -v args_oap || echo MISSING)"
echo "diamond:  $(command -v diamond || echo MISSING)"
echo "bwa:      $(command -v bwa || echo MISSING)"
echo "samtools: $(command -v samtools || echo MISSING)"
echo "blastp:   $(command -v blastp || echo MISSING)"
blastp -version | head -n 2 || true
diamond version || true
echo "============================================================"

command -v args_oap >/dev/null 2>&1 || { echo "[ERROR] args_oap not found. Abort."; exit 10; }
[[ -s "$R1" ]] || { echo "[ERROR] Missing R1: $R1"; exit 20; }
[[ -s "$R2" ]] || { echo "[ERROR] Missing R2: $R2"; exit 21; }

if [[ -s "$DONE" ]]; then
  echo "[SKIP] Found done file: $DONE"
  echo "DONE: ${SAMPLE}   $(date)"
  exit 0
fi

mkdir -p "$INDIR" "$OUTDIR"

ln -sf "$R1" "${INDIR}/${SAMPLE}_1.fq.gz"
ln -sf "$R2" "${INDIR}/${SAMPLE}_2.fq.gz"
ln -sf "$R1" "${INDIR}/${SAMPLE}_clean_1.fastq.gz"
ln -sf "$R2" "${INDIR}/${SAMPLE}_clean_2.fastq.gz"

echo "[INFO] Input links:"
ls -lh "$INDIR"

# 全局 DB 初始化锁（避免并发 Building databases 冲突）
if [[ ! -f "$READY" ]]; then
  if mkdir "$LOCKDIR" 2>/dev/null; then
    echo "[DB-INIT] I got the lock. I will initialize the default DB (may take time)..."
    trap 'rmdir "$LOCKDIR" 2>/dev/null || true' EXIT
  else
    echo "[DB-WAIT] Waiting for DB init to complete..."
    while [[ ! -f "$READY" ]]; do
      sleep 15
    done
    echo "[DB-WAIT] READY detected. Continue."
  fi
fi

echo "[RUN] stage_one"
args_oap stage_one \
  -i "$INDIR" \
  -o "$STAGE1" \
  -t "${SLURM_CPUS_PER_TASK}" \
  -f fq

META="${STAGE1}/metadata.txt"
if [[ ! -s "$META" ]]; then
  echo "[ERROR] stage_one finished but metadata.txt missing/empty: $META"
  echo "[DEBUG] stage_one listing (top 200):"
  find "$STAGE1" -maxdepth 2 -type f -ls | head -n 200 || true
  exit 30
fi

# 若本任务持锁，则标记 DB 已初始化
if [[ -d "$LOCKDIR" && ! -f "$READY" ]]; then
  echo "[DB-INIT] stage_one OK -> mark READY"
  touch "$READY"
  rmdir "$LOCKDIR" 2>/dev/null || true
fi

echo "[INFO] metadata.txt lines:"
wc -l "$META" || true
head -n 5 "$META" || true

echo "[RUN] stage_two"
args_oap stage_two \
  -i "$STAGE1" \
  -o "$STAGE2" \
  -t "${SLURM_CPUS_PER_TASK}"

if [[ ! -s "$DONE" ]]; then
  echo "[ERROR] stage_two finished but done file missing/empty: $DONE"
  echo "[DEBUG] stage_two listing (top 200):"
  find "$STAGE2" -maxdepth 2 -type f -ls | head -n 200 || true
  exit 40
fi

echo "[OK] stage_two output: $DONE"
echo "[INFO] stage_two gene outputs:"
ls -lh "$STAGE2"/*gene*.txt 2>/dev/null || true

echo "DONE: ${SAMPLE}   $(date)"
