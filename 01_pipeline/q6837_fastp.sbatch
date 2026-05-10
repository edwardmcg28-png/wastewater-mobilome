#!/bin/bash
#SBATCH --job-name=q6837_fastp
#SBATCH --array=1-43%10
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/fastp_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/fastp_%A_%a.err

module purge
module load fastp/0.23.2-gcc-11.3.0

TSV=/QRISdata/Q6837/samples.tsv
LINE=$(awk "NR==$((SLURM_ARRAY_TASK_ID+1))" "$TSV")
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE"     | cut -f2)
R2=$(echo "$LINE"     | cut -f3)

OUTDIR=/QRISdata/Q6837/clean_fastp
OUT1="${OUTDIR}/${SAMPLE}_clean_1.fastq.gz"
OUT2="${OUTDIR}/${SAMPLE}_clean_2.fastq.gz"

# 跳过已完成
[[ -s "$OUT1" && -s "$OUT2" ]] && echo "[SKIP] $SAMPLE" && exit 0

echo "START: $SAMPLE  $(date)"
fastp \
    -i "$R1" -I "$R2" \
    -o "$OUT1" -O "$OUT2" \
    -j "${OUTDIR}/${SAMPLE}_fastp.json" \
    -h "${OUTDIR}/${SAMPLE}_fastp.html" \
    --thread 4 \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --detect_adapter_for_pe

echo "DONE: $SAMPLE  $(date)"
