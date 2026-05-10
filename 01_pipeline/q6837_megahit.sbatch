#!/bin/bash
#SBATCH --job-name=q6837_megahit
#SBATCH --array=1-43%8
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/megahit_%A_%a.out
#SBATCH --error=/QRISdata/Q6837/logs/megahit_%A_%a.err

module purge
module load megahit/1.2.9-gcccore-11.3.0

TSV=/QRISdata/Q6837/samples.tsv
SAMPLE=$(awk "NR==$((SLURM_ARRAY_TASK_ID+1))" "$TSV" | cut -f1)
R1=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_1.fastq.gz
R2=/QRISdata/Q6837/clean_fastp/${SAMPLE}_clean_2.fastq.gz
OUTDIR=/QRISdata/Q6837/assemblies/${SAMPLE}
FINAL=/QRISdata/Q6837/assemblies/${SAMPLE}_contigs.fa

[[ -s "$FINAL" ]] && echo "[SKIP] $SAMPLE" && exit 0

echo "START: $SAMPLE  $(date)"
megahit -1 "$R1" -2 "$R2" \
    -o "$OUTDIR" \
    --min-contig-len 1000 \
    -t ${SLURM_CPUS_PER_TASK}

# 移动最终contigs到标准位置
cp "${OUTDIR}/final.contigs.fa" "$FINAL"
echo "DONE: $SAMPLE  $(date)"
