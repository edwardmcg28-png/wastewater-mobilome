#!/bin/bash
#SBATCH --job-name=q6837_analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/analysis_%j.out
#SBATCH --error=/QRISdata/Q6837/logs/analysis_%j.err

set -euo pipefail
module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh

conda run -p /scratch/user/uqnzhai/conda_envs/r_analysis \
    python3 /scratch/user/uqnzhai/scripts/01_mge_risk_classification.py

echo "Analysis complete"
