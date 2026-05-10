#!/bin/bash
#SBATCH --job-name=q6837_nmds
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --account=a_obrien
#SBATCH --output=/QRISdata/Q6837/logs/nmds_%j.out
#SBATCH --error=/QRISdata/Q6837/logs/nmds_%j.err

set -euo pipefail
module purge
module load anaconda3/2023.09-0
source /sw/auto/rocky8c/epyc3/software/Anaconda3/2023.09-0/etc/profile.d/conda.sh

conda run -p /scratch/user/uqnzhai/conda_envs/r_env \
    Rscript /scratch/user/uqnzhai/scripts/03_nmds_permanova.R
