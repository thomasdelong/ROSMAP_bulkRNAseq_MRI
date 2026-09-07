#!/bin/bash
#SBATCH --job-name=10_RRHO_map
#SBATCH --output=/scratch/tdelong/MRI_RNA_Scratch/logs/10_RRHO_map_jobID_%j.out
#SBATCH --error=/scratch/tdelong/MRI_RNA_Scratch/logs/10_RRHO_map_jobID_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00

source /scratch/tdelong/miniforge3/etc/profile.d/conda.sh
conda activate MRI_RNA
Rscript /home/tdelong/MRI_RNA_Code/10_RRHO_map.R
