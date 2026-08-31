#!/bin/bash
#SBATCH --job-name=fstest
#SBATCH --output=/scratch/tdelong/MRI_RNA_Scratch/MRI_Processing/fslogs/test.out
#SBATCH --error=/scratch/tdelong/MRI_RNA_Scratch/MRI_Processing/fslogs/test.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00

module load freesurfer/8.2.0-1

# Run FreeSurfer commands to process MRI

recon-all -i /project/rrg-shreejoy/tdelong/MRI_RNA_Data/BIDS/sub-10042633/ses-0/anat/sub-10042633_ses-0_acq-200902113DMPRAGEa_T1w.nii.gz -sd /scratch/tdelong/MRI_RNA_Scratch/MRI_Processing/fsoutput/sub-10042633/ses-0/output -subjid sub-10042633 -all 