#!/bin/bash
#SBATCH --job-name=dkt_stats
#SBATCH --output=/scratch/tdelong/MRI_RNA_Scratch/MRI_Processing/fslogs/dkt_stats_%j.log
#SBATCH --error=/scratch/tdelong/MRI_RNA_Scratch/MRI_Processing/fslogs/dkt_stats_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00

# Load FreeSurfer environment
module load freesurfer/8.2.0-1

# Define base directories
FS_BASE_DIR="$(yq -r '.paths.scratch' 01_environment.yml)/MRI_Processing/fsoutput" 
OUTPUT_DIR="$(yq -r '.paths.scratch' 01_environment.yml)/MRI_Processing/DKTstats"
mkdir -p "$OUTPUT_DIR"

# Loop through subjects in FS_BASE_DIR
for subject_dir in "$FS_BASE_DIR"/sub-*; do
    subject_id=$(basename "$subject_dir")
    
    # Loop through each session
    for session_dir in "$subject_dir"/ses-*; do
        session_id=$(basename "$session_dir")
        # Set SUBJECTS_DIR to the specific session's output directory
        export SUBJECTS_DIR="$session_dir/output"

        # Prepare output file for the session
        LEFT_STATS_FILE="$OUTPUT_DIR/${subject_id}_${session_id}_lh.txt"
        RIGHT_STATS_FILE="$OUTPUT_DIR/${subject_id}_${session_id}_rh.txt"

        mris_anatomical_stats -a "$SUBJECTS_DIR/$subject_id/label/lh.aparc.DKTatlas.annot" -f "$LEFT_STATS_FILE" "$subject_id" "lh"
        mris_anatomical_stats -a "$SUBJECTS_DIR/$subject_id/label/rh.aparc.DKTatlas.annot" -f "$RIGHT_STATS_FILE" "$subject_id" "rh"
    done
done
