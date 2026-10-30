import os
import csv
import pandas as pd
import yaml
import more_itertools

with open(os.path.join(os.getcwd(), '01_environment.yml')) as f:
    environment = yaml.load(f, Loader = yaml.SafeLoader)
data_path = environment['paths']['data']
scratch_path = environment['paths']['scratch']
code_path = environment['paths']['code']

batch_size = 16

# Set directory containing MRI data and output directory
data_dir = os.path.join(data_path, 'BIDS')
output_dir = os.path.join(scratch_path, "MRI_Processing/fsscripts")
os.makedirs(output_dir, exist_ok = True)
output_file_dir = os.path.join(scratch_path, "MRI_Processing/fsoutput")
os.makedirs(output_file_dir, exist_ok = True)
logs_dir = os.path.join(scratch_path, "MRI_Processing/fslogs")
os.makedirs(logs_dir, exist_ok = True)
# List MRI files matching the conditions
mri_files = []
for root, dirs, files in os.walk(data_dir):
    T1s = []
    print("Current directory:", root)
    for file in files:
        print("Current file:", file)
        if "T1w.nii.gz" in file and "sub-" in root:
            subject_id = root.split("sub-")[-1].split("/")[0]
            print("Subject ID:", subject_id)
            full_path = os.path.join(root, file)
            print("Full path:", full_path)
            T1s.append(full_path)
    if T1s:
        most_recent = max(T1s)  #some people have multiple t1 files since they got redone, they're named alphabetically (IE MPRAGEa is more recent than MPRAGE)
        mri_files.append(most_recent)
print(len(mri_files))
#print(mri_files)

# Generate and submit Slurm job scripts for each MRI
mri_files_batched = list(more_itertools.batched(mri_files, batch_size))
batch_number = 0
for batch in mri_files_batched:
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=recon-all_{batch_number+1}
#SBATCH --output={logs_dir}/recon-all_{batch_number+1}_jobID_%j.out
#SBATCH --error={logs_dir}/recon-all_{batch_number+1}_jobID_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={batch_size}
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00

module load freesurfer/8.2.0-1

# Run FreeSurfer commands to process MRI
"""
    for mri_file in mri_files_batched[batch_number]:
        mri_name = os.path.splitext(os.path.basename(mri_file))[0]
        # Extract subject and session IDs from file path
        path_parts = os.path.dirname(mri_file).split("/")
        subject_id = path_parts[-3]
        session_id = path_parts[-2]
        # Generate output directory for the session
        session_output_dir = os.path.join(output_file_dir, subject_id, session_id, "output")
        # Create session output directory if it doesn't exist
        os.makedirs(session_output_dir, exist_ok=True)

        slurm_script += f'\nrecon-all -i {mri_file} -sd {session_output_dir} -subjid {subject_id} -all -clean &'
    slurm_script += f'\nwait'
    
    script_file = os.path.join(output_dir, f"recon-all_{batch_number+1}.sh")
    with open(script_file, "w") as file:
        file.write(slurm_script)
    
    batch_number += 1


dkt_script = f"""#!/bin/bash
#SBATCH --job-name=dkt_stats
#SBATCH --output={logs_dir}/dkt_stats_%j.log
#SBATCH --error={logs_dir}/dkt_stats_%j.err
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
        LEFT_STATS_FILE="$OUTPUT_DIR/${{subject_id}}_${{session_id}}_lh.txt"
        RIGHT_STATS_FILE="$OUTPUT_DIR/${{subject_id}}_${{session_id}}_rh.txt"

        mris_anatomical_stats -a "$SUBJECTS_DIR/$subject_id/label/lh.aparc.DKTatlas.annot" -f "$LEFT_STATS_FILE" "$subject_id" "lh"
        mris_anatomical_stats -a "$SUBJECTS_DIR/$subject_id/label/rh.aparc.DKTatlas.annot" -f "$RIGHT_STATS_FILE" "$subject_id" "rh"
    done
done
"""

dkt_script_file = os.path.join(code_path, f"05_DKT_stats.sh")
with open(dkt_script_file, "w") as file:
    file.write(dkt_script)