import os
import subprocess
import yaml
import more_itertools
import glob

with open(os.path.join(os.getcwd(), '01_environment.yml')) as f:
    environment = yaml.load(f, Loader=yaml.SafeLoader)
data_path = environment['paths']['data']
scratch_path = environment['paths']['scratch']
code_path = environment['paths']['code']
batch_size = 16

data_dir = os.path.join(data_path, 'BIDS')
output_dir = os.path.join(scratch_path, "MRI_Processing/fsscripts_retry")
os.makedirs(output_dir, exist_ok=True)
output_file_dir = os.path.join(scratch_path, "MRI_Processing/fsoutput")
logs_dir = os.path.join(scratch_path, "MRI_Processing/fslogs")
for f in glob.glob(os.path.join(output_dir, "*.sh")):
    os.remove(f)

# all files
mri_files = []
for root, dirs, files in os.walk(data_dir):
    T1s = []
    for file in files:
        if "T1w.nii.gz" in file and "sub-" in root:
            T1s.append(os.path.join(root, file))
    if T1s:
        mri_files.append(max(T1s))

# check if completed
failed_files = []
completed = 0

for mri_file in mri_files:
    path_parts = os.path.dirname(mri_file).split("/")
    subject_id = path_parts[-3]
    session_id = path_parts[-2]
    session_output_dir = os.path.join(output_file_dir, subject_id, session_id, "output")

    done_marker = os.path.join(session_output_dir, subject_id, "scripts", "recon-all.done")
    log_file = os.path.join(session_output_dir, subject_id, "scripts", "recon-all.log")

    finished_ok = False
    if os.path.exists(done_marker) and os.path.exists(log_file):
        with open(log_file, "r") as f:
            log_content = f.read()
        finished_ok = "finished without error" in log_content

    if finished_ok:
        completed += 1
    else:
        failed_files.append(mri_file)

print(f"Completed: {completed} / {len(mri_files)}")
print(f"Needs rerun: {len(failed_files)}")
for f in failed_files:
    print("  ", f)
if not failed_files:
    print("None failed, run 05_DKT_stats")

# new slurm scripts for the missed ones
failed_batched = list(more_itertools.batched(failed_files, batch_size))
for batch_number, batch in enumerate(failed_batched):
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=recon-all_retry_{batch_number+1}
#SBATCH --output={logs_dir}/recon-all_retry_{batch_number+1}_jobID_%j.out
#SBATCH --error={logs_dir}/recon-all_retry_{batch_number+1}_jobID_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={batch_size}
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
module load freesurfer/8.2.0-1
"""
    for mri_file in batch:
        path_parts = os.path.dirname(mri_file).split("/")
        subject_id = path_parts[-3]
        session_id = path_parts[-2]
        session_output_dir = os.path.join(output_file_dir, subject_id, session_id, "output")
        os.makedirs(session_output_dir, exist_ok=True)
        # recon-all exits when -i is given and the subject folder already exists, so clear the failed run first
        subject_dir = os.path.join(session_output_dir, subject_id)
        slurm_script += f'\n(rm -rf {subject_dir} && recon-all -i {mri_file} -sd {session_output_dir} -subjid {subject_id} -all -clean) &'
    slurm_script += '\nwait'

    script_file = os.path.join(output_dir, f"recon-all_retry_{batch_number+1}.sh")
    with open(script_file, "w") as f:
        f.write(slurm_script)