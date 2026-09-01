import os
import yaml
import itertools
import pandas as pd

with open('01_environment.yml') as f:
    env = yaml.safe_load(f)
scratch_path = env['paths']['scratch']
code_path = env['paths']['code']

dkt = pd.read_csv(f"{scratch_path}/processed_data/DKT_stats_combat_averaged.csv")
region_cols = [c for c in dkt.columns if c not in ('Subject', 'Session')]
tissues = ['dorsolateralprefrontalcortex', 'frontalcortex', 'temporalcortex', 'posteriorcingulatecortex', 'Headofcaudatenucleus']

script_dir = os.path.join(scratch_path, "08_slurm")
jobs_dir = os.path.join(scratch_path, '08_dkt_job_chunks')
logs_dir = os.path.join(scratch_path, 'logs/08_dkt_dream')
os.makedirs(script_dir, exist_ok=True)
os.makedirs(jobs_dir, exist_ok=True)
os.makedirs(logs_dir, exist_ok=True)

jobs = list(itertools.product(tissues, region_cols))
total_jobs = len(jobs)
#joblist_path = os.path.join(scratch_path, 'dkt_main_joblist.txt')
print(f"{total_jobs} total runs of dream")

tasks_per_job = 96
chunks = max(1, (total_jobs + tasks_per_job - 1) // tasks_per_job)

print(f"dream runs per job: {tasks_per_job}")
print(f"jobs: {chunks}")

for i in range(chunks):
    start = i * tasks_per_job
    end = min((i + 1) * tasks_per_job, total_jobs)
    chunk_jobs = jobs[start:end]
    chunk_id = f"{i + 1}"

    chunk_file = os.path.join(jobs_dir, f'joblist_{chunk_id}.txt')
    with open(chunk_file, 'w') as f:
        for tissue, region in chunk_jobs:
            f.write(f"{tissue}\t{region}\n")

    slurm_script = f'''#!/bin/bash
#SBATCH --job-name=dkt_regional_{chunk_id}
#SBATCH --array=1-{len(chunk_jobs)}
#SBATCH --time=04:00:00
#SBATCH --output={logs_dir}/dkt_dream_{chunk_id}_%a_jobID_%j.out
#SBATCH --error={logs_dir}/dkt_dream_{chunk_id}_%a_jobID_%j.err
#SBATCH --ntasks={tasks_per_job}
#SBATCH --nodes=1

CHUNK_ID={chunk_id}
LINE=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {chunk_file})
TISSUE=$(echo "$LINE" | cut -f1)
REGION=$(echo "$LINE" | cut -f2)

source /scratch/tdelong/miniforge3/etc/profile.d/conda.sh

Rscript 08_dkt_regions_dream.R "$TISSUE" "$REGION"'''

    with open(os.path.join(script_dir, f'submit_chunk_{i+1}.sh'), "w") as file:
        file.write(slurm_script)