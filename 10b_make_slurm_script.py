import os
import yaml

with open(os.path.join(os.getcwd(), '01_environment.yml')) as f:
    environment = yaml.load(f, Loader = yaml.SafeLoader)
data_path = environment['paths']['data']
scratch_path = environment['paths']['scratch']
code_path = environment['paths']['code']

script_dir = os.path.join(code_path, "10c_slurm.sh")
logs_dir = os.path.join(scratch_path, 'logs')
os.makedirs(logs_dir, exist_ok = True)

slurm_script = f"""#!/bin/bash
#SBATCH --job-name=10_RRHO_map
#SBATCH --output={logs_dir}/10_RRHO_map_jobID_%j.out
#SBATCH --error={logs_dir}/10_RRHO_map_jobID_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /scratch/tdelong/miniforge3/etc/profile.d/conda.sh
conda activate MRI_RNA
Rscript {code_path}/10_RRHO_map.R
"""

with open(script_dir, "w") as file:
    file.write(slurm_script)