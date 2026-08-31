import os
import yaml

with open(os.path.join(os.getcwd(), '01_environment.yml')) as f:
    environment = yaml.load(f, Loader = yaml.SafeLoader)
data_path = environment['paths']['data']
scratch_path = environment['paths']['scratch']
code_path = environment['paths']['code']

script_dir = os.path.join(code_path, "07c_slurm.sh")
logs_dir = os.path.join(scratch_path, 'logs')
os.makedirs(logs_dir, exist_ok = True)

slurm_script = f"""#!/bin/bash
#SBATCH --job-name=07_demographics_globalfeatures
#SBATCH --output={logs_dir}/07_demographics_globalfeatures_jobID_%j.out
#SBATCH --error={logs_dir}/07_demographics_globalfeatures_jobID_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00


"""