#might need to just copy/paste this into a bash terminal depending on permissions
script_path=$(yq -r '.paths.scratch' 01_environment.yml)/MRI_Processing/fsscripts/*
for file in $script_path; do
    sbatch $file
done