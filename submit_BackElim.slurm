#!/bin/bash -l
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dimml002@umn.edu
#SBATCH --job-name RF.slurm
#SBATCH -o ./jobs/%j.training_RF.out
#SBATCH -e ./jobs/%j.training_RF.err
#SBATCH -p amdsmall,amdlarge

set -e

cd $SLURM_SUBMIT_DIR

micromamba activate random_forest

snakemake -s RanFoG.smk \
    --profile slurm.RF \
    --keep-going \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --configfile configfile.yaml

# Define the number of times to run submit.sh
num_runs=45

# Loop to run submit.sh multiple times
for ((i=1; i<=$num_runs; i++))
do
	echo "Running iteration $i"
	snakemake -s BackElim.smk \
	    --profile slurm.RF \
	    --keep-going \
	    --rerun-incomplete \
	    --rerun-triggers mtime \
	    --configfile configfile.yaml \
	    --restart-times 3
	
	python3 rename.py

	wait
done

python3 rename.py
