#!/bin/bash -l
#SBATCH -t 6:00:00
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

conda activate nugen

snakemake -s RanFoG.smk \
    --profile slurm.RF \
    --keep-going \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --configfile configfile.yaml \
