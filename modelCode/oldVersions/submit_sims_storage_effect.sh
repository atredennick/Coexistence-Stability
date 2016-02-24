#!/bin/bash
#SBATCH --array=1-13
#SBATCH --job-name=R_simsstorage_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --constraint=gonium

# Send mail
#SBATCH --mail-type=ALL
#SBATCH --mail-user=atredenn@gmail.com

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH -$SLURM_ARRAY_TASK_ID $HOME/storage_effect_simulations/wrapper_storage_effect_simulations_hpc.R