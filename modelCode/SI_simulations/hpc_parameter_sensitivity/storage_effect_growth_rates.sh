#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=144
#SBATCH --ntasks=1
#SBATCH --array=1-485
#SBATCH --job-name=storage_growth
#SBATCH --account=adlerp     
#SBATCH --partition=lonepeak  

export FILENAME=sim_storage_effect_growth_rates.R
export SCR_DIR=/scratch/general/lustre/$USER/$SLURM_JOBID
export WORK_DIR=$HOME/coexist_stability/

# Load R (version 3.2.3)
module load R/3.2.3.omp.$USER

# Run an array of serial jobs
export OMP_NUM_THREADS=1

R CMD BATCH -$SLURM_ARRAY_TASK_ID WORK_DIR/FILENAME