#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1

#SBATCH --account=adlerp     
#SBATCH --partition=lonepeak
#SBATCH --job-name=storage_growth

#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=atredenn@gmail.com 
#SBATCH -o out.%j 
#SBATCH -e err.%j

export EXE=./rwrapper.sh
export WORK_DIR=$HOME/coexist_stability
export FILENAME=sim_storage_effect_growth_rates.R
export SCRATCH_DIR=/scratch/local/$SLURM_JOBID
export SCRIPT_DIR=$WORK_DIR/Rfiles
export OUT_DIR=$WORK_DIR/results

# Load R (version 3.2.3)
module load R/3.2.3.omp.$USER

# Run an array of serial jobs
export OMP_NUM_THREADS=1

for (( i=1; i < $SLURM_TASKS_PER_NODE ; i++ )); \
 do echo $i $EXE $WORK_DIR/$FILENAME $i ; \
done > my.config.$UUFSCELL.$SLURM_JOBID

# Run a task on each core
cd $WORK_DIR
srun --multi-prog my.config.$UUFSCELL.$SLURM_JOBID

# Clean-up the root scratch dir
rm -rf $SCRATCH_DIR