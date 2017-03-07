#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --tasks-per-node=8

#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
#SBATCH --job-name=att_test

#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=atredenn@gmail.com 
#SBATCH -o out.%j 
#SBATCH -e err.%j

export EXE=./etas_rwrapper.sh
export WORK_DIR=$HOME/coex_stability
export FILENAME=storageeffect_sims_div+envvar_stability_varyetas.R
export SCRATCH_DIR=/scratch/local/$SLURM_JOBID
export SCRIPT_DIR=$WORK_DIR/Rfiles
export OUT_DIR=$WORK_DIR/results

# Load R (version 3.3.2)
module use ~/MyModules
module load R/3.3.2.$USER

# Run an array of serial jobs
export OMP_NUM_THREADS=1

# initial counter for the jobs to allow multiple jobs
START_COUNTER=1
for (( i=0; i < $SLURM_NTASKS ; i++ )); do  
   ip=$((START_COUNTER+i))
   echo $i $EXE $ip $SCRATCH_DIR/$ip $SCRIPT_DIR $OUT_DIR ; \
done > my.config.$UUFSCELL.$SLURM_JOBID

# Run a task on each core
cd $WORK_DIR
srun --multi-prog my.config.$UUFSCELL.$SLURM_JOBID

# Clean-up the root scratch dir
rm -rf $SCRATCH_DIR
