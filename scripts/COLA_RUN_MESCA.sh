#!/bin/bash

#SBATCH --job-name=cola_calib
#SBATCH --output=cola_%A_%a.out
#SBATCH --error=cola_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --partition=mesca2
#SBATCH -t 1-12:00:00
#SBATCH --mem-per-cpu=46000
#SBATCH --mail-user=joao.reboucas@unesp.br

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Number of tasks is $SLURM_NTASKS
echo Number of cpus per task is $SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR
module load anaconda3/2020.11
conda activate /scratch/decola/joao.reboucas2/.colaenv/
source start_cola

export OMP_PROC_BIND=close
export OMP_NUM_THREADS=1

LUAFILE=/scratch/decola/joao.reboucas2/COLA_projects/calibration/lua_files/test_${SLURM_ARRAY_TASK_ID}_htable.lua
mpirun -n ${SLURM_NTASKS} --report-bindings --mca btl_tcp_if_include ib0 --bind-to core --map-by node:pe=${OMP_NUM_THREADS} ./FML/COLASolver/nbody $LUAFILE

wait

echo Ending run at `date`
