#!/bin/bash

#SBATCH --job-name=tiling

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --exclusive

. /etc/profile
export OMP_NUM_THREADS=28
export MKL_NUM_THREADS=28

cd /home/cp2021cp202107

srun python3 auto.py -t 4 6 10
