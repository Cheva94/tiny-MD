#!/bin/bash

#SBATCH --job-name=qdah

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --exclusive

. /etc/profile
export OMP_NUM_THREADS=28
export MKL_NUM_THREADS=28

srun python3 auto.py -q 4 10 3
