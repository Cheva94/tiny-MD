#!/bin/bash

#SBATCH --job-name=test

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --exclusive

. /etc/profile
export OMP_NUM_THREADS=28
export MKL_NUM_THREADS=28

srun make clean && make && ./Qdah && make clean
# srun make clean && make && cuda-memcheck ./Qdah && make clean
