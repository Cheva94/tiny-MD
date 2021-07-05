#!/bin/bash

#SBATCH --job-name=roof

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --exclusive

. /etc/profile
export OMP_NUM_THREADS=28
export MKL_NUM_THREADS=28

srun make clean && make && TMPDIR=/home/cp2021cp202107 /opt/cuda/11.2.2/nsight-compute-2020.3.1/ncu -c 10 -o 2048 --set full ./Qdah && make clean
