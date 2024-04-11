#!/bin/bash 

#SBATCH -J OptiSwan
#SBATCH -p gm_base
#SBATCH -c 1
#SBATCH -o console_init.out
#SBATCH -e error_init.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH --ntasks-per-node 1

srun python Mini-Optimorph_SWAN.py