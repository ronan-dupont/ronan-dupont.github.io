#!/bin/bash
# Example of running python script with a job array
#SBATCH -J Run_test
#SBATCH -p gm_gladys
#SBATCH --account=shoremotion
#SBATCH -c 1
#SBATCH -o console.out
#SBATCH -e erreur.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH --ntasks-per-core 1
# one CPU core per task
# Run python script with a command line argument
srun python optimorph/Mini-Optimorph.py
