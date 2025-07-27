#!/bin/bash
#SBATCH --job-name=post1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=domenico.zaza@polito.it
#SBATCH --partition=global
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=post1.log
#SBATCH --mem-per-cpu=1024M
module load fftw
module load openblas
mpirun POST.x
