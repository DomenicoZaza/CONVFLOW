#!/bin/bash
#SBATCH --job-name=post
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --output=post.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=domenico.zaza@polito.it
#SBATCH --account=IscrB_ALPCF
#SBATCH --partition=g100_usr_smem
#SBATCH --time=03:00:00

ml autoload fftw/3.3.10--openmpi--4.1.1--gcc--10.2.0-pmi
module load openblas/0.3.18--gcc--10.2.0

srun --cpu-bind=cores -m block:block ./POST.x
