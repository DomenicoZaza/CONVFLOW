module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0-cuda-12.2 
module load fftw/3.3.10--openmpi--4.1.6--gcc--12.2.0-spack0.22
module load openblas/0.3.26--gcc--12.2.0

make -f Makefile clean
make -f Makefile