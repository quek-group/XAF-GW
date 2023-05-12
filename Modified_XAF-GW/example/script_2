#!/bin/bash
#BSUB -q scalable
#BSUB -n 24
#BSUB -R "span[ptile=24]"
#BSUB -R "rusage[mem=250000]"
#BSUB -a openmpi
#BSUB -W 48:00
#BSUB -J "gw"
#BSUB -o lsf%J.o
#BSUB -e lsf%J.e
module load BerkeleyGW/3.0.1-intel-2019b-Python-3.8.6
export LOCAL_SCRATCHDIR=/scratch/tmp_$LSB_JOBID
export WDIR=$(pwd)
export OMP_NUM_THREADS=1
cd $WDIR/3a-epsilon
mpirun -np $LSB_DJOB_NUMPROC epsilon.cplx.x > epsilon.out

