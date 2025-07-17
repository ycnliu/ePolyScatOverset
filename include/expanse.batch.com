#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.o%j
#SBATCH --partition=debug
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --account=nsf106
#                                ##SBATCH --qos=lr_normal

module load intel
module load intel-mkl
module load intel-mpi
module load slurm
module list

date
export pe=${e}
export pt=$pe/tests
export po=$ode
export pd=$pt/$t.d
cd $TMPDIR
pwd

$bd/ePolyScat $pe/tests/${t}.inp >&$ode/$t.out
#
date
exit


