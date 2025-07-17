#!/bin/bash -l

#SBATCH -J test
#SBATCH -p regular
#SBATCH -C haswell
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH -o test.o%j
#SBATCH -L SCRATCH
#

date
export pe=${e}
export pt=$pe/tests
export po=$ode
export pd=$pt/$t.d

$bd/ePolyScat $pe/tests/${tt}.inp >&$ode/$tt.out
#
date
exit


