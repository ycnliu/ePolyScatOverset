#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.o%j
#SBATCH --partition=lr3
#SBATCH --time=24:00:00
#SBATCH --qos=condo_axl
#SBATCH --account=lr_axl
#SBATCH --exclude=n0237.lr3
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8

set

date
export pe=${e}
export pt=$pe/tests
export po=$ode
export pd=$pt/$tt.d

$bd/ePolyScat $pt/${tt}.inp >&$ode/$tt.out2
#
date
exit


