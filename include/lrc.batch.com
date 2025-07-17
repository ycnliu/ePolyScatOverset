#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.o%j
#SBATCH --partition=lr3
#SBATCH --time=4:00:00
#SBATCH --qos=condo_axl
#SBATCH --account=lr_axl

date
export pe=${e}
export pt=$pe/tests
export po=$ode
export pd=$pt/$t.d

$bd/ePolyScat $pe/tests/${tt}.inp >&$ode/$tt.out
#
date
exit


