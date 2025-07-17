#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.o%j
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --partition=csd_lr6_192
#SBATCH --qos=condo_amos
#SBATCH --account=lr_amos


set

date
export pe=${e}
export pt=$pe/tests
export po=$ode
export pd=$pt/$tt.d

$bd/ePolyScat $pt/${tt}.inp >&$ode/$tt.out_lr6_n32
#
date
exit


