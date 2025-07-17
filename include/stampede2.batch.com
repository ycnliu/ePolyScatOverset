#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.o%j
#SBATCH --partition=normal      ##lr3
#SBATCH --time=1:00:00
#SBATCH --account=TG-PHY180023
#                                ##SBATCH --qos=lr_normal

date
export pe=${e}
export pt=$pe/tests
export po=$ode
export pd=$pt/$t.d

$bd/ePolyScat $pe/tests/${t}.inp >&$ode/$t.out
#
date
exit


