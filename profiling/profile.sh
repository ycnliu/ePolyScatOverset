#!/bin/sh

#SBATCH --job-name=profile_epolyscat
#SBATCH --output=test.o%j
#SBATCH --partition=csd_lr6_192
#SBATCH --qos=condo_amos
#SBATCH --account=lr_amos
##SBATCH --nodes=1
##SBATCH --ntasks=16
##SBATCH --perf=vtune
#SBATCH --time=24:00:00

export pe=$HOME/epolyscat
export machcomp=$MACH$COMPILER
export bindir=$pe/bin/$machcomp

mkdir -p $pe/profiling/outdir.$MACH
export od=$pe/profiling/outdir.$MACH

mkdir -p $SCRATCH_DATA/res_dir

# The following command is for running as a batch job and analyze the results in the
# GUI afterwards. At the time of writing this, this batch job has been running for
# two hours and doesn't seem promising and not sure why.  Should continue emailing
# Wei about VTune profiling on Lawrencium.  She is aware of my attemps etc. and
# will be able to help. Only hotspots and threading is available on Lawrencium
# (replace 'hotspots' with 'threading' to profile for parallelism)

srun --nodes=1 --ntasks=16 amplxe-cl -collect hotspots -r $SCRATCH_DATA/res_dir -finalization-mode=none $bindir/ePolyScat.exe $pe/tests/test14.inp >& $od/test14.out

# For running in the GUI, uncomment the line below for profiling on the login
# node inside the GUI.

# $bindir/ePolyScat $pe/tests/test63.inp >& $od/test63.out
