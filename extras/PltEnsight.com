#
cd $TMPDIR
mkdir PltEnsight.d
cd PltEnsight.d
workdir=$PWD
pe=~/Applications/ePolyScat.E3
$pe/bin/ePolyScat <<eoi
#
# input file for test09
#
# Expand HOMO and LUMO of SF6
#
  LMax   15     # maximum l to be used for wave functions
  LMaxI  40     # maximum l value used to determine numerical angular grids
  EMax  50.0    # EMax, maximum asymptotic energy in eV
CnvOrbSel  33 36
Convert '$pe/tests/test09.g03' 'g03'
GetBlms
ExpOrb
FileName 'ViewOrb' '$workdir/ViewOrb.dat' 'REWIND'
FileName 'ViewOrbGeom' '$workdir/ViewOrbGeom.dat' 'REWIND'
ViewOrbGrid
  0.0 0.0 0.0
  0.0 0.0 1.0
  1.0 0.0 0.0
  -2.5 2.5 0.1
  -2.5 2.5 0.1
  -2.5 2.5 0.1
ViewOrb 'ExpOrb' 1 3
ViewOrb 'ExpOrb' 2
eoi

$pe/bin/ViewOrb2ensight.exe <<eoi
'$workdir/ViewOrb.dat' '$workdir/ViewOrbGeom.dat' ' '
-1 -
1.75
eoi



#
# Start EnSight
# Preferences, Data, Default Data Directory change to $TMPDIR/PltEnsight.d, click Save To Preferences File button
# Restart EnSight if you needed to change the default directory
# File -> Command
# Load: Browse -> Open the Vieworb.cmd file
# Click Start Playing COmmands button
# Close Command window
# Change time to see different orbitals
#
