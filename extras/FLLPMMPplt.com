#
rm -f FLLPMMP.dat
ln -s $1 FLLPMMP.dat
matlab -nosplash -nodesktop <FLLPMMPplt.m
rm FLLPMMP.dat
#
