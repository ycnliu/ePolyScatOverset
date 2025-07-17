#
GITSHA1=$(git log -1 --format=%H)
GITDATE=$(git log -1 --format=%cD)
appdir=~/Applications
cd $appdir
rm -r -f ePolyScatDistVer
git clone https://rlucchese@bitbucket.org/rlucchese/epolyscat.git ePolyScatDistVer
cd ePolyScatDistVer
rm -r -f .git .gitignore
mv Makefile OldMakefile
sed -e "s/\$(shell git log -1 --format=%H)/$GITSHA1/" <OldMakefile | sed -e "s/\$(shell git log -1 --format=%cD)/$GITDATE/" >Makefile
rm OldMakefile
cd $appdir
tar cfz ePolyScatDistVer${GITSHA1:0:7}.tgz ePolyScatDistVer
scp ePolyScatDistVer${GITSHA1:0:7}.tgz $xlrc:Applications
#
