#
# Script to seet up bash functions to do comparisons
#
compstring="MaxIter =|alcoef|^RotOrAs|^EPhi|SelcRoots|CalcInt value|CalcIntL value|Normalization Error|^Test2Re|^Tst5I2R|^TsR5I2R"
compbin=${PWD/\/extras/}/bin/$MACH
if [ -n "$COMPILER" ]; then
   compbin=$compbin$COMPILER
fi
echo CompDiff \<dir .out 1\> \<dir .out 2\> \<test name\>
CompDiff () {
if [ $# -lt 3 ] ; then
   echo CompDiff \<dir .out 1\> \<dir .out 2\> \<test name\>
   echo CompDiffStnd \<dir .out\> \<dir .ostnd\> \<test name\>
   echo CompDiffSD \<file extension 1\> \<file extension 2\> \<test name\> \# assuming current directory
   echo CompDiffAll \<dir .out 1\> \<dir .out 2\> \<dir with .ostnd files\>
   echo CompDiffSDAll \<file extension 1\> \<file extension 2\> \# all test[0-9]\* in current directory
   echo CompDiffStndAll \<dir .out\> \<dir .ostnd\> \# all test[0-9]\*
else
   egrep "$compstring" $1/$3.out >grep1$$
   egrep "$compstring" $2/$3.out >grep2$$
   $compbin/CompDiff.exe grep1$$ grep2$$ $3
   rm grep1$$ grep2$$
fi
}
echo CompDiffStnd \<dir .out\> \<dir .ostnd\> \<test name\>
CompDiffStnd () {
egrep "$compstring" $1/$3.out >grep1$$
egrep "$compstring" $2/$3.ostnd >grep2$$
$compbin/CompDiff.exe grep1$$ grep2$$ $3
rm grep1$$ grep2$$
}
echo CompDiffSD \<file extension 1\> \<file extension 2\> \<test name\> \# assuming current directory
CompDiffSD () {
egrep "$compstring" $3.$1 >grep1$$
egrep "$compstring" $3.$2 >grep2$$
$compbin/CompDiff.exe grep1$$ grep2$$ $3
rm grep1$$ grep2$$
}
echo CompDiffAll \<dir .out 1\> \<dir .out 2\> \<dir with .ostnd files\>
CompDiffAll () {
# argument 1 is the first directory containing .out files
# argument 2 is the second directory containing .out files
# argument 3 is the directory containing the .ostnd files
for TEST in $3/test[0-9]?.ostnd ; do
TESTN=${TEST/.ostnd/}
TESTN2=${TESTN##*/}
CompDiff $1 $2 $TESTN2
done
}
echo CompDiffSDAll \<file extension 1\> \<file extension 2\> \# all test[0-9]\* in current directory
CompDiffSDAll () {
for NAME in test[0-9]?.$1
do 
   CompDiffSD $1 $2 ${NAME/.$1/}
done
}
echo CompDiffStndAll \<dir .out\> \<dir .ostnd\> \# all test[0-9]\*
CompDiffStndAll () {
# argument 1 is the directory containing the .out files
# argument 2 is the directory containing the .ostnd files
for TEST in $2/test[0-9]?.ostnd ; do
TESTN=${TEST/.ostnd/}
TESTN2=${TESTN##*/}
CompDiffStnd $1 $2 $TESTN2
done
}
#
