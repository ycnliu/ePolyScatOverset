#
cd include
echo Machine $MACH batch scripts are
for file in ${MACH}.batch*.com ; do
tmp1=${file/${MACH}./}
echo ${tmp1/.com/}
done

