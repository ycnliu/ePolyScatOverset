#!/bin/bash -l
#
#  Master script for parallel execution of ePolyScat
#
#  Raffaele Montuoro, June 2006
#
#SBATCH -p regular
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -L SCRATCH
#SBATCH -J ePolyScat
#SBATCH -e ePolyScat.e%j%a
#SBATCH -o ePolyScat.o%j%a
#

pe="/global/homes/l/loreng/Projects/epolyscat/"
export pe
sexec="srun -n ${SLURM_NTASKS}"
exec="/global/homes/l/loreng/Projects/epolyscat/bin/intelmkl/ePolyScat.exe"
help_msg="Usage: ${0##*/}[-X OPTIONS (direct passed)] [-o OUTFILE] [INPUTFILE]"
opt=
while test $# ; do
case "$1" in
	-o)
	    outfile=$2
	    if [ "${2/#\/}" == "$2" ] ; then
		outfile="`pwd`/$2"
	    else
		outfile="$2"
	    fi
	    shift 
	    shift
	    ;;
	-X)
	    opt="$opt $2"
	    shift 
	    shift 
	    ;;
        -m)
	    memfile=ePolyScat.mem
	    shift
	    ;;
	-*)
	    echo "$help_msg"
	    exit 0
            ;;
	*)
	    break
	    ;;
esac
done

temp=${SCRATCH}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
if [ $memfile ]; then memfile=$temp/$memfile; fi
if [ ${SLURM_ARRAY_TASK_ID} ] ; then
	temp=${SCRATCH}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
fi
mkdir $temp
if [ "$1" != "" ] ; then
# read input from indicated file
    if [ "${1/#\/}" == "$1" ] ; then
	inpfile="`pwd`/$1"
    else
	inpfile="$1"
    fi
    eval "$(echo -e "cat >$temp/inp.dat <<eoi \n$(cat $inpfile)\neoi")"
else
# read from the standard input
    rm -f $temp/inp.dat
    while read NLINE;
      do
      echo "$NLINE" >> $temp/inp.dat ;
    done
fi
if [ $memfile ]; then
   while [ 1 ]; do
        free -g >> $memfile
	sleep 10;
   done &
   mempid=$!
fi 
if [ ! -f $temp/inp.dat ]; then
    echo "No INPUTFILE specified"
    echo "$help_msg"
else
    export MPI_MEMMAP_OFF=1 #??
    export MPI_DSM_DISTRIBUTE=1 #??


    if [ "$outfile" != "" ]; then
        export temp ncpus opt exec outfile
	perl -ne ' /\#\#\*\s+((?:[.a-zA-Z0-9_\$]+\/)*)(\w+\.\w+)/ and system "cp -v $1$2 $ENV{'temp'}/$2" ' $inpfile
        echo Execution on $(hostname) >$outfile
        (cd $temp ; eval "$sexec $exec >> $outfile")
    else
	outfile=$SLURM_SUBMIT_DIR
        export temp ncpus opt exec outfile
	perl -ne ' /\#\#\*\s+((?:[.a-zA-Z0-9_\$]+\/)*)(\w+\.\w+)/ and system "cp -v $1$2 $ENV{'temp'}/$2" ' $inpfile
        echo Execution on $(hostname)
        (cd $temp ; eval "$sexec $exec " )
    fi

    if [ $mempid ]; then kill -s SIGKILL $mempid; fi
 
fi

outdir=$outfile
cp -rv $temp $outdir

exit 0

