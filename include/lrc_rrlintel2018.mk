#
# for use with the ifort version 18.0.5
#
# module load intel/2018.5.274.par
# module load impi/2018.5.274
# export COMPILER=intel2018
#
   SHELL = /bin/bash

   FORTMPI := mpiifort
   FORT := ifort
   FFLAGSMPI := -O2 -heap-arrays 1000 -assume realloc-lhs
   FFLAGSMPIO := -O1 -heap-arrays 1000 -assume realloc-lhs
   FFLAGS := -O2 -heap-arrays 1000 -assume realloc-lhs
   LFLAGS := -mkl
   M4FLAGS :=

   DPRE := OFF
   CRAY := OFF
   DOSF := OFF
   IRIXP := ON
   UAIX := OFF
   UAIXMC := 100
   GEN := OFF
   LAPACKR := D
   LAPACKC := Z
   MPI := OFF
   MPIF90 := ON
   RUNBATCH := sbatch
