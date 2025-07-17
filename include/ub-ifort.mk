   SHELL = /bin/bash

   FORTMPI := mpifort
   FORT := ifort
   FFLAGSMPI := -O2 -heap-arrays 1000 -assume realloc-lhs
   FFLAGSMPIO := -O1 -heap-arrays 1000 -assume realloc-lhs
   FFLAGS := -O2 -heap-arrays 1000 -assume realloc-lhs
   LFLAGS := -llapack -lblas
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
