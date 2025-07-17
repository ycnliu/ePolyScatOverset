   SHELL = /bin/bash

   FORTMPI := mpif90
   FORT := ifort
   FFLAGSMPI := -O2 -heap-arrays 1000
   FFLAGSMPIO := -O1 -heap-arrays 1000
   FFLAGS := -O2 -heap-arrays 1000
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
