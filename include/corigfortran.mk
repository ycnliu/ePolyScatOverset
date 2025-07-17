   SHELL = /bin/bash

   FORTMPI := ftn
   FORT := ftn
   FFLAGSMPI := 
   FFLAGSMPIO := 
   FFLAGS := 
   LFLAGS :=    -L${MKLROOT}/lib/intel64 -lmkl_core 
   CPPLIBFLAGS := -lstdc++
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
   RUNBATCH := squeue
