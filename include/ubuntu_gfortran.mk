   SHELL = /bin/bash

   FORTMPI := /usr/bin/mpif90
   FORT := /usr/bin/mpif90
   FFLAGSMPI := -I. -O2 -m64
   FFLAGSMPIO := -I. -O2 -m64
   FFLAGS := -I. -O2 -m64
   LFLAGS := -L/usr/lib -llapack -L/usr/lib -lblas
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
   RUNBATCH := echo "Not programmed"
