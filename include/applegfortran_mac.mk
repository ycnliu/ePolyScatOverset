   SHELL = /bin/bash

#   FORTMPI := /opt/openmpi-64-gfortran/bin/mpif90
   FORTMPI := /opt/local/libexec/openmpi-mp/mpif90
#   FORT := /opt/openmpi-64-gfortran/bin/mpif90
   FORT := /opt/local/libexec/openmpi-mp/mpif90
   FFLAGSMPI := -I. -O2 -m64 -ffpe-summary=invalid,zero,overflow -fcheck=bounds
   FFLAGSMPIO := -I. -O2 -m64 -ffpe-summary=invalid,zero,overflow -fcheck=bounds
   FFLAGS := -I. -O2 -m64 -ffpe-summary=invalid,zero,overflow -fcheck=bounds
   LFLAGS :=  -Wl,-framework -Wl,Accelerate
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

