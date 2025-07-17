   SHELL = /bin/bash

   FORTMPI := mpiifort
   FORT := ifort
   FFLAGSMPI := -O2 -heap-arrays 1000 -assume realloc-lhs -mkl=sequential
   FFLAGSMPIO := -O1 -heap-arrays 1000 -assume realloc-lhs -mkl=sequential
   FFLAGS := -O2 -heap-arrays 1000 -assume realloc-lhs -mkl=sequential
#   LFLAGS :=  -Wl,--start-group ${INTEL_MKLHOME}/mkl/lib/intel64/libmkl_intel_lp64.a ${INTEL_MKLHOME}/mkl/lib/intel64/libmkl_sequential.a ${INTEL_MKLHOME}/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#    LFLAGS := -L${INTEL_MKLHOME}/mkl/lib/intel64 -lpthread -lm -ldl
#   LFLAGS := -mkl -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#   LFLAGS := -L${INTEL_MKLHOME}/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lpthread -lm -ldl
#    LFLAGS :=  -lmkl_intel_lp64 -lmkl_core -lpthread -lm -ldl
   LFLAGS := -L${INTEL_MKLHOME}/mkl/lib/intel64  -lpthread -lm -ldl
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
