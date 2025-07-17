   SHELL = /bin/bash

   FORTMPI := tau mpif90
   FORT := tau mpif90

   FFLAGSMPI := -O2 -heap-arrays 1000 -assume realloc-lhs
   FFLAGSMPIO := -O1 -heap-arrays 1000 -assume realloc-lhs
   FFLAGS := -O2 -heap-arrays 1000 -assume realloc-lhs

#    FFLAGSMPI := -O2 -heap-arrays 1000 -check bounds -assume realloc-lhs
#    FFLAGSMPIO := -O1 -heap-arrays 1000 -check bounds -assume realloc-lhs
#    FFLAGS := -O2 -heap-arrays 1000 -check bounds -assume realloc-lhs

   LFLAGS := -mkl
   M4FLAGS :=

   CPPLIBFLAGS := -cxxlib
   LIBINT_TOPDIR := /global/home/users/rlucchese/Applications/libint
#   From $(LIBINT_TOPDIR)/libint-2.1.0-beta2/MakeVars
   CXX := tau mpic++
   CXXFLAGS := -O2 -DHAVE_CONFIG_H -std=c++11

   EIGENDIR := $(EIGEN_DIR)/include/eigen3

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
   DEFT := batch
