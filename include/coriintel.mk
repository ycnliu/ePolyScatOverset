   SHELL = /bin/bash

   FORTMPI := ftn
   FORT := ftn
   FFLAGSMPI := -fast -no-ipo -assume realloc_lhs
   FFLAGSMPIO := -fast -no-ipo -assume realloc_lhs
   FFLAGS := -fast -no-ipo -assume realloc_lhs
#    FFLAGSMPI := -O0 -g -check bounds
#    FFLAGSMPIO := -O0 -g -check bounds
#    FFLAGS := -O0 -g -check bounds
#    LFLAGS := -mkl
#    LFLAGS := -mkl -L/global/homes/l/loreng/Projects/gmp/lib
   LFLAGS := 
   M4FLAGS :=

   CPPLIBFLAGS := -cxxlib
   LIBINT_TOPDIR := /global/homes/l/lucchese/Applications/cori/libint
#   From $(LIBINT_TOPDIR)/libint-2.1.0-beta2/MakeVars
   CXX := icpc
   CXXFLAGS := -O2 -DHAVE_CONFIG_H -std=c++11 -Wno-deprecated-declarations 

   EIGENDIR := $(EIGEN3_DIR)/include/eigen3

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
