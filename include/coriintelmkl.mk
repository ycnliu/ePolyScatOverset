   SHELL = /bin/bash

   FORTMPI := ftn
   FORT := ftn
   FFLAGSMPI := -axMIC-AVX512,AVX
   FFLAGSMPIO := -axMIC-AVX512,AVX
   FFLAGS := -axMIC-AVX512,AVX
   LFLAGS := -mkl

   CPPLIBFLAGS := -cxxlib
   LIBINT_TOPDIR := /global/homes/l/loreng/Projects/libint2lib
   CXX := CC
   CXXFLAGS := -DHAVE_CONFIG_H -std=c++11 -Wno-deprecated-declarations

   EIGENDIR := /global/common/cori/software/eigen3/3.3.3/include/eigen3
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
