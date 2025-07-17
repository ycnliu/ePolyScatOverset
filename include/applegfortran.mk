   SHELL = /bin/bash

#   FORTMPI := /opt/openmpi-64-gfortran/bin/mpif90
#   FORT := /opt/openmpi-64-gfortran/bin/mpif90
#   FFLAGSMPI := -I. -O2 -m64 -ffpe-summary=invalid,zero,overflow -fcheck=bounds
#   FFLAGSMPIO := -I. -O2 -m64 -ffpe-summary=invalid,zero,overflow -fcheck=bounds
#   FFLAGS := -I. -O2 -m64 -ffpe-summary=invalid,zero,overflow -fcheck=bounds
   FORTMPI := /usr/local/opt/openmpi/bin/mpif90
   FORT := /usr/local/opt/openmpi/bin/mpif90
   FFLAGSMPI := -I. -O2 -m64
   FFLAGSMPIO := -I. -O2 -m64
   FFLAGS := -I. -O2 -m64

   LFLAGS :=  -Wl,-framework -Wl,Accelerate

   LIBINT_TOPDIR := /Users/yuchen/Applications/libint
#   From $(LIBINT_TOPDIR)/libint-2.1.0-beta2/MakeVars
   CXX := clang++
   CXXFLAGS := -O2 -DHAVE_CONFIG_H -std=c++11 -Wno-deprecated-declarations
   CPPLIBFLAGS := -v -lm -lc++

   EIGENDIR := /opt/local/include/eigen3
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