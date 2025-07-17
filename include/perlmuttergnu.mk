HELL = /bin/bash

# On August 29, 2023
# using modules
#
# module load PrgEnv-gnu
# module unload craype-accel-nvidia80
#
# with the LAPACK library compiled using lapack-3.11.0
# in make.inc:
# CC = cc
# CFLAGS = -O3
# FC=ftn
# FFLAGS = -O2 -frecursive
# FFLAGS_DRV = $(FFLAGS)
# FFLAGS_NOOPT = -O0 -frecursive
# TIMER = INT_ETIME

  FORTMPI := ftn
  FORT := ftn
  FFLAGSMPI := -O2 -w -fallow-argument-mismatch
  FFLAGSMPIO := -O2 -w -fallow-argument-mismatch
  FFLAGS := -O2 -w -fallow-argument-mismatch
  LIBINT_TOPDIR := /global/homes/l/lucchese/Applications/perlmutter/libint
	CXX := CC
  CXXFLAGS := -O2 -DHAVE_CONFIG_H -std=c++11 -Wno-deprecated-declarations
	CPPLIBFLAGS :=   -lm -lstdc++

  EIGENDIR := $(EIGEN_DIR)/include/eigen3
#  FORTMPI := ftn
#  FORT := ftn
#  FFLAGSMPI := -fallow-argument-mismatch -fbounds-check
#  FFLAGSMPIO := -fallow-argument-mismatch -fbounds-check
#  FFLAGS := -fbounds-check

#  LFLAGS :=
   LFLAGS := -L/global/homes/y/ycliu/Applications/lapack-3.12.0 -llapack
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
