   SHELL = /bin/bash

#   export OMPI_MCA_shmem_mmap_enable_nfs_warning=0

   FORTMPI := mpifort
   FORT := ifort
   FFLAGSMPI := -O2 -heap-arrays 1000 -assume realloc-lhs
   FFLAGSMPIO := -O1 -heap-arrays 1000 -assume realloc-lhs
   FFLAGS := -O2 -heap-arrays 1000 -assume realloc-lhs
   LFLAGS := -mkl
   M4FLAGS :=

   CPPLIBFLAGS := -cxxlib
   LIBINT_TOPDIR := /global/home/users/camarante/libint-2.1.0-beta2
#   From $(LIBINT_TOPDIR)/libint-2.1.0-beta2/MakeVars
   CXX := icpc
   CXXFLAGS := -O2 -DHAVE_CONFIG_H -std=c++11 -Wno-deprecated-declarations 

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
