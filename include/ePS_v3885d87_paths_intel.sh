#!/bin/bash
#
# Paths for ePS compilation. Run with 'source' in build terminal.
#
# NOTE: Set for specific builds & versions.
# NOTE: for >E3 may also need to make dirs obj/mach_comp and bin/mach_comp, or get issues with build script.
#
# 13/03/21 version testing openmpi-4.1 + Intel compilers v2021.1.2
#	    mk file modified from original `lrcintelmkl.mk`
#     OpenMPI built using Intel compilers, v4.1.0 (see https://www.open-mpi.org/faq/?category=building#easy-build)

source ~/intel/oneapi/setvars.sh
source ~/intel/oneapi/compiler/2021.1.2/env/vars.sh

export PATH=/opt/openmpi-4.1.0/bin/:$PATH
export LD_LIBRARY_PATH=/opt/openmpi-4.1.0/lib/:$LD_LIBRARY_PATH
export MACH=ub-ifort
export TMPDIR=/tmp
export pe=/opt/ePolyScat.3885d87


# Make dirs otherwise make throws errors
mkdir $pe/obj/$MACH
mkdir $pe/bin/$MACH
