#!/bin/bash
#
# Paths for ePS compilation. Run with 'source' in build terminal.
#
# NOTE: Set for specific builds & versions.
# NOTE: for >E3 may also need to make dirs obj/mach_comp and bin/mach_comp, or get issues with build script.
#
# 13/03/21 version for OpenMPI + gfortran global installation (Ubuntu 20.04LTS, OpenMPI v4.0.3, gfortran v9.3.0)

# May require additional packages as follows (if not already installed):
#  sudo apt-get install openmpi-bin libopenmpi-dev
#  sudo apt-get install libscalapack-openmpi-dev
#  sudo apt-get install libblas-dev
#  sudo apt-get install liblapack-dev
#  sudo apt-get install libatlas-base-dev

export MACH=ubuntu_gfortran
export TMPDIR=/tmp
export pe=/opt/ePolyScat.3885d87


# Make dirs otherwise make throws errors
mkdir $pe/obj/$MACH
mkdir $pe/bin/$MACH
