cmake_minimum_required(VERSION 3.12)
project(EPS)

# Compiler settings
set(CMAKE_CXX_COMPILER "/usr/local/opt/openmpi/bin/mpicxx")
set(CMAKE_Fortran_COMPILER "/usr/local/opt/openmpi/bin/mpif90")
set(CMAKE_CXX_STANDARD 11)

# Compiler flags
set(FFLAGS "-g -O2 -mtune=native -fbounds-check")
set(FFLAGSMPI "-g -O2 -mtune=native -fbounds-check")
set(CMAKE_CXX_FLAGS "-O2 -DHAVE_CONFIG_H -std=c++11 -Wno-deprecated-declarations")

# Directories
set(LIBINT_TOPDIR "/Users/yuchen/Applications/libint")
set(LIBINT_INTERFACE_DIR "/Users/yuchen/AMO/epolyscat/libint_interface")
set(EIGENDIR "/usr/local/Cellar/eigen/3.4.0_1/include/eigen3")

include_directories(${LIBINT_TOPDIR}/include ${EIGENDIR})

# Source files
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Fortran source files
# Fortran source files
set(FORTRAN_SOURCES
    ${SRC_DIR}/Orient.f90
    ${SRC_DIR}/VcpPol.f90
    ${SRC_DIR}/ExpOrb.f90
    ${SRC_DIR}/GaussianCnv.f90
    ${SRC_DIR}/GamessCnv.f90
    ${SRC_DIR}/MesaCnv.f90
    ${SRC_DIR}/SymGen.f90
    ${SRC_DIR}/GenGrid.f90
    ${SRC_DIR}/AngGCt.f90
    ${SRC_DIR}/RotOrb.f90
    ${SRC_DIR}/Den.f90
    ${SRC_DIR}/StPot.f90
    ${SRC_DIR}/AsyPol.f90
    ${SRC_DIR}/Fege.f90
    ${SRC_DIR}/GetRMax.f90
    ${SRC_DIR}/ScatStabA.f90
    ${SRC_DIR}/SymNormMode.f90
    ${SRC_DIR}/SymProd.f90
    ${SRC_DIR}/MatEle.f90
    ${SRC_DIR}/DipoleOp.f90
    ${SRC_DIR}/GetDataRecordDef.f90
    ${SRC_DIR}/TotalCrossSection.f90
    ${SRC_DIR}/MoldenCnv.f90
    ${SRC_DIR}/MFDCS.f90
    ${SRC_DIR}/EDCS.f90
    ${SRC_DIR}/VcpPN.f90
    ${SRC_DIR}/VcpBN.f90
    ${SRC_DIR}/SaveData.f90
    ${SRC_DIR}/DPot.f90
    ${SRC_DIR}/Resonance.f90
    ${SRC_DIR}/ViewOrb.f90
    ${SRC_DIR}/Vposfit.f90
    ${SRC_DIR}/Rotate.f90
    ${SRC_DIR}/OrientN.f90
    ${SRC_DIR}/CalcInt.f90
    ${SRC_DIR}/DumpMesa.f90
    ${SRC_DIR}/RotOrientAsym.f90
    ${SRC_DIR}/SchmidtOrth.f90
    ${SRC_DIR}/MFTimeDelay.f90
    ${SRC_DIR}/LFTimeDelay.f90
    # Add other Fortran source files...
)

# C++ source files
set(CPP_SOURCES
    ${LIBINT_INTERFACE_DIR}/IntegrationTools.cc
    ${LIBINT_INTERFACE_DIR}/ElectrostaticPointIntegrals.cc
    # Add other C++ source files...
)

# Additional source files
set(ADDITIONAL_SOURCES
    ${SRC_DIR}/LMPI.f90
    ${SRC_DIR}/SetUp.for
    ${SRC_DIR}/Vectors.m4
    ${SRC_DIR}/Modules.f90
    ${SRC_DIR}/LPP.f90
    # Add other additional source files...
)

# Executable
add_executable(ePolyScat ${FORTRAN_SOURCES} ${CPP_SOURCES} ${ADDITIONAL_SOURCES})

# Linking
target_link_libraries(ePolyScat ${LIBINT_TOPDIR}/lib/libint2.a)

# Additional configuration for tests, manuals, etc. can be added as needed.
