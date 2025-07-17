/*
 * This header file sets up the translation between libint and ePolyScat
 *  Loren Greenman, July 1, 2015
 *
 */

#ifndef _libint2_ePolyscat_interface_cpp_header_
#define _libint2_ePolyscat_interface_cpp_header_

typedef double REAL_XR; // Note that the XR FORTRAN KIND is just set to double-precision here
// This is a possible source of error, but unlikely on most machines

#include<libint2.hpp>
#include<IntegrationTools.h>


// The following objects mimic FORTRAN objects from ePolyScat
// They are initialized from FORTRAN using the libint2_finterface_initialize function
class NCoord{
  public:
    REAL_XR vec[3];
    int z, zs;
};

struct{
    // Key variables for interface are 
    //  EX: exponents
    //  CSN: coefficients
    //  NOPRIM: # primitives per function
    //  IATMFR/IATMTO: function indices
    REAL_XR **CMatGrp, *EX, *CSN, **CMAT;
    REAL_XR* OrbEng;
    int *OrbOrig, *NoPrim, *IAtmFr, *IAtmTo;
    int NBASIS, NMO, ndggrp;
    int *XEXP, *YEXP, *ZEXP, *neledg;
    int **modgen, **nmondg;
    char* SYMNAM;
    int* LMaxAtomBasis;
    int LMaxBasisOrigin;
    int NFcn;
    int* funset; // ordering functions properly to match LIBINT cartesian shell ordering
    REAL_XR* rnfacs; // renormalization factors for the shell members that aren't copied over
} eps_gbaorb;

struct{
    NCoord* ACoord;
    int NAtom;
    char* gtitle;
    REAL_XR* AMass;
    int NNormMode;
    REAL_XR *FreqNormMode, *RedMassNormMode, *ForceKNormMode;
    REAL_XR** VecNormMod;
    REAL_XR* XCT;
} eps_atgeom;

extern "C"{
    // Loren Greenman
    // These routines form the crux of the interface between LIBINT2 and ePolyScat
    // They are callable from FORTRAN
    //
    // Initialize the relevant variables
    void libint2_finterface_initialize_(int* FNAtom, NCoord* ACoord, int* FIAtmFr, int* FIAtmTo, int* FXEXP,
					int* FYEXP, int* FZEXP, int* FNoPrim, REAL_XR* FCSN, REAL_XR* FEX,
					int* FFirstInShell);
    // Debugging functions
    void libint2_finterface_check_initialization_();
    void libint2_finterface_print_shell_normalizations_();
    void libint2_finterface_check_wfnorm_(REAL_XR* FRHODEN, REAL_XR* FNELEC);
    void libint2_finterface_check_wfkinetic_(REAL_XR* FRHODEN);
    void libint2_finterface_do_na_debugasym_(REAL_XR* FRSX, REAL_XR* FRSY, REAL_XR* FRSZ, REAL_XR* FRHODEN, REAL_XR* FINTOUT);
    // Perform the nuclear attraction integral (charge = -1) at point (FRS)X, Y, Z, w/ density matrix (F)RHODEN, put the result in FINTOUT
    void libint2_finterface_do_nucl_attraction_(REAL_XR* FRSX, REAL_XR* FRSY, REAL_XR* FRSZ, REAL_XR* FRHODEN, REAL_XR* FINTOUT);
    // Finalize the interface
    void libint2_finterface_finalize_();
}

// These functions are internal to the interface: 
//  This function generates the Shell objects from the basis set given in FORTRAN
//  The Shell objects are the key to doing the integrals with LIBINT
std::vector<libint2::Shell> generate_libint2_shells();
//  A vector of Shell objects, basically the basis set is held here
std::vector<libint2::Shell> AllShells;
//  An engine for doing nuclear integrals, declared globally so that it is long-lived and doesn't need to be recreated for every integral
libint2::OneBodyEngine NucEngine;
//  The main routine for computing the integrals, using the longlived objects we create once
Matrix compute_nuclear_ints(const std::vector<libint2::Shell>& shells,
                          libint2::OneBodyEngine &nucengine,
                          const std::vector<Atom>& atoms);

int dfac[20];
int* FirstInShell;

#endif /* header guard */

