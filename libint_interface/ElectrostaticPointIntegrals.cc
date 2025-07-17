/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  Modified June 30, 2015 by Loren Greenman 
 *      in order to interface with grid-based ePolyScat code and 
 *      produce electrostatic point charge integrals on the grid
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>

// Libint Gaussian integrals library
#include <libint2.hpp>

// Interface to ePolyScat
#include <ElectrostaticPointIntegrals.h>


extern "C"{
 void libint2_finterface_initialize_(int* FNAtom, NCoord* FACoord, int* FIAtmFr, int* FIAtmTo, int* FXEXP, int* FYEXP,
                       int* FZEXP, int* FNoPrim, REAL_XR* FCSN, REAL_XR* FEX, int* FFirstInShell){
     // Loren Greenman
     // Initialize the interface by porting the FORTRAN variables to C structs
     // Note that only pointers are imported, the data still exists in the FORTRAN arrays
    eps_atgeom.NAtom = *FNAtom;
    eps_atgeom.ACoord = FACoord;
    eps_gbaorb.IAtmFr = FIAtmFr;
    eps_gbaorb.IAtmTo = FIAtmTo;
    eps_gbaorb.XEXP = FXEXP;
    eps_gbaorb.YEXP = FYEXP;
    eps_gbaorb.ZEXP = FZEXP;
    eps_gbaorb.NoPrim = FNoPrim;
    eps_gbaorb.CSN = FCSN;
    eps_gbaorb.EX = FEX;
    FirstInShell = FFirstInShell;
    // eps_gbaorb.NFcn = FIAtmTo[*FNAtom-1] - FIAtmFr[0] + 1; // Total number of AOs
    eps_gbaorb.NFcn = FIAtmTo[0];
    for (int i=1; i< eps_atgeom.NAtom; i++) eps_gbaorb.NFcn = std::max(eps_gbaorb.NFcn, FIAtmTo[i]);

    eps_gbaorb.funset = new int [eps_gbaorb.NFcn];
    // std::cout<<"DEBUG LIBINT ALLOC RNFACS"<<std::endl;
    eps_gbaorb.rnfacs = new REAL_XR [eps_gbaorb.NFcn];
    // std::cout<<"DEBUG LIBINT DONE ALLOC RNFACS"<<std::endl;
    // std::cout<<"DEBUG LIBINT "<<0<<" "<<dfac[0]<<std::endl;
    dfac[0] = 1;
    dfac[1] = 1;
    for (int i=2; i<20; i++){
        // std::cout<<"DEBUG LIBINT "<<i<<" "<<dfac[i-1]<<std::endl;
        dfac[i] = (2*i-1)*dfac[i-1];
    }
    libint2::init(); // Initialize the library
    AllShells = generate_libint2_shells(); // Generate the library objects
    NucEngine = libint2::OneBodyEngine(libint2::OneBodyEngine::nuclear, max_nprim(AllShells), max_l(AllShells), 0); // Create the long-lived nuclear engine
 }

 void libint2_finterface_check_initialization_(){
     // Loren Greenman
     // Print out the initialized variables, this function serves only debugging purposes
     int NFcn = eps_gbaorb.NFcn;
    std::cout<<"NATOM "<<eps_atgeom.NAtom<<std::endl;
    std::cout<<"NFCN "<<NFcn<<std::endl;
    for (int i=0;i<eps_atgeom.NAtom;i++){
        std::cout<<"INFO FOR ATOM "<<i<<std::endl;
        std::cout<<"ACOORD X "<<eps_atgeom.ACoord[i].vec[0]<<std::endl;
        std::cout<<"ACOORD Y "<<eps_atgeom.ACoord[i].vec[1]<<std::endl;
        std::cout<<"ACOORD Z "<<eps_atgeom.ACoord[i].vec[2]<<std::endl;
        std::cout<<"IATMFRTO "<<eps_gbaorb.IAtmFr[i]<<" "<<eps_gbaorb.IAtmTo[i]<<std::endl;
        for (int j=eps_gbaorb.IAtmFr[i];j<=eps_gbaorb.IAtmTo[i];j++){
            std::cout<<"INFO FOR FUNCTION "<<j<<std::endl;
            std::cout<<"NOPRIM "<<eps_gbaorb.NoPrim[j-1]<<std::endl;
            std::cout<<"XEXP X "<<eps_gbaorb.XEXP[j-1]<<std::endl;
            std::cout<<"YEXP Y "<<eps_gbaorb.YEXP[j-1]<<std::endl;
            std::cout<<"ZEXP Z "<<eps_gbaorb.ZEXP[j-1]<<std::endl;
            for (int k=0;k<eps_gbaorb.NoPrim[j-1];k++){
                std::cout<<"INFO FOR PRIMITIVE "<<j<<" "<<k<<std::endl;
                std::cout<<"CSN "<<eps_gbaorb.CSN[k*NFcn + j-1]<<std::endl;
                std::cout<<"EX "<<eps_gbaorb.EX[k*NFcn + j-1]<<std::endl;
            }
        }
    }
    auto S = compute_1body_ints(AllShells, libint2::OneBodyEngine::overlap);
    std::cout<<"Overlap Integrals \n";
    std::cout<<S<<std::endl;
    std::cout<<"CHECK COMPLETE"<<std::endl;
 }

 void libint2_finterface_check_wfnorm_(REAL_XR* FRHODEN, REAL_XR* FNELEC){
    auto S = compute_1body_ints(AllShells, libint2::OneBodyEngine::overlap);
    int mx = nbasis(AllShells);
    if (mx>10) mx=10;
    std::cout<<"Overlap Integrals \n";
    std::cout<<S.block(0,0,mx,mx)<<std::endl;
    int NFcn = eps_gbaorb.NFcn;
    double wfnorm = 0.0;
    for (int iao = 0; iao < eps_gbaorb.NFcn; ++iao){
        int ilao = eps_gbaorb.funset[iao];
        for (int jao = 0; jao < eps_gbaorb.NFcn; ++jao){
            int jlao = eps_gbaorb.funset[jao];
            /* if (iao == jao){
                printf("DEBUG LIBINT OVERLAP %6i%12.8f%3i(%3i)%3i(%3i)%3i(%3i)%6.2f%12.8f\n",iao,S(ilao,jlao),
                        eps_gbaorb.XEXP[iao], dfac[eps_gbaorb.XEXP[iao]],
                        eps_gbaorb.YEXP[iao], dfac[eps_gbaorb.YEXP[iao]],
                        eps_gbaorb.ZEXP[iao], dfac[eps_gbaorb.ZEXP[iao]],
                        eps_gbaorb.rnfacs[iao],S(ilao,jlao)*sqrt(eps_gbaorb.rnfacs[iao]*eps_gbaorb.rnfacs[jao]));
            } */
            wfnorm += FRHODEN[jao*NFcn+iao] * S(ilao,jlao) * sqrt(eps_gbaorb.rnfacs[iao] * eps_gbaorb.rnfacs[jao]);
        }
    }
    std::cout<<"WFNORM CHECK COMPLETE "<<std::fixed<<std::setprecision(8)<<wfnorm<<std::endl;
    int nelec = (int) *FNELEC;
    *FNELEC = wfnorm;
    // std::cout<<"WFNORM RENORMALIZING TO "<<std::fixed<<std::setprecision(8)<<nelec<<std::endl;
    // if (wfnorm != nelec){
    //     for (int ijao = 0; ijao < NFcn*NFcn; ++ijao){
    //         FRHODEN[ijao] *= nelec/wfnorm;
    //     }
    // }
 }

 void libint2_finterface_check_wfkinetic_(REAL_XR* FRHODEN){
    auto T = compute_1body_ints(AllShells, libint2::OneBodyEngine::kinetic);
    int mx = nbasis(AllShells);
    if (mx>10) mx=10;
    std::cout<<"K. E. Integrals \n";
    std::cout<<T.block(0,0,mx,mx)<<std::endl;
    int NFcn = eps_gbaorb.NFcn;
    double wfkinetic = 0.0;
    for (int iao = 0; iao < eps_gbaorb.NFcn; ++iao){
        int ilao = eps_gbaorb.funset[iao];
        for (int jao = 0; jao < eps_gbaorb.NFcn; ++jao){
            int jlao = eps_gbaorb.funset[jao];
            wfkinetic += FRHODEN[jao*NFcn+iao] * T(ilao,jlao) * sqrt(eps_gbaorb.rnfacs[iao] * eps_gbaorb.rnfacs[jao]);
        }
    }
    std::cout<<"WFNORM KINETIC ENERGY CHECK COMPLETE "<<std::fixed<<std::setprecision(8)<<wfkinetic<<std::endl;
 }


 void libint2_finterface_wfkinetic_(REAL_XR* FRHODEN, REAL_XR* Telement){
    auto T = compute_1body_ints(AllShells, libint2::OneBodyEngine::kinetic);
    int NFcn = eps_gbaorb.NFcn;
    double wfkinetic = 0.0;
    for (int iao = 0; iao < eps_gbaorb.NFcn; ++iao){
        int ilao = eps_gbaorb.funset[iao];
        for (int jao = 0; jao < eps_gbaorb.NFcn; ++jao){
            int jlao = eps_gbaorb.funset[jao];
            wfkinetic += FRHODEN[jao*NFcn+iao] * T(ilao,jlao) * sqrt(eps_gbaorb.rnfacs[iao] * eps_gbaorb.rnfacs[jao]);
        }
    }
    *Telement = wfkinetic;
 }


  void libint2_finterface_fock_(REAL_XR* FRHODEN, REAL_XR* Felement){
    //    auto T = compute_1body_ints(AllShells, libint2::OneBodyEngine::coulomb);
    printf("1111111111111");
    int NFcn = eps_gbaorb.NFcn;
    Matrix D = Matrix::Zero(NFcn,NFcn) ;
    printf("222222222");
    //
    for (int iao = 0; iao < NFcn; ++iao){
      for (int jao = 0; jao < NFcn; ++jao){
	D(iao,jao)= FRHODEN[jao*NFcn+iao];
      }
    }
    printf("3333333333");
    //
        auto F = compute_2body_fock(AllShells,D);
    //    auto F = compute_2body_fock_simple(AllShells,D);
    double wfock = 0.0;
    for (int iao = 0; iao < NFcn; ++iao){
      int ilao = eps_gbaorb.funset[iao];
      for (int jao = 0; jao < NFcn; ++jao){
	int jlao = eps_gbaorb.funset[jao];
	wfock += FRHODEN[jao*NFcn+iao] * F(ilao,jlao) * sqrt(eps_gbaorb.rnfacs[iao] * eps_gbaorb.rnfacs[jao]);
      }
    }
    printf("44444444");
    *Felement = wfock;
  }
  

  void libint2_finterface_nuc_atraction_(REAL_XR* FRHODEN, int* NUC_CHARGE, REAL_XR* NUC_X, REAL_XR* NUC_Y, REAL_XR* NUC_Z, REAL_XR* Nelement){

    int rscharge = *NUC_CHARGE;
    double rsx = *NUC_X;
    double rsy = *NUC_Y;
    double rsz = *NUC_Z;
    
    std::vector<Atom> atomlike;
    { // Initiate an atom with charge -1, this could probably be made more efficient by passing a large number of grid points at once
      Atom ThisAtom = {rscharge,rsx,rsy,rsz};
      atomlike.push_back(ThisAtom);
    }

    // compute nuclear-attraction integrals
    auto N = compute_nuclear_ints(AllShells, NucEngine, atomlike);
    
    // Trace over density matrix to get nuclear-type attraction at the point given
    double epitr = 0.;
    int NFcn = eps_gbaorb.NFcn;
    for (int iao = 0; iao < eps_gbaorb.NFcn; ++iao){
      int ilao = eps_gbaorb.funset[iao]; // LIBINT ilao = FUNSET [iao (polyscat idnex)]
      for (int jao = 0; jao < eps_gbaorb.NFcn; ++jao){
	int jlao = eps_gbaorb.funset[jao];
	epitr += FRHODEN[jao*NFcn+iao] * N(ilao,jlao) * sqrt(eps_gbaorb.rnfacs[iao] * eps_gbaorb.rnfacs[jao]);
      }
    }
    
    *Nelement = epitr;
  }


 void libint2_finterface_wfoverlap_(REAL_XR* FRHODEN, REAL_XR* Selement){
    auto S = compute_1body_ints(AllShells, libint2::OneBodyEngine::overlap);
    int NFcn = eps_gbaorb.NFcn;
    double wfoverlap = 0.0;
    for (int iao = 0; iao < eps_gbaorb.NFcn; ++iao){
        int ilao = eps_gbaorb.funset[iao];
        for (int jao = 0; jao < eps_gbaorb.NFcn; ++jao){
            int jlao = eps_gbaorb.funset[jao];
            wfoverlap += FRHODEN[jao*NFcn+iao] * S(ilao,jlao) * sqrt(eps_gbaorb.rnfacs[iao] * eps_gbaorb.rnfacs[jao]);
        }
    }
    *Selement = wfoverlap;
 }
}


std::vector<libint2::Shell> generate_libint2_shells() {
    // Loren Greenman
    // This function generates LIBINT2 data structures for the basis elements
    // A libint2::Shell object is created for each ePolyScat "function" (basis function)
    // A std::vector of such shells is created and output
     int NFcn = eps_gbaorb.NFcn;
    int NAtom = eps_atgeom.NAtom;
    std::vector<libint2::Shell> AllShells;
    // Shell - exponents, contraction (l, bool, coefs), origin
    for (int iatom = 0; iatom < NAtom; ++iatom){
        for (int kfcn = eps_gbaorb.IAtmFr[iatom] - 1; kfcn < eps_gbaorb.IAtmTo[iatom]; ++kfcn){
            int l = eps_gbaorb.XEXP[kfcn] + eps_gbaorb.YEXP[kfcn] + eps_gbaorb.ZEXP[kfcn];
            // Determine the correspondence between ePolyScat and LIBINT order
            int LIBINTord = nbasis(AllShells);
            if (FirstInShell[kfcn] != 1) LIBINTord-=((l+1)*(l+2))/2; // reset count so that it starts at the first function of this shell
            /* std::cout<<"DEBUG LIBINT Add function "<<kfcn<<" "<<l<<
                " "<<eps_gbaorb.XEXP[kfcn]<<
                " "<<eps_gbaorb.YEXP[kfcn]<<
                " "<<eps_gbaorb.ZEXP[kfcn]<<std::endl; */
            // std::cout<<"DEBUG LIBINT AllShells.size() "<<LIBINTord<<"\n";
            for (int ilxord = 0; ilxord <= l; ilxord++){
                int nx = l - ilxord; // LIBINT ordering: x exponent
                for (int jlxord = 0; jlxord <= ilxord; jlxord++){
                    int ny = ilxord - jlxord;
                    int nz = l - nx - ny; // y and z exponents
                    // std::cout<<"DEBUG LIBINT order "<<LIBINTord<<" "<<nx<<" "<<ny<<" "<<nz<<" ";
                    if ((eps_gbaorb.XEXP[kfcn] == nx)&&
                            (eps_gbaorb.YEXP[kfcn] == ny)&&
                            (eps_gbaorb.ZEXP[kfcn] == nz)){
                        // std::cout<<"YES"<<std::endl;
                        eps_gbaorb.funset[kfcn] = LIBINTord;
                    }
                    // else{
                    //     std::cout<<"NO"<<std::endl;
                    // }
                    ++LIBINTord;
                }
            }
            // Determine the renormalization constant for unity-normalized shell members
            eps_gbaorb.rnfacs[kfcn] = 1.0;
            if (1 != FirstInShell[kfcn]){
                int dfx = dfac[eps_gbaorb.XEXP[kfcn]];
                int dfy = dfac[eps_gbaorb.YEXP[kfcn]];
                int dfz = dfac[eps_gbaorb.ZEXP[kfcn]];
                eps_gbaorb.rnfacs[kfcn] = dfac[l]/(dfx*dfy*dfz);
            }
            if (1 == FirstInShell[kfcn]){
                std::vector<REAL_XR> ContrCoefs;
                std::vector<REAL_XR> Exponents;
                for (int iprim = 0; iprim < eps_gbaorb.NoPrim[kfcn]; ++iprim){
                    ContrCoefs.push_back(eps_gbaorb.CSN[iprim*NFcn+kfcn]);
                    Exponents.push_back(eps_gbaorb.EX[iprim*NFcn+kfcn]);
                    /* int dfac = 1.0;
                    for (int nx = eps_gbaorb.XEXP[kfcn]; nx >= 1; nx++){
                        dfac *= 2 * nx - 1;
                    }
                    double nfac = pow((2./3.14159265359),0.75)*pow(2.,l)*pow(eps_gbaorb.EX[iprim*NFcn+kfcn],0.25*(2.0*l+3.0))/dfac;
                    std::cout<<"DEBUG LIBINT Norm "<<kfcn<<" "<<iprim<<" "<<nfac<<" "<<eps_gbaorb.CSN[iprim*NFcn+kfcn]<<std::endl; */
                }
                libint2::Shell::Contraction ThisContraction = {l, false, ContrCoefs};
                std::vector<libint2::Shell::Contraction> ThisContrV;
                ThisContrV.push_back(ThisContraction);
                std::array<REAL_XR,3> Origin = {eps_atgeom.ACoord[iatom].vec[0],
                                                eps_atgeom.ACoord[iatom].vec[1],
                                                eps_atgeom.ACoord[iatom].vec[2]};
                // This constructor does some weird renormalization, check if there's a problem later
                // libint2::Shell ThisShell(Exponents, ThisContrV, Origin);
                // By building the shells this way, the renormalization is avoided
                libint2::Shell ThisShell;
                ThisShell.alpha = Exponents;
                ThisShell.contr = ThisContrV;
                ThisShell.O = Origin;
                AllShells.push_back(ThisShell);
                // std::cout<<"GENERATED SHELL FOR FCN "<<kfcn<<" "<<nbasis(AllShells)<<std::endl;
            }
        }
    }
    return AllShells;
}

extern "C"{
 void libint2_finterface_print_shell_normalizations_() {
     // Loren Greenman
     // A debugging function, print out sum^2 for the contraction coefficients
    // libint2::Shell::do_enforce_unit_normalization(false);
   // std::cout<<"PSNORM: THE VALUE OF THE NORMALIZATION BOOL IS "<<libint2::Shell::do_enforce_unit_normalization()<<std::endl;
   int NAtom = eps_atgeom.NAtom;
   std::vector<libint2::Shell> AllShells = generate_libint2_shells();
   std::vector<libint2::Shell>::iterator ishell;
   for ( ishell = AllShells.begin(); ishell != AllShells.end(); ishell++ ){
       libint2::Shell ThisShell = *ishell;
       std::vector<libint2::Shell::Contraction> ThisContrs = ThisShell.contr;
       std::vector<libint2::Shell::Contraction>::iterator icontr;
       double nrmshell = 0.0;
       for ( icontr = ThisContrs.begin(); icontr != ThisContrs.end(); icontr++ ){
           double nrmcontr = 0.0;
           libint2::Shell::Contraction ThisContr = *icontr;
           std::vector<REAL_XR> ContrCoefs = ThisContr.coeff;
           for ( std::vector<REAL_XR>::iterator icoef = ContrCoefs.begin(); icoef != ContrCoefs.end(); icoef++){
               nrmcontr += *icoef * *icoef;
               std::cout<<"\t\t\t\t"<<*icoef<<" "<<(*icoef)*(*icoef)<<std::endl;
           }
           std::cout<<"\t\tThe Norm of the contraction is "<<nrmcontr<<std::endl;
           nrmshell += nrmcontr;
       }
       std::cout<<"\t\tThe Norm of the shell is "<<nrmshell<<std::endl;
   }
 }

// Key function to call from FORTRAN to do nuclear-type integrals
void libint2_finterface_do_nucl_attraction_(REAL_XR* FRSX, REAL_XR* FRSY, REAL_XR* FRSZ, REAL_XR* FRHODEN, REAL_XR* FINTOUT){
    // Set the x,y,z point from the FORTRAN input
    double rsx = *FRSX;
    double rsy = *FRSY;
    double rsz = *FRSZ;

    /* -- Uncomment if you want to see how many integrals you're doing
    static int ntimes = 0;

    if (ntimes % 1 == 0){
        std::cout<<"DOING A NUCLEAR-TYPE INTEGRAL FOR THE POINT ("<<rsx<<","<<rsy<<","<<rsz<<") : R = "<<rsx*rsx+rsy*rsy+rsz*rsz<<std::endl;
    }
    */

    std::vector<Atom> atomlike;
    { // Initiate an atom with charge -1, this could probably be made more efficient by passing a large number of grid points at once
        Atom ThisAtom = {-1, rsx,rsy,rsz};
        atomlike.push_back(ThisAtom);
    }

    // compute nuclear-attraction integrals
    auto EPI = compute_nuclear_ints(AllShells, NucEngine, atomlike);

    // Trace over density matrix to get nuclear-type attraction at the point given
    double epitr = 0.;
    int NFcn = eps_gbaorb.NFcn;
    for (int iao = 0; iao < eps_gbaorb.NFcn; ++iao){
        int ilao = eps_gbaorb.funset[iao]; // LIBINT ilao = FUNSET [iao (polyscat idnex)]
        for (int jao = 0; jao < eps_gbaorb.NFcn; ++jao){
            int jlao = eps_gbaorb.funset[jao];
            // printf("%6i%6i%6i%6i%20.10e%20.10e\n",iao,jao,ilao,jlao,FRHODEN[jao*NFcn+iao],EPI(ilao,jlao));
            epitr += FRHODEN[jao*NFcn+iao] * EPI(ilao,jlao) * sqrt(eps_gbaorb.rnfacs[iao] * eps_gbaorb.rnfacs[jao]);
        }
    }

    /*
    if (ntimes % 1 == 0){
        std::cout<<ntimes<<" "<<epitr<<std::endl;
    }
    ++ntimes;
    */

    // Set the return variable
    *FINTOUT = epitr;
 } 

 void libint2_finterface_finalize_(){
     delete [] eps_gbaorb.funset;
     delete [] eps_gbaorb.rnfacs;
     libint2::finalize();
 }
}

Matrix compute_nuclear_ints(const std::vector<libint2::Shell>& shells,
                          libint2::OneBodyEngine &nucengine,
                          const std::vector<Atom>& atoms)
{
  const auto n = nbasis(shells);
  Matrix result(n,n);

  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (true){ // (obtype == libint2::OneBodyEngine::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    nucengine.set_params(q);
  }

  auto shell2bf = map_shell_to_basis_function(shells);

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
  for(auto s1=0; s1!=shells.size(); ++s1) {
    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      const auto* buf = nucengine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
      result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

    }
  }

  return result;

} 
