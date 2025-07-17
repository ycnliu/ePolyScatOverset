#include "StPot_overset_sub.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace OversetGrid {

// Example global/config data (replace with your own in production!)
extern int NATOM;  // Number of atoms
extern std::vector<NCoord> ACOORD;
extern std::vector<int> IATMFR, IATMTO, NOPRIM;
extern std::vector<std::vector<double>> CSN, EX;
extern std::vector<int> XEXP, YEXP, ZEXP;

void CalcAO(const std::vector<Cvec>& XYZ, std::vector<std::vector<double>>& valao)
{
    size_t npe = XYZ.size();
    if (npe == 0) return;

    // KMAX = total number of basis functions (from your config)
    const int KMAX = /* determine from your basis set */ 100;
    valao.assign(npe, std::vector<double>(KMAX, 0.0));

    std::vector<std::vector<double>> scrv(npe, std::vector<double>(4, 0.0));
    std::vector<double> valbfn(npe, 0.0);

    // Loop over atoms
    for (int IN = 0; IN < NATOM; ++IN) {
        // Calculate local frame coordinates
        for (size_t j = 0; j < npe; ++j) {
            scrv[j][0] = XYZ[j].x - ACOORD[IN].vec[0];
            scrv[j][1] = XYZ[j].y - ACOORD[IN].vec[1];
            scrv[j][2] = XYZ[j].z - ACOORD[IN].vec[2];
            scrv[j][3] = scrv[j][0]*scrv[j][0] + scrv[j][1]*scrv[j][1] + scrv[j][2]*scrv[j][2];
        }

        // For each contracted basis function on atom IN
        for (int K = IATMFR[IN]; K <= IATMTO[IN]; ++K) {
            std::fill(valbfn.begin(), valbfn.end(), 0.0);

            // Sum over primitives for contraction
            for (int I = 0; I < NOPRIM[K]; ++I) {
                for (size_t j = 0; j < npe; ++j) {
                    valbfn[j] += CSN[K][I] * std::exp(-EX[K][I] * scrv[j][3]);
                }
            }
            // Cartesian powers
            if (XEXP[K] == 1) {
                for (size_t j = 0; j < npe; ++j)
                    valbfn[j] *= scrv[j][0];
            } else if (XEXP[K] > 1) {
                for (size_t j = 0; j < npe; ++j)
                    valbfn[j] *= std::pow(scrv[j][0], XEXP[K]);
            }
            if (YEXP[K] == 1) {
                for (size_t j = 0; j < npe; ++j)
                    valbfn[j] *= scrv[j][1];
            } else if (YEXP[K] > 1) {
                for (size_t j = 0; j < npe; ++j)
                    valbfn[j] *= std::pow(scrv[j][1], YEXP[K]);
            }
            if (ZEXP[K] == 1) {
                for (size_t j = 0; j < npe; ++j)
                    valbfn[j] *= scrv[j][2];
            } else if (ZEXP[K] > 1) {
                for (size_t j = 0; j < npe; ++j)
                    valbfn[j] *= std::pow(scrv[j][2], ZEXP[K]);
            }

            for (size_t j = 0; j < npe; ++j)
                valao[j][K] = valbfn[j];
        }
    }
}

} // namespace OversetGrid
#include "StPot_overset_sub.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <mpi.h>
#include <iostream>
#include <cassert>

namespace OversetGrid {
extern int LMPI_myid, LMPI_master, LMPI_numprocs;
extern int NBasis, Natom, NumOrbFrm;
extern std::vector<std::vector<double>> cmat;           // cmat[NBasis][max_mo]
extern std::vector<int> OrbDegn;                        // OrbDegn[NumOrbFrm]
extern std::vector<double> OrbOccFrm;                   // OrbOccFrm[NumOrbFrm]
extern std::vector<NCoord> aCoord;                      // aCoord[Natom]
extern std::vector<int> IATMFR, IATMTO, NOPRIM;
extern std::vector<std::vector<double>> CSN, EX;
extern std::vector<int> XEXP, YEXP, ZEXP;
extern std::vector<int> MyFirstPoint, LastPoint;        // per processor
extern std::vector<GridDefType> GridDef;                // GridDef[0]: Ubergrid, GridDef[ig]: subgrids
extern std::vector<AngBIType> AngBI;
extern std::vector<int> WhichAR;
extern int N_SG;
extern std::vector<SubGridListType> List_SG;
extern std::vector<int> MyList_SG;
extern double PClose;
extern int iprnfg;
extern double unit_stdout;
extern double xfauang;
// Helper: Get the current MPI rank and size
static void get_mpi_info(int& myid, int& numprocs) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
}

// Helper: Returns which processor owns a global point index [Fortran: MyFirstPoint/LastPoint arrays]
int get_proc_for_index(const std::vector<int>& MyFirstPoint, const std::vector<int>& LastPoint, int idx) {
    for (size_t p = 0; p < MyFirstPoint.size(); ++p) {
        if (MyFirstPoint[p] <= idx && idx <= LastPoint[p])
            return static_cast<int>(p);
    }
    return -1; // Not found (should not happen)
}

void ContinuityCheck_CD(const Overset_CD_t& Orb_CD)
{
    constexpr double dtol = 1e-12;
    int myid, numprocs;
    get_mpi_info(myid, numprocs);

    // These should be initialized with global point index per processor
    extern std::vector<int> MyFirstPoint, LastPoint;   // Global point index ranges per proc
    extern std::vector<int> WhichAR;                   // Angular region for each grid point
    extern int FEMDef0_ni, FEMDef0_nr;                 // FEM info for Ubergrid
    extern int MAX_ANG_POINTS;                         // Max number of angular points per shell

    // --- Check Ubergrid element boundaries ---
    for (int iel = 1; iel < FEMDef0_ni; ++iel) {
        int ileft = iel * FEMDef0_nr;
        int iright = ileft + 1;
        int procleft = get_proc_for_index(MyFirstPoint, LastPoint, ileft);
        int procright = get_proc_for_index(MyFirstPoint, LastPoint, iright);

        std::vector<std::complex<double>> pleft(MAX_ANG_POINTS, {0.0, 0.0});
        std::vector<std::complex<double>> pright(MAX_ANG_POINTS, {0.0, 0.0});

        if (procleft == procright) {
            if (myid == procleft) {
                int irgp_left = ileft - MyFirstPoint[myid];
                int irgp_right = iright - MyFirstPoint[myid];
                // Copy data (assumes Orb_CD.UG[irgp].m is std::vector<std::complex<double>>)
                size_t n = Orb_CD.UG[irgp_left].m.size();
                for (size_t j = 0; j < n; ++j) pleft[j] = Orb_CD.UG[irgp_left].m[j];
                n = Orb_CD.UG[irgp_right].m.size();
                for (size_t j = 0; j < n; ++j) pright[j] = Orb_CD.UG[irgp_right].m[j];
                if (WhichAR[ileft] == WhichAR[iright]) {
                    for (size_t j = 0; j < n; ++j)
                        if (std::abs(pleft[j] - pright[j]) >= dtol)
                            throw std::runtime_error("CONTINUITY CHECK FAILED");
                }
            }
        } else {
            if (myid == procright) {
                int irgp_right = iright - MyFirstPoint[myid];
                size_t n = Orb_CD.UG[irgp_right].m.size();
                for (size_t j = 0; j < n; ++j) pright[j] = Orb_CD.UG[irgp_right].m[j];
                int LMPI_tag = ileft;
                MPI_Send(pright.data(), 2*n, MPI_DOUBLE, procleft, LMPI_tag, MPI_COMM_WORLD);
            } else if (myid == procleft) {
                int LMPI_tag = ileft;
                int irgp_left = ileft - MyFirstPoint[myid];
                size_t n = Orb_CD.UG[irgp_left].m.size();
                for (size_t j = 0; j < n; ++j) pleft[j] = Orb_CD.UG[irgp_left].m[j];
                MPI_Status status;
                MPI_Recv(pright.data(), 2*n, MPI_DOUBLE, procright, LMPI_tag, MPI_COMM_WORLD, &status);
                if (WhichAR[ileft] == WhichAR[iright]) {
                    for (size_t j = 0; j < n; ++j)
                        if (std::abs(pleft[j] - pright[j]) >= dtol)
                            throw std::runtime_error("CONTINUITY CHECK FAILED");
                }
            }
        }
    }

    // --- Check subgrid element boundaries ---
    extern std::vector<int> MyList_SG; // List of subgrids handled by this proc
    extern std::vector<int> List_SG_ig; // Index in FEMDef for each subgrid
    extern std::vector<std::vector<std::vector<std::complex<double>>>> SGm; // SGm[isg][ang][radial]
    for (int myisg : MyList_SG) {
        int isg = myisg;
        int ig = List_SG_ig[isg];
        int ni = FEMDef0_ni; // TODO: Replace with actual value for subgrid ig
        int nr = FEMDef0_nr; // TODO: Replace with actual value for subgrid ig
        for (int iel = 1; iel < ni; ++iel) {
            int ileft = iel * nr;
            int iright = ileft + 1;
            std::vector<std::complex<double>> pleft(MAX_ANG_POINTS, {0.0, 0.0});
            std::vector<std::complex<double>> pright(MAX_ANG_POINTS, {0.0, 0.0});
            // Fill from Orb_CD.SG[isg].m[ang][radial]
            for (size_t j = 0; j < SGm[isg].size(); ++j) {
                pleft[j] = SGm[isg][j][ileft];
                pright[j] = SGm[isg][j][iright];
            }
            for (size_t j = 0; j < SGm[isg].size(); ++j)
                if (std::abs(pleft[j] - pright[j]) >= dtol)
                    throw std::runtime_error("SUBGRID CONTINUITY CHECK FAILED");
        }
    }
}
void StPot_overset(
    // Optionally, pass references for output arrays
) {
    // === 1. Banner ===
    if (LMPI_myid == LMPI_master)
        std::cout << "StPot_overset - calculate static potential on Overset Grid" << std::endl;

    // === 2. Construct AO Density Matrix ===
    std::vector<std::vector<double>> RhoDen(NBasis, std::vector<double>(NBasis, 0.0));
    int mo = 0;
    for (int iform = 0; iform < NumOrbFrm; ++iform) {
        double OrbOcc = OrbOccFrm[iform];
        for (int modegen = 0; modegen < OrbDegn[iform]; ++modegen) {
            mo += 1;
            for (int iao = 0; iao < NBasis; ++iao)
                for (int jao = 0; jao < NBasis; ++jao)
                    RhoDen[iao][jao] += OrbOcc * cmat[iao][mo] * cmat[jao][mo];
        }
    }

    // === 3. Ubergrid: Loop over radial points on this processor ===
    int iAR = -1;
    std::vector<double> sphv, cphv, EPIIntegrals;
    std::vector<Cvec> CartGridPt;
    std::vector<std::vector<double>> valao;
    int totug = 0;

    for (int irg = MyFirstPoint[LMPI_myid]; irg <= LastPoint[LMPI_myid]; ++irg) {
        int irgp = irg - MyFirstPoint[LMPI_myid];
        int newAR = WhichAR[irg];

        // (Re)allocate arrays for new angular region if necessary
        if (newAR != iAR) {
            iAR = newAR;
            int nphii   = AngBI[iAR].nphi;
            int ntheti  = AngBI[iAR].ntheta;
            int nthphi  = AngBI[iAR].nthph;
            sphv.assign(nphii, 0.0);
            cphv.assign(nphii, 0.0);
            CartGridPt.assign(nthphi, Cvec{});
            valao.assign(nthphi, std::vector<double>(NBasis, 0.0));
            EPIIntegrals.assign(nthphi, 0.0);
            for (int iph = 0; iph < nphii; ++iph) {
                cphv[iph] = std::cos(AngBI[iAR].phi[iph]);
                sphv[iph] = std::sin(AngBI[iAR].phi[iph]);
            }
        }
        double r = GridDef[0].r[irg];

        // Fill CartGridPt with all (x,y,z) for the shell
        int ij = -1;
        for (int ith = 0; ith < AngBI[iAR].ntheta; ++ith) {
            double cost = AngBI[iAR].ctheta[ith];
            double sint = std::sqrt(std::max(0.0, 1.0 - cost*cost));
            for (int iph = 0; iph < AngBI[iAR].nphi; ++iph) {
                ij += 1;
                CartGridPt[ij].x = r * sint * cphv[iph];
                CartGridPt[ij].y = r * sint * sphv[iph];
                CartGridPt[ij].z = r * cost;
            }
        }
        assert(ij + 1 == (int)CartGridPt.size());

        // === Evaluate AO values at all points in shell ===
        CalcAO(CartGridPt, valao);

        // === Compute nuclear attraction part at each point ===
        for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
            double nucpart = 0.0;
            for (int iatom = 0; iatom < Natom; ++iatom) {
                double dx = CartGridPt[ij].x - aCoord[iatom].vec[0];
                double dy = CartGridPt[ij].y - aCoord[iatom].vec[1];
                double dz = CartGridPt[ij].z - aCoord[iatom].vec[2];
                double rnuc = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (rnuc > PClose)
                    nucpart -= aCoord[iatom].z / rnuc;
                else
                    nucpart -= aCoord[iatom].z / (rnuc + PClose);
            }
            EPIIntegrals[ij] = nucpart; // Add more terms if needed (e.g. via external libraries)
        }

        // === Contract AO values to build density at each grid point ===
        std::vector<double> Den_UG(CartGridPt.size(), 0.0);
        for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
            double den = 0.0;
            for (int iao = 0; iao < NBasis; ++iao)
                for (int jao = 0; jao < NBasis; ++jao)
                    den += RhoDen[iao][jao] * valao[ij][iao] * valao[ij][jao];
            Den_UG[ij] = den;
        }

        // === Output for each point (can be stored in Vstat_CD, etc.) ===
        // For simplicity, print a few values
        if (iprnfg > 1 && LMPI_myid == LMPI_master && irgp < 5) {
            for (int ij = 0; ij < (int)CartGridPt.size() && ij < 3; ++ij) {
                std::cout << "Grid point r=" << r*xfauang
                          << " x=" << CartGridPt[ij].x << " y=" << CartGridPt[ij].y << " z=" << CartGridPt[ij].z
                          << " Pot=" << EPIIntegrals[ij] << " Den=" << Den_UG[ij] << std::endl;
            }
        }
        totug += CartGridPt.size();
    }

    // === 4. Subgrids ===
    int totsg = 0;
    for (int myisg : MyList_SG) {
        int isg = myisg;
        int ig = List_SG[isg].ig;
        int center = List_SG[isg].center;
        int iAR = /* logic for iAR as in Fortran */;
        int nphii   = AngBI[iAR].nphi;
        int ntheti  = AngBI[iAR].ntheta;
        int nthphi  = AngBI[iAR].nthph;
        sphv.assign(nphii, 0.0);
        cphv.assign(nphii, 0.0);
        CartGridPt.assign(nthphi, Cvec{});
        valao.assign(nthphi, std::vector<double>(NBasis, 0.0));
        EPIIntegrals.assign(nthphi, 0.0);
        for (int iph = 0; iph < nphii; ++iph) {
            cphv[iph] = std::cos(AngBI[iAR].phi[iph]);
            sphv[iph] = std::sin(AngBI[iAR].phi[iph]);
        }

        for (int irg = 0; irg < (int)GridDef[ig].n; ++irg) {
            double r = GridDef[ig].r[irg];
            int ij = -1;
            for (int ith = 0; ith < ntheti; ++ith) {
                double cost = AngBI[iAR].ctheta[ith];
                double sint = std::sqrt(std::max(0.0, 1.0 - cost*cost));
                for (int iph = 0; iph < nphii; ++iph) {
                    ij += 1;
                    CartGridPt[ij].x = r * sint * cphv[iph] + aCoord[center].vec[0];
                    CartGridPt[ij].y = r * sint * sphv[iph] + aCoord[center].vec[1];
                    CartGridPt[ij].z = r * cost + aCoord[center].vec[2];
                }
            }
            CalcAO(CartGridPt, valao);

            // Repeat same as above: compute potentials, density, output/store
            for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
                double nucpart = 0.0;
                for (int iatom = 0; iatom < Natom; ++iatom) {
                    double dx = CartGridPt[ij].x - aCoord[iatom].vec[0];
                    double dy = CartGridPt[ij].y - aCoord[iatom].vec[1];
                    double dz = CartGridPt[ij].z - aCoord[iatom].vec[2];
                    double rnuc = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (rnuc > PClose)
                        nucpart -= aCoord[iatom].z / rnuc;
                    else
                        nucpart -= aCoord[iatom].z / (rnuc + PClose);
                }
                EPIIntegrals[ij] = nucpart;
            }
            std::vector<double> Den_SG(CartGridPt.size(), 0.0);
            for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
                double den = 0.0;
                for (int iao = 0; iao < NBasis; ++iao)
                    for (int jao = 0; jao < NBasis; ++jao)
                        den += RhoDen[iao][jao] * valao[ij][iao] * valao[ij][jao];
                Den_SG[ij] = den;
            }
            totsg += CartGridPt.size();
        }
    }

    // === 5. Reduce, report, finalize ===
    int total_ug, total_sg;
    MPI_Allreduce(&totug, &total_ug, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&totsg, &total_sg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (LMPI_myid == LMPI_master) {
        std::cout << "Total Ubergrid points: " << total_ug << std::endl;
        std::cout << "Total Subgrid points: " << total_sg << std::endl;
        std::cout << "Grid size would take " << ((total_ug + total_sg)*sizeof(double)*2.0/1024.0/1024.0)
                  << " MB per grid function" << std::endl;
    }
}

void FegeCorrect_overset(double engv, double fegescale) {
    if (LMPI_myid == LMPI_master)
        std::cout << "Fege_overset - calculate fege potential on Overset Grid" << std::endl;

    // === 1. Construct AO Density Matrix ===
    std::vector<std::vector<double>> RhoDen(NBasis, std::vector<double>(NBasis, 0.0));
    int mo = 0;
    for (int iform = 0; iform < NumOrbFrm; ++iform) {
        double OrbOcc = OrbOccFrm[iform];
        for (int modegen = 0; modegen < OrbDegn[iform]; ++modegen) {
            mo += 1;
            for (int iao = 0; iao < NBasis; ++iao)
                for (int jao = 0; jao < NBasis; ++jao)
                    RhoDen[iao][jao] += OrbOcc * cmat[iao][mo] * cmat[jao][mo];
        }
    }

    // === 2. Read ecor (correlation energy parameter for Fege, set as needed) ===
    double ecor = 0.5; // TODO: Replace with your project logic or input reader

    double fscl = fegescale;
    // Fege grid potentials to write (output containers)
    std::vector<std::vector<double>> Vxstat_CD_UG, Den_CD_UG;

    // === 3. Main loop: Ubergrid ===
    int iAR = -1;
    std::vector<double> sphv, cphv;
    std::vector<Cvec> CartGridPt;
    std::vector<std::vector<double>> valao;
    for (int irg = MyFirstPoint[LMPI_myid]; irg <= LastPoint[LMPI_myid]; ++irg) {
        int irgp = irg - MyFirstPoint[LMPI_myid];
        int newAR = WhichAR[irg];
        if (newAR != iAR) {
            iAR = newAR;
            int nphii   = AngBI[iAR].nphi;
            int ntheti  = AngBI[iAR].ntheta;
            int nthphi  = AngBI[iAR].nthph;
            sphv.assign(nphii, 0.0);
            cphv.assign(nphii, 0.0);
            CartGridPt.assign(nthphi, Cvec{});
            valao.assign(nthphi, std::vector<double>(NBasis, 0.0));
            for (int iph = 0; iph < nphii; ++iph) {
                cphv[iph] = std::cos(AngBI[iAR].phi[iph]);
                sphv[iph] = std::sin(AngBI[iAR].phi[iph]);
            }
        }
        double r = GridDef[0].r[irg];

        // Set up all grid points for this shell
        int ij = -1;
        for (int ith = 0; ith < AngBI[iAR].ntheta; ++ith) {
            double cost = AngBI[iAR].ctheta[ith];
            double sint = std::sqrt(std::max(0.0, 1.0 - cost*cost));
            for (int iph = 0; iph < AngBI[iAR].nphi; ++iph) {
                ij += 1;
                CartGridPt[ij].x = r * sint * cphv[iph];
                CartGridPt[ij].y = r * sint * sphv[iph];
                CartGridPt[ij].z = r * cost;
            }
        }

        // === Evaluate AO values at all points ===
        CalcAO(CartGridPt, valao);

        // === Compute density at each point ===
        std::vector<double> Den_UG(CartGridPt.size(), 0.0);
        for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
            double den = 0.0;
            for (int iao = 0; iao < NBasis; ++iao)
                for (int jao = 0; jao < NBasis; ++jao)
                    den += RhoDen[iao][jao] * valao[ij][iao] * valao[ij][jao];
            Den_UG[ij] = den;
        }
        // === Compute Fege potential at each point ===
        std::vector<double> Vxstat_UG(CartGridPt.size(), 0.0);
        for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
            Vxstat_UG[ij] = efege(Den_UG[ij], engv, ecor) * fscl;
        }
        Vxstat_CD_UG.push_back(Vxstat_UG);
        Den_CD_UG.push_back(Den_UG);
    }

    // === 4. Subgrids (repeat, shifted coordinates) ===
    for (int myisg : MyList_SG) {
        int isg = myisg;
        int ig = List_SG[isg].ig;
        int center = List_SG[isg].center;
        int iAR = /* logic for iAR as in Fortran */;
        int nphii   = AngBI[iAR].nphi;
        int ntheti  = AngBI[iAR].ntheta;
        int nthphi  = AngBI[iAR].nthph;
        sphv.assign(nphii, 0.0);
        cphv.assign(nphii, 0.0);
        CartGridPt.assign(nthphi, Cvec{});
        valao.assign(nthphi, std::vector<double>(NBasis, 0.0));
        for (int iph = 0; iph < nphii; ++iph) {
            cphv[iph] = std::cos(AngBI[iAR].phi[iph]);
            sphv[iph] = std::sin(AngBI[iAR].phi[iph]);
        }

        for (int irg = 0; irg < (int)GridDef[ig].n; ++irg) {
            double r = GridDef[ig].r[irg];
            int ij = -1;
            for (int ith = 0; ith < ntheti; ++ith) {
                double cost = AngBI[iAR].ctheta[ith];
                double sint = std::sqrt(std::max(0.0, 1.0 - cost*cost));
                for (int iph = 0; iph < nphii; ++iph) {
                    ij += 1;
                    CartGridPt[ij].x = r * sint * cphv[iph] + aCoord[center].vec[0];
                    CartGridPt[ij].y = r * sint * sphv[iph] + aCoord[center].vec[1];
                    CartGridPt[ij].z = r * cost + aCoord[center].vec[2];
                }
            }
            CalcAO(CartGridPt, valao);

            std::vector<double> Den_SG(CartGridPt.size(), 0.0);
            for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
                double den = 0.0;
                for (int iao = 0; iao < NBasis; ++iao)
                    for (int jao = 0; jao < NBasis; ++jao)
                        den += RhoDen[iao][jao] * valao[ij][iao] * valao[ij][jao];
                Den_SG[ij] = den;
            }
            std::vector<double> Vxstat_SG(CartGridPt.size(), 0.0);
            for (int ij = 0; ij < (int)CartGridPt.size(); ++ij) {
                Vxstat_SG[ij] = efege(Den_SG[ij], engv, ecor) * fscl;
            }
            // Store or output Vxstat_SG as needed
        }
    }

    // --- Optionally, add MPI reduction, output, timing, as in StPot_overset ---

    if (LMPI_myid == LMPI_master)
        std::cout << "FegeCorrect_overset completed." << std::endl;
}
} // namespace OversetGrid

