#include "PatchAlgo.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "linalg_utils.hpp" // for LAPACK-style solvers

namespace PatchAlgo {

// This helper solves a small real-valued linear system using Eigen or LAPACK.
static void solve_linear_system(std::vector<std::vector<double>>& A, std::vector<double>& b) {

    int n = static_cast<int>(b.size());
    for (int i = 0; i < n; ++i) {
        // Pivot
        double maxv = std::abs(A[i][i]);
        int maxrow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > maxv) {
                maxv = std::abs(A[k][i]);
                maxrow = k;
            }
        }
        if (maxrow != i) {
            std::swap(A[i], A[maxrow]);
            std::swap(b[i], b[maxrow]);
        }
        // Eliminate
        for (int k = i + 1; k < n; ++k) {
            double c = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) A[k][j] -= c * A[i][j];
            b[k] -= c * b[i];
        }
    }
    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        b[i] /= A[i][i];
        for (int k = 0; k < i; ++k) b[k] -= A[k][i] * b[i];
    }
}

void LucchesePartitionEG_Re(
    std::vector<std::vector<double>>& PartialWave,
    const AngExp& AExpX,
    int NumRadPt,
    int EG_idx,
    int iType,
    bool useabel
) {
    const auto& Actual_GridDef = useabel ? GridDef_Ab : GridDef;
    const auto& Actual_lvals_EG = useabel ? lvalsb_EG : lvals_EG;
    const auto& Actual_FEMDef = useabel ? FEMDef_Ab : FEMDef;

    int PolyOrder = 2;
    const auto& grid = Actual_GridDef[EG_idx];
    double r_max = grid.r[NumRadPt - 1];
    int nr = Actual_FEMDef[EG_idx].nr;
    int ni = Actual_FEMDef[EG_idx].ni;
    int LeftBound = (ni - 1) * nr;
    int pow_offset = useabel ? 0 : 1;

    std::vector<double> r_powers(PolyOrder);
    std::vector<std::vector<double>> DerivSubMat(nr, std::vector<double>(nr));
    for (int i = 0; i < nr; ++i)
        for (int j = 0; j < nr; ++j)
            DerivSubMat[i][j] = grid.deriv.m[LeftBound + i][LeftBound + j];
    
    for (int i = 0; i < AExpX.nlht; ++i) {
        std::vector<std::vector<double>> LeftHandSide(PolyOrder, std::vector<double>(PolyOrder, 0.0));
        std::vector<double> RightHandSide(PolyOrder, 0.0);
        int l_ang = Actual_lvals_EG[i][iType][EG_idx];

        for (int k = 0; k < PolyOrder; ++k) {
            int r_pow = l_ang + pow_offset + 2 * k;
            r_powers[k] = std::pow(r_max, r_pow);
            int prefactor = 1;
            for (int j = 0; j < PolyOrder; ++j) {
                if (r_pow < 0) break;
                LeftHandSide[j][k] = prefactor * r_powers[k];
                prefactor *= r_pow--;
            }
        }

        std::vector<double> CD_deriv_temp(nr);
        for (int p = 0; p < nr; ++p)
            CD_deriv_temp[p] = PartialWave[i][LeftBound + p];

        for (int j = 0; j < PolyOrder; ++j) {
            RightHandSide[j] = CD_deriv_temp[nr - 1];
            std::fill(CD_deriv_temp.begin(), CD_deriv_temp.end(), 0.0);
            for (int k = 0; k < nr; ++k)
                for (int m = 0; m < nr; ++m)
                    CD_deriv_temp[k] += DerivSubMat[k][m] * CD_deriv_temp[m];
        }

        solve_linear_system(LeftHandSide, RightHandSide);

        for (int j = 0; j < NumRadPt; ++j) {
            double grid_r_j = grid.r[j];
            for (int k = 0; k < PolyOrder; ++k) {
                int r_pow = l_ang + pow_offset + 2 * k;
                PartialWave[i][j] -= RightHandSide[k] * std::pow(grid_r_j, r_pow);
            }
        }
    }
}
void EvalEGOnUG_PW(
    Overset_CD_t& Phi_CD,
    const Overset_PW_t& Psi_PW,
    const InterpolatingMatrices& IM,
    bool useabel
) {
    // Select actual arrays for abelian/non-abelian
    const auto& Actual_List_SG      = useabel ? List_SG_Ab : List_SG;
    const auto& Actual_WhichProc_EG = useabel ? WhichProc_EG_Ab : WhichProc_EG;
    const auto& Actual_FEMDef       = useabel ? FEMDef_Ab : FEMDef;
    const auto& irUGov              = useabel ? irUGov_Ab : irUGov;
    const auto& iwUGov              = useabel ? iwUGov_Ab : iwUGov;

    int nprocs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // Allocate temporary buffers for all-to-all comms
    std::vector<std::vector<int>> opUGindr(nprocs);
    std::vector<std::vector<int>> opUGindw(nprocs);
    std::vector<std::vector<std::complex<double>>> opUGval(nprocs);

    int isgovpt = 0;

    for (size_t isg = 0; isg < N_SG; ++isg) {
        int ig = Actual_List_SG[isg].ig;
        if (Actual_WhichProc_EG[ig] == myid) {
            const auto& irUGs = irUGov[isg].v;
            const auto& iwUGs = iwUGov[isg].v;
            for (size_t iovpt = 0; iovpt < irUGs.size(); ++iovpt) {
                int irug = irUGs[iovpt];
                int iwug = iwUGs[iovpt];
                const auto& imatsgR = IM.SGIMatR[isgovpt + iovpt];
                const auto& imatsgBLM = IM.SGIMatBLM[isgovpt + iovpt];
                int jel = IM.SGIElInd[isgovpt + iovpt];
                int nr = Actual_FEMDef[ig].nr;

                for (int iproc = 0; iproc < nprocs; ++iproc) {
                    if (MyFirstPoint[iproc] <= irug && irug <= LastPoint[iproc]) {
                        std::complex<double> val = 0.0;
                        for (size_t ilh = 0; ilh < imatsgBLM.m.size(); ++ilh) {
                            auto it_beg = Psi_PW.EG[ig].m[ilh].begin() + nr * (jel - 1);
                            auto it_end = it_beg + nr;
                            val += std::inner_product(imatsgR.m.begin(), imatsgR.m.end(), it_beg, 0.0) * imatsgBLM.m[ilh];
                        }
                        if (iproc == myid) {
                            int irgp = irug - MyFirstPoint[myid];
                            Phi_CD.UG[irgp].m[iwug] -= val;
                        } else {
                            opUGindr[iproc].push_back(irug);
                            opUGindw[iproc].push_back(iwug);
                            opUGval[iproc].push_back(val);
                        }
                    }
                }
            }
        }
        isgovpt += irUGov[isg].n;
    }

    // All-to-all communication: send updates to Ubergrid points on other ranks
    for (int iproc = 0; iproc < nprocs; ++iproc) {
        for (int iprocto = 0; iprocto < nprocs; ++iprocto) {
            if (iproc == iprocto) continue;
            if (iprocto == myid) {
                int nrecv = 0;
                MPI_Status status;
                MPI_Recv(&nrecv, 1, MPI_INT, iproc, 0, MPI_COMM_WORLD, &status);
                std::vector<int> irrecv(nrecv), iwrecv(nrecv);
                std::vector<std::complex<double>> valrecv(nrecv);
                MPI_Recv(irrecv.data(), nrecv, MPI_INT, iproc, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(iwrecv.data(), nrecv, MPI_INT, iproc, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(valrecv.data(), nrecv * 2, MPI_DOUBLE, iproc, 3, MPI_COMM_WORLD, &status);
                for (int k = 0; k < nrecv; ++k) {
                    int irgp = irrecv[k] - MyFirstPoint[myid];
                    Phi_CD.UG[irgp].m[iwrecv[k]] -= valrecv[k];
                }
            } else if (iproc == myid) {
                int nsend = opUGval[iprocto].size();
                MPI_Send(&nsend, 1, MPI_INT, iprocto, 0, MPI_COMM_WORLD);
                MPI_Send(opUGindr[iprocto].data(), nsend, MPI_INT, iprocto, 1, MPI_COMM_WORLD);
                MPI_Send(opUGindw[iprocto].data(), nsend, MPI_INT, iprocto, 2, MPI_COMM_WORLD);
                MPI_Send(opUGval[iprocto].data(), nsend * 2, MPI_DOUBLE, iprocto, 3, MPI_COMM_WORLD);
            }
        }
    }
}
void Dump_CD_UG(
    const Overset_CD_t& Orb_CD_in,
    const std::string& debug_filename,
    int datlabel
) {
    using namespace detail;
    const double tol = 0.02;
    const double CosThetaTarget = CosThetaTargetDefault;
    const double PhiTarget = PhiTargetDefault;
    const double rTarget = rTargetDefault;

    std::ofstream ofs = openDebugOut(debug_filename);

    ofs << " TranCDtoPW Debug " << std::endl;
    ofs << "DataTag  ProcID   ithph  ctheta      phi      ir     r         Phi_CD%UG" << std::endl;

    for (int irg = MyFirstPoint[LMPI_myid]; irg <= LastPoint[LMPI_myid]; irg += 5) {
        int iAR = WhichAR(irg);
        int irgp = irg - MyFirstPoint[LMPI_myid];
        int nphii = AngBI[iAR].nphi;
        int nthetai = AngBI[iAR].ntheta;
        // Find closest theta/phi to the target
        double minThetaDiff = 999.0, minPhiDiff = 999.0;
        int iThetaTarget = 0, iPhiTarget = 0;
        for (int itheta = 0; itheta < nthetai; ++itheta) {
            double diff = std::abs(AngBI[iAR].ctheta[itheta] - CosThetaTarget);
            if (diff < minThetaDiff) {
                minThetaDiff = diff;
                iThetaTarget = itheta;
            }
        }
        for (int iphi = 0; iphi < nphii; ++iphi) {
            double diff = std::abs(AngBI[iAR].phi[iphi] - PhiTarget);
            if (diff < minPhiDiff) {
                minPhiDiff = diff;
                iPhiTarget = iphi;
            }
        }
        ofs << "DUMP_CD_UG "
            << std::setw(7) << datlabel
            << std::setw(5) << 1 // ithph always 1
            << std::fixed << std::setw(10) << std::setprecision(5) << AngBI[iAR].ctheta[iThetaTarget]
            << std::setw(10) << AngBI[iAR].phi[iPhiTarget]
            << std::setw(10) << (GridDef[0].r[irg] * XFAuAng)
            << std::setw(17) << std::scientific << std::setprecision(8)
            << std::real(Orb_CD_in.UG[irgp].m[0])
            << std::endl;
    }
}

// Dump a 2D slice of the Ubergrid over (theta, phi) for diagnostic/debugging
void Dump_CD_UG_2D(
    const Overset_CD_t& Orb_CD_in,
    const std::string& debug_filename,
    int datlabel
) {
    using namespace detail;
    std::ofstream ofs = openDebugOut(debug_filename);

    for (int irg = MyFirstPoint[LMPI_myid]; irg <= LastPoint[LMPI_myid]; ++irg) {
        int iAR = WhichAR(irg);
        int irgp = irg - MyFirstPoint[LMPI_myid];
        int nphii = AngBI[iAR].nphi;
        int nthetai = AngBI[iAR].ntheta;
        for (int itheta = 0; itheta < nthetai; itheta += 3) {
            for (int iphi = 0; iphi < nphii; iphi += 3) {
                int idx = itheta * nphii + iphi;
                ofs << "DUMP_CD_UG_2D "
                    << std::setw(6) << datlabel
                    << std::fixed << std::setw(10) << std::setprecision(5) << AngBI[iAR].ctheta[itheta]
                    << std::setw(10) << AngBI[iAR].phi[iphi]
                    << std::setw(10) << (GridDef[0].r[irg] * XFAuAng)
                    << std::setw(17) << std::scientific << std::setprecision(8)
                    << std::real(Orb_CD_in.UG[irgp].m[idx])
                    << std::endl;
            }
        }
    }
}

// Print out the Ubergrid gridpoint coordinates for debugging or plotting
void print_UG_grid(int datlabel) {
    for (int irg = MyFirstPoint[LMPI_myid]; irg <= LastPoint[LMPI_myid]; ++irg) {
        int iAR = WhichAR(irg);
        int irgp = irg - MyFirstPoint[LMPI_myid];
        int nphii = AngBI[iAR].nphi;
        int nthetai = AngBI[iAR].ntheta;
        int ithph = 1;
        for (int itheta = 0; itheta < nthetai; ++itheta) {
            for (int iphi = 0; iphi < nphii; ++iphi) {
                std::cout << "Print UG points DataTag " << datlabel
                          << "  ctheta " << std::fixed << std::setprecision(5) << AngBI[iAR].ctheta[itheta]
                          << "  phi " << AngBI[iAR].phi[iphi]
                          << "  itheta " << itheta
                          << " iphi " << iphi
                          << " irg " << irg
                          << "  r(AU) " << GridDef[0].r[irg]
                          << std::endl;
                ++ithph;
            }
        }
    }
}
} // namespace PatchAlgo
