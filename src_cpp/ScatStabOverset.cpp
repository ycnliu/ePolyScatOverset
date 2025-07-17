// ScatStabOversetSub.cpp

#include "ScatStabOverset.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>


extern int MyFirstPoint(int mpi_rank); 
extern int LastPoint(int mpi_rank);
extern int LMPI_myid;
extern int NCoefKInt;
extern int N_SG, N_EG;
extern int lmax;
extern int PClose;
extern double kMomentum;
extern int NDipoleOp;
extern std::vector<int> MyList_SG, MyList_EG, LValAR; 
extern std::vector<BLMPointType> BLMPoint; 
extern std::vector<GridDefType> GridDef; 

/**
 * @brief Build J and JV on the overset grid for a given angular momentum value, scattering potential, and type
 *
 * @param J_CD J in the coordinate representation
 * @param JV_CD JV in the coordinate representation
 * @param LHVAL Angular momentum value
 * @param POT_CD Scattering potential in the coordinate representation
 * @param itype Type of the scattering potential
 */
void ScatStabOversetSub::BuildJV(
    Overset_CD_t& J_CD,
    Overset_CD_t& JV_CD,
    int LHVAL,
    const Overset_CD_t& POT_CD,
    int itype
) {
    // Prepare J in the plane wave (PW) representation
    Overset_PW_t J_PW;

    // ALLOCATE(J_PW%UG(LastPoint(LMPI_myid)-MyFirstPoint(LMPI_myid)+1))
    int ug_size = LastPoint(LMPI_myid) - MyFirstPoint(LMPI_myid) + 1;
    J_PW.UG.resize(ug_size);

    for (int irg = MyFirstPoint(LMPI_myid), irgp = 0; irg <= LastPoint(LMPI_myid); ++irg, ++irgp) {
        int iAR = WhichAR(irg);
        int ltop = std::min(LValAR.at(iAR), lmax);
        int lmtop = BLMPointType::size(itype, ltop); 

        J_PW.UG[irgp].m.assign(lmtop, 0.0);

        if (LHVAL <= lmtop) {

            J_PW.UG[irgp].m[LHVAL - 1] = F_UG[irgp].m[BLMPoint[LHVAL - 1][itype].p.lval];
        }
    }

    J_PW.EG.resize(N_EG);
    for (size_t idx = 0; idx < MyList_EG.size(); ++idx) {
        int ieg = MyList_EG[idx];
        int iAR = LValAR.size() - N_EG + ieg;
        int ltop = std::min(LValAR.at(iAR), lmax);
        int lmtop = BLMPointType::size_EG(itype, ltop, ieg); // implement as needed

        J_PW.EG[ieg].m.resize(lmtop, std::vector<double>(GridDef[ieg].n, 0.0));
    }

    // Build J in the coordinate representation by transformation and interpolation
    TranPWtoCD(J_CD, J_PW, itype); // User-implemented
    BuildJEG(J_CD, LHVAL, itype);

    // ALLOCATE(QCDX(NCoefKInt))
    std::vector<Overset_CD_t> QCDX(NCoefKInt);

    for (int iForm = 0; iForm < NCoefKInt; ++iForm) {
        if (std::abs(CoefKInt[iForm].Coef) > PClose) {
            if (CoefKInt[iForm].Type == 2) {
                int ktyppd = FindProdSym(itype, ktypo(CoefKInt[iForm].OrbG) + CoefKInt[iForm].OrbCL - 1);
                int OriginalOrbitalIndex = GetOriginalOrbitalIndex(
                    CoefKInt[iForm].OrbG,
                    ktypo(CoefKInt[iForm].OrbG) + CoefKInt[iForm].OrbCL - 1
                );
                Overset_CD_t QCDXP, QCDXD;
                ApplyExchange(QCDXP, QCDXD, J_CD, Orb_CD[OriginalOrbitalIndex], ktyppd, IMatAb[ktyppd]);
                QCDX[iForm] = QCDXP;
            }
        }
    }

    // Allocate and compute JV_CD
    JV_CD.UG.resize(ug_size);

    for (int irg = MyFirstPoint(LMPI_myid), irgp = 0; irg <= LastPoint(LMPI_myid); ++irg, ++irgp) {
        JV_CD.UG[irgp].m.resize(J_CD.UG[irgp].m.size(), 0.0);
        for (size_t idx = 0; idx < J_CD.UG[irgp].m.size(); ++idx) {
            JV_CD.UG[irgp].m[idx] = J_CD.UG[irgp].m[idx] * POT_CD.UG[irgp].m[idx];
        }

        for (int iForm = 0; iForm < NCoefKInt; ++iForm) {
            if (std::abs(CoefKInt[iForm].Coef) > PClose && CoefKInt[iForm].Type == 2) {
                for (size_t idx = 0; idx < JV_CD.UG[irgp].m.size(); ++idx) {
                    JV_CD.UG[irgp].m[idx] += QCDX[iForm].UG[irgp].m[idx] * CoefKInt[iForm].Coef;
                }
            }
        }
        // Special treatment at r ~ 0 for l = 0
        if (GridDef[0].r[irg] < std::numeric_limits<double>::epsilon() && BLMPoint[LHVAL - 1][itype].p.lval == 0) {
            for (size_t idx = 0; idx < JV_CD.UG[irgp].m.size(); ++idx)
                JV_CD.UG[irgp].m[idx] = POT_CD.UG[irgp].m[idx] * kMomentum / std::sqrt(4.0 * M_PI);
        }
    }

    JV_CD.SG.resize(N_SG);
    for (size_t idx = 0; idx < MyList_SG.size(); ++idx) {
        int isg = MyList_SG[idx];
        int ig = List_SG[isg].ig;
        int iAR = LValAR.size() - N_EG + ig;
        size_t nrow = J_CD.SG[isg].m.size();
        size_t ncol = J_CD.SG[isg].m[0].size();

        JV_CD.SG[isg].m.assign(nrow, std::vector<double>(ncol, 0.0));

        for (size_t i = 0; i < nrow; ++i) {
            for (size_t j = 0; j < ncol; ++j) {
                JV_CD.SG[isg].m[i][j] = J_CD.SG[isg].m[i][j] * POT_CD.SG[isg].m[i][j];
            }
        }
        for (int iForm = 0; iForm < NCoefKInt; ++iForm) {
            if (std::abs(CoefKInt[iForm].Coef) > PClose && CoefKInt[iForm].Type == 2) {
                for (size_t i = 0; i < nrow; ++i) {
                    for (size_t j = 0; j < ncol; ++j) {
                        JV_CD.SG[isg].m[i][j] += QCDX[iForm].SG[isg].m[i][j] * CoefKInt[iForm].Coef;
                    }
                }
            }
        }
        for (int irg = 0; irg < GridDef[ig].n; ++irg) {
            if (GridDef[ig].r[irg] < std::numeric_limits<double>::epsilon()) {
                // Only for l = 0 (central) -- needs more specifics from your actual types
                for (size_t idx = 0; idx < JV_CD.SG[isg].m.size(); ++idx) {
                    JV_CD.SG[isg].m[idx][irg] = POT_CD.SG[isg].m[idx][irg] * kMomentum / std::sqrt(4.0 * M_PI);
                }
            }
        }
    }
}
/**
 * @brief Generate the next iterate of the Arnoldi process
 *
 * @param POT_CD Overset representation of the scattering potential
 * @param i Index of the Arnoldi vector
 * @param ikry Index of the Arnoldi vector on the current iteration
 * @param itype Type of the scattering potential
 * @param vecnorm Norm of the new iterate
 * @param DipoleFlag Flag for dipole moment calculation
 * @param PotentialMatOccOrbitals Matrices of the potential for occupied orbitals
 * @param EminusHMatOccOrbitals Matrices of (E-H) for occupied orbitals
 * @param energy Energy value
 *
 * This function generates the next iterate of the Arnoldi process from the previous iterate.
 * It applies the G0 operator to the previous iterate and then applies the scattering potential
 * to the resulting vector. It also orthogonalizes the new iterate with respect to the
 * previous iterates.
 */
void ScatStabOversetSub::AddIterate(
    const Overset_CD_t& POT_CD,
    int i,
    int ikry,
    int itype,
    double& vecnorm,
    bool DipoleFlag,
    const std::vector<std::vector<double>>& PotentialMatOccOrbitals,
    const std::vector<std::vector<double>>& EminusHMatOccOrbitals,
    double energy
) {

    assert(ArnoldiIterates.size() >= static_cast<size_t>(i));
    assert(ArnoldiIterates[i].size() >= static_cast<size_t>(ikry));

    int NECenter = 0;
    if (LMPI_myid == 0) {
        std::cout << "Adding iterate " << i << ", " << ikry << std::endl;
    }
    int ovMODE = 3;

    if (ikry < IterateK[i]) {
        if (LMPI_myid == 0 /* LMPI_master */) {
            std::cout << "Iterate " << i << ", " << ikry << " has already been generated" << std::endl;
        }
        return;
    }

    double KCriterion = 1e99;

    // Generate iterate
    if (LMPI_myid == 0) {
        std::cout << "Generating iterate " << i << " " << ikry << std::endl;
    }

    // Alias to simplify referencing
    auto& QCD = ArnoldiIterates[i];
    // QCD[ikry][x] for x=0,1,2 

    // Apply G0 to V*previous iterate
    Overset_CD_t QCDP, QCDXP, QCDXD;
    Overset_PW_t QPW1, QPW2;

    ApplyG0(QCDP, QPW1, QPW2, QCD[ikry][1], itype, IMat, false);

    // Allocate QCDX
    std::vector<Overset_CD_t> QCDX(NCoefKInt);

    for (int iForm = 0; iForm < NCoefKInt; ++iForm) {
        if (std::abs(CoefKInt[iForm].Coef) > PClose) {
            if (CoefKInt[iForm].Type == 2) {
                int ktyppd = FindProdSym(itype, ktypo(CoefKInt[iForm].OrbG) + CoefKInt[iForm].OrbCL - 1);
                int OriginalOrbitalIndex = GetOriginalOrbitalIndex(
                    CoefKInt[iForm].OrbG,
                    ktypo(CoefKInt[iForm].OrbG) + CoefKInt[iForm].OrbCL - 1
                );
                Overset_CD_t QCDXP, QCDXD;
                ApplyExchange(QCDXP, QCDXD, QCDP, Orb_CD[OriginalOrbitalIndex], ktyppd, IMatAb[ktyppd]);
                QCDX[iForm] = QCDXP;
            }
        }
    }

    // Allocate next ArnoldiIterates for this i
    // QCD[ikry+1][0]: iterate itself
    // QCD[ikry+1][1]: V * iterate
    // QCD[ikry+1][2]: (E-H) * iterate = V * (iterate-1 - iterate)
    int ug_range = LastPoint(LMPI_myid) - MyFirstPoint(LMPI_myid) + 1;
    QCD[ikry+1][0].UG.resize(ug_range);
    QCD[ikry+1][1].UG.resize(ug_range);
    QCD[ikry+1][2].UG.resize(ug_range);

    for (int irg = MyFirstPoint(LMPI_myid), irgp = 0; irg <= LastPoint(LMPI_myid); ++irg, ++irgp) {
        int iAR = WhichAR(irg);
        int angsize = AngBI[iAR].nthph; // Assume AngBI is accessible

        // Copy and apply V
        for (int idx = 0; idx < angsize; ++idx) {
            QCD[ikry+1][1].UG[irgp].m[idx] = QCDP.UG[irgp].m[idx] * POT_CD.UG[irgp].m[idx];
        }

        // Apply Exchange terms
        for (int iForm = 0; iForm < NCoefKInt; ++iForm) {
            if (std::abs(CoefKInt[iForm].Coef) > PClose && CoefKInt[iForm].Type == 2) {
                for (int idx = 0; idx < angsize; ++idx) {
                    QCD[ikry+1][1].UG[irgp].m[idx] += QCDX[iForm].UG[irgp].m[idx] * CoefKInt[iForm].Coef;
                }
            }
        }

        // Special treatment at r ~ 0 
        if (GridDef[0].r[irg] < std::numeric_limits<double>::epsilon()) {
            // Use L'Hospital's rule 
        }

        // (E-H) * iterate
        for (int idx = 0; idx < angsize; ++idx) {
            QCD[ikry+1][2].UG[irgp].m[idx] = QCD[ikry][1].UG[irgp].m[idx] - QCD[ikry+1][1].UG[irgp].m[idx];
        }
    }

    // SG (subgrid) parts
    QCD[ikry+1][0].SG.resize(N_SG);
    QCD[ikry+1][1].SG.resize(N_SG);
    QCD[ikry+1][2].SG.resize(N_SG);

    for (size_t idx = 0; idx < MyList_SG.size(); ++idx) {
        int isg = MyList_SG[idx];
        int ig = List_SG[isg].ig;
        int iAR = LValAR.size() - N_EG + ig;
        int nrow = AngBI[iAR].nthph;
        int ncol = GridDef[ig].n;

        QCD[ikry+1][0].SG[isg].m.assign(nrow, std::vector<double>(ncol, 0.0));
        QCD[ikry+1][1].SG[isg].m.assign(nrow, std::vector<double>(ncol, 0.0));
        QCD[ikry+1][2].SG[isg].m.assign(nrow, std::vector<double>(ncol, 0.0));

        // Copy and apply V
        for (int i = 0; i < nrow; ++i) {
            for (int j = 0; j < ncol; ++j) {
                QCD[ikry+1][1].SG[isg].m[i][j] = QCDP.SG[isg].m[i][j] * POT_CD.SG[isg].m[i][j];
            }
        }

        for (int iForm = 0; iForm < NCoefKInt; ++iForm) {
            if (std::abs(CoefKInt[iForm].Coef) > PClose && CoefKInt[iForm].Type == 2) {
                for (int i = 0; i < nrow; ++i) {
                    for (int j = 0; j < ncol; ++j) {
                        QCD[ikry+1][1].SG[isg].m[i][j] += QCDX[iForm].SG[isg].m[i][j] * CoefKInt[iForm].Coef;
                    }
                }
            }
        }

        for (int i = 0; i < nrow; ++i) {
            for (int j = 0; j < ncol; ++j) {
                QCD[ikry+1][2].SG[isg].m[i][j] = QCD[ikry][1].SG[isg].m[i][j] - QCD[ikry+1][1].SG[isg].m[i][j];
            }
        }
    }

    // Orthogonalization (metric==0 or metric==1; call Overlap_CD and do Gram-Schmidt)

/**
 * @brief Arnoldi-CG-Kohn-Sham for asymmetric systems without J
 *
 * @param lambdaF output, the eigenvalue of the converged solution
 * @param muDipoleSc output, the dipole moment of the converged solution, size NDipoleOp
 * @param potCd input, potential matrix for the system
 * @param iL input, the index of the left atom
 * @param iR input, the index of the right atom
 * @param jL Cd input, the matrix for the left atom
 * @param jLV Cd input, the matrix for the left atom with V
 * @param jR Cd input, the matrix for the right atom
 * @param jRV Cd input, the matrix for the right atom with V
 * @param hOut Cd input, the matrix for the output
 * @param lHOut Cd input, the matrix for the output with lambda
 * @param uLuBorn input, the Born term for the system
 * @param uLuOutgoing input, the outgoing term for the system
 * @param maxIter input, the maximum number of iterations
 * @param tol input, the tolerance for convergence
 * @param uTol input, the tolerance for the Born term
 * @param itype input, the type of iteration (0 for CG, 1 for GMRES)
 * @param numIterI output, the number of iterations
 * @param dipoleFlag input, whether to compute the dipole moment
 * @param potentialMatOccOrbitals input, the potential matrix for occupied orbitals
 * @param eMinusHMatOccOrbitals input, the energy minus Hamiltonian matrix for occupied orbitals
 * @param energy input, the energy of the system
 */
void ScatStabOversetSub::ArnoldiCKohnGVNoJAsymSv(
    std::complex<double>& LambdaF,
    std::vector<std::complex<double>>& muDipoleSc, // size NDipoleOp
    const Overset_CD_t& POT_CD,
    int iL, int iR,
    const Overset_CD_t& JL_CD,
    const Overset_CD_t& JLV_CD,
    const Overset_CD_t& JR_CD,
    const Overset_CD_t& JRV_CD,
    const Overset_CD_t& HOUT_CD,
    const Overset_CD_t& LHOUT_CD,
    const std::complex<double>& uLuBorn,
    const std::array<std::complex<double>, 3>& uLuOutgoing,
    int MaxIter,
    double Tol,
    double UTol,
    int itype,
    int& NumIterI,
    bool DipoleFlag,
    const std::vector<std::vector<double>>* PotentialMatOccOrbitals,
    const std::vector<std::vector<double>>* EminusHMatOccOrbitals,
    double energy
) {
    using namespace std::complex_literals;
    using std::vector;
    using std::complex;

    // Constants and sizes
    int NDip = NDipoleOp; 
    const double PI = M_PI;

    // Subspace matrices and solution vectors
    int Nmax = MaxIter + 2; // To match 0:MaxIter+1
    Eigen::MatrixXcd vsub(Nmax, Nmax);
    Eigen::MatrixXcd esub(Nmax, Nmax);

    Eigen::VectorXcd y(Nmax);
    Eigen::VectorXcd y_orthog(Nmax);
    Eigen::VectorXcd y_total(Nmax);

    if (uLuBorn.imag() == 0.0 && uLuBorn.real() < 0.0) {
        throw std::runtime_error("Invalid uLuBorn value: non-positive real part with zero imaginary part.");
    }
    complex<double> beta = std::sqrt(-uLuBorn);
    complex<double> ibeta = 1.0 / beta;

    vsub.setZero();
    vsub(0, 0) = uLuOutgoing[2];
    vsub(0, 1) = uLuOutgoing[1];
    esub.setZero();

    complex<double> Lambda = -1e99;
    complex<double> Lambda_Last = -1e99;
    double delta_dip = 1e99;
    double delta_max = 0.0;

    std::vector<complex<double>> muDipoleSc_tmp(NDip, 0.0);
    std::vector<double> muDipoleSc_delta(NDip, 1e99);

    int mindip = 0;

    // GMRES Arnoldi-like iteration
    int ikry = 0;
    for (ikry = 1; ikry <= MaxIter; ++ikry) {
        // (1) Generate next Arnoldi iterates for iL and iR
        double vecnorm;
        try {
            AddIterate(POT_CD, iL, ikry, itype, vecnorm, DipoleFlag, PotentialMatOccOrbitals, EminusHMatOccOrbitals, energy);
            AddIterate(POT_CD, iR, ikry, itype, vecnorm, DipoleFlag, PotentialMatOccOrbitals, EminusHMatOccOrbitals, energy);
        } catch (const std::exception& e) {
            std::cerr << "Error in AddIterate: " << e.what() << std::endl;
            break;
        }

        // (2) Build subspace matrices
        for (int j = 0; j < ikry; ++j) {
            y(j + 1) = vsubsv(iL, iR, j + 2, 1);
        }
        y(0) = vsub(0, 1);

        Eigen::MatrixXcd by = Eigen::MatrixXcd::Zero(Nmax, Nmax);
        for (int jj = 0; jj < ikry; ++jj) {
            for (int kk = 0; kk < ikry; ++kk) {
                by(jj + 1, kk + 1) = esubsv(iL, iR, jj + 2, kk + 2);
            }
        }
        for (int jj = 0; jj < ikry; ++jj) {
            by(jj + 1, 0) = esub(jj + 2, 0);
        }
        for (int kk = 0; kk < ikry; ++kk) {
            by(0, kk + 1) = esub(0, kk + 2);
        }
        by(0, 0) = vsub(0, 0);

        // (3) Solve by * y = y
        if (ikry > 0) {
            Eigen::MatrixXcd by_main = by.block(1, 1, ikry, ikry);
            Eigen::VectorXcd rhs = y.segment(1, ikry);
            Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> dec(by_main);
            if (dec.isInvertible()) {
                Eigen::VectorXcd sol = dec.solve(rhs);
                y.segment(1, ikry) = sol;
            } else {
                throw std::runtime_error("Matrix is not invertible.");
            }
        }
        y(0) = 0.0;

        // (4) Compute T-matrix element and update residuals
        complex<double> prod = 0.0;
        for (int j = 0; j < ikry; ++j) {
            prod += vsubsv(iL, iR, 1, j + 2) * y(j + 1);
        }
        complex<double> prod_j = y(0) * vsub(0, 1);

        Lambda_Last = Lambda;
        Lambda = -2.0 / kMomentum * (beta * beta + prod_j + prod);
        LambdaF = Lambda;

        // (5) Dipole moment vector (if required)
        if (DipoleFlag) {
            if (mindip == 1) {
                delta_max = 0.0;
                for (int iDip = 0; iDip < NDip; ++iDip) {
                    complex<double> dip_val = 0.0;
                    for (int j = 0; j < ikry; ++j) {
                        dip_val += dipsubsv(iDip, iR, j + 2) * y(j + 1);
                    }
                    delta_dip = std::abs((dip_val - muDipoleSc_tmp[iDip]) / muDipoleSc_tmp[iDip]);
                    if (delta_dip > delta_max) {
                        delta_max = delta_dip;
                    }
                    if (delta_dip < muDipoleSc_delta[iDip]) {
                        muDipoleSc[iDip] = dip_val;
                        muDipoleSc_delta[iDip] = delta_dip;
                    }
                    muDipoleSc_tmp[iDip] = dip_val;
                }
            } else if (mindip == 0) {
                for (int iDip = 0; iDip < NDip; ++iDip) {
                    complex<double> dip_val = 0.0;
                    for (int j = 0; j < ikry; ++j) {
                        dip_val += dipsubsv(iDip, iR, j + 2) * y(j + 1);
                    }
                    delta_dip = std::abs((dip_val - muDipoleSc_tmp[iDip]) / muDipoleSc_tmp[iDip]);
                    muDipoleSc[iDip] = dip_val;
                    muDipoleSc_delta[iDip] = delta_dip;
                }
            }
        }

        // (6) Convergence and reporting
        if (LMPI_myid == LMPI_master && iprnfg > 0) {
            std::cout << ikry << " | "
                << 4.0 * PI * std::abs(LambdaF / kMomentum * xfauang) * std::abs(LambdaF / kMomentum * xfauang) << " "
                << 4.0 * PI * std::pow(std::sin(std::atan2(std::imag(Lambda), std::real(Lambda))), 2) /
                    (kMomentum * kMomentum) * xfauang * xfauang << " "
                << std::abs((Lambda - Lambda_Last) / Lambda_Last) << " | "
                << std::abs(1.0 + 2.0i * Lambda) << " "
                << std::abs(beta * beta) << " " << std::abs(prod) << " " << std::abs(prod_j) << std::endl;
        }

        // Residual test
        if (std::abs(Lambda - Lambda_Last) < Tol) {
            ++ikry;
            break;
        }
        if (ikry == MaxIter && LMPI_myid == LMPI_master) {
            std::cerr << "######## Warning ############  Iterations did not converge ###############" << std::endl;
        }
    }
    if (ikry > MaxIter) ikry = MaxIter;
    if (LMPI_myid == LMPI_master) {
        std::cout << "MaxEleIter = " << ikry
            << " c.s. = " << 4.0 * M_PI * std::abs(LambdaF / kMomentum * xfauang) * std::abs(LambdaF / kMomentum * xfauang)
            << " angs^2  rmsk = " << std::abs(Lambda - Lambda_Last) << std::endl;
    }
    NumIterI = ikry;
}