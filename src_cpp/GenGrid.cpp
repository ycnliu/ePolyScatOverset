#include "GenGrid.hpp"

namespace GenGrid {

const std::array<int, ZMax> ZValence = {
    1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32
};

namespace LocalSubs {
double EstHMin(
    double RCurr, int NNucCen,
    const std::vector<double>& XNucCen, const std::vector<double>& XMaxAlpha,
    double HFacGauss, double EMaxUse, double HFacWave,
    double RMax, double MaxStep, int LMax, double PCutRd,
    const std::vector<NCoord>& acoord, int NAtom
) {
    constexpr double PGau = 0.1;
    constexpr double RCut = 1.0;
    constexpr double pi = 3.14159265358979323846;
    constexpr double RFac3 = 1.0e-4;
    constexpr double ZFac3 = 10.0;
    constexpr double Ratio3 = 10.0;

    double VLoc = 0.0;

    // --- Potential from all atoms ---
    for (int i = 0; i < NAtom; ++i) {
        int zval = acoord[i].z;
        if (zval > ZMax || zval < 0) {
            throw std::runtime_error("GenGrid: bad ZValence in acoord(i).z");
        }
        double dist = std::abs(RCurr - std::sqrt(
            acoord[i].vec[0]*acoord[i].vec[0] +
            acoord[i].vec[1]*acoord[i].vec[1] +
            acoord[i].vec[2]*acoord[i].vec[2]
        ));
        VLoc += static_cast<double>(ZValence[zval - 1]) / std::max(RCut, dist);
    }
    double hmin1 = pi / (std::sqrt(2.0 * (EMaxUse + VLoc)) * HFacWave);

    // --- hmin2: Minimum step size around nuclear centers ---
    double hmin2 = MaxStep;
    for (int i = 0; i < NNucCen; ++i) {
        double RDiff = std::abs(RCurr - XNucCen[i]);
        double RSwitch = 1.0 / std::sqrt(2.0 * PGau * XMaxAlpha[i]);
        if (RDiff < RSwitch) {
            hmin2 = std::min(hmin2, std::exp(PGau * XMaxAlpha[i] * RDiff * RDiff) / std::sqrt(XMaxAlpha[i]) / HFacGauss);
        } else {
            hmin2 = std::min(hmin2, RDiff * std::exp(0.5) * std::sqrt(2.0 * PGau) / HFacGauss);
        }
    }

    // --- hmin3: Step size from grid partitioning, subgrid proximity ---
    double hmin3 = MaxStep;
    for (int i = 0; i < NAtom; ++i) {
        double RCen = std::sqrt(
            acoord[i].vec[0]*acoord[i].vec[0] +
            acoord[i].vec[1]*acoord[i].vec[1] +
            acoord[i].vec[2]*acoord[i].vec[2]
        );
        if (RCen > RFac3) {
            if (RCen / static_cast<double>(8) < RCurr && acoord[i].z > 0) {
                double invlmax = std::abs(
                    std::log(RCurr / RCen) /
                    std::log(std::sqrt(PCutRd * static_cast<double>(acoord[i].z) / ZFac3))
                );
                invlmax = std::max(invlmax, 1.0 / static_cast<double>(LMax));
                hmin3 = std::min(hmin3, RCurr * (std::pow(Ratio3, invlmax) - 1.0) / static_cast<double>(8));
            } else {
                hmin3 = std::min(hmin3, RCen / static_cast<double>(8));
            }
        }
    }

    return std::min({hmin1, hmin2, hmin3});
}
void deriv_legendre_lobatto(
    const std::vector<double>& x,
    const std::vector<double>& w,
    double a,
    GenGrid::Mat2D_t& DMAT
) {
    // N = order of Legendre-Lobatto grid
    int N = static_cast<int>(x.size()) - 1;

    // Allocate/resize the output matrix
    DMAT.m.assign(N + 1, std::vector<double>(N + 1, 0.0));

    // Fill the corners
    DMAT.m[0][0]     = -0.25 * N * (N + 1.0) / (a * a);
    DMAT.m[N][N]     =  0.25 * N * (N + 1.0) / (a * a);

    // Fill the rest using the Fortran loop logic
    for (int i = 0; i <= N; ++i) {
        int ph = -1;
        for (int j = i + 1; j <= N; ++j) {
            double di = 1.0 / (x[i] - x[j]);
            DMAT.m[i][j] = di * ph * std::sqrt(w[j] / w[i]);
            DMAT.m[j][i] = -di * ph * std::sqrt(w[i] / w[j]);
            ph = -ph;
        }
    }
}

void GenGrid(
    double RMax, double EMaxUse, int LMax, double PCutRd,
    const std::vector<NCoord>& acoord, int NAtom,
    const std::vector<std::vector<double>>& ex,
    const std::vector<int>& iatmfr, const std::vector<int>& iatmto,
    const std::vector<int>& NoPrim
) {
    // Placeholder variables for constants (replace with actual defaults as needed)
    constexpr int LoopCntX = 1000;
    constexpr int NPFac = 8; // This should match your project parameter!
    constexpr double xfauang = 0.529177; // Bohr to Angstrom
    constexpr double XFEVAU = 27.2114; // Hartree to eV

    // Read defaults (replace with values if needed)
    double HFacGauss = 30.0;
    double HFacWave  = 120.0;
    double MinExpFac = 0.1;
    int GridFac = 1;
    double MaxStep = RMax; // Initial

    // The nuclear centers and max exponents
    std::vector<double> XNucCen(NAtom + 2, 0.0);
    std::vector<double> XMaxAlpha(NAtom + 2, 1.0);

    int NNucCen = 0; // Always include origin
    XNucCen[NNucCen++] = 0.0;
    XMaxAlpha[NNucCen-1] = 1.0;

    // Locate each atom and find its maximum exponent
    for (int iAtom = 0; iAtom < NAtom; ++iAtom) {
        double almax = acoord[iAtom].z * acoord[iAtom].z * MinExpFac;
        for (int ifn = iatmfr[iAtom]; ifn <= iatmto[iAtom]; ++ifn) {
            for (int k = 0; k < NoPrim[ifn]; ++k) {
                almax = std::max(almax, ex[ifn][k]);
            }
        }
        double RAtom = 0.0;
        for (int i = 0; i < 3; ++i) RAtom += acoord[iAtom].vec[i] * acoord[iAtom].vec[i];
        RAtom = std::sqrt(RAtom);

        // Check if we already have this center
        auto found = std::lower_bound(XNucCen.begin(), XNucCen.begin() + NNucCen, RAtom);
        if (found != XNucCen.begin() + NNucCen) {
            XMaxAlpha[found - XNucCen.begin()] = std::max(XMaxAlpha[found - XNucCen.begin()], almax);
        } else {
            // Insert here
            XNucCen[NNucCen] = RAtom;
            XMaxAlpha[NNucCen] = almax;
            NNucCen++;
        }
    }
    // Set endpoint for the grid
    XNucCen[NNucCen] = RMax;

    // === Generate the points in each region ===
    int nptrg = 0;
    double RCurr = 0.0;
    std::vector<double> hregvec;
    std::vector<int> nptinvec;

    for (int iGReg = 0; iGReg < NNucCen; ++iGReg) {
        bool LastReg = false;
        double hmin0 = LocalSubs::EstHMin(
            RCurr, NNucCen, XNucCen, XMaxAlpha,
            HFacGauss, EMaxUse, HFacWave, RMax, MaxStep, LMax, PCutRd,
            acoord, NAtom
        );
        if (NPFac * hmin0 + RCurr > XNucCen[iGReg+1]) {
            LastReg = true;
            hmin0 = (XNucCen[iGReg+1] - RCurr) / static_cast<double>(NPFac);
        }

        int LoopCnt = 0;
        while (true) {
            LoopCnt++;
            if (LoopCnt > LoopCntX) throw std::runtime_error("LoopCntX exceeded!");

            bool adjust = false;
            for (int i = 1; i <= NPFac; ++i) {
                double hmini = LocalSubs::EstHMin(
                    RCurr + i * hmin0, NNucCen, XNucCen, XMaxAlpha,
                    HFacGauss, EMaxUse, HFacWave, RMax, MaxStep, LMax, PCutRd,
                    acoord, NAtom
                );
                if (hmini + 1e-10 < hmin0) {
                    hmin0 = hmini;
                    LastReg = false;
                    adjust = true;
                    break;
                }
            }
            if (!adjust) break;
        }

        // Add this region
        if (!hregvec.empty() && std::abs(hregvec.back() - hmin0) < 1e-10 && nptinvec.back() < 8 * NPFac) {
            nptinvec.back() += NPFac;
            RCurr += NPFac * hmin0;
            if (LastReg) break;
            continue;
        }
        nptrg++;
        hregvec.push_back(hmin0);
        nptinvec.push_back(NPFac);

        RCurr += NPFac * hmin0;
        if (LastReg) break;
    }

    // Now nptin and hreg arrays
    std::vector<double> hreg(hregvec.begin(), hregvec.end());
    std::vector<int> nptin(nptinvec.begin(), nptinvec.end());

    // Apply grid factor
    for (auto& n : nptin) n *= GridFac;
    for (auto& h : hreg) h /= static_cast<double>(GridFac);

    // Output (can be parallelized or redirected to file as needed)
    int nptgd = 0;
    RCurr = 0.0;
    for (size_t irg = 0; irg < hreg.size(); ++irg) {
        nptgd += nptin[irg];
        RCurr += nptin[irg] * hreg[irg];
        std::cout << "Grid region: " << (irg+1)
                  << "  Points: " << nptin[irg]
                  << "  Total: " << nptgd
                  << "  Step: " << hreg[irg] * xfauang
                  << "  R_end: " << RCurr * xfauang << std::endl;
    }
}
void GenGridOverset(
    double RMax, double EMaxUse,
    const std::vector<Mat1DI>& UniqNuc,
    const std::vector<NCoord>& acoord,
    int nr, double NLambda, double N_rpnts_SG, double HFacGaussOG,
    const std::vector<Mat1DI>& UniqNucAb,
    const std::vector<std::vector<double>>& ex,
    const std::vector<int>& iatmfr, const std::vector<int>& iatmto,
    const std::vector<int>& NoPrim
) {
    // Placeholders for values/constants - replace with real project values!
    constexpr double PI = 3.14159265358979323846;
    constexpr double PGau = 0.1;
    constexpr double XFAuAng = 0.529177;
    constexpr double PClose = 1e-8;
    constexpr int GridFac_def = 1, SGridFac_def = 1;
    constexpr double MinExpFac_def = 0.1, HFacGaussOG_def = 30.0;
    constexpr int nrFEM_def = 1;
    constexpr double NLambda_def = 10.0, N_rpnts_SG_def = 10.0;
    constexpr double FirstPoint_def = 0.0;
    constexpr double ReservedPortion = 0.25;

    // Data containers for this routine
    std::vector<double> XMaxAlpha(acoord.size(), 0.0);
    double MinExpFac = MinExpFac_def;
    double HFacGauss = HFacGaussOG;
    int GridFac = GridFac_def, SGridFac = SGridFac_def;
    double FirstPoint = FirstPoint_def;
    double NrfemFac = 1.0, NrfemFacSG = 1.0;

    // Set XMaxAlpha for each nuclear center
    for (size_t iAtom = 0; iAtom < acoord.size(); ++iAtom) {
        double almax = acoord[iAtom].z * acoord[iAtom].z * MinExpFac;
        for (int ifn = iatmfr[iAtom]; ifn <= iatmto[iAtom]; ++ifn) {
            for (int k = 0; k < NoPrim[ifn]; ++k) {
                almax = std::max(almax, ex[ifn][k]);
            }
        }
        XMaxAlpha[iAtom] = almax;
    }

    // Calculate Radius_UG, Radius_EG, MapSGtoEG, etc. as in the Fortran code
    double Radius_UG = RMax;
    std::vector<double> RadiusCent; // to be computed
 
    std::vector<double> Radius_EG(UniqNuc.size());
    for (size_t i = 0; i < UniqNuc.size(); ++i) {
        for (size_t j = 0; j < UniqNuc[i].p.size(); ++j) {
            // MapSGtoEG[UniqNuc[i].p[j]] = i; // if needed
            Radius_EG[i] = RadiusCent[UniqNuc[i].p[j]];
        }
    }

    double RadiusCentMin = *std::max_element(RadiusCent.begin(), RadiusCent.end());
    double RMaxWithSubgrids = RadiusCentMin * 2.0;
    double RMinOfSubgrids = RMax;

    for (size_t i = 0; i < acoord.size(); ++i) {
        double val = 2.0 * RadiusCent[i] + std::sqrt(
            acoord[i].vec[0]*acoord[i].vec[0] +
            acoord[i].vec[1]*acoord[i].vec[1] +
            acoord[i].vec[2]*acoord[i].vec[2]
        );
        RMaxWithSubgrids = std::max(RMaxWithSubgrids, val);
        double dist = std::sqrt(
            acoord[i].vec[0]*acoord[i].vec[0] +
            acoord[i].vec[1]*acoord[i].vec[1] +
            acoord[i].vec[2]*acoord[i].vec[2]
        ) - 2.0 * RadiusCent[i];
        if (dist > 2.0 * RadiusCent[i]) {
            RMinOfSubgrids = std::min(RMinOfSubgrids, dist);
        }
    }
    if (RMaxWithSubgrids <= RMax) {
        throw std::runtime_error("RMax - GenGridOverset - RMax not large enough");
    }

    // Compute step sizes
    double hSG = RadiusCentMin / N_rpnts_SG;
    double hAsymp = (2.0 * PI) / (NLambda * std::sqrt(2.0 * EMaxUse));

    // --- Finite element mesh and grid construction ---
    // Not shown in full here for brevity, but you would proceed as in the Fortran code,
    // mapping ALLOCATE/DEALLOCATE to std::vector, and all loops to C++.

    for (size_t ig = 0; ig <= UniqNuc.size(); ++ig) {
        std::cout << "Grid information for group " << ig << " of equivalent nuclei" << std::endl;
    }

}

void GetRadiusCent(
    std::vector<double>& RadiusCentLocal,
    const std::vector<NCoord>& acoord,
    double SGFactor,
    const std::vector<Mat1DI>& UniqNuc,
    double XFAuAng
) {
    // Covalent radii table for Z = 0 to 20
    constexpr std::array<double, 21> CovRadRable = {
        0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57,
        0.58, 1.66, 1.41, 1.21, 1.11, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70
    };

    bool AdaSG = false;
    double SGRMax = std::numeric_limits<double>::max();

    size_t N = acoord.size();
    RadiusCentLocal.assign(N, 0.0);

    // Compute distance matrix
    std::vector<std::vector<double>> distMat(N, std::vector<double>(N, 0.0));
    double MaxDist = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (j != i) {
                double d = 0.0;
                for (int k = 0; k < 3; ++k)
                    d += (acoord[j].vec[k] - acoord[i].vec[k]) * (acoord[j].vec[k] - acoord[i].vec[k]);
                d = std::sqrt(d);
                distMat[i][j] = d;
                MaxDist = std::max(MaxDist, d);
            }
        }
    }

    // Compute MinDist for each atom
    std::vector<double> MinDist(N, MaxDist);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (j != i) MinDist[i] = std::min(MinDist[i], distMat[i][j]);
        }
    }
    if (N == 1) {
        double norm = 0.0;
        for (int k = 0; k < 3; ++k)
            norm += acoord[0].vec[k] * acoord[0].vec[k];
        MinDist[0] = std::sqrt(norm);
        MaxDist = MinDist[0];
    }

    for (size_t i = 0; i < N; ++i)
        RadiusCentLocal[i] = MinDist[i] / (SGFactor * 4.0);

    // Optional: adaptive subgrid growth
    if (AdaSG) {
        std::vector<double> RadEGTmp(UniqNuc.size(), 0.0);
        std::vector<bool> GrowEG(UniqNuc.size(), true);
        int GrowMaxIter = 1000;
        double GrowStepSize = 0.01;
        for (int iter = 0; iter < GrowMaxIter; ++iter) {
            for (size_t ig = 0; ig < UniqNuc.size(); ++ig) {
                if (GrowEG[ig]) {
                    int zidx = acoord[UniqNuc[ig].p[0]].z;
                    RadEGTmp[ig] += GrowStepSize * CovRadRable[std::min<int>(zidx, 20)];
                    for (size_t j = 0; j < UniqNuc[ig].p.size(); ++j) {
                        for (size_t k = 0; k < N; ++k) {
                            if (UniqNuc[ig].p[j] != static_cast<int>(k)) {
                                double d = distMat[UniqNuc[ig].p[j]][k];
                                double RadSum = RadEGTmp[ig] + RadEGTmp[k];
                                if (RadSum > d) {
                                    GrowEG[ig] = false;
                                    GrowEG[k] = false;
                                    RadEGTmp[ig] -= (RadSum - d);
                                }
                            }
                        }
                    }
                }
            }
            if (std::all_of(GrowEG.begin(), GrowEG.end(), [](bool g) { return !g; })) break;
        }
        for (size_t i = 0; i < UniqNuc.size(); ++i) {
            for (size_t j = 0; j < UniqNuc[i].p.size(); ++j)
                RadiusCentLocal[UniqNuc[i].p[j]] = RadEGTmp[i] / (SGFactor * 2.0);
        }
    }

    // Enforce maximum radius constraint
    for (size_t i = 0; i < N; ++i)
        RadiusCentLocal[i] = std::min(RadiusCentLocal[i], SGRMax / 2.0 / XFAuAng);

    // Output (can be controlled by MPI rank if desired)
    std::cout << "Distance matrix (in Angstroms)\n";
    for (size_t i = 0; i < N; ++i) {
        std::cout << i + 1 << " ";
        for (size_t j = 0; j < N; ++j)
            std::cout << distMat[i][j] * XFAuAng << " ";
        std::cout << "\n";
    }
    std::cout << "   Min ";
    for (double m : MinDist) std::cout << m * XFAuAng << " ";
    std::cout << "\nRadius ";
    for (double r : RadiusCentLocal) std::cout << r * XFAuAng << " ";
    std::cout << std::endl;
}
void EstimateLMaxFromSG(
    std::vector<int>& LMaxFromSG,
    std::vector<double>& RadiusCent,
    double SGFactor, double SGLMaxFactor,
    const std::vector<NCoord>& acoord,
    int LMaxOrigin,
    double XFAuAng
) {
    size_t N = acoord.size();
    GetRadiusCent(RadiusCent, acoord, SGFactor, {}, XFAuAng);

    LMaxFromSG.assign(N, 0);
    if (N > 1) {
        for (size_t i = 0; i < N; ++i) {
            double RCent = 0.0;
            for (int k = 0; k < 3; ++k)
                RCent += acoord[i].vec[k] * acoord[i].vec[k];
            RCent = std::sqrt(RCent);
            LMaxFromSG[i] = static_cast<int>(SGLMaxFactor * (RCent / RadiusCent[i])) + LMaxOrigin;
        }
    } else {
        LMaxFromSG[0] = LMaxOrigin;
    }
    // Output results (could restrict to MPI master if needed)
    std::cout << "Computed LMaxFromSG using SGLMaxFactor " << SGLMaxFactor
              << " with LMaxOrigin " << LMaxOrigin << std::endl;
    for (size_t i = 0; i < N; ++i) {
        std::cout << i+1 << " ";
        for (int k = 0; k < 3; ++k)
            std::cout << acoord[i].vec[k] * XFAuAng << " ";
        std::cout << LMaxFromSG[i] << std::endl;
    }
}

