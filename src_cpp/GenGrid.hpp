#ifndef GENGRID_HPP
#define GENGRID_HPP

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace GenGrid {

constexpr int ZMax = 86;
extern const std::array<int, ZMax> ZValence;

// Utility structs (these may need extension as implementation proceeds)
struct NCoord {
    std::array<double, 3> vec;
    int z;
};
struct Mat1DI {
    std::vector<int> p;
};
struct ChVecILL {
    std::vector<std::string> v;
};
struct Mat2D_t {
    std::vector<std::vector<double>> m;
};
struct RVec {
    std::vector<double> v;
};
struct IVec {
    std::vector<int> v;
};

// ===== Module GenGrid_localsubs =====

namespace LocalSubs {
    double EstHMin(
        double RCurr, int NNucCen,
        const std::vector<double>& XNucCen, const std::vector<double>& XMaxAlpha,
        double HFacGauss, double EMaxUse, double HFacWave,
        double RMax, double MaxStep, int LMax, double PCutRd,
        const std::vector<NCoord>& acoord, int NAtom
    );

    void deriv_legendre_lobatto(
        const std::vector<double>& x,
        const std::vector<double>& w,
        double a,
        Mat2D_t& DMAT
    );
}

// ===== Module GenGrid_sub =====

void GenGrid(
    double RMax, double EMaxUse, int LMax, double PCutRd,
    const std::vector<NCoord>& acoord, int NAtom,
    const std::vector<std::vector<double>>& ex, // For basis
    const std::vector<int>& iatmfr, const std::vector<int>& iatmto,
    const std::vector<int>& NoPrim
);

void GenGridOverset(
    double RMax, double EMaxUse,
    const std::vector<Mat1DI>& UniqNuc,
    const std::vector<NCoord>& acoord,
    int nr, double NLambda, double N_rpnts_SG, double HFacGaussOG,
    const std::vector<Mat1DI>& UniqNucAb,
    // Global variables for OversetGrid_cmn etc. will be in implementation
    const std::vector<std::vector<double>>& ex, // For basis
    const std::vector<int>& iatmfr, const std::vector<int>& iatmto,
    const std::vector<int>& NoPrim
);

void GetRadiusCent(
    std::vector<double>& RadiusCentLocal,
    const std::vector<NCoord>& acoord,
    double SGFactor,
    const std::vector<Mat1DI>& UniqNuc,
    // MapSGtoEG will be managed as an internal static/global in .cpp
    double XFAuAng
);

void EstimateLMaxFromSG(
    std::vector<int>& LMaxFromSG,
    std::vector<double>& RadiusCent,
    double SGFactor, double SGLMaxFactor,
    const std::vector<NCoord>& acoord,
    int LMaxOrigin,
    double XFAuAng
);

} // namespace GenGrid

#endif // GENGRID_HPP
