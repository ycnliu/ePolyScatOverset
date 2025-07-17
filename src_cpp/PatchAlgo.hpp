// PatchAlgo.hpp
#pragma once
#include <vector>
#include <complex>
#include <string>
#include "OversetTypes.hpp"
#include "InterpolatingMatrices.hpp"

// Forward declarations for required Fortran-derived types (replace with your C++ definitions)
struct Overset_CD_t;
struct Overset_PW_t;
struct InterpolatingMatrices;
struct AngExp;
struct Grid_t;
struct ChVecILL;

namespace PatchAlgo {

// Utility output tags
constexpr int BEFORE_TAG = 0;
constexpr int DIFF_TAG   = 1;
constexpr int AFTER_TAG  = 2;
constexpr int SUBGRID_TAG = 3;
constexpr int ERROR_TAG  = 4;

// Default plotting targets
constexpr double CosThetaTargetDefault = 0.8165;
constexpr double PhiTargetDefault = 0.0;
constexpr double PhiTargetDefault2 = 3.1415926535;
constexpr double rTargetDefault = 1.27802;

void TranCDtoPW_alt(
    Overset_PW_t& Orb_PW_out,
    const Overset_CD_t& Orb_CD_in,
    int itype,
    const InterpolatingMatrices& IMatInput,
    bool useabel = false,
    bool dump_tmp_cd = false
);

// Real-valued Lucchese partitioning for the Partial Wave
void LucchesePartitionEG_Re(
    std::vector<std::vector<double>>& PartialWave, // [AExpX.nlht][NumRadPt]
    const AngExp& AExpX,
    int NumRadPt,
    int EG_idx,
    int iType,
    bool useabel = false
);

// Interpolate EG on UG using patching matrices and subtract
void EvalEGOnUG_PW(
    Overset_CD_t& Phi_CD,
    const Overset_PW_t& Psi_PW,
    const InterpolatingMatrices& IM,
    bool useabel = false
);

// Dump coordinate representation of UG (debug output)
void Dump_CD_UG(const Overset_CD_t& Orb_CD_in, int datlabel = 0);
void Dump_CD_UG_2D(const Overset_CD_t& Orb_CD_in, int datlabel = 0);
void print_UG_grid(int datlabel = 0);

} // namespace PatchAlgo
