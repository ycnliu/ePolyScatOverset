// ScatStabOversetSub.hpp

#pragma once

#include "OversetGridTypes.hpp"
#include "InterpolationOverset.hpp"
#include "StPotOversetSub.hpp"
#include <vector>
#include <complex>
#include <string>
#include <optional>

/**
 * @brief Class encapsulating all data and routines for overset grid complex Kohn scattering calculations.
 *        This is the modern C++ equivalent of the Fortran module ScatStab_overset_sub.
 */
class ScatStabOversetSub {
public:
    using real = double;
    using complex = std::complex<double>;

    // --- Module data (shared state, originally module variables) ---

    // Interpolating matrices for overset grid
    InterpolatingMatrices IMat;
    std::vector<InterpolatingMatrices> IMatAb;

    // Iterated dipole operator in coordinate representation
    // 3D: [Gauge][Iteration][VectorType]
    std::vector<std::vector<std::vector<Overset_CD_t>>> Iterated_DipoleOp_CD;

    // Arnoldi iterates: [PartialWave][Iteration][VectorType]
    std::vector<std::vector<std::vector<Overset_CD_t>>> ArnoldiIterates;

    // Which Arnoldi iterate each partial wave is on
    std::vector<int> IterateK;

    // Subspace matrices for Arnoldi (G0*V)
    std::vector<std::vector<std::vector<complex>>> hsv;

    // Subspace matrices for V and E (vsubsv, esubsv)
    std::vector<std::vector<std::vector<std::vector<complex>>>> vsubsv;
    std::vector<std::vector<std::vector<std::vector<complex>>>> esubsv;

    // Subspace matrix for the dipole operator
    std::vector<std::vector<std::vector<complex>>> dipsubsv;

    // Test output flag
    bool FoundTestOut = false;

    // Orthogonalization/normalization intermediates
    complex OrthogComp = 0.0, OrthogComp2 = 0.0, OrthogComp3 = 0.0;
    complex NormFactor = 1.0, NormFactor2 = 1.0;

    // Parameters
    static constexpr bool DoOccOrbitalOrthogonalization = true;
    static constexpr bool UseOrthogEffectivePot = true;
    static constexpr bool UseHartreeEnergies = false;
    static constexpr const char* DefaultOrientationLabel = "1";

    // Orthogonal subspace matrices
    std::vector<std::vector<std::vector<std::vector<complex>>>> vsubsv_orthog;
    std::vector<std::vector<std::vector<std::vector<complex>>>> esubsv_orthog;

    // Orbital vectors and effective potentials
    std::vector<Overset_CD_t> VxOrbitalVec, HminusE_x_OrbitalVec;

    // Ubergrid (central grid) basis and overlap matrix
    std::vector<Overset_PW_t> UG_Basis;
    std::vector<std::vector<real>> UGUG_Overlap;

    // Algorithmic parameters
    int ItAlg = 0;
    int Metric = 0;

    ScatStabOversetSub() = default;
    ~ScatStabOversetSub() = default;


    /**
     * @brief Main scattering calculation routine (Complex Kohn, overset grid).
     */
    void ScatStab_overset(
        double engv,
        bool DipoleFlag,
        int iukmat,
        const std::string& FileNameMatrixElements,
        const std::optional<std::vector<std::vector<double>>>& KineticMatOccOrbitals = std::nullopt,
        const std::optional<std::vector<std::vector<double>>>& NucAtractionMatOccOrbitals = std::nullopt
    );

    /**
     * @brief Born-Arnoldi process for T-matrix calculation.
     */
    void ArnoldiCKohnGVNoJAsymSv(
        complex& LambdaF,
        std::vector<complex>& muDipoleSc,
        const Overset_CD_t& POT_CD,
        int iL, int iR,
        const Overset_CD_t& JL_CD,
        const Overset_CD_t& JLV_CD,
        const Overset_CD_t& JR_CD,
        const Overset_CD_t& JRV_CD,
        const Overset_CD_t& HOUT_CD,
        const Overset_CD_t& LHOUT_CD,
        complex uLuBorn,
        const std::vector<complex>& uLuOutgoing,
        int MaxIter,
        double Tol,
        double UTol,
        int itype,
        int& NumIterI,
        bool DipoleFlag = false,
        const std::vector<std::vector<double>>* PotentialMatOccOrbitals = nullptr,
        const std::vector<std::vector<double>>* EminusHMatOccOrbitals = nullptr,
        double energy = 0.0
    );

    /**
     * @brief Add a Born-Arnoldi iterate to the partial wave i.
     */
    void AddIterate(
        const Overset_CD_t& POT_CD,
        int i,
        int ikry,
        int itype,
        double& vecnorm,
        bool DipoleFlag = false,
        const std::vector<std::vector<double>>* PotentialMatOccOrbitals = nullptr,
        const std::vector<std::vector<double>>* EminusHMatOccOrbitals = nullptr,
        double energy = 0.0
    );

    /**
     * @brief Build J and JV on the overset grid.
     */
    void BuildJV(
        Overset_CD_t& J_CD,
        Overset_CD_t& JV_CD,
        int LHVAL,
        const Overset_CD_t& POT_CD,
        int itype
    );

    /**
     * @brief Build the j_l(kr_ug)Y_lm(w_ug) function on equivalent grids.
     */
    void BuildJEG(
        Overset_CD_t& JEG_CD,
        int LHVAL,
        int itype
    );

    /**
     * @brief Match a short-range T-Matrix to the Coulomb T-Matrix.
     */
    void CoulombMatching(
        std::vector<std::vector<complex>>& TMatrix,
        const std::vector<std::vector<complex>>& TMatrixCutoff,
        double RMatch,
        int LMaxC,
        int itype,
        bool DipoleFlag = false
    );



};
