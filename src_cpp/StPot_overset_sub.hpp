#ifndef STPOT_OVERSET_SUB_HPP
#define STPOT_OVERSET_SUB_HPP

#include <vector>
#include <array>
#include <string>
#include <complex>
#include <memory>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace OversetGrid {

    // === Data structures mapping Fortran types ===
    struct Cvec {
        double x, y, z;
    };

    struct ChVecILL {
        std::vector<std::string> v;
    };

    struct Overset_CD_t {
        // Maps to Fortran: type(Overset_CD_t)
        // Add actual members as needed
        // Example:
        struct UGridPoint {
            std::vector<double> m;
        };
        std::vector<UGridPoint> UG;    // Ubergrid
        std::vector<UGridPoint> SG;    // Subgrids
    };

    struct Overset_PW_t {
        struct UGridPoint {
            std::vector<double> m;
        };
        struct EGridPoint {
            std::vector<std::vector<double>> m;
        };
        std::vector<UGridPoint> UG; // Ubergrid
        std::vector<EGridPoint> EG; // Equivalent groups
    };

    // ========== Main Routines ==========
    // Calculate static potential
    void StPot_overset(
        /* Provide all the output/optional matrices as std::vector<std::vector<double>>& or similar,
           optional arguments can be overloaded or defaulted,
           use smart pointers or std::optional for true optional support if desired */
    );

    // Add Fege potential for local exchange calculation
    void FegeCorrect_overset(
        double engv,
        double fegescale = 1.0
    );

    // Check for continuity of the wavefunction
    void ContinuityCheck_CD(
        const Overset_CD_t& Orb_CD
    );

    // Calculate AO basis set in coordinate representation
    void CalcAO(
        const std::vector<Cvec>& XYZ,
        std::vector<std::vector<double>>& valao
    );

} // namespace OversetGrid

#endif // STPOT_OVERSET_SUB_HPP
