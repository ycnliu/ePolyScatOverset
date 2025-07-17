
// ApplyG0.cu
// Factored CUDA implementation of ApplyG0, ApplyG0UG, and ApplyG0EG in one file

#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

using Real = double;
using Complex = thrust::complex<double>;

// ---------------------------
// Utility Indexing Functions
// ---------------------------
#define IDX3D(ieg, ilh, irg, wmax, lp) (((ieg) * (wmax) * (lp)) + ((ilh) * (lp)) + (irg))

// ---------------------------
// Data Structures
// ---------------------------
namespace G0Overset {

struct GridParams {
    int N;        // Number of EG grids
    int wmax;     // Number of partial waves per grid
    int lp;       // Number of radial points per grid
};

struct GridUGParams {
    int N;        // Number of UG points
    int wmax;     // Number of partial waves per UG point
};

struct OperatorArraysEG {
    const Real* F;    // F_EG[N][wmax][lp]
    const Real* G;    // G_EG[N][wmax][lp]
    const Real* W;    // W_EG[N][wmax]
    const Real* Psi;  // Psi_EG[N][wmax][lp]
    const Real* grid_w; // grid weights per [N][lp]
    Real* Phi;        // Output: Phi_EG[N][wmax][lp]
};

struct OperatorArraysUG {
    const Real* F;
    const Real* G;
    const Real* W;
    const Real* Psi;
    const Real* grid_w;
    Real* Phi;
};

// ---------------------------
// CUDA Kernels
// ---------------------------

// EG kernel: full parallelization in (ieg, ilh, irg)
/**
 * @brief CUDA kernel to apply G0 operation on EG grid
 *
 * This kernel performs the G0 operation on an EG grid by
 * calculating contributions from F and G arrays and storing 
 * the result in the Phi array.
 *
 * @param N The number of grids.
 * @param wmax The maximum number of angular components per grid.
 * @param lp The number of radial points per grid.
 * @param F Pointer to the F array.
 * @param G Pointer to the G array.
 * @param W Pointer to the W array.
 * @param Psi Pointer to the Psi array.
 * @param grid_w Pointer to the grid weights array.
 * @param Phi Pointer to the output Phi array.
 */
__global__
void ApplyG0EG_kernel(int N, int wmax, int lp,
                      const Real* __restrict__ F,
                      const Real* __restrict__ G,
                      const Real* __restrict__ W,
                      const Real* __restrict__ Psi,
                      const Real* __restrict__ grid_w,
                      Real* __restrict__ Phi)
{
    // Calculate global thread indices for each dimension
    int ieg = blockIdx.z * blockDim.z + threadIdx.z;
    int ilh = blockIdx.y * blockDim.y + threadIdx.y;
    int irg = blockIdx.x * blockDim.x + threadIdx.x;

    // Ensure indices are within bounds
    if (ieg < N && ilh < wmax && irg < lp) {
        Real WF = 0.0; // Accumulator for contributions from F
        Real WG = 0.0; // Accumulator for contributions from G

        // Accumulate contributions over all radial points
        for (int jrg = 0; jrg < lp; ++jrg) {
            WF += F[IDX3D(ieg, ilh, jrg, wmax, lp)] * Psi[IDX3D(ieg, ilh, jrg, wmax, lp)];
            WG += G[IDX3D(ieg, ilh, jrg, wmax, lp)] * Psi[IDX3D(ieg, ilh, jrg, wmax, lp)];
        }

        // Calculate the wronskian factor for the current radial point
        Real wf = sqrt(grid_w[ieg * lp + irg]);

        // Retrieve the weight value from W
        Real wval = W[ieg * wmax + ilh];

        // Compute the result for the current grid point and store in Phi
        Phi[IDX3D(ieg, ilh, irg, wmax, lp)] = wf * wval * (WF * G[IDX3D(ieg, ilh, irg, wmax, lp)] + WG * F[IDX3D(ieg, ilh, irg, wmax, lp)]);
    }
}

// UG kernel: full parallelization in (irg, ilh)
/**
 * @brief CUDA kernel to apply G0 operation on UG grid
 *
 * This kernel performs the G0 operation on an UG grid by
 * calculating the sum of contributions from F and G arrays
 * and storing the result in the Phi array.
 *
 * @param N The number of radial grid points.
 * @param wmax The maximum number of angular components.
 * @param F Pointer to the F array.
 * @param G Pointer to the G array.
 * @param W Pointer to the W array.
 * @param Psi Pointer to the Psi array.
 * @param grid_w Pointer to the grid_w array.
 * @param Phi Pointer to the output Phi array.
 */
__global__
void ApplyG0UG_kernel(int N, int wmax,
                      const Real* __restrict__ F,
                      const Real* __restrict__ G,
                      const Real* __restrict__ W,
                      const Real* __restrict__ Psi,
                      const Real* __restrict__ grid_w,
                      Real* __restrict__ Phi)
{
    int irg = blockIdx.x * blockDim.x + threadIdx.x;
    int ilh = blockIdx.y * blockDim.y + threadIdx.y;

    if (irg < N && ilh < wmax) {
        Real sum = 0;

        // Accumulate contributions from F before current radial point
        for (int jrg = 0; jrg <= irg; ++jrg)
            sum += F[jrg * wmax + ilh] * Psi[jrg * wmax + ilh];

        // Accumulate contributions from G after current radial point
        for (int jrg = irg; jrg < N; ++jrg)
            sum += G[jrg * wmax + ilh] * Psi[jrg * wmax + ilh];

        // Apply weights and store the result in Phi
        Phi[irg * wmax + ilh] = sqrt(grid_w[irg]) * W[ilh] * sum;
    }
}

// ---------------------------
// Host-side Wrappers
// ---------------------------

// EG grid host-side wrapper
void ApplyG0EG_GPU(const GridParams& params, const OperatorArraysEG& arrays)
{
    int x = std::min(256, params.lp);
    int y = std::min(16, params.wmax);
    dim3 block(x, y);
    dim3 grid((params.lp + block.x - 1) / block.x,
              (params.wmax + block.y - 1) / block.y,
              params.N);

    ApplyG0EG_kernel<<<grid, block, 0, 0>>>(
        params.N, params.wmax, params.lp,
        arrays.F, arrays.G, arrays.W, arrays.Psi, arrays.grid_w, arrays.Phi
    );
    cudaDeviceSynchronize();
}

// UG grid host-side wrapper
void ApplyG0UG_GPU(const GridUGParams& params, const OperatorArraysUG& arrays)
{
    int x = std::min(1024, params.N);
    int y = std::min(16, params.wmax);
    dim3 block(x, y);
    dim3 grid((params.N + block.x - 1) / block.x,
              (params.wmax + block.y - 1) / block.y);

    ApplyG0UG_kernel<<<grid, block>>>(
        params.N, params.wmax,
        arrays.F, arrays.G, arrays.W, arrays.Psi, arrays.grid_w, arrays.Phi
    );
    cudaDeviceSynchronize();
}

// Main dispatcher: ApplyG0 (applies either or both, then could combine results)
void ApplyG0(const GridParams& eg_params, const OperatorArraysEG& eg_arrays,
             const GridUGParams& ug_params, const OperatorArraysUG& ug_arrays,
             bool use_eg, bool use_ug)
{
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    if (use_eg)
        ApplyG0EG_GPU(eg_params, eg_arrays, stream);
    if (use_ug)
        ApplyG0UG_GPU(ug_params, ug_arrays, stream);

    cudaStreamSynchronize(stream);
    cudaStreamDestroy(stream);
}

// EG grid host-side wrapper
void ApplyG0EG_GPU(const GridParams& params, const OperatorArraysEG& arrays, cudaStream_t stream)
{
    dim3 block(16, 8);
    dim3 grid((params.lp + block.x - 1) / block.x,
              (params.wmax + block.y - 1) / block.y,
              params.N);

    ApplyG0EG_kernel<<<grid, block, 0, stream>>>(
        params.N, params.wmax, params.lp,
        arrays.F, arrays.G, arrays.W, arrays.Psi, arrays.grid_w, arrays.Phi
    );
}

// UG grid host-side wrapper
void ApplyG0UG_GPU(const GridUGParams& params, const OperatorArraysUG& arrays, cudaStream_t stream)
{
    dim3 block(32, 8);
    dim3 grid((params.N + block.x - 1) / block.x,
              (params.wmax + block.y - 1) / block.y);

    ApplyG0UG_kernel<<<grid, block, 0, stream>>>(
        params.N, params.wmax,
        arrays.F, arrays.G, arrays.W, arrays.Psi, arrays.grid_w, arrays.Phi
    );
}

} // namespace G0Overset

// ---------------------------
// Example Main/Test Harness
// ---------------------------
void test_comm()
{
    const int N = 2;        // EG grids
    const int wmax = 4;     // Partial waves
    const int lp = 8;       // Radial points
    const int N_UG = 8;     // Ubergrid points

    std::vector<Real> F_EG(N * wmax * lp), G_EG(N * wmax * lp), W_EG(N * wmax), Psi_EG(N * wmax * lp), grid_w(N * lp), Phi_EG(N * wmax * lp);
    std::vector<Real> F_UG(N_UG * wmax), G_UG(N_UG * wmax), W_UG(wmax), Psi_UG(N_UG * wmax), grid_w_UG(N_UG), Phi_UG(N_UG * wmax);

    Real *d_F_EG, *d_G_EG, *d_W_EG, *d_Psi_EG, *d_grid_w, *d_Phi_EG;
    cudaMalloc(&d_F_EG, F_EG.size() * sizeof(Real));
    cudaMalloc(&d_G_EG, G_EG.size() * sizeof(Real));
    cudaMalloc(&d_W_EG, W_EG.size() * sizeof(Real));
    cudaMalloc(&d_Psi_EG, Psi_EG.size() * sizeof(Real));
    cudaMalloc(&d_grid_w, grid_w.size() * sizeof(Real));
    cudaMalloc(&d_Phi_EG, Phi_EG.size() * sizeof(Real));

    Real *d_F_UG, *d_G_UG, *d_W_UG, *d_Psi_UG, *d_grid_w_UG, *d_Phi_UG;
    cudaMalloc(&d_F_UG, F_UG.size() * sizeof(Real));
    cudaMalloc(&d_G_UG, G_UG.size() * sizeof(Real));
    cudaMalloc(&d_W_UG, W_UG.size() * sizeof(Real));
    cudaMalloc(&d_Psi_UG, Psi_UG.size() * sizeof(Real));
    cudaMalloc(&d_grid_w_UG, grid_w_UG.size() * sizeof(Real));
    cudaMalloc(&d_Phi_UG, Phi_UG.size() * sizeof(Real));

    cudaMemcpy(d_F_EG, F_EG.data(), F_EG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_G_EG, G_EG.data(), G_EG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_W_EG, W_EG.data(), W_EG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Psi_EG, Psi_EG.data(), Psi_EG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_grid_w, grid_w.data(), grid_w.size() * sizeof(Real), cudaMemcpyHostToDevice);

    cudaMemcpy(d_F_UG, F_UG.data(), F_UG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_G_UG, G_UG.data(), G_UG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_W_UG, W_UG.data(), W_UG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Psi_UG, Psi_UG.data(), Psi_UG.size() * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(d_grid_w_UG, grid_w_UG.data(), grid_w_UG.size() * sizeof(Real), cudaMemcpyHostToDevice);

    GridParams eg_params{N, wmax, lp};
    GridUGParams ug_params{N_UG, wmax};
    OperatorArraysEG eg_arrays{d_F_EG, d_G_EG, d_W_EG, d_Psi_EG, d_grid_w, d_Phi_EG};
    OperatorArraysUG ug_arrays{d_F_UG, d_G_UG, d_W_UG, d_Psi_UG, d_grid_w_UG, d_Phi_UG};

    ApplyG0(eg_params, eg_arrays, ug_params, ug_arrays, /*use_eg=*/true, /*use_ug=*/true);

    cudaMemcpy(Phi_EG.data(), d_Phi_EG, Phi_EG.size() * sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(Phi_UG.data(), d_Phi_UG, Phi_UG.size() * sizeof(Real), cudaMemcpyDeviceToHost);

    std::cout << "Phi_EG result: ";
    for (int i = 0; i < std::min(10, int(Phi_EG.size())); ++i)
        std::cout << Phi_EG[i] << " ";
    std::cout << std::endl;

    std::cout << "Phi_UG result: ";
    for (int i = 0; i < std::min(10, int(Phi_UG.size())); ++i)
        std::cout << Phi_UG[i] << " ";
    std::cout << std::endl;

    cudaFree(d_F_EG); cudaFree(d_G_EG); cudaFree(d_W_EG); cudaFree(d_Psi_EG); cudaFree(d_grid_w); cudaFree(d_Phi_EG);
    cudaFree(d_F_UG); cudaFree(d_G_UG); cudaFree(d_W_UG); cudaFree(d_Psi_UG); cudaFree(d_grid_w_UG); cudaFree(d_Phi_UG);
}
