#include<cmath>
#include "core-sim-functions.hpp"
void setup_inital_state(float* ψ, float* ω, float* x, float* u, const int nx, const int ny, const float* dims, const float u0) {
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float dξ = 1.0f/nx;
    const float dη = 1.0f/ny;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < nx * ny; i++) {
        ψ[i] = 0.0f;
        ω[i] = 0.0f;
        u[2*i] = 0.0f;
        u[2*i + 1] = 0.0f;
        const float ξ = (i % nx) * dξ;
        const float η = (i / nx) * dη;
        x[2*i] = a*ξ + b*η * std::cos(θ);
        x[2*i + 1] = b*η * std::sin(θ);
        u[2*i] = 0.0f;
        u[2*i + 1] = 0.0f;
    }
    for(int i = 0; i < nx; i++) {
        u[nx*(ny-1)*2 + 2*i] = u0;
    }
    solve_boundary_vorticity_values(ω, u0, ψ, nx, ny,dims);
}

