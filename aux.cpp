#include <iostream>
#include <string>
#include <cmath>
#include "constants.hpp"

void print(const std::string& message) {
    std::cout << message << std::endl;
}

void find_velocity_at_point(float& u_x, float& u_y, const float px, const float py, const float* u, const float u0, const int nx, const int ny, const float* dims) {
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float pξ = (px - py/std::tan(θ)) / a;
    const float pη = py/(b * std::sin(θ));
    const float dξ = 1.0f / nx;
    const float dη = 1.0f / ny;
    const int i = static_cast<int>(std::floor(pξ / dξ));
    const int j = static_cast<int>(std::floor(pη / dη));
    // const int i = static_cast<int>(px / dx);
    // const int j = static_cast<int>(py / dy);
    if(i < 0 || i >= nx-1 || j < 0 || j >= ny-1) {
        u_x = u0 * (j >= ny-1); // top boundary has velocity u0, others have velocity 0
        u_y = 0.0f;
        return;
    }
    const float i_frac = (pξ - i * dξ) / dξ;
    const float j_frac = (pη - j * dη) / dη;
    float ux[4] = {0, 0, 0, 0}; // bottom left, top left, top right, bottom right
    float uy[4] = {0, 0, 0, 0};
    ux[0] = u[2*(j*nx + i)];
    uy[0] = u[2*(j*nx + i) + 1];

    ux[1] = u[2*((j+1)*nx + i)];
    uy[1] = u[2*((j+1)*nx + i) + 1];

    ux[2] = u[2*((j+1)*nx + (i+1))];
    uy[2] = u[2*((j+1)*nx + (i+1)) + 1];

    ux[3] = u[2*(j*nx + (i+1))];
    uy[3] = u[2*(j*nx + (i+1)) + 1];


    u_x = (ux[0] * (1-i_frac) + ux[3] * i_frac) * (1-j_frac) + (ux[1] * (1-i_frac) + ux[2] * i_frac) * j_frac;
    u_y = (uy[0] * (1-i_frac) + uy[3] * i_frac) * (1-j_frac) + (uy[1] * (1-i_frac) + uy[2] * i_frac) * j_frac;
}

void find_vorticity_at_point(float& ω_val, const float px, const float py, const float* ω, const int nx, const int ny, const float* dims) {
    const float a = dims[0];
    const float b = dims[1];
    const float θ = dims[2];
    const float pξ = (px - py/std::tan(θ)) / a;
    const float pη = py/(b * std::sin(θ));
    const float dξ = 1.0f / nx;
    const float dη = 1.0f / ny;
    const int i = static_cast<int>(std::floor(pξ / dξ));
    const int j = static_cast<int>(std::floor(pη / dη));
    const float i_frac = (pξ - i * dξ) / dξ;
    const float j_frac = (pη - j * dη) / dη;
    float ω_local[4] = {0, 0, 0, 0}; // bottom left, top left, top right, bottom right
    // if((i < 0 && (1-i_frac) > EPSILON1) || (i >= nx-1 && i_frac > EPSILON1) || (j < 0 && (1-j_frac) > EPSILON1) || (j >= ny-1 && j_frac > EPSILON1)) {
    //     ω_val = 0.0f;
    //     return;
    // }
    if((i < 0) || (i >= nx-1) || (j < 0) || (j >= ny-1)) {
        ω_val = 0.0f;
        return;
    }
    ω_local[0] = ω[j*nx + i];

    ω_local[1] = ω[(j+1)*nx + i];

    ω_local[2] = ω[(j+1)*nx + (i+1)];
    ω_local[3] = ω[j*nx + (i+1)];
    ω_val = (ω_local[0] * (1-i_frac) + ω_local[3] * i_frac) * (1-j_frac) + (ω_local[1] * (1-i_frac) + ω_local[2] * i_frac) * j_frac;
}