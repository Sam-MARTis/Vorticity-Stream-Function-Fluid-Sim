#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <chrono>

void setup_inital_state(float* psi, float* omega, float* x, float* u, const int nx, const int ny, const float* dims, const float u0);
void set_solver_parallelization(bool enable_parallelization, int max_threads);
void solve_vorticity_transport(float* omega, const float* x, const float* u, const float u0, const float* psi, const int nx, const int ny, const float nu, const float dt, const float* dims, const bool use_operator_splitting, const bool enable_apply_viscosity, const bool enable_advect_vorticity);
void solve_stream_function_update(float* psi, const float* omega, const int nx, const int ny, const float* dims, const int max_iterations, const float tolerance);
void solve_velocity_update(float* u, const float u0, const float* psi, const int nx, const int ny, const float* dims);
bool export_velocity_centerlines(const float* x, const float* u, int nx, int ny, const float* dims, const char* filename);
bool export_state_cache(const float* x, const float* omega, const float* psi, const float* u, int nx, int ny, const float* dims, const char* filename);
void compute_physics_centroid(const float* x, const int nx, const int ny, float centroid[2]);

int main() {
    const int NX = 100;
    const int NY = 100;
    const float a = 1.0f;
    const float b = 1.0f;
    const float theta = static_cast<float>(M_PI/3.0);
    const float dims[3] = {a, b, theta};

    float reynolds = 100.0f;
    float dt = 5e-4f;
    int iterations = 10000;
    int stream_iters = 10000;
    float stream_tolerance = 1e-4f;
    float u0 = 1.0f;

    float line_width = 0.915f;
    float line_alpha = 1.0f;
    int contour_lines = 60;
    int contour_levels = 25;

    std::cout << "Headless Fluid Solver (CLI)\n";
    std::cout << "Press Enter to accept defaults shown in [brackets].\n";

    std::string line;
    std::cout << "Reynolds number [" << reynolds << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) reynolds = std::stof(line);

    std::cout << "Time step dt [" << dt << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) dt = std::stof(line);

    std::cout << "Number of iterations [" << iterations << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) iterations = std::stoi(line);

    std::cout << "Stream-function max iterations [" << stream_iters << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) stream_iters = std::stoi(line);

    std::cout << "Stream-function tolerance [" << stream_tolerance << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) stream_tolerance = std::stof(line);

    std::cout << "Lid speed u0 [" << u0 << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) u0 = std::stof(line);

    std::cout << "Line width [" << line_width << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) line_width = std::stof(line);

    std::cout << "Line alpha [" << line_alpha << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) line_alpha = std::stof(line);

    std::cout << "Contour lines [" << contour_lines << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) contour_lines = std::stoi(line);

    std::cout << "Contour levels [" << contour_levels << "]: ";
    std::getline(std::cin, line);
    if(!line.empty()) contour_levels = std::stoi(line);

    std::cout << "Plot after run? (Y/n): ";
    std::getline(std::cin, line);
    bool do_plot = (line.empty() || line[0] == 'y' || line[0] == 'Y');

    std::cout << "Save images when plotting? (y/N): ";
    std::getline(std::cin, line);
    bool save_images = (!line.empty() && (line[0] == 'y' || line[0] == 'Y'));

    const float nu = u0 * std::sqrt(a * b) / reynolds;

    std::cout << "Allocating arrays (" << NX << "x" << NY << ")...\n";
    float* psi = new float[NX * NY];
    float* omega = new float[NX * NY];
    float* u = new float[2 * NX * NY];
    float* x = new float[2 * NX * NY];

    std::cout << "Initializing fields...\n";
    setup_inital_state(psi, omega, x, u, NX, NY, dims, u0);

    set_solver_parallelization(false, 1);

    std::cout << "Running " << iterations << " iterations...\n";
    auto t0 = std::chrono::steady_clock::now();
    int report_interval = std::max(1, iterations / 100);
    for(int it = 0; it < iterations; ++it) {
        solve_vorticity_transport(omega, x, u, u0, psi, NX, NY, nu, dt, dims, true, true, true);
        solve_stream_function_update(psi, omega, NX, NY, dims, stream_iters, stream_tolerance);
        solve_velocity_update(u, u0, psi, NX, NY, dims);
        if(((it+1) % report_interval) == 0 || it == iterations-1) {
            std::cout << "Iteration " << (it+1) << " / " << iterations << "\n";
        }
    }
    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;
    std::cout << "Sim completed in " << elapsed.count() << " s\n";

    const char* state_fname = "state_cache.txt";
    const char* center_fname = "plotting_values.txt";
    std::cout << "Exporting state cache to " << state_fname << "...\n";
    bool ok_state = export_state_cache(x, omega, psi, u, NX, NY, dims, state_fname);
    std::cout << (ok_state ? "OK\n" : "FAILED\n");

    std::cout << "Exporting centerlines to " << center_fname << "...\n";
    bool ok_center = export_velocity_centerlines(x, u, NX, NY, dims, center_fname);
    std::cout << (ok_center ? "OK\n" : "FAILED\n");

    if(do_plot) {
        std::cout << "Calling Python plotters...\n";

        const std::string base_cmd = std::string("python3 plot_state_cache.py ") + state_fname;
        const std::string state_cmd_common = " --levels " + std::to_string(contour_levels) +
            " --lines " + std::to_string(contour_lines) +
            " --linewidth " + std::to_string(line_width) +
            " --alpha " + std::to_string(line_alpha) +
            " --save-images " + (save_images ? "1" : "0") +
            " --re " + std::to_string(reynolds) +
            " --nx " + std::to_string(NX) +
            " --ny " + std::to_string(NY) +
            " --u0 " + std::to_string(u0);

        const char* props[] = {"vorticity", "psi", "u", "v"};
        for(const char* prop : props) {
            std::string cmd = base_cmd + " " + prop + state_cmd_common;
            std::cout << "Running: " << cmd << "\n";
            std::system(cmd.c_str());
        }

        std::string center_cmd = std::string("python3 plot_velocity_centerlines.py ") + center_fname +
            " --linewidth " + std::to_string(line_width) +
            " --alpha " + std::to_string(line_alpha) +
            " --save-images " + (save_images ? "1" : "0") +
            " --re " + std::to_string(reynolds) +
            " --nx " + std::to_string(NX) +
            " --ny " + std::to_string(NY) +
            " --u0 " + std::to_string(u0);
        std::cout << "Running: " << center_cmd << "\n";
        std::system(center_cmd.c_str());
    }

    delete[] psi;
    delete[] omega;
    delete[] u;
    delete[] x;

    std::cout << "Done.\n";
    return 0;
}
