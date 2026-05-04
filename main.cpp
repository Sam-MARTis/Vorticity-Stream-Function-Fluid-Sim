#include<iostream>
#include<cmath>
#include<vector>
#include<string>
#include<cstdlib>
#include<filesystem>
#include<random>
#include<SFML/Graphics.hpp>
#include "aux.hpp"
#include "constants.hpp"
#include "core-sim-functions.hpp"
#include "display-functions.hpp"
#include "initializations.hpp"
#include "imgui.h"
#include "imgui-SFML.h"

sf::RenderWindow window;

// Grid sizing
const int NX = 100;
const int NY = 100;
const int SCREEN_WIDTH = SCREEN_WIDTH_default;
const int SCREEN_HEIGHT = SCREEN_HEIGHT_default;
const int SCREEN_OFFSET_X = SCREEN_OFFSET_X_default;
const int SCREEN_OFFSET_Y = SCREEN_OFFSET_Y_default;
const int SCREEN_END_X_PADDING = SCREEN_END_X_PADDING_default;
const int SCREEN_END_Y_PADDING = SCREEN_END_Y_PADDING_default;



int main() {
    const float a = 1.0f;  // Domain width
    const float b = 1.0f;  // Domain height
    const float θ = M_PI/3.0f;  // Domain skew angle

    const float dt = 0.0001f;  // Temporal step size
    const float u0 = 1.0f;  // Lid velocity magnitude
    float reynolds_number = 100.0f;  // Reynolds number
    float nu = u0 * std::sqrt(a * b) / reynolds_number;  // Kinematic viscosity
    const float dims[3] = {a, b, θ};
    float* ψ = new float[NX*NY];  // Streamfunction
    float* ω = new float[NX*NY];  // Vorticity
    float* u = new float[2*NX*NY];  // Velocity field (u, v)
    float* x = new float[2*NX*NY];  // Coordinate field


    window.create(sf::VideoMode({SCREEN_WIDTH + (SCREEN_OFFSET_X + SCREEN_END_X_PADDING), SCREEN_HEIGHT + (SCREEN_OFFSET_Y + SCREEN_END_Y_PADDING)}, 10), "Fluid Simulation");
    window.setFramerateLimit(FRAME_RATE_LIMIT);
    ImGui::SFML::Init(window);
    float render_center[2] = {280.0f, 500.0f};
    float render_magnification = 1.4f;
    const float base_scaling[2] = {SCREEN_WIDTH/(a + b*std::cos(θ))*0.5f, SCREEN_WIDTH/(a + b*std::cos(θ))*0.5f};
    float scaling[2] = {base_scaling[0] * render_magnification, base_scaling[1] * render_magnification};

    int render_mode = 0;
    bool render_velocities_enabled = true;
    float normalization_constant = 1.0f;
    int velocity_thickness = 2;
    float velocity_head_fraction = 0.2f;
    int arrow_density_x = 24;
    int arrow_density_y = 24;
    bool standardize_velocity_sizes = false;
    float standardized_vector_size = 40.0f;
    float low_colour[3] = {0.0f, 0.0f, 1.0f};
    float high_colour[3] = {1.0f, 0.0f, 0.0f};
    int iter = 0;
    bool is_running = false;
    bool use_vorticity_operator_splitting = false;  // Use split advection and diffusion
    bool enable_advect_vorticity = true;  // Enable vorticity advection
    bool enable_apply_viscosity = true;  // Enable viscous diffusion
    bool has_selected_point = false;
    int stream_solver_iter_exponent = 4;  // Poisson solver max iterations = 10^exp
    float stream_solver_tolerance_exponent = -4.0f;  // Convergence tolerance = 10^exp
    bool enable_solver_parallelization = true;  // Use OpenMP
    int solver_max_threads = 19;
    int iterations_per_render = 100;  // Solver steps between render frames
    int plot_property_index = 0;
    int plot_num_levels = 20;
    int plot_contour_lines = 10;
    bool plot_enable_contour_areas = true;
    bool plot_enable_contour_lines = true;
    float plot_linewidth = 0.5f;
    float plot_alpha = 0.3f;
    bool plot_save_images = false;
    float stream_convergence_tolerance = 1e-5f;  // Poisson convergence check tolerance
    float stream_max_residual = 0.0f;
    bool stream_is_converged = false;
    float selected_world_x = 0.0f;
    float selected_world_y = 0.0f;
    float selected_u_x = 0.0f;
    float selected_u_y = 0.0f;
    float selected_omega = 0.0f;
    float physics_centroid[2] = {0.0f, 0.0f};
    int streamline_count = 100;
    float streamline_dt = dt * 100.0f;
    float streamline_max_time_exponent = 0.0f;
    float streamline_max_time = 1.0f;
    bool use_uniform_streamline_seeds = false;
    unsigned int streamline_seed = 1337u;
    std::vector<std::vector<sf::Vector2f>> cached_streamlines;
    std::vector<sf::Vector2f> clicked_streamline;
    bool has_clicked_streamline = false;
    sf::Clock deltaClock;

    auto build_streamline_path = [&](const float start_x, const float start_y) {
        const float safe_streamline_dt = std::max(1e-6f, streamline_dt);
        const int history_length = std::max(2, static_cast<int>(std::floor(streamline_max_time / safe_streamline_dt)));

        std::vector<float> pos_history(static_cast<std::size_t>(history_length) * 2u);
        obtain_streamline_path(x, u, u0, NX, NY, start_x, start_y, pos_history.data(), history_length, safe_streamline_dt, dims);

        std::vector<sf::Vector2f> positions(static_cast<std::size_t>(history_length));
        for(int i = 0; i < history_length; i++) {
            positions[static_cast<std::size_t>(i)] = {pos_history[2 * i], pos_history[2 * i + 1]};
        }

        return positions;
    };

    auto set_reynolds_number = [&](const float new_reynolds_number) {
        reynolds_number = std::max(1.0f, new_reynolds_number);
        nu = u0 * std::sqrt(a * b) / reynolds_number;
    };

    auto reset_simulation = [&]() {
        setup_inital_state(ψ, ω, x, u, NX, NY, dims, u0);  // Initialize vorticity and streamfunction
        compute_physics_centroid(x, NX, NY, physics_centroid);
        stream_is_converged = check_stream_function_convergence(ψ, ω, NX, NY, dims, stream_convergence_tolerance, stream_max_residual);
        has_selected_point = false;
        has_clicked_streamline = false;
        clicked_streamline.clear();
        cached_streamlines.clear();
        iter = 0;
        is_running = false;
        selected_world_x = 0.0f;
        selected_world_y = 0.0f;
        selected_u_x = 0.0f;
        selected_u_y = 0.0f;
        selected_omega = 0.0f;
        streamline_seed = 1337u;
    };



    setup_inital_state(ψ, ω, x, u, NX, NY, dims, u0);  // Initialize vorticity and streamfunction
    compute_physics_centroid(x, NX, NY, physics_centroid);
    stream_is_converged = check_stream_function_convergence(ψ, ω, NX, NY, dims, stream_convergence_tolerance, stream_max_residual);
    
    while(window.isOpen()) {
        while(const auto event = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *event);
            if(event->is<sf::Event::Closed>()) {
                window.close();
            }
            if(const auto* mouse_button = event->getIf<sf::Event::MouseButtonPressed>()) {
                if(mouse_button->button == sf::Mouse::Button::Left && !ImGui::GetIO().WantCaptureMouse) {
                    const sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
                    float world_x = 0.0f;
                    float world_y = 0.0f;
                    const bool is_inside_domain = map_screen_to_world(mouse_pos, world_x, world_y, physics_centroid, render_center, scaling, dims);
                    if(is_inside_domain) {
                        has_selected_point = true;
                        selected_world_x = world_x;
                        selected_world_y = world_y;
                        find_velocity_at_point(selected_u_x, selected_u_y, selected_world_x, selected_world_y, u, u0, NX, NY, dims);
                        find_vorticity_at_point(selected_omega, selected_world_x, selected_world_y, ω, NX, NY, dims);
                        clicked_streamline = build_streamline_path(selected_world_x, selected_world_y);
                        has_clicked_streamline = true;
                    } else {
                        has_selected_point = false;
                    }
                }
            }
        }

        ImGui::SFML::Update(window, deltaClock.restart());

        ImGui::Begin("Render Controls");
        bool do_single_step = false;
        ImGui::Text("Iteration: %d", iter);
        if(!is_running) {
            if(ImGui::Button("Start")) {
                is_running = true;
            }
            ImGui::SameLine();
            if(ImGui::Button("Step")) {
                do_single_step = true;
            }
        } else {
            if(ImGui::Button("Stop")) {
                is_running = false;
            }
            ImGui::SameLine();
            if(ImGui::Button("Step")) {
                do_single_step = true;
            }
        }
        ImGui::Separator();
        ImGui::Text("Vorticity");
        ImGui::Checkbox("Vorticity Operator Splitting", &use_vorticity_operator_splitting);
        if(use_vorticity_operator_splitting) {
            ImGui::Checkbox("Apply Viscosity", &enable_apply_viscosity);
            ImGui::Checkbox("Advect Vorticity", &enable_advect_vorticity);
        }
        if(ImGui::Button("Re = 100")) {
            set_reynolds_number(100.0f);
        }
        ImGui::SameLine();
        if(ImGui::Button("Re = 200")) {
            set_reynolds_number(200.0f);
        }
        ImGui::SameLine();
        if(ImGui::Button("Re = 500")) {
            set_reynolds_number(500.0f);
        }
        ImGui::SameLine();
        if(ImGui::Button("Re = 1000")) {
            set_reynolds_number(1000.0f);
        }
        const bool reynolds_changed = ImGui::SliderFloat("Reynolds Number", &reynolds_number, 1.0f, 10000.0f, "%.2f");
        if(reynolds_changed) {
            nu = u0 * std::sqrt(a * b) / std::max(1.0f, reynolds_number);
        }
        ImGui::SliderFloat("Viscosity", &nu, 0.0f, 50.0f, "%.4f");
        if(ImGui::IsItemEdited()) {
            reynolds_number = u0 * std::sqrt(a * b) / std::max(1e-6f, nu);
        }
        if(ImGui::Button("Reset Simulation")) {
            reset_simulation();
        }
        ImGui::Separator();
        ImGui::Text("Stream Function Solver");
        ImGui::SliderFloat("Domain Magnification", &render_magnification, 0.1f, 10.0f, "%.3f");
        ImGui::SliderFloat("Centroid Screen X", &render_center[0], 0.0f, static_cast<float>(window.getSize().x), "%.1f");
        ImGui::SliderFloat("Centroid Screen Y", &render_center[1], 0.0f, static_cast<float>(window.getSize().y), "%.1f");
        scaling[0] = base_scaling[0] * render_magnification;
        scaling[1] = base_scaling[1] * render_magnification;
        if(ImGui::Checkbox("Enable Parallelization", &enable_solver_parallelization)) {
            if(enable_solver_parallelization && solver_max_threads < 19) {
                solver_max_threads = 19;
            }
        }
        if(enable_solver_parallelization) {
            ImGui::SliderInt("Max Threads", &solver_max_threads, 1, 64);
        }
        ImGui::SliderInt("Iterations Per Render", &iterations_per_render, 1, 100);
        ImGui::SliderInt("Iterations Exponent", &stream_solver_iter_exponent, 2, 7);
        ImGui::SliderFloat("Tolerance Exponent", &stream_solver_tolerance_exponent, -10.0f, -2.0f, "%.1f");
        const int stream_solver_max_iterations = static_cast<int>(std::pow(10.0f, stream_solver_iter_exponent));
        const float stream_solver_tolerance = std::pow(10.0f, stream_solver_tolerance_exponent);
        stream_convergence_tolerance = stream_solver_tolerance;
        ImGui::Text("Iterations: %d", stream_solver_max_iterations);
        ImGui::Text("Tolerance: %.2e", stream_solver_tolerance);
        const char* render_modes[] = {"None", "Vorticity", "Stream Function"};
        ImGui::Combo("Field", &render_mode, render_modes, IM_ARRAYSIZE(render_modes));
        ImGui::ColorEdit3("Low Colour", low_colour);
        ImGui::ColorEdit3("High Colour", high_colour);
        ImGui::Checkbox("Render Velocities", &render_velocities_enabled);
        ImGui::SliderFloat("Normalization", &normalization_constant, 0.01f, 20.0f, "%.3f");
        ImGui::SliderInt("Thickness", &velocity_thickness, 1, 10);
        ImGui::SliderFloat("Arrow Head Fraction", &velocity_head_fraction, 0.01f, 0.5f, "%.2f");
        ImGui::SliderInt("Arrow Density X", &arrow_density_x, 1, 100);
        ImGui::SliderInt("Arrow Density Y", &arrow_density_y, 1, 100);
        ImGui::Checkbox("Standardize sizes", &standardize_velocity_sizes);
        if(standardize_velocity_sizes) {
            ImGui::SliderFloat("Vector Size", &standardized_vector_size, 1.0f, 200.0f, "%.2f");
        }
        ImGui::Separator();
        ImGui::Text("Streamlines");
        ImGui::Checkbox("Uniform Seed Spacing", &use_uniform_streamline_seeds);
        ImGui::SliderInt("Num Seed Points", &streamline_count, 2, 500);
        ImGui::SliderFloat("Streamline dt", &streamline_dt, 0.0001f, 0.5f, "%.6f");
        ImGui::SliderFloat("Streamline Max Time Exponent", &streamline_max_time_exponent, -2.0f, 4.5f, "%.2f");
        streamline_max_time = std::pow(10.0f, streamline_max_time_exponent);
        ImGui::Text("Streamline Max Time: %.3e", streamline_max_time);
        const bool update_streamline_pressed = ImGui::Button("update-streamline");
        ImGui::Separator();
        ImGui::Text("Stream Function Convergence");
        ImGui::Text("Tolerance: %.2e", stream_convergence_tolerance);
        ImGui::Text("Max residual: %.3e", stream_max_residual);
        ImGui::Text("Status: %s", stream_is_converged ? "Converged" : "Not converged");
        ImGui::Separator();
        ImGui::Text("Point Characteristics");
        if(has_selected_point) {
            ImGui::Text("x: %.5f", selected_world_x);
            ImGui::Text("y: %.5f", selected_world_y);
            ImGui::Text("vx: %.5f", selected_u_x);
            ImGui::Text("vy: %.5f", selected_u_y);
            ImGui::Text("vorticity: %.5f", selected_omega);
        } else {
            ImGui::Text("Click inside the fluid domain to sample.");
        }
        ImGui::Separator();
        ImGui::Text("Plot Options");
        ImGui::SliderFloat("Line Width", &plot_linewidth, 0.1f, 3.0f);
        ImGui::SliderFloat("Line Alpha", &plot_alpha, 0.0f, 1.0f);
        ImGui::Checkbox("Save Images to Images/", &plot_save_images);
        ImGui::Separator();
        const char* plot_properties[] = {"vorticity", "psi", "u", "v"};
        if(ImGui::Button("Export Centerlines & Plot")) {
            const char* outfname = "plotting_values.txt";
            const bool ok = export_velocity_centerlines(x, u, NX, NY, dims, outfname);
            if(ok) {
                try {
                    std::filesystem::create_directories("images");
                } catch(...) {
                }

                std::string cmd = std::string("python3 plot_velocity_centerlines.py ") + outfname +
                    " --linewidth " + std::to_string(plot_linewidth) +
                    " --alpha " + std::to_string(plot_alpha) +
                    " --save-images " + (plot_save_images ? "1" : "0") +
                    " --re " + std::to_string(reynolds_number) +
                    " --nx " + std::to_string(NX) +
                    " --ny " + std::to_string(NY) +
                    " --u0 " + std::to_string(u0);
                const int rc = std::system(cmd.c_str());
                (void)rc;
            } else {
                std::cerr << "export_velocity_centerlines failed\n";
            }
        }
        ImGui::Separator();
        ImGui::Text("State Plotting");
        ImGui::Combo("Property", &plot_property_index, plot_properties, IM_ARRAYSIZE(plot_properties));
        ImGui::SliderInt("Contour Levels", &plot_num_levels, 5, 100);
        ImGui::SliderInt("Contour Lines", &plot_contour_lines, 2, 100);
        ImGui::Checkbox("Plot Contour Areas", &plot_enable_contour_areas);
        ImGui::Checkbox("Plot Contour Lines", &plot_enable_contour_lines);
        if(ImGui::Button("Plot State")) {
            const char* outfname = "state_cache.txt";
            const bool ok = export_state_cache(x, ω, ψ, u, NX, NY, dims, outfname);
            if(ok) {
                try {
                    std::filesystem::create_directories("cache");
                } catch(...) {
                }
                const std::string prop_name = plot_properties[plot_property_index];
                std::string cmd = std::string("python3 plot_state_cache.py ") + outfname + " " + prop_name + 
                    " --levels " + std::to_string(plot_num_levels) +
                    " --lines " + std::to_string(plot_contour_lines) +
                    " --areas " + (plot_enable_contour_areas ? "1" : "0") +
                    " --plot-lines " + (plot_enable_contour_lines ? "1" : "0") +
                    " --linewidth " + std::to_string(plot_linewidth) +
                    " --alpha " + std::to_string(plot_alpha) +
                    " --save-images " + (plot_save_images ? "1" : "0") +
                    " --re " + std::to_string(reynolds_number) +
                    " --nx " + std::to_string(NX) +
                    " --ny " + std::to_string(NY) +
                    " --u0 " + std::to_string(u0);
                const int rc = std::system(cmd.c_str());
                (void)rc;
            } else {
                std::cerr << "export_state_cache failed\n";
            }
        }
        ImGui::End();

        set_solver_parallelization(enable_solver_parallelization, solver_max_threads);

        if(is_running || do_single_step) {
            const int stream_solver_max_iterations = static_cast<int>(std::pow(10.0f, stream_solver_iter_exponent));
            const float stream_solver_tolerance = std::pow(10.0f, stream_solver_tolerance_exponent);
            stream_convergence_tolerance = stream_solver_tolerance;
            const int solver_steps = is_running ? iterations_per_render : 1;
            for(int step = 0; step < solver_steps; step++) {
                solve_vorticity_transport(ω, x, u, u0, ψ, NX, NY, nu, dt, dims, use_vorticity_operator_splitting, enable_apply_viscosity, enable_advect_vorticity);  // Advect and diffuse vorticity
                solve_stream_function_update(ψ, ω, NX, NY, dims, stream_solver_max_iterations, stream_solver_tolerance);  // Solve Poisson for streamfunction
                stream_is_converged = check_stream_function_convergence(ψ, ω, NX, NY, dims, stream_convergence_tolerance, stream_max_residual);
                solve_velocity_update(u, u0, ψ, NX, NY, dims);  // Compute velocity from streamfunction
                iter++;

                if(has_selected_point) {
                    find_velocity_at_point(selected_u_x, selected_u_y, selected_world_x, selected_world_y, u, u0, NX, NY, dims);
                    find_vorticity_at_point(selected_omega, selected_world_x, selected_world_y, ω, NX, NY, dims);
                }
            }
        }

        window.clear(sf::Color::Black);

        const sf::Color low_sf(static_cast<unsigned char>(low_colour[0] * 255.0f),
                               static_cast<unsigned char>(low_colour[1] * 255.0f),
                               static_cast<unsigned char>(low_colour[2] * 255.0f));
        const sf::Color high_sf(static_cast<unsigned char>(high_colour[0] * 255.0f),
                                static_cast<unsigned char>(high_colour[1] * 255.0f),
                                static_cast<unsigned char>(high_colour[2] * 255.0f));

        if(render_mode == 1) {
            render_scalar_field(x, ω, NX, NY, physics_centroid, render_center, scaling, low_sf, high_sf, window);
        } else if(render_mode == 2) {
            render_scalar_field(x, ψ, NX, NY, physics_centroid, render_center, scaling, low_sf, high_sf, window);
        }

        if(render_velocities_enabled) {
            render_velocities(x, u, NX, NY, normalization_constant, velocity_thickness, velocity_head_fraction, standardize_velocity_sizes, standardized_vector_size, arrow_density_x, arrow_density_y, physics_centroid, render_center, scaling, window);
        }

        if(update_streamline_pressed) {
            streamline_seed += 1u;
            cached_streamlines.clear();
            cached_streamlines.reserve(streamline_count);

            const int grid_cols = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<float>(streamline_count)))));
            const int grid_rows = std::max(1, static_cast<int>(std::floor(std::sqrt(static_cast<float>(streamline_count)))));
            // std::cout << "Grid: " << grid_cols << " cols x " << grid_rows << " rows." << std::endl;
            // std::cout << "Streamline count: " << streamline_count << ", effective history length: " << history_length << std::endl;
            cached_streamlines.resize(static_cast<std::size_t>(streamline_count));

            #pragma omp parallel for schedule(dynamic) if(enable_solver_parallelization)
            for(int streamline_idx = 0; streamline_idx < streamline_count; streamline_idx++) {
                const int row = streamline_idx / grid_cols;
                const int col = streamline_idx % grid_cols;

                float xi = 0.0f;
                float eta = 0.0f;
                if(use_uniform_streamline_seeds) {
                    xi = (static_cast<float>(col) + 0.5f) / static_cast<float>(grid_cols);
                    eta = (static_cast<float>(row) + 0.5f) / static_cast<float>(grid_rows);
                } else {
                    std::mt19937 local_rng(streamline_seed + static_cast<unsigned int>(streamline_idx));
                    std::uniform_real_distribution<float> local_unit_dist(0.0f, 1.0f);
                    xi = local_unit_dist(local_rng);
                    eta = local_unit_dist(local_rng);
                }

                const float start_x = a * xi + b * eta * std::cos(θ);
                const float start_y = b * eta * std::sin(θ);
                cached_streamlines[static_cast<std::size_t>(streamline_idx)] = build_streamline_path(start_x, start_y);
            }
        }

        for(const auto& streamline_positions : cached_streamlines) {
            render_streamline_path(streamline_positions, physics_centroid, render_center, scaling, window, sf::Color::Cyan);
        }

        if(has_clicked_streamline) {
            const sf::Color orange_colour(255, 165, 0);
            render_streamline_path(clicked_streamline, physics_centroid, render_center, scaling, window, orange_colour);
        }

        ImGui::SFML::Render(window);
        window.display();

    }

    ImGui::SFML::Shutdown();
    delete[] ψ;
    delete[] ω;
    delete[] u;
    delete[] x;
    return 0;
}