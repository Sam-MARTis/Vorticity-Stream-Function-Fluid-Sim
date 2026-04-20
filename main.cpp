#include<iostream>
#include<cmath>
#include<vector>
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
const int NX = 50;
const int NY = 50;
const int SCREEN_WIDTH = SCREEN_WIDTH_default;
const int SCREEN_HEIGHT = SCREEN_HEIGHT_default;
const int SCREEN_OFFSET_X = SCREEN_OFFSET_X_default;
const int SCREEN_OFFSET_Y = SCREEN_OFFSET_Y_default;
const int SCREEN_END_X_PADDING = SCREEN_END_X_PADDING_default;
const int SCREEN_END_Y_PADDING = SCREEN_END_Y_PADDING_default;



int main() {
    const float a = 1.0f;
    const float b = 1.0f;
    const float θ = M_PI/2.0f;

    // const int NX = 100;
    // const int NY = 100;
    // const float dx = ;
    // const float dy = 0.01f;

    const float dt = 0.001f;
    const float u0 = 1.0f;
    float reynolds_number = 100.0f;
    float nu = u0 * std::sqrt(a * b) / reynolds_number;
    const float dims[3] = {a, b, θ};
    float* ψ = new float[NX*NY];
    float* ω = new float[NX*NY];
    float* u = new float[2*NX*NY];
    float* x = new float[2*NX*NY];


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
    float low_colour[3] = {0.0f, 0.0f, 1.0f};
    float high_colour[3] = {1.0f, 0.0f, 0.0f};
    int iter = 0;
    bool is_running = false;
    bool use_vorticity_operator_splitting = false;
    bool enable_advect_vorticity = true;
    bool enable_apply_viscosity = true;
    bool has_selected_point = false;
    int stream_solver_iter_exponent = 4;
    float stream_solver_tolerance_exponent = -4.0f;
    bool enable_solver_parallelization = false;
    int solver_max_threads = 4;
    int iterations_per_render = 1;
    float stream_convergence_tolerance = 1e-5f;
    float stream_max_residual = 0.0f;
    bool stream_is_converged = false;
    float selected_world_x = 0.0f;
    float selected_world_y = 0.0f;
    float selected_u_x = 0.0f;
    float selected_u_y = 0.0f;
    float selected_omega = 0.0f;
    float physics_centroid[2] = {0.0f, 0.0f};
    const int streamline_count = 24;
    int streamline_history_length = 120;
    float streamline_dt = dt;
    unsigned int streamline_seed = 1337u;
    std::vector<std::vector<sf::Vector2f>> cached_streamlines;
    sf::Clock deltaClock;



    setup_inital_state(ψ, ω, x, u, NX, NY, dims, u0);
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
        const bool reynolds_changed = ImGui::SliderFloat("Reynolds Number", &reynolds_number, 1.0f, 10000.0f, "%.2f");
        if(reynolds_changed) {
            nu = u0 * std::sqrt(a * b) / std::max(1.0f, reynolds_number);
        }
        ImGui::SliderFloat("Viscosity", &nu, 0.0f, 50.0f, "%.4f");
        if(ImGui::IsItemEdited()) {
            reynolds_number = u0 * std::sqrt(a * b) / std::max(1e-6f, nu);
        }
        ImGui::Separator();
        ImGui::Text("Stream Function Solver");
        ImGui::SliderFloat("Domain Magnification", &render_magnification, 0.1f, 10.0f, "%.3f");
        ImGui::SliderFloat("Centroid Screen X", &render_center[0], 0.0f, static_cast<float>(window.getSize().x), "%.1f");
        ImGui::SliderFloat("Centroid Screen Y", &render_center[1], 0.0f, static_cast<float>(window.getSize().y), "%.1f");
        scaling[0] = base_scaling[0] * render_magnification;
        scaling[1] = base_scaling[1] * render_magnification;
        ImGui::Checkbox("Enable Parallelization", &enable_solver_parallelization);
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
        ImGui::SliderInt("Streamline Points", &streamline_history_length, 2, 1000);
        ImGui::SliderFloat("Streamline dt", &streamline_dt, 0.0001f, 0.05f, "%.5f");
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
        ImGui::End();

        set_solver_parallelization(enable_solver_parallelization, solver_max_threads);

        if(is_running || do_single_step) {
            const int stream_solver_max_iterations = static_cast<int>(std::pow(10.0f, stream_solver_iter_exponent));
            const float stream_solver_tolerance = std::pow(10.0f, stream_solver_tolerance_exponent);
            stream_convergence_tolerance = stream_solver_tolerance;
            const int solver_steps = is_running ? iterations_per_render : 1;
            for(int step = 0; step < solver_steps; step++) {
                solve_vorticity_transport(ω, x, u, u0, ψ, NX, NY, nu, dt, dims, use_vorticity_operator_splitting, enable_apply_viscosity, enable_advect_vorticity);
                solve_stream_function_update(ψ, ω, NX, NY, dims, stream_solver_max_iterations, stream_solver_tolerance);
                stream_is_converged = check_stream_function_convergence(ψ, ω, NX, NY, dims, stream_convergence_tolerance, stream_max_residual);
                solve_velocity_update(u, u0, ψ, NX, NY, dims);
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
            render_velocities(x, u, NX, NY, normalization_constant, velocity_thickness, physics_centroid, render_center, scaling, window);
        }

        if(update_streamline_pressed) {
            streamline_seed += 1u;
            std::mt19937 rng(streamline_seed);
            std::uniform_real_distribution<float> unit_dist(0.0f, 1.0f);
            cached_streamlines.clear();
            cached_streamlines.reserve(streamline_count);

            for(int streamline_idx = 0; streamline_idx < streamline_count; streamline_idx++) {
                const float xi = unit_dist(rng);
                const float eta = unit_dist(rng);
                const float start_x = a * xi + b * eta * std::cos(θ);
                const float start_y = b * eta * std::sin(θ);

                std::vector<float> pos_history(static_cast<std::size_t>(streamline_history_length) * 2u);
                obtain_streamline_path(x, u, u0, NX, NY, start_x, start_y, pos_history.data(), streamline_history_length, streamline_dt, dims);

                std::vector<sf::Vector2f> positions;
                positions.reserve(static_cast<std::size_t>(streamline_history_length));
                for(int i = 0; i < streamline_history_length; i++) {
                    positions.emplace_back(pos_history[2 * i], pos_history[2 * i + 1]);
                }

                cached_streamlines.push_back(std::move(positions));
            }
        }

        for(const auto& streamline_positions : cached_streamlines) {
            render_streamline_path(streamline_positions, physics_centroid, render_center, scaling, window, sf::Color::Cyan);
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