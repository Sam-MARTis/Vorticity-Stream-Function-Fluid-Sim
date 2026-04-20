#pragma once
#include <string>
#include "SFML/Graphics.hpp"

void print(const std::string &message);
void find_velocity_at_point(float &u_x, float &u_y, const float px, const float py, const float *u, const float u0, const int nx, const int ny, const float *dims);
void find_vorticity_at_point(float &ω_val, const float px, const float py, const float *ω, const int nx, const int ny, const float *dims);
bool map_screen_to_world(const sf::Vector2i &screen_point, float &world_x, float &world_y, const float centroid[2], const float render_center[2], const float scaling[2], const float *dims);
void compute_physics_centroid(const float* x, const int nx, const int ny, float centroid[2]);