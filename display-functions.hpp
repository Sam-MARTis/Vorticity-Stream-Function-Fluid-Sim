#pragma once 
#include<SFML/Graphics.hpp>
#include <vector>

void render_scalar_field(const float* x, const float* values, int nx, int ny, const float centroid[2], const float render_center[2], const float scaling[2], sf::Color low_colour, sf::Color high_colour, sf::RenderWindow& window);
void render_velocities(const float* x, const float* u, int nx, int ny, float normalization_factor, int thickness, const float centroid[2], const float render_center[2], const float scaling[2], sf::RenderWindow& window);
void render_streamline_path(const std::vector<sf::Vector2f>& positions, const float centroid[2], const float render_center[2], const float scaling[2], sf::RenderWindow& window, sf::Color colour = sf::Color::White);
void render_random_streamlines(const float* x, const float* u, const float u0, const int nx, const int ny, const int streamline_count, const int history_length, const float dt, const float* dims, const float centroid[2], const float render_center[2], const float scaling[2], sf::RenderWindow& window, sf::Color colour = sf::Color::White, unsigned int seed = 1337u);