#pragma once 
#include<SFML/Graphics.hpp>

void render_scalar_field(const float* x,
						 const float* values,
						 int nx,
						 int ny,
						 const float centroid[2],
						 const float render_center[2],
						 const float scaling[2],
						 sf::Color low_colour,
						 sf::Color high_colour,
						 sf::RenderWindow& window);

void render_velocities(const float* x,
					   const float* u,
					   int nx,
					   int ny,
					   float normalization_factor,
					   int thickness,
				   const float centroid[2],
				   const float render_center[2],
					   const float scaling[2],
					   sf::RenderWindow& window);