#include<SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include "constants.hpp"
#include "aux.hpp"

void draw_arrow(const float px, const float py, const float dPx, const float dPy, sf::RenderWindow& window, const int thickness = 2, const float head_fraction = 0.2f, sf::Color colour = sf::Color::Red);
void obtain_streamline_path(const float* x, const float* u, const float u0, const int nx, const int ny, const float posx, const float posy, float* pos_history, const int history_length, const float dt, const float* dims);

sf::Color interpolate_colour(const sf::Color low, const sf::Color high, const float t) {
    const float clamped_t = std::clamp(t, 0.0f, 1.0f);
    const auto lerp_channel = [clamped_t](const unsigned char a, const unsigned char b) {
        return static_cast<unsigned char>(a + (b - a) * clamped_t);
    };
    return sf::Color(lerp_channel(low.r, high.r),
                     lerp_channel(low.g, high.g),
                     lerp_channel(low.b, high.b),
                     lerp_channel(low.a, high.a));
}


void render_scalar_field(const float* x, const float* values, const int nx, const int ny, const float centroid[2], const float render_center[2], const float scaling[2], const sf::Color low_colour, const sf::Color high_colour, sf::RenderWindow& window) {
    const int count = nx * ny;
    if(count <= 0 || nx < 2 || ny < 2) {
        return;
    }

    float min_val = values[0];
    float max_val = values[0];
    for(int idx = 1; idx < count; idx++) {
        min_val = std::min(min_val, values[idx]);
        max_val = std::max(max_val, values[idx]);
    }

    const float range = std::max(max_val - min_val, 1e-12f);
    sf::VertexArray cells(sf::PrimitiveType::Triangles);

    for(int j = 0; j < ny - 1; j++) {
        for(int i = 0; i < nx - 1; i++) {
            const int idx_bl = j * nx + i;
            const int idx_br = j * nx + (i + 1);
            const int idx_tl = (j + 1) * nx + i;
            const int idx_tr = (j + 1) * nx + (i + 1);

            const float cell_value = 0.25f * (values[idx_bl] + values[idx_br] + values[idx_tl] + values[idx_tr]);
            const float normalized = (cell_value - min_val) / range;
            const sf::Color colour = interpolate_colour(low_colour, high_colour, normalized);

            auto to_screen = [&](const int idx) {
                return sf::Vector2f(render_center[0] + (x[2 * idx] - centroid[0]) * scaling[0],
                                    render_center[1] - (x[2 * idx + 1] - centroid[1]) * scaling[1]);
            };

            const sf::Vector2f p_bl = to_screen(idx_bl);
            const sf::Vector2f p_br = to_screen(idx_br);
            const sf::Vector2f p_tl = to_screen(idx_tl);
            const sf::Vector2f p_tr = to_screen(idx_tr);

            const sf::Vertex vertices[6] = {
                {p_bl, colour},
                {p_br, colour},
                {p_tr, colour},
                {p_bl, colour},
                {p_tr, colour},
                {p_tl, colour}
            };

       
            for(const auto& vertex : vertices) {
                cells.append(vertex);
            }
        }
    }

    window.draw(cells);
}

void render_velocities(const float* x, const float* u, const int nx, const int ny, const float normalization_factor, const int thickness, const float head_fraction, const bool standardize_sizes, const float standardized_vector_size, const int density_x, const int density_y, const float centroid[2], const float render_center[2], const float scaling[2], sf::RenderWindow& window) {
    const float safe_norm = std::max(1e-6f, normalization_factor);
    const float inv_normalization_factor = BASE_ARROW_SIZE / safe_norm;
    const int stride_x = std::max(1, nx / std::max(1, density_x));
    const int stride_y = std::max(1, ny / std::max(1, density_y));

    for(int j = 0; j < ny; j += stride_y) {
        for(int i = 0; i < nx; i += stride_x) {
            const int idx = j * nx + i;
            const float px = render_center[0] + (x[2 * idx] - centroid[0]) * scaling[0];
            const float py = render_center[1] - (x[2 * idx + 1] - centroid[1]) * scaling[1];
            float u_x = 0.0f;
            float u_y = 0.0f;
            if(standardize_sizes) {
                const float velocity_x = u[2 * idx];
                const float velocity_y = -u[2 * idx + 1];
                const float velocity_length = std::sqrt(velocity_x * velocity_x + velocity_y * velocity_y);
                if(velocity_length > 1e-6f) {
                    const float scale = standardized_vector_size / velocity_length;
                    u_x = velocity_x * scale;
                    u_y = velocity_y * scale;
                }
            } else {
                u_x = u[2 * idx] * inv_normalization_factor;
                u_y = -u[2 * idx + 1] * inv_normalization_factor;
            }
            draw_arrow(px, py, u_x, u_y, window, thickness, head_fraction);
        }
    }
}
void draw_arrow(const float px, const float py, const float dPx, const float dPy, sf::RenderWindow& window, const int thickness, const float head_fraction, sf::Color colour) {
    sf::Vector2f start(px, py);
    sf::Vector2f direction(dPx, dPy);
    sf::Vector2f end = start + direction;
    const float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
    const sf::Angle angleDegrees = std::atan2(direction.y, direction.x) * sf::radians(1);
    const float head_length = length * head_fraction;
    sf::RectangleShape line(sf::Vector2f(length - head_length, (float)thickness));
    line.setOrigin({0, (float)thickness * 0.5f});
    line.setPosition(start);
    line.setRotation(angleDegrees);
    line.setFillColor(colour);

    sf::ConvexShape head;
    head.setPointCount(3);
    head.setPoint(0, {0, 0});
    head.setPoint(1, {-head_length, head_length * 0.5f});
    head.setPoint(2, {-head_length, -head_length * 0.5f});
    head.setFillColor(colour);

    head.setPosition(end);
    head.setRotation(angleDegrees);

    window.draw(line);
    window.draw(head);
}

void render_streamline_path(const std::vector<sf::Vector2f>& positions,
								 const float centroid[2],
								 const float render_center[2],
								 const float scaling[2],
								 sf::RenderWindow& window,
								 sf::Color colour) {
    if(positions.size() < 2) {
        return;
    }

    sf::VertexArray line_strip(sf::PrimitiveType::LineStrip);
    for(const auto& position : positions) {
        sf::Vertex vertex;
        vertex.position = {
            render_center[0] + (position.x - centroid[0]) * scaling[0],
            render_center[1] - (position.y - centroid[1]) * scaling[1]
        };
        vertex.color = colour;
        line_strip.append(vertex);
    }

    window.draw(line_strip);
}

void render_random_streamlines(const float* x, const float* u, const float u0, const int nx, const int ny, const int streamline_count, const int history_length, const float dt, const float* dims, const float centroid[2], const float render_center[2], const float scaling[2], sf::RenderWindow& window, sf::Color colour, unsigned int seed) {

    const float a = dims[0];
    const float b = dims[1];
    const float theta = dims[2];

    std::mt19937 rng(seed); // Please kill me
    std::uniform_real_distribution<float> unit_dist(0.0f, 1.0f); 

    for(int streamline_idx = 0; streamline_idx < streamline_count; streamline_idx++) {
        const float xi = unit_dist(rng);
        const float eta = unit_dist(rng);
        const float start_x = a * xi + b * eta * std::cos(theta);
        const float start_y = b * eta * std::sin(theta);

        std::vector<float> pos_history(static_cast<std::size_t>(history_length) * 2u);
        obtain_streamline_path(x, u, u0, nx, ny, start_x, start_y, pos_history.data(), history_length, dt, dims);

        std::vector<sf::Vector2f> positions;
        positions.reserve(static_cast<std::size_t>(history_length));
        for(int i = 0; i < history_length; i++) {
            positions.emplace_back(pos_history[2 * i], pos_history[2 * i + 1]);
        }

        render_streamline_path(positions, centroid, render_center, scaling, window, colour);
    }
}



