#include<SFML/Graphics.hpp>
#include <algorithm>
#include <cmath>
#include "constants.hpp"

void draw_arrow(const float px, const float py, const float dPx, const float dPy, sf::RenderWindow& window, const int thickness = 2, const float head_fraction = 0.2f, sf::Color colour = sf::Color::Red);

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

void render_scalar_field(const float* x,
                         const float* values,
                         const int nx,
                         const int ny,
                         const float origin[2],
                         const float scaling[2],
                         const sf::Color low_colour,
                         const sf::Color high_colour,
                         sf::RenderWindow& window) {
    const int count = nx * ny;
    if(count <= 0) {
        return;
    }

    float min_val = values[0];
    float max_val = values[0];
    for(int idx = 1; idx < count; idx++) {
        min_val = std::min(min_val, values[idx]);
        max_val = std::max(max_val, values[idx]);
    }

    const float range = std::max(max_val - min_val, 1e-12f);
    const sf::Vector2f cell_size(std::max(1.0f, scaling[0] * 0.9f), std::max(1.0f, scaling[1] * 0.9f));

    sf::RectangleShape cell(cell_size);
    for(int idx = 0; idx < count; idx++) {
        const float normalized = (values[idx] - min_val) / range;
        cell.setFillColor(interpolate_colour(low_colour, high_colour, normalized));
        cell.setPosition(sf::Vector2f(origin[0] + x[2 * idx] * scaling[0], origin[1] - x[2 * idx + 1] * scaling[1]));
        window.draw(cell);
    }
}

void render_velocities(const float* x,
                       const float* u,
                       const int nx,
                       const int ny,
                       const float normalization_factor,
                       const int thickness,
                       const float origin[2],
                       const float scaling[2],
                       sf::RenderWindow& window) {
    const float safe_norm = std::max(1e-6f, normalization_factor);
    const float inv_normalization_factor = BASE_ARROW_SIZE / safe_norm;
    const int stride = std::max(1, nx / 24);

    for(int j = 0; j < ny; j += stride) {
        for(int i = 0; i < nx; i += stride) {
            const int idx = j * nx + i;
            const float px = origin[0] + x[2 * idx] * scaling[0];
            const float py = origin[1] - x[2 * idx + 1] * scaling[1];
            const float u_x = u[2 * idx] * inv_normalization_factor;
            const float u_y = -u[2 * idx + 1] * inv_normalization_factor;
            draw_arrow(px, py, u_x, u_y, window, thickness);
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

// void drawArrow(sf::RenderWindow &window, sf::Vector2f start, sf::Vector2f end,
//                int thickness = 2, float head_fraction = 0.2f,
//                sf::Color colour = sf::Color::Red)
// {
//     sf::Vector2f direction = end - start;
//     float length = std::sqrt(direction.x * direction.x + direction.y * direction.y);
//     sf::Angle angleDegrees = std::atan2(direction.y, direction.x) * sf::radians(1);

//     float head_length = length * head_fraction;
//     sf::RectangleShape line(sf::Vector2f(length - head_length, (float)thickness));
//     line.setOrigin({0, (float)thickness * 0.5f});
//     line.setPosition(start);
//     line.setRotation(angleDegrees);
//     line.setFillColor(colour);

//     sf::ConvexShape head;
//     head.setPointCount(3);
//     head.setPoint(0, {0, 0});
//     head.setPoint(1, {-head_length, head_length * 0.5f});
//     head.setPoint(2, {-head_length, -head_length * 0.5f});
//     head.setFillColor(colour);

//     head.setPosition(end);
//     head.setRotation(angleDegrees);

//     window.draw(line);
//     window.draw(head);
// }
