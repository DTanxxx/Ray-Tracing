#ifndef COLOR_H
#define COLOR_H

#include "vec3.hpp"
#include <iostream>

// handles multi-sampled color computation (for antialiasing)
void write_color(std::ostream& out, color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // replace NaN components with zero, to remove surface acne
    // here we are using the fact that a NaN does not equal itself, so we can check whether
    // a component is not equal to itself and if so, we know it is a NaN and we will replace 
    // it with 0
    if (r != r) {
        r = 0.0;
    }
    if (g != g) {
        g = 0.0;
    }
    if (b != b) {
        b = 0.0;
    }

    // divide the color by the number of samples and gamma-correct for gamma = 2.0
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(r * scale);  // apply gamma correction (gamma = 2.0, so raise the color to the power 1/gamma = 1/2 = square root)
    g = sqrt(g * scale);
    b = sqrt(b * scale);

    // write the translated [0,255] value of each color component (pixel_color's values are from [0,1])
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' ' 
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' ' 
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

#endif
