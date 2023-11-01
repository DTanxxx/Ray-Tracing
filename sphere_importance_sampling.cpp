#include "constants.hpp"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

// in ray tracers we pick random directions instead of points, we can represent
// directions by using points on the unit sphere, so now we need to have a PDF 
// defined over 2D - the PDF of these random points will be 1/area of the sphere
// or 1/(4 * pi) since it is a density on the unit sphere

inline double pdf(const vec3& p) {
    return 1 / (4 * pi);
}

int main() {
    // generate importance-sampled points on the unit sphere
    int N = 1000000;
    auto sum = 0.0;
    for (int i = 0; i < N; ++i) {
        vec3 d = random_unit_vector();
        // assume that the integrand is cos^2(theta), and theta is the angle with the z axis
        auto cosine_squared = d.z() * d.z();
        sum += cosine_squared / pdf(d);
    }

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "I = " << sum / N << '\n';
}
