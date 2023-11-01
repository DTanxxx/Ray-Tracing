#include "constants.hpp"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

// Monte Carlo ray tracer
// 1. You have an integral of f(x) over some domain [a,b]
// 2. You pick a PDF p that is non-zero over [a,b]
// 3. You average a whole ton of f(r)/p(r) where r is a random number with PDF p
// note: any choice of PDF p will always converge to the right answer, but the 
// closer that p approximates f, the faster that it will converge

inline double accurate_pdf(double x) {
    return 3 * x * x / 8;
}

void perfect_importance_sampling() {
    // we only need one sample to see the effect
    int N = 1;
    auto sum = 0.0;
    for (int i = 0; i < N; ++i) {
        auto x = pow(random_double(0, 8), 1.0 / 3.0);
        sum += x * x / accurate_pdf(x);
    }

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "I = " << sum / N << '\n';  // this always returns the exact answer
}

inline double constant_pdf(double x) {
    return 0.5;
}

void importance_sampling() {
    int N = 1000000;
    auto sum = 0.0;
    for (int i = 0; i < N; ++i) {
        auto x = random_double(0, 2);
        sum += x * x / constant_pdf(x);
    }

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "I = " << sum / N << '\n';
}

inline double pdf(double x) {
    return 0.5 * x;
}

void integrate_with_pdf() {
    int N = 1000000;
    auto sum = 0.0;
    for (int i = 0; i < N; ++i) {
        auto x = sqrt(random_double(0, 4));
        sum += x * x / pdf(x);
    }

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "I = " << sum / N << '\n';
}

// let's say that we want to integrate x^2 from 0 to 2, we will write this as area(x^2, 0, 2)
// we can also write it as 2 * average(x^2, 0, 2)

int main() {
    int N = 1000000;
    auto sum = 0.0;
    for (int i = 0; i < N; ++i) {
        auto x = random_double(0, 2);
        sum += x * x;
    }

    std::cout << std::fixed << std::setprecision(12);

    // apply the formula 2 * average(x^2, 0, 2) where the average(x^2, 0, 2) term is equal
    // to our "sum" divided by "N" - Monte Carlo approximation
    std::cout << "I = " << 2 * sum / N << '\n';

    // we can integrate x^2 with PDF (the maths is too complicated to understand lmao)
    //integrate_with_pdf();

    // importance sampling - by using a non-uniform PDF, we are steering our samples
    // toward the parts of the distribution that are more important (we are sampling more
    // where the integrand is big, so we can expect less noise and thus faster convergence)
    //importance_sampling();

    // perfect importance sampling - testing our importance sampling code by using a known
    // accurate PDF for our distribution
    //perfect_importance_sampling();

    return 0;
}
