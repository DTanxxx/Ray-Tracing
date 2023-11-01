#include "constants.hpp"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

// although the stratified method converges with a better asymptotic rate,
// this advantage decreases with the dimension of the problem - this is the
// Curse of Dimensionality
// since we will be very high dimensional (each reflection adds two dimensions),
// we won't be stratifying here, although if you want to do single-reflection
// or shadowing or some strictly 2D problem, stratifying is essential
void stratified_sampling() {
    // take a hundred million samples
    int inside_cirle = 0;
    int inside_circle_stratified = 0;
    int sqrt_N = 10000;
    for (int i = 0; i < sqrt_N; ++i) {
        for (int j = 0; j < sqrt_N; ++j) {
            auto x = random_double(-1, 1);
            auto y = random_double(-1, 1);
            if (x * x + y * y < 1) {
                inside_cirle += 1;
            }

            x = 2 * ((i + random_double()) / sqrt_N) - 1;
            y = 2 * ((j + random_double()) / sqrt_N) - 1;
            if (x * x + y * y < 1) {
                inside_circle_stratified += 1;
            }
        }
    }

    auto N = static_cast<double>(sqrt_N) * sqrt_N;
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "Regular Estimate of Pi = "
              << 4 * double(inside_cirle) / (sqrt_N * sqrt_N) << '\n'
              << "Stratified Estimate of Pi = "
              << 4 * double(inside_circle_stratified) / (sqrt_N * sqrt_N) << '\n';
}

void running_estimate() {
    int inside_circle = 0;
    int runs = 0;
    std::cout << std::fixed << std::setprecision(12);
    while (true) {
        runs += 1;
        auto x = random_double(-1, 1);
        auto y = random_double(-1, 1);
        if (x * x + y * y < 1) {
            inside_circle += 1;
        }

        if (runs % 100000 == 0) {
            std::cout << "Estimate of Pi = " << 4 * double(inside_circle) / runs << '\n';
        }
    }
}

int main() {
    // we will estimate pi here, suppose we have a circle inscribed inside a square
    // and we pick random points inside the square. The fraction of those random
    // points that end up inside the circle should be proportional to the area of 
    // the circle
    // in fact, this fraction value is equal to (pi * r^2)/(2 * r)^2 = pi / 4
    // since r cancels out, we will choose arbitrarily a circle of radius 1 centered
    // at the origin
    int N = 1000;
    int inside_circle = 0;
    for (int i = 0; i < N; ++i) {
        auto x = random_double(-1, 1);
        auto y = random_double(-1, 1);
        if (x * x + y * y < 1) {
            // Pythagoras, check if the random point is inside a circle
            // of radius 1, given that the square has side length = 2
            inside_circle += 1;
        }
    }

    std::cout << std::fixed << std::setprecision(12);
    
    // rearrange the equation to estimate pi: pi = 4 * fraction
    std::cout << "Estimate of Pi = " << 4 * double(inside_circle) / N << '\n';

    // we can also just run the iteration forever and print out a running estimate
    //running_estimate();

    // using previous approaches, we are limited by the Law of Diminishing Returns,
    // where each sample helps less than the last in estimating pi
    // we can mitigate this diminishing return by stratifying the samples (aka jittering),
    // where instead of taking random samples, we take a grid and take one sample within
    // each - we will need to know how many samples we are taking in advance because we
    // need to know the grid
    //stratified_sampling();

    return 0;
}
