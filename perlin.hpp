#ifndef PERLIN_H
#define PERLIN_H

#include "constants.hpp"

// perlin noise is used to add some randomness to solid textures
// perlin noise is repeatable: it takes a 3D point as input and always returns the same
// randomish number. Nearby points return similar numbers
// perlin noise is simple and fast to execute
class perlin {
    public:
        perlin() {
            // tile all of space with a 3D array of random numbers and use them in blocks
            //ranfloat = new double[point_count];

            // use random unit vectors on the lattice points to remove blocky-looking noises
            // the blocky-looking noises are due to the min and max of the pattern always
            // landing exactly on the integer x/y/z
            // using random unit vectors instead of floats and dot-producting them to move
            // the min and max off the lattice will give us smoother noise effects
            ranvec = new vec3[point_count];
            for (int i = 0; i < point_count; ++i) {
                // these random unit vectors can be in irregular directions
                ranvec[i] = unit_vector(vec3::random(-1, 1));
            }

            perm_x = perlin_generate_perm();
            perm_y = perlin_generate_perm();
            perm_z = perlin_generate_perm();
        }

        ~perlin() {
            delete[] ranvec;
            delete[] perm_x;
            delete[] perm_y;
            delete[] perm_z;
        }

        double noise(const point3& p) const {
            auto u = p.x() - floor(p.x());
            auto v = p.y() - floor(p.y());
            auto w = p.z() - floor(p.z());

            auto i = static_cast<int>(floor(p.x()));
            auto j = static_cast<int>(floor(p.y()));
            auto k = static_cast<int>(floor(p.z()));
            vec3 c[2][2][2];

            for (int di = 0; di < 2; ++di) {
                for (int dj = 0; dj < 2; ++dj) {
                    for (int dk = 0; dk < 2; ++dk) {
                        c[di][dj][dk] = ranvec[
                            perm_x[(i + di) & 255] ^
                            perm_y[(j + dj) & 255] ^
                            perm_z[(k + dk) & 255]
                        ];
                    }
                }
            }

            //return trilinear_interp(c, u, v, w);
            return perlin_interp(c, u, v, w);
        }

        // a turbulence is a composite noise that has multiple summed frequencies, 
        // in other words it is a sum of repeated calls to noise()
        double turbulence(const point3& p, int depth = 7) const {
            auto accum = 0.0;
            auto temp_p = p;
            auto weight = 1.0;

            for (int i = 0; i < depth; ++i) {
                accum += weight * noise(temp_p);
                weight *= 0.5;
                temp_p *= 2;
            }

            return fabs(accum);
        }

    private:
        static const int point_count = 256;
        vec3* ranvec;
        int* perm_x;
        int* perm_y;
        int* perm_z;

        static int* perlin_generate_perm() {
            auto p = new int[point_count];

            for (int i = 0; i < perlin::point_count; ++i) {
                p[i] = i;
            }

            permute(p, point_count);

            return p;
        }

        // permute/shuffle a list of integers
        static void permute(int* p, int n) {
            for (int i = n - 1; i > 0; --i) {
                // choose a random index to swap with
                int target = random_int(0, i);

                // perform the swap
                int tmp = p[i];
                p[i] = p[target];
                p[target] = tmp;
            }
        }

        // apply trilinear interpolation to smooth out the noise effect
        static double trilinear_interp(double c[2][2][2], double u, double v, double w) {
            auto accum = 0.0;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        accum += (i * u + (1 - i) * (1 - u)) *
                                (j * v + (1 - j) * (1 - v)) *
                                (k * w + (1 - k) * (1 - w)) * c[i][j][k];
                    }
                }
            }

            return accum;
        }

        // interpolation but for vectors
        static double perlin_interp(vec3 c[2][2][2], double u, double v, double w) {
            // apply Hermite cubic to round off the interpolation to remove Mach bands
            auto uu = u * u * (3 - 2 * u);
            auto vv = v * v * (3 - 2 * v);
            auto ww = w * w * (3 - 2 * w);
            auto accum = 0.0;

            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        vec3 weight_v(u - i, v - j, w - k);
                        accum += (i * uu + (1 - i) * (1 - uu)) *
                                (j * vv + (1 - j) * (1 - vv)) *
                                (k * ww + (1 - k) * (1 - ww)) *
                                dot(c[i][j][k], weight_v);
                    }
                }
            }

            return accum;
        }
};

#endif
