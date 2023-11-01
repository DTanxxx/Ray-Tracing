#ifndef ONB_H
#define ONB_H

// onb stands for Orthonormal Basis

#include "constants.hpp"

// this class builds an orthonormal basis of s, t, n
class onb {
    public:
        onb() {}

        inline vec3 operator[](int i) const {
            return axis[i];
        }

        vec3 u() const {
            return axis[0];
        }

        vec3 v() const {
            return axis[1];
        }

        vec3 w() const {
            return axis[2];
        }

        vec3 local(double a, double b, double c) const {
            // calculate the local coordinates (a, b, c) relative to orthonormal basis (axes) u, v, w
            return a * u() + b * v() + c * w();
        }

        vec3 local(const vec3& a) const {
            return a.x() * u() + a.y() * v() + a.z() * w();
        }

        void build_from_w(const vec3&);

    public:
        vec3 axis[3];  // an array of vec3 representing a set of 3 orthonormal bases s, t, n
};

// this method creates an orthonormal basis set given a single vector n
void onb::build_from_w(const vec3& n) {
    axis[2] = unit_vector(n);  // n is one of the bases
    vec3 a = (fabs(w().x()) > 0.9) ? vec3(0, 1, 0) : vec3(1, 0, 0);  // pick an arbitrary vector a such that it is not parallel to n - we choose an axis that is not parallel to n
    axis[1] = unit_vector(cross(w(), a));  // find cross product of n and a, to get another basis t
    axis[0] = cross(w(), v());  // find cross product of n and t, to get another basis s
}

#endif
