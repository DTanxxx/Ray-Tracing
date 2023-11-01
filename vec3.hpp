#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vec3 {
    public:
        // constructors
        vec3() : e{0, 0, 0} {}
        vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

        // getters
        double x() const { return e[0]; }
        double y() const { return e[1]; }
        double z() const { return e[2]; }

        // operator overloading
        vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }  // negating
        double operator[](int i) const { return e[i]; }  // indexing to get value
        double& operator[](int i) { return e[i]; }  // indexing to get reference
        vec3& operator+=(const vec3 &v) {
            e[0] += v.e[0];
            e[1] += v.e[1];
            e[2] += v.e[2];
            return *this;
        }
        vec3& operator*=(const double t) {
            e[0] *= t;
            e[1] *= t;
            e[2] *= t;
            return *this;
        }
        vec3& operator/=(const double t) {
            return *this *= 1/t;
        }
        
        // member functions
        double length() const {
            return sqrt(length_squared());
        }
        double length_squared() const {
            return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
        }
        bool near_zero() const {
            // return true if the vector is close to zero in all dimensions
            const auto s = 1e-8;
            return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
        }
        
        // static functions
        // create a random unit vector
        inline static vec3 random() {
            return vec3(random_double(), random_double(), random_double());
        }
        inline static vec3 random(double min, double max) {
            return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
        }
    
    public:
        double e[3];  // an array representing vec3
};

// Type aliases for vec3
using point3 = vec3;  // 3D point
using color = vec3;  // RGB color

// vec3 Utility Functions (inlined)
inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];  // overload cout operator
}
inline vec3 operator+(const vec3 &u, const vec3& v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
} 
inline vec3 operator-(const vec3& u, const vec3& v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}
inline vec3 operator*(const vec3& u, const vec3& v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}
inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}
inline vec3 operator*(const vec3& v, double t) {
    return t * v;
}
inline vec3 operator/(vec3 v, double t) {
    return 1/t * v;
}
inline double dot(const vec3& u, const vec3& v) {
    return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2];
}
inline vec3 cross(const vec3& u, const vec3& v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}
inline vec3 unit_vector(vec3 v) {
    return v / v.length();
}

// picks a random point in an unit radius sphere, starting from the sphere center
vec3 random_in_unit_sphere() {
    // keep sampling points within an unit cube until the point is also inside the unit radius sphere
    while (true) {
        auto p = vec3::random(-1, 1);
        if (p.length_squared() >= 1) {
            // point outside of unit radius sphere
            continue;
        }
        return p;
    }
}

// normalize the result from random_in_unit_sphere so the resulting point is on the sphere surface instead of in it
vec3 random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}

// generate a random point in a hemisphere region
vec3 random_in_hemisphere(const vec3& normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) {
        // in the same hemisphere as the normal
        return in_unit_sphere;
    }
    else {
        return -in_unit_sphere;
    }
}

// generate a random point in an unit disk region (in xy plane)
vec3 random_in_unit_disk() {
    // keep sampling points within an unit square until the point is also inside the unit disk
    while (true) {
        auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
        if (p.length_squared() >= 1) {
            // point outside of unit disk
            continue;
        }
        return p;
    }
}

// calculate the reflected ray direction given an incident ray v and normal vector n
vec3 reflect(const vec3& v, const vec3& n) {
    // reflected ray is v + 2 * b from the hit point on surface, where the length of b is dot(v, n) and the direction is in the normal
    return v - 2 * dot(v, n) * n;  // since v points "inwards", we need "v - 2 * b" instead of "v + 2 * b"
}

// calculate the refracted ray direction given an incident ray uv, a normal vector n, and a refractive indices ratio etai_over_etat
// Snell's law: sin(theta prime) = etai_over_etat * sin(theta) where theta prime is the angle between refracted ray and normal, and theta is the angle between incident ray and normal
vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = fmin(dot(-uv, n), 1.0);

    // find refracted ray by adding R_perpendicular and R_parallel, where R_perpendicular is the horizontal component of refracted ray 
    // (perpendicular to normal) and R_parallel is the vertical component of refracted ray (parallel to normal)
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);  // formula for computing R_perpendicular
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;  // formula for computing R_parallel
    return r_out_perp + r_out_parallel;
}

#endif
