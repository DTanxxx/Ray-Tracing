#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.hpp"
#include "vec3.hpp"
#include "pdf.hpp"

// sphere inherits from hittable
class sphere : public hittable {
    public:
        sphere() {}
        sphere(point3 cen, double r, shared_ptr<material> m) : center(cen), radius(r), mat_ptr(m) {};

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;
        virtual double pdf_value(const point3& o, const vec3& v) const override;
        virtual vec3 random(const point3& o) const override;

    private:
        static void get_sphere_uv(const point3& p, double& u, double& v) {
            // p: a given point on the sphere of radius one, centered at the origin
            // u: returned value [0,1] of angle around the Y axis from X = -1
            // v: returned value [0,1] of angle from Y = -1 to Y = +1
            // p = <1 0 0> yields u,v = <0.50 0.50>
            // p = <0 1 0> yields u,v = <0.50 1.00>
            // p = <0 0 1> yields u,v = <0.25 0.50>
            // p = <-1 0 0> yields u,v = <0.00 0.50>
            // p = <0 -1 0> yields u,v = <0.50 0.00>
            // p = <0 0 -1> yields u,v = <0.75 0.50>

            // calculate the spherical coordinates theta and phi for point p
            // theta is the angle up from the bottom pole of sphere (up from -Y)
            // phi is the angle around the Y-axis (from -X to +Z to +X to -Z back to -X)

            // start with equations for the corresponding Cartesian coordinates:
            // y = -cos(theta)
            // x = -cos(phi)sin(theta)
            // z = sin(phi)sin(theta)

            // find phi by phi = atan2(z, -x), however atan2 returns a value in the range [-pi,pi]
            // whereas we want a value in the range [0,2pi], so we will transform our result using
            // the formula atan2(a, b) = atan2(-a, -b) + pi, leading to phi = atan2(-z, x) + pi which
            // gives a value in the range [0,2pi]
            // find theta by theta = acos(-y)
            auto theta = acos(-p.y());
            auto phi = atan2(-p.z(), p.x()) + pi;

            // map theta and phi to texture coordinates u and v, each in [0,1], where 
            // (u = 0, v = 0) maps to the bottom-left corner of the texture
            u = phi / (2 * pi);  // normalize phi
            v = theta / pi;  // normalize theta
        }

    public:
        point3 center;
        double radius;
        shared_ptr<material> mat_ptr;
};

// implementation of hit
// ray intersection with a sphere
bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // ray r hits the specified sphere if there are 2 real solutions for t to this
    // quadratic equation: t^2*b*b + 2t*b*(A-C) + (A-C)*(A-C)-r^2 = 0 where b is 
    // ray direction, A is ray origin, C is sphere center, and r is sphere radius
    // we check if the determinant is positive; if so, then we have 2 real solutions for t
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();  // equivalent to dot(r.direction(), r.direction())
    auto half_b = dot(oc, r.direction());  // we use b/2 here for a simplified quadratic formula
    auto c = oc.length_squared() - radius * radius;

    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) {
        // no solutions for t
        return false;
    }
    
    auto sqrtd = sqrt(discriminant);

    // find the nearest root (closest to the center of camera) that lies in the 
    // acceptable range [t_min,t_max]
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max) {
        // this root is outside the range, try the further root
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max) {
            // no valid root in range
            return false;
        }
    }

    // populate the passed-in hit_record reference
    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;  // unit normal vector (from sphere center to ray hit position)
    rec.set_face_normal(r, outward_normal);
    get_sphere_uv(outward_normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;
    
    return true;
}

bool sphere::bounding_box(double time0, double time1, aabb& output_box) const {
    output_box = aabb(center - vec3(radius, radius, radius),
                      center + vec3(radius, radius, radius));
    return true;
}

// PDF for sampling a sphere object
double sphere::pdf_value(const point3& o, const vec3& v) const {
    hit_record rec;
    if (!this->hit(ray(o, v), 0.001, infinity, rec)) {
        return 0;
    }

    auto cos_theta_max = sqrt(1 - radius * radius / (center - o).length_squared());
    auto solid_angle = 2 * pi * (1 - cos_theta_max);

    return 1 / solid_angle;
}

// retrieve a random value based on the PDF
vec3 sphere::random(const point3& o) const {
    vec3 direction = center - o;
    auto distance_squared = direction.length_squared();
    onb uvw;
    uvw.build_from_w(direction);
    return uvw.local(random_to_sphere(radius, distance_squared));
}

#endif
