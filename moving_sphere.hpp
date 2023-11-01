#ifndef MOVING_SPHERE_H
#define MOVING_SPHERE_H

#include "constants.hpp"
#include "hittable.hpp"
#include "aabb.hpp"

// this sphere class will have its center move linearly from center0 at time0 to center1 at time1
class moving_sphere : public hittable {
    public:
        moving_sphere() {}
        moving_sphere(point3 cen0, point3 cen1, double _time0, double _time1, double r, shared_ptr<material> m)
            : center0(cen0), center1(cen1), time0(_time0), time1(_time1), radius(r), mat_ptr(m)
        {};

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

        point3 center(double time) const;

    public:
        point3 center0, center1;
        double time0, time1;
        double radius;
        shared_ptr<material> mat_ptr;
};

// calculates the position of sphere center given a time between time0 and time1
point3 moving_sphere::center(double time) const {
    return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
}

// hit method is almost identical to stationery sphere's hit, except the sphere center data will now come from the center(time) function
bool moving_sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center(r.time());
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius * radius;

    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) {
        return false;
    }

    auto sqrtd = sqrt(discriminant);

    // find the nearest root that lies in the acceptable range
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max) {
            return false;
        }
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    auto outward_normal = (rec.p - center(r.time())) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    return true;
}

// for a moving sphere, we can take the box of the sphere at time0 and the box of the
// sphere at time1, and compute the box of those two boxes
bool moving_sphere::bounding_box(double time0, double time1, aabb& output_box) const {
    // box of the sphere at time0
    aabb box0(
        center(time0) - vec3(radius, radius, radius),
        center(time0) + vec3(radius, radius, radius)
    );

    // box of the sphere at time 1
    aabb box1(
        center(time1) - vec3(radius, radius, radius),
        center(time1) + vec3(radius, radius, radius)
    );

    // compute the bounding box of box0 and box1
    output_box = surrounding_box(box0, box1);
    return true;
}

#endif
