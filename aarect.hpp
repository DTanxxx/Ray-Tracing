#ifndef AARECT_H
#define AARECT_H

#include "constants.hpp"
#include "hittable.hpp"

// this class represents a rectangle in the xy plane, defined by z = k and x = x0, x = x1, y = y0, y = y1 (four edges)
class xy_rect : public hittable {
    public:
        xy_rect() {}
        xy_rect(double _x0, double _x1, double _y0, double _y1, double _k, shared_ptr<material> mat)
            : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {};

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
            // the bounding box must have non-zero width in each dimension, so pad the Z dimension
            // a small amount
            output_box = aabb(point3(x0, y0, k - 0.0001), point3(x1, y1, k + 0.0001));
            return true;
        }

    public:
        shared_ptr<material> mp;
        double x0, x1, y0, y1, k;
};

bool xy_rect::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // to find out whether a ray hits the rectangle, we need to first determine where the
    // ray hits the plane
    // looking at the Z coordinate only, we know that P_z(t) = A_z + t * b_z must be true
    // if a ray hits the rectangle, therefore we can calculate t where z = k with this
    // rearranged version of the equation: t = (k - A_z) / b_z
    auto t = (k - r.origin().z()) / r.direction().z();

    // ensure that our ray is actually able to capture this t value
    if (t < t_min || t > t_max) {
        return false;
    }

    // plug t into the equivalent ray equations for x and y coordinates: P_x(t) = A_x + t * b_x and P_y(t) = A_y + t * b_y
    auto x = r.origin().x() + t * r.direction().x();
    auto y = r.origin().y() + t * r.direction().y();

    // make sure both x and y values lie within the bounding axes intervals
    if (x < x0 || x > x1 || y < y0 || y > y1) {
        return false;
    }

    // at this point, we know the ray has hit our rectangle; populate rec accordingly
    rec.u = (x - x0) / (x1 - x0);  // fractional x position
    rec.v = (y - y0) / (y1 - y0);  // fractional y position
    rec.t = t;

    auto outward_normal = vec3(0, 0, 1);  // right hand rule is used
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mp;
    rec.p = r.at(t);
    return true;
}

// this class represents a rectangle in the xz plane, defined by y = k and x = x0, x = x1, z = z0, z = z1 (four edges)
class xz_rect : public hittable {
    public:
        xz_rect() {}
        xz_rect(double _x0, double _x1, double _z0, double _z1, double _k, shared_ptr<material> mat)
            : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
            // the bounding box must have non-zero width in each dimension, so pad the Y dimension
            // a small amount
            output_box = aabb(point3(x0, k - 0.0001, z0), point3(x1, k + 0.0001, z1));
            return true;
        }

        virtual double pdf_value(const point3& origin, const vec3& v) const override {
            hit_record rec;
            if (!this->hit(ray(origin, v), 0.001, infinity, rec)) {
                return 0;
            }

            auto area = (x1 - x0) * (z1 - z0);
            auto distance_squared = rec.t * rec.t * v.length_squared();
            auto cosine = fabs(dot(v, rec.normal) / v.length());

            return distance_squared / (cosine * area);
        }

        virtual vec3 random(const point3& origin) const override {
            auto random_point = point3(random_double(x0, x1), k, random_double(z0, z1));
            return random_point - origin;
        }

    public:
        shared_ptr<material> mp;
        double x0, x1, z0, z1, k;
};

bool xz_rect::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // similar logic to xy_rect's hit()
    auto t = (k - r.origin().y()) / r.direction().y();

    if (t < t_min || t > t_max) {
        return false;
    }

    auto x = r.origin().x() + t * r.direction().x();
    auto z = r.origin().z() + t * r.direction().z();

    if (x < x0 || x > x1 || z < z0 || z > z1) {
        return false;
    }

    rec.u = (x - x0) / (x1 - x0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;

    auto outward_normal = vec3(0, 1, 0);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mp;
    rec.p = r.at(t);
    return true;
}

// this class represents a rectangle in the yz plane, defined by x = k and y = y0, y = y1, z = z0, z = z1 (four edges)
class yz_rect : public hittable {
    public:
        yz_rect() {}
        yz_rect(double _y0, double _y1, double _z0, double _z1, double _k, shared_ptr<material> mat)
            : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
            // the bounding box must have non-zero width in each dimension, so pad the X dimension
            // a small amount
            output_box = aabb(point3(k - 0.0001, y0, z0), point3(k + 0.0001, y1, z1));
            return true;
        }

    public:
        shared_ptr<material> mp;
        double y0, y1, z0, z1, k;
};

bool yz_rect::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // similar logic to xy_rect's hit()
    auto t = (k - r.origin().x()) / r.direction().x();

    if (t < t_min || t > t_max) {
        return false;
    }

    auto y = r.origin().y() + t * r.direction().y();
    auto z = r.origin().z() + t * r.direction().z();

    if (y < y0 || y > y1 || z < z0 || z > z1) {
        return false;
    }

    rec.u = (y - y0) / (y1 - y0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;

    auto outward_normal = vec3(1, 0, 0);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mp;
    rec.p = r.at(t);
    return true;
}

#endif 
