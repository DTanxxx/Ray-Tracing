#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.hpp"
#include "constants.hpp"
#include "aabb.hpp"

class material;  // forward declaration

struct hit_record {
    point3 p;  // ray hit point
    vec3 normal;
    shared_ptr<material> mat_ptr;
    double t;
    double u;  // for texture coordinates of the ray-object hit point
    double v;  // for texture coordinates of the ray-object hit point
    bool front_face;  // whether the ray comes from the outside (true) or inside (false) the surface

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        // if the ray is in the opposite direction as outward_normal, it is from outside
        front_face = dot(r.direction(), outward_normal) < 0;

        // determine the normal vector relative to this ray; if ray comes from outside
        // then we keep outward_normal as it is, otherwise ray is from inside and 
        // the normal to it is the negation of outward_normal
        normal = front_face ? outward_normal : -outward_normal;
    }
};

// an abstract class for anything a ray might hit
class hittable {
    public:
        // we want to specify an interval along the ray that are valid for hits [t_min,t_max], so
        // hits only count if t_min < t < t_max
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;

        // compute the bounding box for a hittable object; return a bool indicating whether a bounding box is 
        // constructed successfully for a primitive geometry
        // moving objects will have a bounding box that encloses the object for the entire time interval [time0,time1]
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const = 0;

        virtual double pdf_value(const point3& o, const vec3& v) const {
            return 0.0;
        }

        virtual vec3 random(const vec3& o) const {
            return vec3(1, 0, 0);
        }
};

// a translate instance: a geometric primitive that has been moved
// in ray tracing, instead of moving the primitive itself during translation, we move
// the rays in the opposite direction
class translate : public hittable {
    public:
        translate(shared_ptr<hittable> p, const vec3& displacement) : ptr(p), offset(displacement) {}

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

    public:
        shared_ptr<hittable> ptr;  // the object to be translated
        vec3 offset;  // translation offset
};

bool translate::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // move rays in opposite direction to translation offset
    ray moved_r(r.origin() - offset, r.direction(), r.time());
    if (!ptr->hit(moved_r, t_min, t_max, rec)) {
        return false;
    }

    // restore the hit point 
    rec.p += offset;
    rec.set_face_normal(moved_r, rec.normal);

    return true;
}

bool translate::bounding_box(double time0, double time1, aabb& output_box) const {
    if (!ptr->bounding_box(time0, time1, output_box)) {
        return false;
    }

    // translate the bounding box by the offset
    output_box = aabb(output_box.min() + offset, output_box.max() + offset);

    return true;
}

// a rotate Y instance: a geometric primitive that has been rotated counterclockwise around the Y axis
class rotate_y : public hittable {
    public:
        rotate_y(shared_ptr<hittable> p, double angle);

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
            output_box = bbox;
            return hasBox;
        }

    public:
        shared_ptr<hittable> ptr;
        double sin_theta;
        double cos_theta;
        bool hasBox;
        aabb bbox;
};

// create the bounding box for the rotated primitive object
rotate_y::rotate_y(shared_ptr<hittable> p, double angle) : ptr(p) {
    auto radians = degrees_to_radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasBox = ptr->bounding_box(0, 1, bbox);

    point3 min(infinity, infinity, infinity);
    point3 max(-infinity, -infinity, -infinity);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                auto x = i * bbox.max().x() + (1 - i) * bbox.min().x();
                auto y = j * bbox.max().y() + (1 - j) * bbox.min().y();
                auto z = k * bbox.max().z() + (1 - k) * bbox.min().z();

                // the new X and Z coordinates after rotating counterclockwise by angle theta about Y is:
                // x' = cos(theta) * x + sin(theta) * z
                // z' = -sin(theta) * x + cos(theta) * z
                auto newX = cos_theta * x + sin_theta * z;
                auto newZ = -sin_theta * x + cos_theta * z;

                vec3 tester(newX, y, newZ);

                for (int c = 0; c < 3; ++c) {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
    }
    
    bbox = aabb(min, max);
}

bool rotate_y::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    auto origin = r.origin();
    auto direction = r.direction();

    origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
    origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];

    direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
    direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];

    ray rotated_r(origin, direction, r.time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec)) {
        return false;
    }

    auto p = rec.p;
    auto normal = rec.normal;

    p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
    p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];

    // the surface normal vector changes with rotation, so we need to transform its direction using 
    // the same formula as before when calculating the bounding box
    normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
    normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];

    rec.p = p;
    rec.set_face_normal(rotated_r, normal);

    return true;
}

// flip a given hittable's front face (in our case, we want to flip the light rectangle's front face downwards since it only emits light
// on the front face side, to support unidirectional light emission; originally, the front face is in the +y direction, but we want 
// the light to emit downwards so we flip front face to -y direction)
class flip_face : public hittable {
    public:
        flip_face(shared_ptr<hittable> p) : ptr(p) {}

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
            if (!ptr->hit(r, t_min, t_max, rec)) {
                return false;
            }

            rec.front_face = !rec.front_face;
            return true;
        }

        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
            return ptr->bounding_box(time0, time1, output_box);
        }

    public:
        shared_ptr<hittable> ptr;
};

#endif
