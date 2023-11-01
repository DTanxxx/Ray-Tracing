#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.hpp"
#include <memory>
#include <vector>
#include "aabb.hpp"

using std::shared_ptr;
using std::make_shared;

class hittable_list : public hittable {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) {
            add(object);
        }

        void clear() {
            objects.clear();
        }

        void add(shared_ptr<hittable> object) {
            objects.push_back(object);
        }

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;
        virtual double pdf_value(const point3& o, const vec3& v) const override;
        virtual vec3 random(const vec3& o) const override;

    public:
        // a list of hittable objects where each hittable object is referenced by a 
        // shared_ptr which has reference-counting semantics -> helps freeing pointers automatically
        std::vector<shared_ptr<hittable>> objects;
};

// ray intersection with the hittable that is as close as possible to the camera
bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    // (assume object 1 and object 2 overlap at where the ray points)
    // for object 1 that is closer to camera than object 2, after the loop reaches object 1,
    // the temp_rec.t will be small enough such that when the loop reaches object 2 later on, the
    // closest_so_far value will prohibit the ray from hitting object 2, hence rendering
    // only the closest hittable object (object 1) for this given ray
    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

bool hittable_list::bounding_box(double time0, double time1, aabb& output_box) const {
    if (objects.empty()) {
        // no hittable objects, no bounding box
        return false;
    }

    aabb temp_box;
    bool first_box = true;

    for (const auto& object : objects) {
        // try to compute a bounding box for the given hittable object
        if (!object->bounding_box(time0, time1, temp_box)) {
            // unable to compute bounding box, return false
            return false;
        }

        // bounding box successfully created, compute a bigger surrounding bounding box
        // enclosing the existing output_box and the new temp_box
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }

    return true;
}

// create a mixture of PDF densities from the objects list
double hittable_list::pdf_value(const point3& o, const vec3& v) const {
    auto weight = 1.0 / objects.size();
    auto sum = 0.0;

    for (const auto& object : objects) {
        sum += weight * object->pdf_value(o, v);
    }

    return sum;
}

// retrieve a random value based on the PDF mixture
vec3 hittable_list::random(const vec3& o) const {
    auto int_size = static_cast<int>(objects.size());
    return objects[random_int(0, int_size - 1)]->random(o);
}

#endif
