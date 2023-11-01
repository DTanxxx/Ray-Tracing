#ifndef CONSTANT_MEDIUM_H
#define CONSTANT_MEDIUM_H

#include "constants.hpp"
#include "hittable.hpp"
#include "material.hpp"
#include "texture.hpp"

// create a region of fog (aka volume) using subsurface scattering
// a volume can be created using a random surface, where at every point in the volume
// the surface probabilistically might or might not be there
// we will start with a volume of constant density - a ray going through there can either
// scatter inside the volume, or it can make it all the way through
// more thin transparent volumes (eg light fog) are more likely to have rays going all the
// way through, so the denser the volume, the more likely a ray will scatter; the probability
// that the ray scatters in any small distance delta L is C * delta L where C is proportional
// to the optical density of the volume
// if the distance needed for scattering is outside the volume, then there is no "hit"
class constant_medium : public hittable {
    public:
        constant_medium(shared_ptr<hittable> b, double d, shared_ptr<texture> a) 
            : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<isotropic>(a)) {}
        constant_medium(shared_ptr<hittable> b, double d, color c)
            : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<isotropic>(c)) {}

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
            return boundary->bounding_box(time0, time1, output_box);
        }
    
    public:
        shared_ptr<hittable> boundary;
        double neg_inv_density;
        shared_ptr<material> phase_function;
};

bool constant_medium::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    // print occasional samples when debugging - to enable, set enableDebug true
    const bool enableDebug = false;
    const bool debugging = enableDebug && random_double() < 0.00001;

    hit_record rec1, rec2;

    // we need to be careful about the logic around the boundary because we need to 
    // make sure this works for ray origins inside the volume (since rays would bounce
    // around a lot)
    if (!boundary->hit(r, -infinity, infinity, rec1)) {
        return false;
    }

    if (!boundary->hit(r, rec1.t + 0.0001, infinity, rec2)) {
        return false;
    }

    if (debugging) {
        std::cerr << "\nt_min=" << rec1.t << ", t_max=" << rec2.t << '\n';
    }

    if (rec1.t < t_min) {
        rec1.t = t_min;
    }

    if (rec2.t > t_max) {
        rec2.t = t_max;
    }

    if (rec1.t >= rec2.t) {
        return false;
    }

    if (rec1.t < 0) {
        rec1.t = 0;
    }

    // we assume that once a ray exits the constant medium boundary, it will continue
    // forever outside the boundary - it assumes that the boundary shape is convex (so
    // this implementation will work for boundaries like boxes or spheres, but will not
    // work with toruses or shapes that contain voids)
    const auto ray_length = r.direction().length();
    const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
    const auto hit_distance = neg_inv_density * log(random_double());

    if (hit_distance > distance_inside_boundary) {
        return false;
    }

    rec.t = rec1.t + hit_distance / ray_length;
    rec.p = r.at(rec.t);

    if (debugging) {
        std::cerr << "hit_distance = " << hit_distance << '\n'
                  << "rec.t = " << rec.t << '\n'
                  << "rec.p = " << rec.p << '\n';
    }

    rec.normal = vec3(1, 0, 0);  // arbitrary
    rec.front_face = true;  // also arbitrary
    rec.mat_ptr = phase_function;

    return true;
}

#endif
