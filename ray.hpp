#ifndef RAY_H
#define RAY_H

#include "vec3.hpp"

class ray {
    public:
        // constructors
        ray() {}
        ray(const point3& origin, const vec3& direction, double time = 0.0)
            : orig(origin), dir(direction), tm(time) 
        {}  // a ray has a starting point, direction, and the time it exists at (for motion blur effect)

        // getters
        point3 origin() const { return orig; }
        vec3 direction() const { return dir; }
        double time() const { return tm; }

        // retrieves the position at the end of ray
        point3 at(double t) const {
            return orig + dir * t;
        }

    public:
        point3 orig;
        vec3 dir;
        double tm;  // for motion blur - general idea: generate rays at random times while the shutter is open
                    // and intersect the model at that one time. This is done by having the camera move and the 
                    // objects move, but have each ray exist at exactly one time
};

#endif
