#ifndef CAMERA_H
#define CAMERA_H

#include "constants.hpp"

class camera {
    public:
        camera(point3 lookfrom,  // the position where we place the camera
               point3 lookat,  // the point we look at
               vec3 vup,  // an arbitrary world "up" vector, typically (0, 1, 0)
               double vfov,  // vertical field-of-view in degrees
               double aspect_ratio,
               double aperture,
               double focus_dist,  // focus distance is the distance between the projection point and the plane where everything is in perfect focus
               double _time0 = 0,  // camera will generate rays at a random time between _time0 and _time1 (for motion blur)
               double _time1 = 0) {  
            auto theta = degrees_to_radians(vfov);
            auto h = tan(theta / 2);  // the h value depends on the distance between the projection point and the projection plane
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;
            
            w = unit_vector(lookfrom - lookat);  // vector perpendicular to the camera plane; acts as the "z direction" for the camera
            u = unit_vector(cross(vup, w));  // a basis vector for the camera plane
            v = cross(w, u);  // another basis vector for the camera plane

            origin = lookfrom;  // "eye" of the camera
            horizontal = focus_dist * viewport_width * u;  // positive x is rightwards
            vertical = focus_dist * viewport_height * v;  // positive y is upwards
            
            // negative z is "into the screen"
            // here we find the coordinate of bottom left corner of the viewport plane
            lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w;

            lens_radius = aperture / 2;
            time0 = _time0;
            time1 = _time1;
        }
        
        // create a ray with origin (0,0,0) and direction that goes to each pixel position on viewport plane
        ray get_ray(double s, double t) const {
            // to achieve defocus blur effect: generate random scene rays originating from inside
            // a disk centered at the lookfrom point - the larger the radius of the disk, the greater
            // the defocus blur
            vec3 rd = lens_radius * random_in_unit_disk();  // scale disk radius by the aperture size
            vec3 offset = u * rd.x() + v * rd.y();

            return ray(origin + offset,
                       lower_left_corner + s * horizontal + t * vertical - origin - offset,
                       random_double(time0, time1));
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        double lens_radius;
        double time0, time1;  // shutter open/close times
};

#endif
