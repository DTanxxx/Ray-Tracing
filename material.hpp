#ifndef MATERIAL_H
#define MATERIAL_H

#include "constants.hpp"
#include "hittable.hpp"
#include "texture.hpp"
#include "onb.hpp"
#include "pdf.hpp"

struct scatter_record {
    ray specular_ray;
    bool is_specular;
    color attenuation;
    shared_ptr<pdf> pdf_ptr;
};
 
// materials will tell us how rays interact with the surface
class material {
    public:
        virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec) const {
            // taking in pdf for importance sampling
            return false;
        }
        virtual double scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
            // for importance sampling
            return 0;
        }
        virtual color emitted(const ray& r_in, const hit_record& rec, double u, double v, const point3& p) const {
            // non-emitting materials will by default return black for its light color
            return color(0, 0, 0);
        }

    protected:
        // generate random directions relative to the Z-axis with cos(theta) / pi as the density function
        vec3 random_cosine_direction() const {
            auto r1 = random_double();
            auto r2 = random_double();
            auto z = sqrt(1 - r2);

            auto phi = 2 * pi * r1;
            auto x = cos(phi) * sqrt(r2);
            auto y = sin(phi) * sqrt(r2);

            return vec3(x, y, z);
        }
};

// Lambertian material - scatters incident ray to produce a diffuse effect
class lambertian : public material {
    public:
        lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {}
        lambertian(shared_ptr<texture> a) : albedo(a) {}

        virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec) const override {
            /*
            // calculate the scattered ray's direction
            // scattered ray's direction should be towards a point on the unit radius sphere surface adjacent to rec.p,
            // to achieve Lambertian effect
            // Lambertian: more uniform scattering of light rays, so fewer rays are scattering towards the normal and diffuse objects will appear lighter
            auto scatter_direction = rec.normal + random_unit_vector();

            // ALTERNATIVE to Lambertian diffuse
            // have a uniform scatter direction for all angles away from the hit pooint, with no
            // dependence on the angle from the normal
            //auto scatter_direction = random_in_hemisphere(rec.normal);

            // catch degenerate scatter direction - if the random unit vector is exactly opposite
            // to the normal vector, scatter_direction will be a zero vector, which can be problematic
            if (scatter_direction.near_zero()) {
                // make scatter_direction the same as the normal
                scatter_direction = rec.normal;
            }
            */
            /*
            //auto scattered_direction = random_in_hemisphere(rec.normal);  // random hemisphere sampling - choose randomly from the hemisphere above the surface
            onb uvw;
            uvw.build_from_w(rec.normal);
            auto scattered_direction = uvw.local(random_cosine_direction());  // scattering using orthonormal basis
            
            scattered = ray(rec.p, unit_vector(scattered_direction), r_in.time());  // create a scattered ray, taking into account the time of ray intersection
            alb = albedo->value(rec.u, rec.v, rec.p);
            //pdf = dot(rec.normal, scattered.direction()) / pi;  // for importance sampling
            //pdf = 0.5 / pi;  // PDF for random hemisphere sampling
            pdf = dot(uvw.w(), scattered.direction()) / pi;  // for orthonormal basis
            */
            srec.is_specular = false;
            srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
            srec.pdf_ptr = make_shared<cosine_pdf>(rec.normal);
            return true;  // diffuse materials can always scatter rays
        }

        double scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
            // for importance sampling, calculates the relevant PDF
            auto cosine = dot(rec.normal, unit_vector(scattered.direction()));
            return cosine < 0 ? 0 : cosine / pi;
        }
    
    public:
        shared_ptr<texture> albedo;
};

// metal material - reflects incident ray
class metal : public material {
    public:
        metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

        virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec) const override {
            vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);

            // randomize the reflected ray's direction by using a small "utility" sphere and choosing a new endpoint for the reflected ray
            // to produce a "fuzzy" reflection effect
            // the larger the utility sphere, the fuzzier the reflections will be
            srec.specular_ray = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
            srec.attenuation = albedo;
            srec.is_specular = true;
            srec.pdf_ptr = 0;
            //return (dot(scattered.direction(), rec.normal) > 0);  // return true iif successful reflection
            return true;
        }
    
    public:
        color albedo;
        double fuzz;  // fuzziness parameter
};

// dielectric material - refracts incident ray
class dielectric : public material {
    public:
        dielectric(double index_of_refraction) : ir(index_of_refraction) {}

        virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec) const override {
            srec.attenuation = color(1.0, 1.0, 1.0);  // glass surface does not absorb any colors
            srec.is_specular = true;
            srec.pdf_ptr = nullptr;

            // if ray is coming from outside the surface (rec.front_face is true), refraction
            // ratio is 1/ir, where 1 is air's refraction index and ir is the passed-in refraction
            // index of the material; otherwise the refraction ratio is ir/1 = ir
            double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

            vec3 unit_direction = unit_vector(r_in.direction());

            // determine if total internal reflection or refraction occurs
            double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);  // cos(theta) = dot(R, n)
            double sin_theta = sqrt(1.0 - cos_theta * cos_theta);  // sin(theta) = sqrt(1 - cos(theta)^2)

            // if refraction_ratio * sin(theta) > 1, there is no way sin(theta prime) can equate to this
            // because the largest value a sin function can take is 1 -> in this case, total internal reflection occurs
            bool cannot_refract = refraction_ratio * sin_theta > 1.0;
            vec3 direction;

            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double()) {
                // reflection, happens either when the incident ray comes from a medium with a higher refractive index
                // (total internal reflection) or when reflectance value is high
                direction = reflect(unit_direction, rec.normal);
            }
            else {
                // refraction
                direction = refract(unit_direction, rec.normal, refraction_ratio);
            }

            srec.specular_ray = ray(rec.p, direction, r_in.time());
            return true;
        }   

    private:
        static double reflectance(double cosine, double ref_idx) {
            // use Schlick's approximation for reflectance
            // reflectance is used to determine if the material still gives reflection even
            // if total internal reflection does not occur (i.e. Snell's law is still valid)
            // this occurs because real glass has reflectivity that varies with angle
            auto r0 = (1 - ref_idx) / (1 + ref_idx);
            r0 = r0 * r0;
            return r0 + (1 - r0) * pow((1 - cosine), 5);
        }

    public:
        double ir;  // index of refraction
};

// light emitting material
class diffuse_light : public material {
    public:
        diffuse_light(shared_ptr<texture> a) : emit(a) {}
        diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}

        virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec) const override {
            return false;
        }

        virtual color emitted(const ray& r_in, const hit_record& rec, double u, double v, const point3& p) const override {
            // light emission is unidirectional, so that the light at the ceiling only emits downwards (to remove noise at the ceiling)
            if (rec.front_face) {
                return emit->value(u, v, p);
            }
            else {
                return color(0, 0, 0); 
            }
        }

    public:
        shared_ptr<texture> emit;
};

// material used for a constant density volume (eg fog)
class isotropic : public material {
    public:
        isotropic(color c) : albedo(make_shared<solid_color>(c)) {}
        isotropic(shared_ptr<texture> a) : albedo(a) {}

        virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec) const override {
            // pick an uniform random direction for our scattered ray
            srec.specular_ray = ray(rec.p, random_in_unit_sphere(), r_in.time());
            srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
            return true;
        }

    public:
        shared_ptr<texture> albedo;
};

#endif
