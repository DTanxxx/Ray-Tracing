/*
This project is from Raytracing In One Weekend Series : Book 1-3
*/

#include "camera.hpp"
#include "hittable_list.hpp"
#include "material.hpp"
#include "sphere.hpp"
#include "constants.hpp"
#include "color.hpp"
#include "moving_sphere.hpp"
#include "bvh.hpp"
#include "aarect.hpp"
#include "box.hpp"
#include "constant_medium.hpp"

#include <iostream>

// returns a color to render the pixel for given a ray and a world containing hittable objects
// if no hit: linearly blends white and blue depending on the height of the y coordinate after scaling the ray direction to unit length
// if hit: recursively generate reflected rays and factor in all their resulting colors
color ray_color(const ray& r, const color& background, const hittable& world, const shared_ptr<hittable>& lights, int depth) {
    hit_record rec;

    // if we have exceeded the ray bounce limit, no more light is gathered
    if (depth <= 0) {
        // return black color
        return color(0, 0, 0);
    }

    // OLD RAY HIT LOGIC FOR SKYBOX BACKGROUND
    /*
    // give our ray an allowed range of [0.001,infinity] so that hits very close to the ray origin are ignored
    if (world.hit(r, 0.001, infinity, rec)) {
        // if the ray hits the hittable world (t's range is [0.001,infinity]), create a scattered
        // ray and factor in its color for scatter effect
        ray scattered;  // out parameter
        color attenuation;  // out parameter
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            // if ray hits a hittable that can scatter rays, recursively calculate and apply new ray colors with an attenuation
            return attenuation * ray_color(scattered, world, depth - 1);
        }

        // hittable does not have a material that scatters rays, return black color (fully absorbed)
        return color(0, 0, 0);

        // LEGACY (uses normal vector to calculate a color)
        //return 0.5 * (rec.normal + color(1, 1, 1));
    }

    vec3 unit_direction = unit_vector(r.direction());  // unit_direction.y() is in range [-1,1]

    // graphics trick: scale unit_direction.y() from [-1,1] to [0,1]
    auto t = 0.5 * (unit_direction.y() + 1.0);

    // when t = 1, want blue; when t = 0, want white
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
    */

    // NEW RAY HIT LOGIC FOR SOLID COLOR BACKGROUND
    // we would like the only light in the scene to come from the emitters

    // if the ray hits nothing, return the background color
    if (!world.hit(r, 0.001, infinity, rec)) {
        return background;
    }

    scatter_record srec;
    color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, srec)) {
        // hittable does not have a material that scatters rays, so return its emitted
        // color (usually just black i.e. non-light-emitting, but for light emitting 
        // materials this may be a color RGB value greater than 1 - the higher this 
        // value, the brighter the emitted light will be)
        return emitted;
    }

    if (srec.is_specular) {
        // generate an implicitly sampled ray
        return srec.attenuation * ray_color(srec.specular_ray, background, world, lights, depth - 1);
    }

    /*
    // add light sampling - we want to send more rays to the light by picking a random point on the 
    // light and sending a ray in that direction, because otherwise just by sampling uniformly over
    // all directions the lights are not sampled any more than unimportant directions (which leads to noise)
    auto on_light = point3(random_double(213, 343), 554, random_double(227, 332));  // random point on light
    auto to_light = on_light - rec.p;  // direction to that point
    auto distance_squared = to_light.length_squared();
    to_light = unit_vector(to_light);

    if (dot(to_light, rec.normal) < 0) {
        // scattered ray is unable to hit the light because the angle between the direction to light and surface normal is obtuse
        return emitted;
    }

    double light_area = (343 - 213) * (332 - 227);
    auto light_cosine = fabs(to_light.y());
    if (light_cosine < 0.000001) {
        // light is almost at the same horizontal level as the incident ray hit point/scattered ray origin, so scattered ray 
        // is unable to hit the light
        return emitted;
    }

    pdf_val = distance_squared / (light_cosine * light_area);  // PDF(direction) = distance(incident ray hit point, light)^2 / (light_area * light_cosine)
    scattered = ray(rec.p, to_light, r.time());
    */
    /*
    // cosine density PDF
    cosine_pdf p(rec.normal);
    scattered = ray(rec.p, p.generate(), r.time());
    pdf_val = p.value(scattered.direction());
    */
    /*
    // hittable PDF (for light sampling)
    hittable_pdf light_pdf(lights, rec.p);
    scattered = ray(rec.p, light_pdf.generate(), r.time());
    pdf_val = light_pdf.value(scattered.direction());
    */
    // combine the cosine density PDF and hittable PDF using a mixture PDF
    auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
    //auto p1 = make_shared<cosine_pdf>(rec.normal);
    mixture_pdf mixed_pdf(light_ptr, srec.pdf_ptr);

    ray scattered = ray(rec.p, mixed_pdf.generate(), r.time());
    auto pdf_val = mixed_pdf.value(scattered.direction());

    // note that the "emitted" value is directly added onto the scattered ray's color,
    // to build up "brightness"
    return emitted + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered) *
           ray_color(scattered, background, world, lights, depth - 1) / pdf_val;  // for importance sampling
}
/*
// generate a scene with random spheres
hittable_list random_scene() {
    hittable_list world;

    // render a sphere that acts as the ground (larger and lower on the screen)
    // note that geometrically, the smaller spheres are closer to the camera than the ground sphere is
    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

    // render a lot of radius = 0.2 spheres
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();  // RNG for selecting which material to apply
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());  // random spawn point

            // decide whether the generated spawn position is appropriate
            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse material
                    auto albedo = color::random() * color::random();  // random albedo
                    sphere_material = make_shared<lambertian>(albedo);

                    // create moving spheres that move from its center C at time t=0 to C+(0,r/2,0) at time t=1, where r is a random number in [0,1)
                    auto center2 = center + vec3(0, random_double(0, 0.5), 0);
                    world.add(make_shared<moving_sphere>(center, center2, 0.0, 1.0, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal material
                    auto albedo = color::random(0.5, 1);  // random bounded albedo
                    auto fuzz = random_double(0, 0.5);  // random bounded fuzz
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass material
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    // render three radius = 1.0 spheres
    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));
    
    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

// a scene with two checkered spheres
hittable_list two_spheres() {
    hittable_list objects;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));

    objects.add(make_shared<sphere>(point3(0, -10, 0), 10, make_shared<lambertian>(checker)));
    objects.add(make_shared<sphere>(point3(0, 10, 0), 10, make_shared<lambertian>(checker)));

    return objects;
}

// a scene with two perlin-noise spheres
hittable_list two_perlin_spheres() {
    hittable_list objects;

    auto perlin_texture = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(perlin_texture)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(perlin_texture)));

    return objects;
}

// a scene with a sphere that uses an image texture
hittable_list earth() {
    auto earth_texture = make_shared<image_texture>("earthmap.jpg");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);

    return hittable_list(globe);
}

// a scene with an emissive rectangle
hittable_list simple_light() {
    hittable_list objects;

    auto perlin_texture = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(perlin_texture)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(perlin_texture)));

    auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
    objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));

    return objects;
}
*/
// a Cornell Box scene containing 5 walls and a light on the ceiling
hittable_list cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(0.65, 0.05, 0.05));
    auto white = make_shared<lambertian>(color(0.73, 0.73, 0.73));
    auto green = make_shared<lambertian>(color(0.12, 0.45, 0.15));

    auto light = make_shared<diffuse_light>(color(15, 15, 15));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    /*
    // add two rotated boxes into the scene
    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);  // rotate by 15 degrees counterclockwise around Y axis
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));  // translate by (265, 0, 295)
    objects.add(box1);

    shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));
    objects.add(box2);
    */
    shared_ptr<material> aluminium = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), aluminium);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    objects.add(box1);

    auto glass = make_shared<dielectric>(1.5);
    objects.add(make_shared<sphere>(point3(190, 90, 190), 90, glass));

    return objects;
}
/*
// a Cornell Box scene similar to above, but with a block of smoke and a block of fog
hittable_list cornell_smoke() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(0.65, 0.05, 0.05));
    auto white = make_shared<lambertian>(color(0.73, 0.73, 0.73));
    auto green = make_shared<lambertian>(color(0.12, 0.45, 0.15));

    auto light = make_shared<diffuse_light>(color(7, 7, 7));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xz_rect>(113, 443, 127, 432, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));

    shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));

    // box1 acts as the boundary for smoke volume (dark albedo)
    objects.add(make_shared<constant_medium>(box1, 0.01, color(0, 0, 0)));

    // box2 acts as the boundary for fog volume (light albedo)
    objects.add(make_shared<constant_medium>(box2, 0.01, color(1, 1, 1)));

    return objects;
}

// a scene containing a big thin mist covering everything, and a blue subsurface reflection
// sphere (which is essentially a dielectric sphere with a volume inside)
hittable_list final_scene() {
    shared_ptr<hittable_list> boxes1 = make_shared<hittable_list>();
    auto ground = make_shared<lambertian>(color(0.48, 0.83, 0.53));

    // create a bunch of box geometries for an uneven ground
    const int boxes_per_side = 20;
    for (int i = 0; i < boxes_per_side; ++i) {
        for (int j = 0; j < boxes_per_side; ++j) {
            auto w = 100.0;
            auto x0 = -1000.0 + i * w;
            auto z0 = -1000.0 + j * w;
            auto y0 = 0.0;
            auto x1 = x0 + w;
            auto y1 = random_double(1, 101);
            auto z1 = z0 + w;

            boxes1->add(make_shared<box>(point3(x0, y0, z0), point3(x1, y1, z1), ground));
        }
    }

    hittable_list objects;

    objects.add(boxes1);

    // light-emitting rectangle
    auto light = make_shared<diffuse_light>(color(7, 7, 7));
    objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));

    // moving sphere
    auto center1 = point3(400, 400, 200);
    auto center2 = center1 + vec3(30, 0, 0);
    auto moving_sphere_material = make_shared<lambertian>(color(0.7, 0.3, 0.1));
    objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50, moving_sphere_material));

    objects.add(make_shared<sphere>(point3(260, 150, 45), 50, make_shared<dielectric>(1.5)));
    objects.add(make_shared<sphere>(point3(0, 150, 145), 50, make_shared<metal>(color(0.8, 0.8, 0.9), 1.0)));

    // construct the blue subsurface reflection sphere
    // start by creating the dielectric sphere boundary
    auto boundary = make_shared<sphere>(point3(360, 150, 145), 70, make_shared<dielectric>(1.5));
    objects.add(boundary);

    // now add a constant density volume inside that boundary
    objects.add(make_shared<constant_medium>(boundary, 0.2, color(0.2, 0.4, 0.9)));

    // construct the big thin mist that occupies the entire scene
    // start again by creating the dielectric sphere boundary, with a large radius
    boundary = make_shared<sphere>(point3(0, 0, 0), 5000, make_shared<dielectric>(1.5));

    // add a constant density volume representing the mist, inside that boundary
    objects.add(make_shared<constant_medium>(boundary, 0.0001, color(1, 1, 1)));

    // sphere with Earth image texture 
    auto emat = make_shared<lambertian>(make_shared<image_texture>("earthmap.jpg"));
    objects.add(make_shared<sphere>(point3(400, 200, 400), 100, emat));

    // perlin noise sphere
    auto perlin_texture = make_shared<noise_texture>(0.1);
    objects.add(make_shared<sphere>(point3(220, 280, 300), 80, make_shared<lambertian>(perlin_texture)));

    shared_ptr<hittable_list> boxes2 = make_shared<hittable_list>();
    auto white = make_shared<lambertian>(color(0.73, 0.73, 0.73));
    int ns = 1000;

    // create a clustered bunch of sphere geometries 
    for (int j = 0; j < ns; ++j) {
        boxes2->add(make_shared<sphere>(point3::random(0, 165), 10, white));
    }

    // translate these sphere geometries by (-100, 270, 395) and rotate them by 15 degrees counterclockwise about Y
    objects.add(make_shared<translate>(make_shared<rotate_y>(boxes2, 15), vec3(-100, 270, 395)));

    return objects;
}

int main() {
    // Image
    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 400;
    int samples_per_pixel = 100;  // for antialiasing
    const int max_depth = 50;  // for ray reflection

    // World
    hittable_list world_list;

    point3 lookfrom;
    point3 lookat;
    auto vfov = 40.0;
    auto aperture = 0.0;
    color background(0, 0, 0);  // by default, background color is black

    // scene selection
    switch (0) {
        case 1:
            // create a random scene with lots of spheres
            world_list = random_scene();

            // set a background color
            background = color(0.70, 0.80, 1.00);

            // the camera resides at (13, 2, 3) facing towards (0, 0, 0)
            lookfrom = point3(13, 2, 3);
            lookat = point3(0, 0, 0);

            // the camera has a 20 degrees vertical field-of-view
            vfov = 20.0;

            // the aperture size is 0.1
            aperture = 0.1;
            break;
        case 2:
            // create a scene with two checkered spheres
            world_list = two_spheres();
            background = color(0.70, 0.80, 1.00);
            lookfrom = point3(13, 2, 3);
            lookat = point3(0, 0, 0);
            vfov = 20.0;
            break;
        case 3:
            world_list = two_perlin_spheres();
            background = color(0.70, 0.80, 1.00);
            lookfrom = point3(13, 2, 3);
            lookat = point3(0, 0, 0);
            vfov = 20.0;
            break;
        case 4:
            world_list = earth();
            background = color(0.70, 0.80, 1.00);
            lookfrom = point3(13, 2, 3);
            lookat = point3(0, 0, 0);
            vfov = 20.0;
            break;
        case 5:
            world_list = simple_light();
            samples_per_pixel = 400;
            background = color(0.0, 0.0, 0.0);
            lookfrom = point3(26, 3, 6);
            lookat = point3(0, 2, 0);
            vfov = 20.0;
            break;
        case 6:
            world_list = cornell_box();
            aspect_ratio = 1.0;
            image_width = 600;
            samples_per_pixel = 200;
            background = color(0, 0, 0);
            lookfrom = point3(278, 278, -800);
            lookat = point3(278, 278, 0);
            vfov = 40.0;
            break;
        case 7:
            world_list = cornell_smoke();
            aspect_ratio = 1.0;
            image_width = 600;
            samples_per_pixel = 200;
            lookfrom = point3(278, 278, -800);
            lookat = point3(278, 278, 0);
            vfov = 40.0;
            break;
        default:
        case 8:
            world_list = final_scene();
            aspect_ratio = 1.0;
            image_width = 800;
            samples_per_pixel = 10000;
            background = color(0, 0, 0);
            lookfrom = point3(478, 278, -600);
            lookat = point3(278, 278, 0);
            vfov = 40.0;
            break;
    }

    bvh_node world = bvh_node(world_list, 0.0, 1.0);

    // initialise some materials that we are going to use
    //auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    //auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    //auto material_left = make_shared<dielectric>(1.5);
    //auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 0.0);

    // render a sphere located on the center of viewport plane with radius 0.5
    //world.add(make_shared<sphere>(point3(0.0, 0.0, -1.0), 0.5, material_center));

    // render a dielectric sphere with a negative radius - the geometry is unaffected but the surface normal points inward, 
    // which "hollows-out" the spherical region; if there's a larger sphere at the same position as this one then
    // the resulting render will produce a hollow sphere
    //world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), -0.45, material_left));

    //world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
    //world.add(make_shared<sphere>(point3(1.0, 0.0, -1.0), 0.5, material_right));
    
    // Camera
    // the camera has vup vector (0, 1, 0) and the focus distance is 10.0
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // camera will open its shutter at time 0 and close at time 1 (for motion blur)
    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

    // Render
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j) {
        // display progress
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            
            // send multiple rays to a given pixel position
            // remember that each pixel is a distance of "1" away from each other, therefore
            // by appending a decimal from [0,1) to the original position, the original pixel is 
            // essentially "broken up" into smaller pixels
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);

                // accumulate sample colors
                // the rays are essentially hitting "smaller" broken-up pixels from the original
                // pixel, allowing different hit results to occur on the original pixel; namely,
                // sometimes the rays will hit a smaller pixel that is outside of the sphere
                // geometry, resulting in a "no-hit"; other times, the rays will hit a smaller 
                // pixel that is within the sphere boundary, thus resuling in a "hit".
                // To be concrete, the hit_record.t values may be different between each
                // ray sample - when the smaller pixels are outside the sphere boundary, we won't
                // have any roots and we will receive a background color sample; when the smaller
                // pixels are within the boundary, we will have some values for t which means that
                // our sample color will be the sphere color (depending on the unit normal vector)
                pixel_color += ray_color(r, background, world, max_depth);
            }

            // average the pixel_color from multiple samples
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
    
    return 0;
}
*/

int main() {
    const auto aspect_ratio = 1.0 / 1.0;
    const int image_width = 600;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 1000;
    const int max_depth = 50;

    auto world = cornell_box();

    // use a mixture of glass and light PDFs in our scene
    shared_ptr<hittable_list> lights = make_shared<hittable_list>();
    lights->add(make_shared<sphere>(point3(190, 90, 190), 90, make_shared<material>()));
    lights->add(make_shared<xz_rect>(213, 343, 227, 332, 554, make_shared<material>()));

    color background(0, 0, 0);

    point3 lookfrom(278, 278, -800);
    point3 lookat(278, 278, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.0;
    auto vfov = 40.0;
    auto time0 = 0.0;
    auto time1 = 1.0;

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, time0, time1);

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, lights, max_depth);
            }

            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
    
    return 0;
}