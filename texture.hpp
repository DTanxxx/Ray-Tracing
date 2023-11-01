#ifndef TEXTURE_H
#define TEXTURE_H

#include "constants.hpp"
#include "perlin.hpp"
#include "stb_image_wrapper.hpp"

#include <iostream>

// a texture is a function that makes the colors on a surface procedural - this procedural
// can be synthesis code or an image lookup (or both)
class texture {
    public:
        virtual color value(double u, double v, const point3& p) const = 0;
};

// make colors a texture, called a "constant texture"
class solid_color : public texture {
    public:
        solid_color() {}
        solid_color(color c) : color_value(c) {}
        solid_color(double red, double green, double blue) : solid_color(color(red, green, blue)) {}

        virtual color value(double u, double v, const vec3& p) const override {
            return color_value;
        }

    private:
        color color_value;
};

// checker texture: since the sign of sine and cosine alternates in a regular way, if we
// multiply trig functions in all three dimensions, the sign of that product forms a 3D
// checker pattern
class checker_texture : public texture {
    public:
        checker_texture() {}
        checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd)
            : even(_even), odd(_odd) {}
        checker_texture(color c1, color c2) 
            : even(make_shared<solid_color>(c1)), odd(make_shared<solid_color>(c2)) {}

        virtual color value(double u, double v, const point3& p) const override {
            auto sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());
            if (sines < 0) {
                return odd->value(u, v, p);
            }
            else {
                return even->value(u, v, p);
            }
        }
    
    public:
        shared_ptr<texture> odd;
        shared_ptr<texture> even;
};

// solid texture with perlin noise applied
class noise_texture : public texture {
    public:
        noise_texture() {}
        noise_texture(double sc) : scale(sc) {}

        virtual color value(double u, double v, const point3& p) const override {
            // since the output of the perlin interpolation can return negative values,
            // these negative values will be passed to the sqrt() function of our gamma
            // function and get turned into NaNs
            // we will need to cast the perlin output back to between 0 and 1
            //return color(1, 1, 1) * 0.5 * (1.0 + noise.noise(scale * p));

            // turbulence gives a sort of camouflage netting appearance when used directly
            //return color(1, 1, 1) * noise.turbulence(scale * p);

            // to produce a marble-like texture, we can use turbulence indirectly
            // in order to create this effect, we need to make color proportional to something
            // like a sine function, and use turbulence to adjust the phase (so it shifts x
            // in sin(x)) which makes the stripes undulate
            return color(1, 1, 1) * 0.5 * (1 + sin(scale * p.z() + 10 * noise.turbulence(p)));
        }

    public:
        perlin noise;
        double scale;  // noise frequency scale - higher scale corresponds to more frequent noise pattern
};

// texture that holds an image, using stb_image which reads in an image into a big array
// of unsigned char (each char is a RGB value in the range [0,255], and every 3 chars make
// up the data needed to specify a pixel: R, G, and B)
class image_texture : public texture {
    public: 
        const static int bytes_per_pixel = 3;  // since each pixel has R, G, and B data and each of those requires 1 byte

        image_texture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
        image_texture(const char* filename) {
            auto components_per_pixel = bytes_per_pixel;

            // load the image data into "data" array
            data = stbi_load(filename, &width, &height, &components_per_pixel, components_per_pixel);

            if (!data) {
                std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
                width = height = 0;
            }

            bytes_per_scanline = bytes_per_pixel * width;
        }

        ~image_texture() {
            delete data;
        }

        virtual color value(double u, double v, const vec3& p) const override {
            // if we have no texture data, then return solid cyan as a debugging aid
            if (data == nullptr) {
                return color(0, 1, 1);
            }

            // clamp input texture coordinates to [0,1] x [1,0]
            u = clamp(u, 0.0, 1.0);
            v = 1.0 - clamp(v, 0.0, 1.0);  // flip v to image coordinates

            // since u, v are values in the range [0,1], we need to scale them up to the
            // full width and height dimensions in order to map them to image coordinates in 
            // pixels - think of u, v as fractional positions along each dimension
            auto i = static_cast<int>(u * width);
            auto j = static_cast<int>(v * height);

            // clamp integer mapping, since actual coordinates should be less than 1.0
            if (i >= width) {
                i = width - 1;
            }
            if (j >= height) {
                j = height - 1;
            }

            const auto color_scale = 1.0 / 255.0;

            // use the mapped i, j image coordinates to retrieve the appropriate pixel data in the "data" array
            auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;

            // normalize the RGB values to be used in our data structure
            return color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
        }
    
    private:
        unsigned char* data;
        int width, height;
        int bytes_per_scanline;
};

#endif
