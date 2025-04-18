#include "healpix_base.h"
#include "healpix_map.h"
#include "pointing.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"

#include <cmath>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

struct pixeldata {
    double theta;
    double phi;
    double value;
};

int register_pixels(double* img_data, int height, int width, int channels, 
    pointing* fov, pointing* off, pixeldata** vecs) 
{
    if (channels != 1) {
        return -1; // Only support Grayscale images
    }

    // empty array to put pointing vectors into
    *vecs = (pixeldata*)malloc(height * width * sizeof(pixeldata));

    // add the pointing vectors to the `vecs` array
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            *vecs[(y * width + x) * 3] = pixeldata{
                (y - height/2) * fov->theta / height,
                (x - width/2) * fov->phi / width,
                (double)img_data[y * width + x]
            };
        }
    }

    return 0;
}

vec3 ang2vec(pointing* angle) {
    return vec3(
        cos(angle->theta) * cos(angle->phi),
        cos(angle->theta) * sin(angle->phi),
        sin(angle->theta)
    );
}

int find_overlapping_pixels(pixeldata* angle, rangeset<int>* pixels, int width, 
    int height, int n_side, pointing* fov, pointing* off, T_Healpix_Base<int>* hp) 
{
    double pixel_size_y = fov->theta / height;
    double pixel_size_x = fov->phi / width;

    std::vector<pointing> vecs(4);

    double a = std::fmod((angle->theta - pixel_size_y/2), M_PI);

    vecs[0] = ang2vec(
        &pointing(
            std::fmod((angle->theta - pixel_size_y/2), M_PI),
            std::fmod((angle->phi - pixel_size_x/2), 2*M_PI)
        )
    );

    vecs[1] = ang2vec(
        &pointing(
            std::fmod((angle->theta + pixel_size_y/2), M_PI),
            std::fmod((angle->phi - pixel_size_x/2), 2*M_PI)
        )
    );

    vecs[2] = ang2vec(
        &pointing(
            std::fmod((angle->theta + pixel_size_y/2), M_PI),
            std::fmod((angle->phi + pixel_size_x/2), 2*M_PI)
        )
    );

    vecs[3] = ang2vec(
        &pointing(
            std::fmod((angle->theta - pixel_size_y/2), M_PI),
            std::fmod((angle->phi + pixel_size_x/2), 2*M_PI)
        )
    );

    *pixels = hp->query_polygon(vecs);
    
    return 0;
}

int main() {
    int nside = 128;

    // HEALPix Base engine (integer indexing)
    T_Healpix_Base<int> hp(nside, RING);

    // camera FOV and rotational offset (theta, phi)
    pointing fov(M_PI/4, M_PI/4);
    pointing off(M_PI/2, 0);

    int width, height, channels;

    // vectors with pixel data from image 
    pixeldata* vecs = nullptr;

    // HEALPix Map (double valuation)
    Healpix_Map<double> map = Healpix_Map<double>(nside, RING);

    double* img_data = (double*)stbi_load(
        "image.png", &width, &height, &channels, 1
    );

    int res = register_pixels(
        img_data, height, width, 
        channels, &pointing(0.1, 0.2), &pointing(0.3, 0.4), &vecs
    );

    if (res != 0) {
        printf("Error registering pixels\n");
        return res;
    }

    rangeset<double>* all_pixels;

    // find the overlapping spherical pixels for each pixel in the image
    for (int i = 0; i < height * width; ++i) {
        rangeset<int> pixels;
        find_overlapping_pixels(
            &vecs[i], &pixels, width, height, 
            nside, &fov, &off, &hp
        );

        std::vector<int> pd = pixels.data();

        for (int i = 0; i < pixels.size(); ++i) {
            map[pd[i]] = vecs[i].value;
        }

        // find the union of the two pixel arrays
        // *all_pixels = all_pixels->op_xor(pixels);
    }

    // save the map to a file
    write_Healpix_map_to_fits(
        "output.fits", map, PLANCK_FLOAT64
    );


    // free the allocated memory
    free(vecs);

    return 0;
}