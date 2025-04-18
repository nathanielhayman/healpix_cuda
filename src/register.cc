#include <healpix_base.h>
#include <pointing.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <cmath>
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

struct pixeldata {
    double theta;
    double phi;
    double value;
};

/*
 * Determine projection pointing vectors for pixels in a 2D image
 */
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
            pixeldata pd;
            pd.theta = (y - height/2) * fov->theta / height;
            pd.phi = (x - width/2) * fov->phi / width;
            pd.value = (double)img_data[y * width + x];

            *vecs[(y * width + x) * 3] = pd;
        }
    }

    return 0;
}

/*
 * Convert an angle (theta, phi) into a normalized vec3 on the unit sphere
 */
vec3 ang2vec(pointing* angle) {
    return vec3(
        cos(angle->theta) * cos(angle->phi),
        cos(angle->theta) * sin(angle->phi),
        sin(angle->theta)
    );
}

vec3 ang2vec(pointing angle) {
    return vec3(
        cos(angle.theta) * cos(angle.phi),
        cos(angle.theta) * sin(angle.phi),
        sin(angle.theta)
    );
}

/*
 * Determine which spherical pixels lie within the bounds of a projected 2D pixel
 */
int find_overlapping_pixels(pixeldata* angle, rangeset<int>* pixels, int width, 
    int height, int n_side, pointing* fov, pointing* off, T_Healpix_Base<int>* hp) 
{
    double pixel_size_y = fov->theta / height;
    double pixel_size_x = fov->phi / width;

    std::vector<pointing> vecs(4);

    double a = std::fmod((angle->theta - pixel_size_y/2), M_PI);

    vecs[0] = ang2vec(pointing(
        std::fmod((angle->theta - pixel_size_y/2), M_PI),
        std::fmod((angle->phi - pixel_size_x/2), 2*M_PI)
    ));

    vecs[1] = ang2vec(pointing(
        std::fmod((angle->theta + pixel_size_y/2), M_PI),
        std::fmod((angle->phi - pixel_size_x/2), 2*M_PI)
    ));

    vecs[2] = ang2vec(pointing(
        std::fmod((angle->theta + pixel_size_y/2), M_PI),
        std::fmod((angle->phi + pixel_size_x/2), 2*M_PI)
    ));

    vecs[3] = ang2vec(pointing(
        std::fmod((angle->theta - pixel_size_y/2), M_PI),
        std::fmod((angle->phi + pixel_size_x/2), 2*M_PI)
    ));

    *pixels = hp->query_polygon(vecs);
    
    return 0;
}

/*
 * Stacking kernel for map additions
 */
void __stack(Healpix_Map<double>* map, int pix, double value) {
    // simple averaged stacking
    (*map)[pix] = ((*map)[pix] + value)/2;
}

/*
 * Add one HEALPix map onto another
 */
int stack_hp(Healpix_Map<double>* map1, Healpix_Map<double>* map2) 
{
    for (int i = 0; i < map1->Npix(); ++i) {
        int loc = i;

        __stack(map1, loc, (*map2)[loc]);
    }

    return 0;
}

/*
 * Add one HEALPix map onto another, inside of region \a region
 */
int stack_hp(Healpix_Map<double>* map1, Healpix_Map<double>* map2,
    std::vector<int>* region) 
{
    for (int i = 0; i < region->size(); ++i) {
        int loc = (*region)[i];

        __stack(map1, loc, (*map2)[loc]);
    }

    return 0;
}

/*
 * Add a single element onto an existing HEALPix map
 */
int stack_hp(Healpix_Map<double>* map, int pix, double value) 
{
    __stack(map, pix, value);

    return 0;
}

/*
 * Add an image from the filesystem to a HEALPix map with 2D projection
 */
int add_image_to_map(Healpix_Map<double>* map, const char* file_loc, 
    pointing* fov, pointing* off, int nside)
{
    int width, height, channels;

    // vectors with pixel data from image 
    pixeldata* vecs = nullptr;

    // HEALPix Base engine (integer indexing)
    T_Healpix_Base<int> hp(nside, RING);

    double* img_data = (double*)stbi_load(
        file_loc, &width, &height, &channels, 1
    );

    int res = register_pixels(
        img_data, height, width, 
        channels, fov, off, &vecs
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
            nside, fov, off, &hp
        );

        std::vector<int> pd = pixels.data();

        // update spherical pixel with corresponding image value
        for (int i = 0; i < pixels.size(); ++i) {
            stack_hp(map, pd[i], vecs[i].value);
            // (*map)[pd[i]] = vecs[i].value;
        }

        // find the union of the two pixel arrays
        // *all_pixels = all_pixels->op_xor(pixels);
    }

    // free the allocated memory
    free(vecs);

    return 0;
}

int main() {
    int nside = 128;

    // camera FOV and rotational offset (theta, phi)
    pointing fov(M_PI/4, M_PI/4);
    pointing off(M_PI/2, 0);

    // HEALPix Map (double valuation)
    Healpix_Map<double> map = Healpix_Map<double>(nside, RING);

    // add the first image
    int res = add_image_to_map(&map, "image.jpg", &fov, &off, nside);

    if (res != 0) {
        printf("Error while adding image to map!\n");
        return res;
    }

    // rotate the next image by pi/4
    off.theta += M_PI/4;

    // add the second image
    res = add_image_to_map(&map, "image2.jpg", &fov, &off, nside);

    if (res != 0) {
        printf("Error while adding image to map!\n");
        return res;
    }

    // save the map to a file
    write_Healpix_map_to_fits(
        "output.fits", map, PLANCK_FLOAT64
    );

    return 0;
}