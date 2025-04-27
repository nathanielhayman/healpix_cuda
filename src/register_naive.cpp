#include <healpix_base.h>
#include <pointing.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <healpix_tables.h>

#include "register.h"

#include <cmath>
#include <string>

#include "stb_image.h"

#define NSIDE 1024
#define MAX_NSIDE 536870912

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
 * Returns a modulus which "wraps around" for negative numbers (as in Python)
 */
template<typename T> T wrapping_mod(T v, T m) {
    return std::fmod(m + std::fmod(v, m), m);
}

/*
 * Determine projection pointing vectors for pixels in a 2D image
 */
int register_pixels(unsigned char* img_data, int height, int width, int channels, 
    pointing* fov, pointing* off, pixeldata* vecs) 
{
    bool gray = false;

    if (channels != 3) {
        if (channels != 1)
            return -1; // only support color or grayscale images
        
        gray = true;
    }

    rgb_t cur_pix;

    pixeldata pd;

    // add the pointing vectors to the `vecs` array
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // printf("image data: (%d, %d) -> %d\n", y, x, img_data[(y * width + x) * 3]);
            
            pd.theta = (y - height/2) * fov->theta / height + off->theta; 
            pd.phi = (x - width/2) * fov->phi / width + off->phi;

            // if image is 3-channel, convert value to single channel via averaging
            if (!gray) {
                cur_pix = ((rgb_t*)img_data)[y * width + x];

                pd.value = (cur_pix.r + cur_pix.g + cur_pix.b)/(3);
            } else {
                pd.value = img_data[(y * width + x) * 3];
            }

            vecs[y * width + x] = pd;
        }
    }

    return 0;
}

/*
 * Determine which spherical pixels lie within the bounds of a projected 2D pixel
 */
int find_overlapping_pixels(pixeldata* angle, rangeset<int>* pixels, int width, 
    int height, int n_side, pointing* fov, T_Healpix_Base<int>* hp) 
{
    double pixel_size_y = fov->theta / height;
    double pixel_size_x = fov->phi / width;

    std::vector<pointing> vecs(4);

    // double a = std::fmod((angle->theta - pixel_size_y/2), M_PI);

    vecs[0] = pointing(
        wrapping_mod((angle->theta - pixel_size_y/2), M_PI),
        wrapping_mod((angle->phi - pixel_size_x/2), 2*M_PI)
    );

    vecs[1] = pointing(
        wrapping_mod((angle->theta + pixel_size_y/2), M_PI),
        wrapping_mod((angle->phi - pixel_size_x/2), 2*M_PI)
    );

    vecs[2] = pointing(
        wrapping_mod((angle->theta + pixel_size_y/2), M_PI),
        wrapping_mod((angle->phi + pixel_size_x/2), 2*M_PI)
    );

    vecs[3] = pointing(
        wrapping_mod((angle->theta - pixel_size_y/2), M_PI),
        wrapping_mod((angle->phi + pixel_size_x/2), 2*M_PI)
    );

    // for (int i = 0; i < vecs.size(); i++) {
    //     printf("index: %d, value: (%f, %f)\n", i, vecs[i].theta, vecs[i].phi);
    // }

    hp->query_polygon(vecs, *pixels);
    
    return 0;
}

/*
 * Stacking kernel for map additions
 */
void __stack(Healpix_Map<int>* map, int pix, int value) {
    // simple averaged stacking, add full value if it is empty
    if ((*map)[pix] == 0)
        (*map)[pix] = value;
    else
        (*map)[pix] = ((*map)[pix] + value)/2;
}

/*
 * Add one HEALPix map onto another
 */
int stack_hp(Healpix_Map<int>* map1, Healpix_Map<int>* map2) 
{
    for (int i = 0; i < map1->Npix(); ++i) {
        int loc = i;

        __stack(map1, loc, (*map2)[loc]);
    }

    return 0;
}

/*
 * Add one HEALPix map onto another, iorder of region \a region
 */
int stack_hp(Healpix_Map<int>* map1, Healpix_Map<int>* map2,
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
int stack_hp(Healpix_Map<int>* map, int pix, int value) 
{
    __stack(map, pix, value);

    return 0;
}

/*
 * Add an image from the filesystem to a HEALPix map with 2D projection
 */
int add_image_to_map(Healpix_Map<int>* map, const char* file_loc, 
    pointing* fov, pointing* off, int order)
{
    int width, height, channels;

    // HEALPix Base engine (integer indexing)
    T_Healpix_Base<int> hp(order, Healpix_Ordering_Scheme::RING);

    unsigned char* img_data = stbi_load(
        file_loc, &width, &height, &channels, 3
    );

    // vectors with pixel data from image 
    pixeldata* vecs = new pixeldata[width * height];

    printf("Number of channels: %d\n", channels);

    int res = register_pixels(
        img_data, height, width, 
        channels, fov, off, vecs
    );

    if (res != 0) {
        printf("Error registering pixels\n");
        return res;
    }

    // rangeset<int>* all_pixels;

    // find the overlapping spherical pixels for each pixel in the image
    for (int i = 0; i < height * width; ++i) {
        printf("pixel attitude at %d: (%f, %f) -> %d\n", i, vecs[i].theta, vecs[i].phi, vecs[i].value);
        rangeset<int>* pixels = new rangeset<int>();
        find_overlapping_pixels(
            &(vecs[i]), pixels, width, height, 
            order, fov, &hp
        );

        std::vector<int> pd = pixels->data();

        // update spherical pixel with corresponding image value
        for (int j = 0; j < pixels->size(); ++j) {
            for (int k = pixels->ivbegin(j); k < pixels->ivend(j); ++k) {
                // printf("processing pixel #%d with value %f\n", k, vecs[i].value);
                stack_hp(map, k, vecs[i].value);
            }
            // (*map)[pd[i]] = vecs[i].value;
        }

        free(pixels);

        // find the union of the two pixel arrays
        // *all_pixels = all_pixels->op_xor(pixels);
    }

    // free the allocated memory
    delete[] vecs;

    return 0;
}

int main(int argc, char** argv) {
    int nside = NSIDE;

    char** p;

    if (argc > 1) {
        int nside = std::strtol(argv[1], p, 10);
    }

    nside = nside < MAX_NSIDE ? nside : NSIDE;

    int order = T_Healpix_Base<int>::nside2order(nside);

    // camera FOV and rotational offset (theta, phi)
    pointing fov(M_PI/4, M_PI/4);
    pointing off(M_PI/2, 0);

    // HEALPix Map (double valuation)
    Healpix_Map<int>* map = new Healpix_Map<int>(order, RING);

    printf("Adding the first image...\n");

    // add the first image
    int res = add_image_to_map(map, (const char*)"image.jpg", &fov, &off, order);

    if (res != 0) {
        printf("Error while adding image to map: %d!\n", res);
        return res;
    }

    // rotate the next image by pi/4
    off.theta += M_PI/4;

    // printf("Adding the second image...\n");

    // add the second image
    // res = add_image_to_map(map, (const char*)"image2.jpg", &fov, &off, order);

    // if (res != 0) {
    //     printf("Error while adding image to map!\n");
    //     return res;
    // }

    // printf("Writing map to file...\n");

    // save the map to a file
    write_Healpix_map_to_fits(
        "output.fits", *map, PLANCK_FLOAT64
    );

    delete map;

    return 0;
}
