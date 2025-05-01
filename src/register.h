#include "healpix_map.h"

struct pixeldata {
    double theta;
    double phi;
    int value;
};

struct rgb_t {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

struct angle_t {
    double theta;
    double phi;
};

struct camera_t {
    angle_t fov;
    angle_t off;
};

struct range_t {
    long start;
    long end;
};

struct hpbound_t {
    int rstart;
    int rend;
    int data_size;
    range_t* pixels;
};

/*
 * Determine projection pointing vectors for pixels in a 2D image
 */
int register_pixels(double* img_data, int height, int width, int channels, 
    pointing* fov, pointing* off, pixeldata** vecs);

/*
 * Convert an angle (theta, phi) into a normalized vec3 on the unit sphere
 */
vec3 ang2vec(pointing* angle);

vec3 ang2vec(pointing angle);

/*
 * Determine which spherical pixels lie within the bounds of a projected 2D pixel
 */
int find_overlapping_pixels(pixeldata* angle, rangeset<int>* pixels, int width, 
    int height, int n_side, pointing* fov, pointing* off, T_Healpix_Base<int>* hp);

/*
 * Stacking kernel for map additions
 */
void __stack(Healpix_Map<double>* map, int pix, double value);

/*
 * Add one HEALPix map onto another
 */
int stack_hp(Healpix_Map<double>* map1, Healpix_Map<double>* map2);

/*
 * Add one HEALPix map onto another, inside of region \a region
 */
int stack_hp(Healpix_Map<double>* map1, Healpix_Map<double>* map2,
    std::vector<int>* region);

/*
 * Add a single element onto an existing HEALPix map
 */
int stack_hp(Healpix_Map<double>* map, int pix, double value);

/*
 * Add an image from the filesystem to a HEALPix map with 2D projection
 */
int add_image_to_map(Healpix_Map<double>* map, const char* file_loc, 
    pointing* fov, pointing* off, int nside);