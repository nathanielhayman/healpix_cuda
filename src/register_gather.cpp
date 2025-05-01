#include <healpix_base.h>
#include <pointing.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <healpix_tables.h>
#include <chrono>
#include <filesystem>

#include "register.h"

#include <cmath>
#include <string>
#include <algorithm>

#include "stb_image.h"

#define ORDER 12
#define MAX_NSIDE 536870912

/*
 * Returns a modulus which "wraps around" for negative numbers (as in Python)
 */
template<typename T> T wrapping_mod(T v, T m) {
    return std::fmod(m + std::fmod(v, m), m);
}

/*
 * Determine which spherical pixels lie within the bounds of a projected 2D pixel
 */
int find_pix(pointing* angle, rangeset<int>* pixels, int width, 
    int height, pointing* fov, T_Healpix_Base<int>* hp) 
{
    double pixel_size_y = fov->theta;
    double pixel_size_x = fov->phi;

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

inline double sigmoid(double x) {
    return x;
    // double k = 0.05;
    // return 1 - 1/(1 + exp(-(abs(x)-0.5)/k));
}

inline int data_at(int y, int x, int channels, const unsigned char* imageData, int width) {
    if (channels == 3) {
        const rgb_t data = ((rgb_t*)imageData)[y * width + x];
        return (data.r + data.g + data.b) / 3;
    } else {
        return imageData[(y * width + x) * 3];
    }
}

inline double bilinear_interpolate(double x, double y, const unsigned char* imageData, int channels, int width, int height) {
    int x_c, x_f, y_c, y_f;
    int v1, v2, v3, v4;
    int q1, q2;
    
    x_c = std::min(width - 1, (int)ceil(x));
    x_f = std::max((int)floor(x), 0);
    y_c = std::min(width - 1, (int)ceil(y));
    y_f = std::max((int)floor(y), 0);

    v1 = data_at(y_f, x_f, channels, imageData, width);
    v2 = data_at(y_f, x_c, channels, imageData, width);
    v3 = data_at(y_c, x_f, channels, imageData, width);
    v4 = data_at(y_c, x_c, channels, imageData, width);

    if (x_c == x_f && y_c == y_f)
        return v1;
    if (x_c == x_f)
        return v1 * sigmoid(y_c - y) + v3 * sigmoid(y - y_f);
    if (y_c == y_f)
        return v1 * sigmoid(x_c - x) + v2 * sigmoid(x - x_f);

    q1 = v1 * sigmoid(x_c - x) + v2 * sigmoid(x - x_f);
    q2 = v3 * sigmoid(x_c - x) + v4 * sigmoid(x - x_f);
    return q1 * sigmoid(y_c - y) + q2 * sigmoid(y - y_f);
}

void populate_spherical(int pix, pointing* fov, pointing* off, const unsigned char* imageData, int channels, int width, int height, Healpix_Map<int>* map, T_Healpix_Base<int>* hp) {
    double x, y;
    int value;
    pointing angle;

    angle = hp->pix2ang(pix);

    x = wrapping_mod(((angle.phi - off->phi) * (double)width) / fov->phi + width/2, (double)(width-1));
    y = wrapping_mod(((angle.theta - off->theta) * (double)height) / fov->theta + height/2, (double)(height-1));

    if (nearbyintf(x) < 0)
        throw -1;
    
    if (nearbyintf(y) < 0)
        throw -1;
    
    double val = bilinear_interpolate(x, y, imageData, channels, width, height);

    (*map)[pix] = (int)val; // imageData[((int)nearbyintf(y) * width + (int)nearbyintf(x)) * 3];
}

/*
 * Add an image from the filesystem to a HEALPix map with 2D projection
 */
int add_image_to_map(Healpix_Map<int>* map, const char* file_loc, 
    pointing* fov, pointing* off)
{
    int width, height, channels;

    // HEALPix Base engine (integer indexing)
    T_Healpix_Base<int> hp(ORDER, RING);

    unsigned char* img_data = stbi_load(
        file_loc, &width, &height, &channels, 3
    );

    printf("%s, %d, %d, %d, %d\n", file_loc, width, height, channels, 3);

    // vectors with pixel data from image 
    pixeldata* vecs = new pixeldata[width * height];

    printf("Number of channels: %d\n", channels);

    auto t1 = std::chrono::high_resolution_clock::now();
    rangeset<int>* pixels = new rangeset<int>();
    find_pix(off, pixels, width, height, fov, &hp);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    printf("find_pix() took %lldms\n", ms_int.count());

    // rangeset<int>* all_pixels;

    // find the overlapping spherical pixels for each pixel in the image
    t1 = std::chrono::high_resolution_clock::now();
    // update spherical pixel with corresponding image value
    for (int j = 0; j < pixels->size(); ++j) {
        for (int k = pixels->ivbegin(j); k < pixels->ivend(j); ++k) {
            // printf("processing pixel #%d with value %f\n", k, vecs[i].value);
            populate_spherical(k, fov, off, img_data, channels, width, height, map, &hp);
        }
        // (*map)[pd[i]] = vecs[i].value;
    }
    t2 = std::chrono::high_resolution_clock::now();

    ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    printf("populate_spherical() took %lldms\n", ms_int.count());

    // free the allocated memory
    delete[] vecs;

    return 0;
}

int main(int argc, char** argv) {
    std::string file_loc = "image.jpg";

    if (argc > 1) {
        file_loc = argv[1];
        // nside = std::strtol(argv[1], p, 10);
    }

    // camera FOV and rotational offset (theta, phi)
    pointing fov(M_PI/4, M_PI/4);
    pointing off(M_PI/2, 0);

    // HEALPix Map (double valuation)
    Healpix_Map<int>* map = new Healpix_Map<int>(ORDER, RING);

    printf("Adding the first image...\n");

    // add the first image
    auto t1 = std::chrono::high_resolution_clock::now();
    int res = add_image_to_map(map, file_loc.c_str(), &fov, &off);
    auto t2 = std::chrono::high_resolution_clock::now();

    if (res != 0) {
        printf("Error while adding image to map: %d!\n", res);
        return res;
    }

    std::chrono::milliseconds ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    printf("add_image_to_map() took %lldms\n", ms_int.count());

    try {
        if (std::__fs::filesystem::remove("output.fits")) {
            std::cout << "Outfile deleted successfully." << std::endl;
        } else {
            std::cout << "Outfile not found." << std::endl;
        }
    } catch (const std::__fs::filesystem::filesystem_error& err) {
        std::cerr << "Filesystem error: " << err.what() << std::endl;
    }

    write_Healpix_map_to_fits(
        "output.fits", *map, PLANCK_FLOAT64
    );

    delete map;

    return 0;
}