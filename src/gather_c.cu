#include "register_alternate.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "stb_image.h"

#define BLOCK_DIM 16
#define ORDER 12
#define MAX_NSIDE 536870912
#define MAX_ORDER

static const int order_ = ORDER;
static const long nside_ = 1 << order_;
static const long npface_ = nside_ << order_;
static const long ncap_ = (npface_ - nside_) << 1;
static const long npix_ = 12 * npface_;
static const double fact2_ = 4. / npix_;
static const double fact1_ = (nside_ << 1) * fact2_;
static const long nring_ = 4 * nside_ - 1;

static const double twothird = 2.0 / 3.0;
static const double pi = 3.141592653589793238462643383279502884197;
static const double twopi = 6.283185307179586476925286766559005768394;
static const double halfpi = 1.570796326794896619231321691639751442099;
static const double inv_halfpi = 0.6366197723675813430755350534900574;
static const double inv_twopi = 1.0 / twopi;
static const double fourpi = 12.56637061435917295385057353311801153679;
static const double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;
static const double ln2 = 0.6931471805599453094172321214581766;
static const double inv_ln2 = 1.4426950408889634073599246810018921;
static const double ln10 = 2.3025850929940456840179914546843642;
static const double onethird = 1.0 / 3.0;
static const double fourthird = 4.0 / 3.0;
static const double degr2rad = pi / 180.0;
static const double arcmin2rad = degr2rad / 60;
static const double rad2degr = 180.0 / pi;

//! Ratio between FWHM and sigma of a Gauss curve (\f$\sqrt{8\ln2}\f$).
static const double sigma2fwhm = 2.3548200450309493; // sqrt(8*log(2.))
static const double fwhm2sigma = 1 / sigma2fwhm;

using Eigen::Vector2i;
using Eigen::Vector3d;

template <typename I>
using stdvec = std::vector<I>;

/*! Writes diagnostic output and exits with an error status. */
#define planck_fail(msg)                     \
    do                                       \
    {                                        \
        printf("Planck failure: %s\n", msg); \
        throw -1;                            \
    } while (0)

/*! Writes diagnostic output and exits with an error status if \a testval
    is \a false. */
#define planck_assert(testval, msg) \
    do                              \
    {                               \
        if (testval)                \
            ;                       \
        else                        \
            planck_fail(msg);       \
    } while (0)

/*
 * Returns a modulus which "wraps around" for negative numbers (as in Python)
 */
template <typename T>
T wrapping_mod(T v, T m)
{
    return std::fmod(m + std::fmod(v, m), m);
}

/*
 * Convert an angle (theta, phi) into a normalized vec3 on the unit sphere
 */
__forceinline__ inline void ang2vec(Vector3d *vec, const angle_t *angle)
{
    double st = sin(angle->theta);
    (*vec)[0] = st * cos(angle->phi);
    (*vec)[1] = st * sin(angle->phi);
    (*vec)[2] = cos(angle->theta);
}

__forceinline__ inline void vec2ang(angle_t *angle, const Vector3d *vec)
{
    angle->theta = acos(vec->z());
    angle->phi = ((0 < vec->y()) - (vec->y() < 0)) * acos(vec->x() / (sqrt(pow(vec->x(), 2) + pow(vec->y(), 2))));
}

// get the ring above the specified
__forceinline__ inline int ring_above(double z)
{
    double az = abs(z);
    if (az <= twothird) // equatorial region
        return nside_ * (2 - 1.5 * z);
    int iring = nside_ * sqrt(3 * (1 - az));
    return (z > 0) ? iring : 4 * nside_ - iring - 1;
}

__forceinline__ inline double ring2z(long ring)
{
    if (ring < nside_)
        return 1 - ring * ring * fact2_;
    if (ring <= 3 * nside_)
        return (2 * nside_ - ring) * fact1_;
    ring = 4 * nside_ - ring;
    return ring * ring * fact2_ - 1;
}

inline void get_ring_info_small(long ring, long &startpix,
                                long &ringpix, bool &shifted)
{
    if (ring < nside_)
    {
        shifted = true;
        ringpix = 4 * ring;
        startpix = 2 * ring * (ring - 1);
    }
    else if (ring < 3 * nside_)
    {
        shifted = ((ring - nside_) & 1) == 0;
        ringpix = 4 * nside_;
        startpix = ncap_ + (ring - nside_) * ringpix;
    }
    else
    {
        shifted = true;
        int nr = 4 * nside_ - ring;
        ringpix = 4 * nr;
        startpix = npix_ - 2 * nr * (nr + 1);
    }
}

__forceinline__ inline bool is_in(int x, int l, int h, int nr)
{
    if (l <= h)
        return l <= x && x < h;
    return x >= l || x < h;
}

/*!
 *  Identifies the spherical pixels which lie within the bounds of a multi-disc query
 *  This query is an intersection, not a union
 *
 *  \a norm : array of disc centers
 *  \a rad : array of disc radii
 *  \a pixset : output pixel ranges
 *
 *  We have explicitly removed all code which relates to inclusive idenitification
 *  since it is not particularly relevant to our use case
 */
__device__ void query_multidisc(const Vector3d *norm, int norm_l,
                     const double *rad, int rad_l, hpbound_t &pixset)
{
    int nv = norm_l; // number of vertices

    planck_assert(nv == rad_l, "inconsistent input arrays");

    // artifacts of inclusive scan code
    int fct = 1;
    double rpsmall, rpbig;
    rpsmall = rpbig = 0;

    int irmin = 1, irmax = nring_;

    // z0 contains the relative heights of each disc (cos of the colatitude)
    double z0[nv], xa[nv], cosrsmall[nv], cosrbig[nv];

    angle_t ptg[nv];

    int cpix[nv];

    /*
     *  Iterates over the normalized vertices of the identified disc (in our case,
     *  derived from the polygon edges of the registered pixels). This should not
     *  be a source of control divergence due to its uniformity and independence
     *  across the dataset
     *
     *  Identifies the rings
     */
    using namespace std;
    for (long i = 0; i < nv; ++i)
    {
        double rsmall = rad[i] + rpsmall;
        if (rsmall < pi)
        {
            double rbig = min(pi, rad[i] + rpbig);

            angle_t pnt;
            vec2ang(&pnt, &norm[i]);

            cosrsmall[i] = (cos(rsmall));

            // cos of the radius
            cosrbig[i] = (cos(rbig));

            // find z-coordinate of center
            double cth = cos(pnt.theta);
            z0[i] = (cth);

            // 1 / sin(theta)
            xa[i] = (1. / sqrt((1 - cth) * (1 + cth)));
            ptg[i] = (pnt);

            // find the min row corresponding to the disc
            double rlat1 = pnt.theta - rsmall;
            double zmax = cos(rlat1);
            int irmin_t = (rlat1 <= 0) ? 1 : ring_above(zmax) + 1;

            // find the max row corresponding to the disc
            double rlat2 = pnt.theta + rsmall;
            double zmin = cos(rlat2);
            long irmax_t = (rlat2 >= pi) ? nring_ : ring_above(zmin);

            if (irmax_t < irmax)
                irmax = irmax_t;
            if (irmin_t > irmin)
                irmin = irmin_t;
        }
    }

    pixset.rstart = irmin;
    pixset.rend = irmax;

    int dr = 0;      // degenerate rows which do not contain any intersected pixels
    bool ff = false; // was the first non-degenerate row found

    /*
     *  Iterates over the identified rings to extract individual pixels. Also
     *  not a source of control divergence because each registered pixel boundary
     *  should encompass a uniform quantity of rings (with error of ~1)
     */
    for (long iz = irmin; iz <= irmax; ++iz)
    {
        double z = ring2z(iz);
        long ipix1, nr;
        long l, h;
        bool shifted;

        // ipix1 is the starting pixel, nr is the number of pixels in the row
        get_ring_info_small(iz, ipix1, nr, shifted);
        double shift = shifted ? 0.5 : 0.;

        hprange_t tr;

        // add range of pixels from the start pixel to the last (+nr)
        tr.start = -1;
        tr.end = -1;
        tr.ring = iz;

        bool nysq = false;

        // for each disc in the query, check which pixels in the current row actually
        // intersect a disc
        for (long j = 0; j < nv; ++j)
        {
            // solve for the azimuthal delta (dphi) of the disc on the current ring,
            // and convert that to a HEALPix coordinate range
            double x = (cosrbig[j] - z * z0[j]) * xa[j];
            double ysq = 1. - z * z - x * x;

            if (ysq > 0) // avoid divide by zero
            {
                nysq = true;

                double dphi = atan2(sqrt(ysq), x);

                long ip_lo = floor(nr * inv_twopi * (ptg[j].phi - dphi) - shift) + 1;
                long ip_hi = floor(nr * inv_twopi * (ptg[j].phi + dphi) - shift);

                if (ip_hi >= nr)
                {
                    ip_lo -= nr;
                    ip_hi -= nr;
                }

                l = ip_lo < 0 ? ipix1 + ip_lo + nr : ipix1 + ip_lo;
                h = ipix1 + ip_hi + 1;

                // ip_lo = ((ip_lo % nr) + nr) % nr;
                // ip_hi = ((ip_hi % nr) + nr) % nr;

                // l = ((ipix1 + ip_lo) % nr + nr) % nr;
                // h = ((ipix1 + ip_hi + 1) % nr + nr) % nr;

                if (tr.start == -1 && tr.end == -1)
                {
                    tr.start = l;
                    tr.end = h;
                    continue;
                }

                bool w1 = l >= h;
                bool w2 = tr.start >= tr.end;

                if (!w1 && !w2)
                {
                    tr.start = max(tr.start, l);
                    tr.end = min(tr.end, h);
                }
                else if (w1 && w2)
                {
                    if (!is_in(tr.start, l, h, nr))
                        tr.start = l;

                    if (!is_in(tr.end - 1, l, h, nr))
                        tr.end = h;
                }
                else
                {
                    int l2 = tr.start;
                    int h2 = tr.end;

                    if (!w1 && w2)
                    {
                        int tmp;
                        tmp = l;
                        l = l2;
                        l2 = tmp;
                        tmp = h;
                        h = h2;
                        h2 = tmp;
                        tmp = w1;
                        w1 = w2;
                        w2 = tmp;
                    }

                    if (is_in(l2, l, h, nr))
                    {
                        tr.start = l2;
                    }
                    else if (is_in(l, l2, h2, nr))
                    {
                        tr.start = l;
                    }
                    else
                    {
                        // planck_assert(false, "Invalid intersection of edges!\n");
                        nysq = false;
                        break;
                    }

                    int end = (h < h2) ? h : h2;
                    if (!is_in(tr.start, l2, h2, nr) || !is_in(tr.start, l, h, nr))
                    {
                        // planck_assert(false, "Invalid intersection of edges!\n");
                        nysq = false;
                        break;
                    }

                    tr.end = end;
                }
            }
        }

        if (tr.start == tr.end || !nysq)
            ++dr;
        else
        {
            if (!ff)
            {
                pixset.rstart = iz;
                ff = true;
            }

            pixset.pixels[iz - dr - irmin] = tr;

            pixset.rend = iz;
        }
    }

    pixset.data_size = irmax - irmin - dr;
}

__device__ void get_circle(const Vector3d *point, long q1, long q2, Vector3d &center,
                double &cosrad)
{
    center = (point[q1] + point[q2]).normalized();
    cosrad = point[q1].dot(center);
    for (long i = 0; i < q1; ++i)
        if (point[i].dot(center) < cosrad) // point outside the current circle
        {
            center = (point[q1] - point[i]).cross(point[q2] - point[i]).normalized();
            cosrad = point[i].dot(center);
            if (cosrad < 0)
            {
                center *= -1;
                cosrad = -cosrad;
            }
        }
}

__device__ void get_circle(const Vector3d *point, long q, Vector3d &center,
                double &cosrad)
{
    center = (point[0] + point[q]).normalized();
    cosrad = point[0].dot(center);
    for (long i = 1; i < q; ++i)
        if (point[i].dot(center) < cosrad) // point outside the current circle
            get_circle(point, i, q, center, cosrad);
}

/*
 *  Assembles the discs whose interesection creates the spherical-space analog
 *  to a four-pixel "square" and performs a query_multidisc to find the corresponding
 *  HEALPix pixels
 */
__device__ void query_square(const angle_t *vertex, hpbound_t &pixset)
{
    // convert all pointing vectors to vec3
    Vector3d vv[4];
    for (long i = 0; i < 4; ++i)
        ang2vec(&vv[i], &vertex[i]);
    Vector3d normal[4];

    int flip = 0;

    for (long i = 0; i < 4; ++i)
    {
        normal[i] = vv[i].cross(vv[(i + 1) % 4]).normalized();

        double hnd = normal[i].dot(vv[(i + 2) % 4]);

        planck_assert(abs(hnd) > 1e-10, "degenerate corner");
        if (i == 0)
            flip = (hnd < 0.) ? -1 : 1;
        else
            planck_assert(flip * hnd > 0, "polygon is not convex");
        normal[i] *= flip;
    }

    double rad[4];

    for (int i = 0; i < 4; i++)
        rad[i] = halfpi; // every disc has radius pi/2

    query_multidisc(normal, 4, rad, 4, pixset);
}

/*
 * Stacking kernel for map additions
 */
__forceinline__ inline void __stack(int *map, int pix, int value)
{
    // simple averaged stacking, add full value if it is empty
    if (map[pix] == 0)
        map[pix] = value;
    else
        map[pix] = (map[pix] + value) / 2;
}

/*
 * Add a single element onto an existing HEALPix map
 */
__forceinline__ inline int stack_hp(int *map, int pix, int value)
{
    __stack(map, pix, value);

    return 0;
}

__global__ void register_pixel(const unsigned char *imageData, Vector2i *imageSize,
                    camera_t *deviceCamProps, int *deviceMap)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    pixeldata pd;

    // Register pixel center to polar coordinates
    pd.theta = (y - imageSize->y() / 2) * deviceCamProps->fov.theta / imageSize->y() + deviceCamProps->off.theta;
    pd.phi = (x - imageSize->x() / 2) * deviceCamProps->fov.phi / imageSize->x() + deviceCamProps->off.phi;

    if (deviceCamProps->n_channels == 3)
    {
        rgb_t cur_pix = ((rgb_t *)imageData)[y * imageSize->x() + x];
        pd.value = (cur_pix.r + cur_pix.g + cur_pix.b) / (3);
    }
    else
    {
        pd.value = imageData[(y * imageSize->x() + x) * 3];
    }

    double pixel_size_y = deviceCamProps->fov.theta / imageSize->y();
    double pixel_size_x = deviceCamProps->fov.phi / imageSize->x();

    angle_t vecs[4];

    // find the boundaries of the pixel on the sphere
    vecs[0].theta = wrapping_mod((pd.theta - pixel_size_y / 2), pi);
    vecs[0].phi = wrapping_mod((pd.phi - pixel_size_x / 2), 2 * pi);

    vecs[1].theta = wrapping_mod((pd.theta + pixel_size_y / 2), pi);
    vecs[1].phi = wrapping_mod((pd.phi - pixel_size_x / 2), 2 * pi);

    vecs[2].theta = wrapping_mod((pd.theta + pixel_size_y / 2), pi);
    vecs[2].phi = wrapping_mod((pd.phi + pixel_size_x / 2), 2 * pi);

    vecs[3].theta = wrapping_mod((pd.theta - pixel_size_y / 2), pi);
    vecs[3].phi = wrapping_mod((pd.phi + pixel_size_x / 2), 2 * pi);

    hpbound_t pixels;

    pixels.pixels = (hprange_t *)malloc(sizeof(hprange_t) * nring_);

    // find corresponding HEALPix pixel indices
    query_square(vecs, pixels);

    hprange_t range;

    long minpix;
    long numpix;
    bool shifted;

    // update spherical pixels with corresponding image value
    for (int row = 0; row < pixels.data_size; row++)
    { // row
        range = pixels.pixels[row];
        if (range.start > range.end)
        {
            get_ring_info_small(range.ring, minpix, numpix, shifted);

            for (int i = range.start; i < minpix + numpix; i++)
                stack_hp(deviceMap, i, pd.value);

            for (int i = minpix; i < range.end; i++)
                stack_hp(deviceMap, i, pd.value);
        }
        else
        {
            for (int i = range.start; i < range.end; i++) // pixel
                stack_hp(deviceMap, i, pd.value);
        }
    }
}

int main(int argc, char *argv[])
{
    int nside, mapSize;
    int width, height, channels;
    int imageSize;

    camera_t *camProps;

    std::string file_loc = "image.jpg";

    if (argc > 1)
    {
        file_loc = argv[1];
    }

    mapSize = 12 * pow((long)nside, 2);

    camProps = (camera_t *)calloc(sizeof(camera_t), 1);

    // camera FOV and rotational offset (theta, phi)
    camProps->fov.theta = halfpi / 2;
    camProps->fov.phi = halfpi / 2;
    camProps->off.theta = halfpi;

    unsigned char *img_data = stbi_load(
        file_loc, &width, &height, &channels, 3);

    // set size of image data for transfer
    imageSize = width * height * sizeof(char);

    cudaMalloc((void **)&deviceInput, imageSize);
    cudaMalloc((void **)&deviceCamProps, sizeof(camera_t));
    cudaMalloc((void **)&deviceMap, mapSize);

    cudaDeviceSynchronize();

    // copy image to device
    cudaMemcpy(deviceInput, imageData, imageSize, cudaMemcpyHostToDevice);
    cudaMemcpy(deviceCamProps, camProps, sizeof(camera_t), cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();

    dim3 dimBlock(BLOCK_DIM, BLOCK_DIM, 1);
    dim3 dimGrid((inputLength + dimBlock.x - 1) / dimBlock.x, 1, 1);

    register_pixel<<<dimGrid, dimBlock>>>(
        deviceInput,
        imageSize,
        deviceCamProps,
        deviceMap);

    cudaMemcpy(hostMap, deviceMap, mapSize, cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();

    cudaFree(deviceInput);
    cudaFree(deviceCamProps);
    cudaFree(deviceMap);

    FILE *fptr;

    fptr = fopen("map_out", "w");

    if (fptr == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    fwrite(&order_, sizeof(int), 1, fptr);

    fwrite(&npix_, sizeof(long), 1, fptr);

    // write data to file
    fwrite(mapData, sizeof(int), npix_, fptr);

    fclose(fptr);

    free(camProps);
    free(imageData);

    return 0;
}