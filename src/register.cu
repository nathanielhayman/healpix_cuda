#include "register.h"

#include <vector>
#include <cmath>
#include <algorithm>

#include <Eigen/Core>

#include "stb_image.h"

#define BLOCK_DIM 16
#define ORDER 4
#define MAX_NSIDE 536870912
#define MAX_ORDER

static const int order_ = ORDER;
static const long nside_  = 1<<order_;
static const long npface_ = nside_<<order_;
static const long ncap_   = (npface_-nside_)<<1;
static const long npix_   = 12*npface_;
static const long fact2_  = 4./npix_;
static const long fact1_  = (nside_<<1)*fact2_;

static const double twothird=2.0/3.0;
static const double pi=3.141592653589793238462643383279502884197;
static const double twopi=6.283185307179586476925286766559005768394;
static const double halfpi=1.570796326794896619231321691639751442099;
static const double inv_halfpi=0.6366197723675813430755350534900574;

using Eigen::Vector3d;
using Eigen::Vector2i;

template<typename I> using stdvec = std::vector<I>;

/*! Writes diagnostic output and exits with an error status. */
#define planck_fail(msg)                                       s    \
  do                                                               \
  {                                                                \
    planck_failure__(__FILE__, __LINE__, PLANCK_FUNC_NAME__, msg); \
    throw PlanckError(msg);                                        \
  } while (0)

/*! Throws a PlanckError without diagnostic message. */
#define planck_fail_quietly(msg) \
  do                             \
  {                              \
    throw PlanckError(msg);      \
  } while (0)

/*! Writes diagnostic output and exits with an error status if \a testval
    is \a false. */
#define planck_assert(testval, msg) \
  do                                \
  {                                 \
    if (testval)                    \
      ;                             \
    else                            \
      planck_fail(msg);             \
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
void ang2vec(Vector3d* vec, const angle_t* angle) {
    (*vec)[0] = cos(angle->theta) * cos(angle->phi);
    (*vec)[1] = cos(angle->theta) * sin(angle->phi);
    (*vec)[2] = sin(angle->theta);
}


void vec2ang(angle_t *angle, const Vector3d *vec) {
  angle->theta = acos(vec->z());
  angle->phi = ((0 < vec->y()) - (vec->y() < 0)) * acos(vec->x()/(sqrt(pow(vec->x(), 2) + pow(vec->y(), 2))));
}

long loc2pix (double z, double phi,
    double sth, bool have_sth)
{
    double za = abs(z);
    double tt = fmodulo(phi*inv_halfpi,4.0); // in [0,4)

    if (za<=twothird) // Equatorial region
    {
        long nl4 = 4*nside_;
        double temp1 = nside_*(0.5+tt);
        double temp2 = nside_*z*0.75;
        long jp = long(temp1-temp2); // index of  ascending edge line
        long jm = long(temp1+temp2); // index of descending edge line

        // ring number counted from z=2/3
        long ir = nside_ + 1 + jp - jm; // in {1,2n+1}
        long kshift = 1-(ir&1); // kshift=1 if ir even, 0 otherwise

        long t1 = jp+jm-nside_+kshift+1+nl4+nl4;
        long ip = (order_>0) ?
            (t1>>1)&(nl4-1) : ((t1>>1)%nl4); // in {0,4n-1}

        return ncap_ + (ir-1)*nl4 + ip;
        }
        else  // North & South polar caps
        {
        double tp = tt-long(tt);
        double tmp = ((za<0.99)||(!have_sth)) ?
                        nside_*sqrt(3*(1-za)) :
                        nside_*sth/sqrt((1.+za)/3.);

        long jp = long(tp*tmp); // increasing edge line index
        long jm = long((1.0-tp)*tmp); // decreasing edge line index

        long ir = jp+jm+1; // ring number counted from the closest pole
        long ip = long(tt*ir); // in {0,4*ir-1}
        planck_assert((ip>=0)&&(ip<4*ir),"must not happen");
        //ip %= 4*ir;

        return (z>0)  ?  2*ir*(ir-1) + ip  :  npix_ - 2*ir*(ir+1) + ip;
        }
}

inline long zphi2pix (double z, double phi)
      { return loc2pix(z,phi,0.,false); }

// get the ring above the specified
inline int ring_above (double z)
{
  double az=abs(z);
  if (az<=twothird) // equatorial region
    return nside_*(2-1.5*z);
  int iring = nside_*sqrt(3*(1-az));
  return (z>0) ? iring : 4*nside_-iring-1;
}

inline double ring2z (long ring)
{
  if (ring<nside_)
    return 1 - ring*ring*fact2_;
  if (ring <=3*nside_)
    return (2*nside_-ring)*fact1_;
  ring=4*nside_ - ring;
  return ring*ring*fact2_ - 1;
}

inline void get_ring_info_small(long ring, long &startpix,
  long &ringpix, bool &shifted)
{
  if (ring < nside_)
  {
    shifted = true;
    ringpix = 4*ring;
    startpix = 2*ring*(ring-1);
  }
  else if (ring < 3*nside_)
  {
    shifted = ((ring-nside_) & 1) == 0;
    ringpix = 4*nside_;
    startpix = ncap_ + (ring-nside_)*ringpix;
  }
  else
  {
    shifted = true;
    int nr= 4*nside_-ring;
    ringpix = 4*nr;
    startpix = npix_-2*nr*(nr+1);
  }
}

void query_multidisc(const stdvec<Vector3d> &norm,
                     const stdvec<double> &rad, stdvec<long> &pixset)
{
  int nv = norm.size(); // number of vertices

  planck_assert(nv == rad.size(), "inconsistent input arrays");

  int fct = 1;

  double rpsmall, rpbig;

  rpsmall = rpbig = 0;

  int irmin = 1, irmax = 4 * nside_ - 1;
  stdvec<double> z0, xa, cosrsmall, cosrbig;
  stdvec<angle_t> ptg;
  stdvec<int> cpix;

  /*
   *  Iterates over the normalized vertices of the identified disc (in our case,
   *  derived from the polygon edges of the registered pixels). This should not
   *  be a source of control divergence due to its uniformity and independence
   *  across the dataset
   */
  using namespace std;
  for (tsize i = 0; i < nv; ++i)
  {
    double rsmall = rad[i] + rpsmall;
    if (rsmall < pi)
    {
      double rbig = min(pi, rad[i] + rpbig);

      angle_t pnt;
      vec2ang(&pnt, &norm[i]);

      cosrsmall.push_back(cos(rsmall));
      cosrbig.push_back(cos(rbig));

      double cth = cos(pnt.theta);

      z0.push_back(cth);

      if (fct > 1)
        cpix.push_back(zphi2pix(cth, pnt.phi));

      xa.push_back(1. / sqrt((1 - cth) * (1 + cth)));
      ptg.push_back(pnt);

      double rlat1 = pnt.theta - rsmall;
      double zmax = cos(rlat1);
      int irmin_t = (rlat1 <= 0) ? 1 : ring_above(zmax) + 1;

      if ((fct > 1) && (rlat1 > 0))
        irmin_t = max(1, irmin_t - 1);

      double rlat2 = pnt.theta + rsmall;
      double zmin = cos(rlat2);
      long irmax_t = (rlat2 >= pi) ? 4 * nside_ - 1 : ring_above(zmin);

      if ((fct > 1) && (rlat2 < pi))
        irmax_t = min(4 * nside_ - 1, irmax_t + 1);

      if (irmax_t < irmax)
        irmax = irmax_t;
      if (irmin_t > irmin)
        irmin = irmin_t;
    }
  }

  /*
   *  Iterates over the identified rings to extract individual pixels. Also
   *  not a source of control divergence because each registered pixel boundary
   *  should encompass a uniform quantity of rings (with error of ~1)
   */
  for (long iz = irmin; iz <= irmax; ++iz)
  {
    double z = ring2z(iz);
    long ipix1, nr;
    bool shifted;
    get_ring_info_small(iz, ipix1, nr, shifted);
    double shift = shifted ? 0.5 : 0.;
    rangeset<I2> tr;
    tr.append(ipix1, ipix1 + nr);

    /*
     *  Granular ring iteration
     */
    for (tsize j = 0; j < z0.size(); ++j)
    {
      double x = (cosrbig[j] - z * z0[j]) * xa[j];
      double ysq = 1. - z * z - x * x;
      if (ysq > 0)
      {
        double dphi = atan2(sqrt(ysq), x);
        I ip_lo = ifloor<I>(nr * inv_twopi * (ptg[j].phi - dphi) - shift) + 1;
        I ip_hi = ifloor<I>(nr * inv_twopi * (ptg[j].phi + dphi) - shift);
        if (ip_hi >= nr)
        {
          ip_lo -= nr;
          ip_hi -= nr;
        }
        if (ip_lo < 0)
          tr.remove(ipix1 + ip_hi + 1, ipix1 + ip_lo + nr);
        else
          tr.intersect(ipix1 + ip_lo, ipix1 + ip_hi + 1);
      }
    }
    pixset.append(tr);
  }
}

void get_circle(const stdvec<Vector3d> &point, tsize q1, tsize q2, Vector3d &center,
                double &cosrad)
{
  center = (point[q1] + point[q2]).normalized();
  cosrad = point[q1].dot(center);
  for (tsize i = 0; i < q1; ++i)
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

void get_circle(const stdvec<Vector3d> &point, tsize q, Vector3d &center,
                double &cosrad)
{
  center = (point[0] + point[q]).normalized();
  cosrad = point[0].dot(center);
  for (tsize i = 1; i < q; ++i)
    if (point[i].dot(center) < cosrad) // point outside the current circle
      get_circle(point, i, q, center, cosrad);
}

void query_polygon(const stdvec<angle_t> &vertex, rangeset<I2> &pixset)
{
  tsize nv = vertex.size();
  tsize ncirc = nv;

  planck_assert(nv >= 3, "not enough vertices in polygon");

  // convert all pointing vectors to vec3
  stdvec<Vector3d> vv(nv);
  for (tsize i = 0; i < nv; ++i)
    ang2vec(&vv[i], &vertex[i]);
  stdvec<Vector3d> normal(ncirc);

  int flip = 0;

  for (tsize i = 0; i < nv; ++i)
  {
    normal[i] = vv[i].cross(vv[(i + 1) % nv]).normalized();

    double hnd = normal[i].dot(vv[(i + 2) % nv]);

    planck_assert(abs(hnd) > 1e-10, "degenerate corner");
    if (i == 0)
      flip = (hnd < 0.) ? -1 : 1;
    else
      planck_assert(flip * hnd > 0, "polygon is not convex");
    normal[i] *= flip;
  }

  stdvec<double> rad(ncirc, halfpi);

  query_multidisc(normal, rad, pixset);
}

void register_pixel(const char *imageData, Vector2i *imageSize, int x, int y,
                               camera_t *deviceCamProps, char *deviceMap, int nside)
{
  rgb_t cur_pix;

  pixeldata pd;

  /*
   *  Register pixel to polar coordinate space
   */
  pd.theta = (y - imageSize->y() / 2) * deviceCamProps->fov.theta / imageSize->y() + deviceCamProps->off.theta;
  pd.phi = (x - imageSize->x() / 2) * deviceCamProps->fov.phi / imageSize->x() + deviceCamProps->off.phi;

  cur_pix = ((rgb_t *)imageData)[y * imageSize->x() + x];
  pd.value = (cur_pix.r + cur_pix.g + cur_pix.b) / (3);

  double pixel_size_y = deviceCamProps->fov.theta / imageSize->y();
  double pixel_size_x = deviceCamProps->fov.phi / imageSize->x();

  // double a = std::fmod((angle->theta - pixel_size_y/2), M_PI);

  stdvec<angle_t> vecs(4);

  /*
   *  Find the boundaries of the pixel in spherical space
   */
  
   vecs[0].theta = wrapping_mod((pd.theta - pixel_size_y / 2), M_PI);
   vecs[0].phi = wrapping_mod((pd.phi - pixel_size_x / 2), 2 * M_PI);

   vecs[1].theta = wrapping_mod((pd.theta + pixel_size_y / 2), M_PI);
   vecs[1].phi = wrapping_mod((pd.phi - pixel_size_x / 2), 2 * M_PI);

   vecs[2].theta = wrapping_mod((pd.theta + pixel_size_y / 2), M_PI);
   vecs[2].phi = wrapping_mod((pd.phi + pixel_size_x / 2), 2 * M_PI);

   vecs[3].theta = wrapping_mod((pd.theta - pixel_size_y / 2), M_PI);
   vecs[3].phi = wrapping_mod((pd.phi + pixel_size_x / 2), 2 * M_PI);

   long *pixels;

  // find corresponding HEALPix pixel indices
  query_polygon(vecs, pixels);

  // update spherical pixel with corresponding image value
  for (int j = 0; j < pixels->size(); ++j)
  {
    for (int k = pixels->ivbegin(j); k < pixels->ivend(j); ++k)
    {
      // printf("processing pixel #%d with value %f\n", k, vecs[i].value);
      stack_hp(map, k, vecs[i].value);
    }
    // (*map)[pd[i]] = vecs[i].value;
  }
}

int main(int argc, char *argv[])
{
  int nside, mapSize;
  int width, height, channels;
  int imageSize;
  camera_t cam_props;
  Healpix_Map<int> *map;

  std::string file_loc = "image.jpg";

  if (argc > 1)
  {
    file_loc = argv[1];
  }

  nside = NSIDE;

  mapSize = 12 * pow((long)nside, 2);

  camProps = (camera_t *)malloc(sizeof(camera_t));

  // camera FOV and rotational offset (theta, phi)
  camProps->fov = pointing(M_PI / 4, M_PI / 4);
  camProps->off = pointing off(M_PI / 2, 0);

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

  register_image<<<dimGrid, dimBlock>>>(
      deviceInput,
      imageSize,
      deviceCamProps,
      deviceMap,
      nside);

  cudaMemcpy(hostMap, deviceMap, mapSize, cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();

  cudaFree(deviceInput);
  cudaFree(deviceCamProps);
  cudaFree(deviceMap);

  free(hostMap);
  free(imageData);

  return 0;
}