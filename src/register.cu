#include "register.h"

#include <thrust/device_vector.h>

#include "stb_image.h"

#define BLOCK_DIM 16
#define ORDER 4
#define MAX_NSIDE 536870912
#define MAX_ORDER

static const order_ = ORDER;
static const nside_  = I(1)<<order_;
static const npface_ = nside_<<order_;
static const ncap_   = (npface_-nside_)<<1;
static const npix_   = 12*npface_;
static const fact2_  = 4./npix_;
static const fact1_  = (nside_<<1)*fact2_;

static const double twothird=2.0/3.0;
static const double pi=3.141592653589793238462643383279502884197;
static const double twopi=6.283185307179586476925286766559005768394;
static const double halfpi=1.570796326794896619231321691639751442099;
static const double inv_halfpi=0.6366197723675813430755350534900574;

/*! Writes diagnostic output and exits with an error status. */
#define planck_fail(msg)                                           \
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


__device__ void vec2ang(vec3 *vec, angle_t *angle) {
  angle->theta = acos(vec->z);
  angle->phi = ((0 < val) - (val < 0)) * acos(vec->x/(sqrt(pow(vec->x, 2), pow(vec->y, 2))));
}

// get the ring above the specified
__forceinline__ inline int ring_above (double z, int nside) const
{
  double az=abs(z);
  if (az<=twothird) // equatorial region
    return I(nside*(2-1.5*z));
  int iring = I(nside*sqrt(3*(1-az)));
  return (z>0) ? iring : 4*nside-iring-1;
}

__forceinline__ inline double ring2z (long ring) const
{
  if (ring<nside_)
    return 1 - ring*ring*fact2_;
  if (ring <=3*nside_)
    return (2*nside_-ring)*fact1_;
  ring=4*nside_ - ring;
  return ring*ring*fact2_ - 1;
}

__forceinline__ inline void get_ring_info_small(long ring, long &startpix,
  long &ringpix, bool &shifted) const
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

__device__ void query_multidisc(int nside, const thrust::device_vector<vec3> &norm,
                     const thrust::device_vector<double> &rad, thrust::device_vector<long> &pixset)
{
  int nv = norm.size(); // number of vertices

  planck_assert(nv == rad_l, "inconsistent input arrays");

  int fct = 1;

  double rpsmall, rpbig;

  rpsmall = rpbig = 0;

  int irmin = 1, irmax = 4 * nside - 1;
  thrust::device_vector<double> z0, xa, cosrsmall, cosrbig;
  thrust::device_vector<angle_t> ptg;
  thrust::device_vector<int> cpix;

  /*
   *  Iterates over the normalized vertices of the identified disc (in our case,
   *  derived from the polygon edges of the registered pixels). This should not
   *  be a source of control divergence due to its uniformity and independence
   *  across the dataset
   */
  for (tsize i = 0; i < nv; ++i)
  {
    double rsmall = rad[i] + rpsmall;
    if (rsmall < pi)
    {
      double rbig = min(pi, rad[i] + rpbig);

      angle_t pnt;
      vec2ang(&norm[i], &pnt);

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
        irmin_t = max(I(1), irmin_t - 1);

      double rlat2 = pnt.theta + rsmall;
      double zmin = cos(rlat2);
      I irmax_t = (rlat2 >= pi) ? 4 * nside_ - 1 : ring_above(zmin);

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
  for (I iz = irmin; iz <= irmax; ++iz)
  {
    double z = ring2z(iz);
    int ipix1, nr;
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

void get_circle(const vector<vec3> &point, tsize q1, tsize q2, vec3 &center,
                double &cosrad)
{
  center = (point[q1] + point[q2]).Norm();
  cosrad = dotprod(point[q1], center);
  for (tsize i = 0; i < q1; ++i)
    if (dotprod(point[i], center) < cosrad) // point outside the current circle
    {
      center = crossprod(point[q1] - point[i], point[q2] - point[i]).Norm();
      cosrad = dotprod(point[i], center);
      if (cosrad < 0)
      {
        center.Flip();
        cosrad = -cosrad;
      }
    }
}

void get_circle(const vector<vec3> &point, tsize q, vec3 &center,
                double &cosrad)
{
  center = (point[0] + point[q]).Norm();
  cosrad = dotprod(point[0], center);
  for (tsize i = 1; i < q; ++i)
    if (dotprod(point[i], center) < cosrad) // point outside the current circle
      get_circle(point, i, q, center, cosrad);
}

template <typename I, typename I2>
void query_polygon(const std::vector<pointing> &vertex, rangeset<I2> &pixset)
{
  tsize nv = vertex.size();
  tsize ncirc = nv;

  planck_assert(nv >= 3, "not enough vertices in polygon");

  // convert all pointing vectors to vec3
  vector<vec3> vv(nv);
  for (tsize i = 0; i < nv; ++i)
    vv[i] = vertex[i].to_vec3();
  arr<vec3> normal(ncirc);

  int flip = 0;

  for (tsize i = 0; i < nv; ++i)
  {
    normal[i] = crossprod(vv[i], vv[(i + 1) % nv]).Norm();

    double hnd = dotprod(normal[i], vv[(i + 2) % nv]);

    planck_assert(abs(hnd) > 1e-10, "degenerate corner");
    if (i == 0)
      flip = (hnd < 0.) ? -1 : 1;
    else
      planck_assert(flip * hnd > 0, "polygon is not convex");
    normal[i] *= flip;
  }

  arr<double> rad(ncirc, halfpi);

  if (inclusive)
  {
    double cosrad;

    tsize np = point.size();
    planck_assert(np >= 2, "too few points");
    center = (point[0] + point[1]).Norm();
    cosrad = dotprod(point[0], center);
    for (tsize i = 2; i < np; ++i)
      if (dotprod(point[i], center) < cosrad) // point outside the current circle
        get_circle(point, i, center, cosrad);

    rad[nv] = acos(cosrad);
  }

  query_multidisc(normal, rad, fact, pixset);
}

__global__ void register_pixel(const char *imageData, int imageSize,
                               camera_t *deviceCamProps, char *deviceMap, int nside)
{
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;

  rgb_t cur_pix;

  pixeldata pd;

  /*
   *  Register pixel to polar coordinate space
   */
  pd.theta = (y - image_size->y / 2) * fov->theta / image_size->y + off->theta;
  pd.phi = (x - image_size->x / 2) * fov->phi / image_size->x + off->phi;

  cur_pix = ((rgb_t *)image_data)[y * image_size->x + x];
  pd.value = (cur_pix.r + cur_pix.g + cur_pix.b) / (3);

  double pixel_size_y = fov->theta / image_size->y;
  double pixel_size_x = fov->phi / image_size->x;

  // double a = std::fmod((angle->theta - pixel_size_y/2), M_PI);

  /*
   *  Find the boundaries of the pixel in spherical space
   */
  vecs[0] = pointing(
      wrapping_mod((angle->theta - pixel_size_y / 2), M_PI),
      wrapping_mod((angle->phi - pixel_size_x / 2), 2 * M_PI));

  vecs[1] = pointing(
      wrapping_mod((angle->theta + pixel_size_y / 2), M_PI),
      wrapping_mod((angle->phi - pixel_size_x / 2), 2 * M_PI));

  vecs[2] = pointing(
      wrapping_mod((angle->theta + pixel_size_y / 2), M_PI),
      wrapping_mod((angle->phi + pixel_size_x / 2), 2 * M_PI));

  vecs[3] = pointing(
      wrapping_mod((angle->theta - pixel_size_y / 2), M_PI),
      wrapping_mod((angle->phi + pixel_size_x / 2), 2 * M_PI));

  // find corresponding HEALPix pixel indices
  // query_polygon(vecs, *pixels);

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

  free(pixels);

  // free the allocated memory
  delete[] vecs;

  return 0;
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