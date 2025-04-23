#include <gputk.h>

#define BLOCK_DIM 16

/*! Writes diagnostic output and exits with an error status. */
#define planck_fail(msg) \
do { planck_failure__(__FILE__,__LINE__,PLANCK_FUNC_NAME__,msg); \
throw PlanckError(msg); } while(0)

/*! Throws a PlanckError without diagnostic message. */
#define planck_fail_quietly(msg) \
do { throw PlanckError(msg); } while(0)

/*! Writes diagnostic output and exits with an error status if \a testval
    is \a false. */
#define planck_assert(testval,msg) \
do { if (testval); else planck_fail(msg); } while(0)

template<typename I> template<typename I2>
  void T_Healpix_Base<I>::query_multidisc (const arr<vec3> &norm,
  const arr<double> &rad, int fact, rangeset<I2> &pixset) const
  {
  bool inclusive = (fact!=0);
  tsize nv=norm.size();
  planck_assert(nv==rad.size(),"inconsistent input arrays");
  pixset.clear();

  if (scheme_==RING)
    {
    I fct=1;
    if (inclusive)
      {
      planck_assert (((I(1)<<order_max)/nside_)>=fact,
        "invalid oversampling factor");
      fct = fact;
      }
    T_Healpix_Base b2;
    double rpsmall, rpbig;
    if (fct>1)
      {
      b2.SetNside(fct*nside_,RING);
      rpsmall = b2.max_pixrad();
      rpbig = max_pixrad();
      }
    else
      rpsmall = rpbig = inclusive ? max_pixrad() : 0;

    I irmin=1, irmax=4*nside_-1;
    vector<double> z0,xa,cosrsmall,cosrbig;
    vector<pointing> ptg;
    vector<I> cpix;



    /*
     *  Iterates over the normalized vertices of the identified disc (in our case,
     *  derived from the polygon edges of the registered pixels). This should not
     *  be a source of control divergence due to its uniformity and independence
     *  across the dataset
     */
    for (tsize i=0; i<nv; ++i)
      {
      double rsmall=rad[i]+rpsmall;
      if (rsmall<pi)
        {
        double rbig=min(pi,rad[i]+rpbig);
        pointing pnt=pointing(norm[i]);
        cosrsmall.push_back(cos(rsmall));
        cosrbig.push_back(cos(rbig));
        double cth=cos(pnt.theta);
        z0.push_back(cth);
        if (fct>1) cpix.push_back(zphi2pix(cth,pnt.phi));
        xa.push_back(1./sqrt((1-cth)*(1+cth)));
        ptg.push_back(pnt);

        double rlat1 = pnt.theta - rsmall;
        double zmax = cos(rlat1);
        I irmin_t = (rlat1<=0) ? 1 : ring_above (zmax)+1;

        if ((fct>1) && (rlat1>0)) irmin_t=max(I(1),irmin_t-1);

        double rlat2 = pnt.theta + rsmall;
        double zmin = cos(rlat2);
        I irmax_t = (rlat2>=pi) ? 4*nside_-1 : ring_above (zmin);

        if ((fct>1) && (rlat2<pi)) irmax_t=min(4*nside_-1,irmax_t+1);

        if (irmax_t < irmax) irmax=irmax_t;
        if (irmin_t > irmin) irmin=irmin_t;
        }
      }

    /*
     *  Iterates over the identified rings to extract individual pixels. Also
     *  not a source of control divergence because each registered pixel boundary
     *  should encompass a uniform quantity of rings (with error of ~1)
     */
    for (I iz=irmin; iz<=irmax; ++iz)
      {
      double z=ring2z(iz);
      I ipix1,nr;
      bool shifted;
      get_ring_info_small(iz,ipix1,nr,shifted);
      double shift = shifted ? 0.5 : 0.;
      rangeset<I2> tr;
      tr.append(ipix1,ipix1+nr);


      /*
       *  Granular ring iteration
       */
      for (tsize j=0; j<z0.size(); ++j)
        {
        double x = (cosrbig[j]-z*z0[j])*xa[j];
        double ysq = 1.-z*z-x*x;
        if (ysq>0)
          {
          double dphi = atan2(sqrt(ysq),x);
          I ip_lo = ifloor<I>(nr*inv_twopi*(ptg[j].phi-dphi) - shift)+1;
          I ip_hi = ifloor<I>(nr*inv_twopi*(ptg[j].phi+dphi) - shift);
          if (fct>1)
            {
            while ((ip_lo<=ip_hi) && check_pixel_ring
              (*this,b2,ip_lo,nr,ipix1,fct,z0[j],ptg[j].phi,cosrsmall[j],cpix[j]))
              ++ip_lo;
            while ((ip_hi>ip_lo) && check_pixel_ring
              (*this,b2,ip_hi,nr,ipix1,fct,z0[j],ptg[j].phi,cosrsmall[j],cpix[j]))
              --ip_hi;
            }
          if (ip_hi>=nr)
            { ip_lo-=nr; ip_hi-=nr;}
          if (ip_lo<0)
            tr.remove(ipix1+ip_hi+1,ipix1+ip_lo+nr);
          else
            tr.intersect(ipix1+ip_lo,ipix1+ip_hi+1);
          }
        }
      pixset.append(tr);
      }
    }
  else // scheme_ == NEST
    {
    int oplus = 0;
    if (inclusive)
      {
      planck_assert ((I(1)<<(order_max-order_))>=fact,
        "invalid oversampling factor");
      planck_assert ((fact&(fact-1))==0,
        "oversampling factor must be a power of 2");
      oplus=ilog2(fact);
      }
    int omax=order_+oplus; // the order up to which we test

    // TODO: ignore all disks with radius>=pi

    arr<T_Healpix_Base<I> > base(omax+1);
    arr3<double> crlimit(omax+1,nv,3);
    for (int o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      for (tsize i=0; i<nv; ++i)
        {
        crlimit(o,i,0) = (rad[i]+dr>pi) ? -1. : cos(rad[i]+dr);
        crlimit(o,i,1) = (o==0) ? cos(rad[i]) : crlimit(0,i,1);
        crlimit(o,i,2) = (rad[i]-dr<0.) ?  1. : cos(rad[i]-dr);
        }
      }

    vector<pair<I,int> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(I(11-i),0));

    int stacktop=0; // a place to save a stack position

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      I pix=stk.back().first;
      int o=stk.back().second;
      stk.pop_back();

      vec3 pv(base[o].pix2vec(pix));

      tsize zone=3;
      for (tsize i=0; i<nv; ++i)
        {
        double crad=dotprod(pv,norm[i]);
        for (tsize iz=0; iz<zone; ++iz)
          if (crad<crlimit(o,i,iz))
            if ((zone=iz)==0) goto bailout;
        }

      check_pixel (o, order_, omax, zone, pixset, pix, stk, inclusive,
        stacktop);
      bailout:;
      }
    }
  }

template<typename I> template<typename I2>
    void T_Healpix_Base<I>::query_polygon
    (const vector<pointing> &vertex, rangeset<I2> &pixset) const
{
    tsize nv=vertex.size();
    tsize ncirc = nv;

    planck_assert(nv>=3,"not enough vertices in polygon");
    
    vector<vec3> vv(nv);
    for (tsize i=0; i<nv; ++i)
    vv[i]=vertex[i].to_vec3();
    arr<vec3> normal(ncirc);

    int flip=0;
    
    /*
     *  As described above, the expected magnitude of \a nv is uniform
     */
    for (tsize i=0; i<nv; ++i) {
        normal[i]=crossprod(vv[i],vv[(i+1)%nv]).Norm();
        
        double hnd=dotprod(normal[i],vv[(i+2)%nv]);

        planck_assert(abs(hnd)>1e-10,"degenerate corner");
        if (i==0)
            flip = (hnd<0.) ? -1 : 1;
        else
            planck_assert(flip*hnd>0,"polygon is not convex");
        normal[i]*=flip;
    }

    arr<double> rad(ncirc,halfpi);

    if (inclusive) {
        double cosrad;
        find_enclosing_circle (vv, normal[nv], cosrad);
        rad[nv]=acos(cosrad);
    }

    query_multidisc(normal,rad,fact,pixset);
}

/*
 * Determine projection pointing vectors for pixels in a 2D image
 */
__global__ int register_pixels(unsigned char* img_data, int height, int width, int channels, 
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
__global__ int find_overlapping_pixels(pixeldata* angle, rangeset<int>* pixels, int width, 
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

__global__ void register_pixel(const char* image_data, const pointing* fov, 
    const pointing* off, const vec2* image_size)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    rgb_t cur_pix;

    pixeldata pd;

    /*
     *  Register pixel to polar coordinate space
     */
    pd.theta = (y - image_size->y/2) * fov->theta / image_size->y + off->theta; 
    pd.phi = (x - image_size->x/2) * fov->phi / image_size->x + off->phi;

    // if image is 3-channel, convert value to single channel via averaging
    if (!gray) {
        cur_pix = ((rgb_t*)img_data)[y * width + x];

        pd.value = (cur_pix.r + cur_pix.g + cur_pix.b)/(3);
    } else {
        pd.value = img_data[(y * width + x) * 3];
    }


    double pixel_size_y = fov->theta / height;
    double pixel_size_x = fov->phi / width;

    std::vector<pointing> vecs(4);

    // double a = std::fmod((angle->theta - pixel_size_y/2), M_PI);

    /*
     *  Find the boundaries of the pixel in spherical space
     */
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
    
    // find corresponding HEALPix pixel indices
    query_polygon(vecs, *pixels);

    // update spherical pixel with corresponding image value
    for (int j = 0; j < pixels->size(); ++j) {
        for (int k = pixels->ivbegin(j); k < pixels->ivend(j); ++k) {
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

int main(int argc, char *argv[]) {
    int nside, order;
    int width, height, channels;
    int image_size;

    const char* file_loc = "image.jpg";

    // read image from file system
    unsigned char* img_data = stbi_load(
        file_loc, &width, &height, &channels, 3
    );
    
    Healpix_Map<int>* map;

    // gpuTKArg_t args;

    // args = gpuTKArg_read(argc, argv);

    nside = NSIDE;

    char** p;

    if (argc > 1) {
        int nside = std::strtol(argv[1], p, 10);
    }

    nside = nside < MAX_NSIDE ? nside : NSIDE;

    order = T_Healpix_Base<int>::nside2order(nside);

    // camera FOV and rotational offset (theta, phi)
    pointing fov(M_PI/4, M_PI/4);
    pointing off(M_PI/2, 0);

    // HEALPix Base engine (integer indexing)
    T_Healpix_Base<int> hp(order, RING);

    // HEALPix Map (integer valuation)
    map = new Healpix_Map<int>(order, RING);

    // set size of image data for transfer
    image_size = width * height * sizeof(char);

    cudaMalloc((void **)&device_input, image_size);
    cudaMalloc((void **)&device_input, image_size);

    cudaDeviceSynchronize();

    // copy image to device
    cudaMemcpy(device_input, image_data, image_size, cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();

    dim3 dimBlock(BLOCK_DIM, BLOCK_DIM, 1);
    dim3 dimGrid((inputLength + dimBlock.x - 1)/dimBlock.x, 1, 1);

    register_image<<<dimGrid, dimBlock>>>(
        device_input,
        image_size,
        hp,
        map->
    );

    cudaMemcpy(hostBins, deviceBins, , cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
    //@@ Free the GPU memory here

    cudaFree(deviceInput);
    cudaFree(deviceBins);

    free(hostBins);
    free(hostInput);
    return 0;
}