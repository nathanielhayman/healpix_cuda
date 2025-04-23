## Compilation Solutions

To properly identify cfitsio and libsharp when running `./configure` for `healpix_cxx`, include the following:

```
CFITSIO_CFLAGS="-I/Users/nathanielhayman/work/healpix_cuda/src/cfitsio/include"
CFITSIO_LIBS="-L/Users/nathanielhayman/work/healpix_cuda/src/cfitsio/lib"

SHARP_CFLAGS="-I/Users/nathanielhayman/work/healpix_cuda/src/libsharp/auto/include"
SHARP_LIBS="-L/Users/nathanielhayman/work/healpix_cuda/src/libsharp/auto/lib -lfftpack -lc_utils"
```

https://gitlab.mpcdf.mpg.de/mtr/levels/-/tree/master

g++ -I./Healpix_3.83/include/healpix_cxx -I./cfitsio/include -L./Healpix_3.83/lib -L./cfitsio/lib -lhealpix_cxx -lsharp -lcfitsio register.cpp stb_impl.cpp -o reg
