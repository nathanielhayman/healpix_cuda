# HEALPix imCUDA
Provides tooling for projecting 2D images onto a HEALPix spherical map.

## Compilation

### Custom Compilation
To build the file via the command line (without the Makefile), use the following format:

g++ -I./Healpix_3.83/include/healpix_cxx -I./cfitsio/include -L./Healpix_3.83/lib -L./cfitsio/lib -lhealpix_cxx -lsharp -lcfitsio register.cpp stb_impl.cpp -o reg


### Dependencies
The cfitsio library can be found here: 
```
https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.6.2.tar.gz
```

The HEALPix library can be found here:
```
https://sourceforge.net/projects/healpix/files/Healpix_3.83/Healpix_3.83_2024Nov13.tar.gz/download
```

### Building Dependencies 

Run the following commands to build the proper dependencies:

```
wget <cfitsio source>

wget <healpix source>

cd cfitsio-x-x-x

./configure

make

make install

cd ..

cd Healpix-xxx
```

Use the following as the path to `fitsio.h`:
```
/<dir>/healpix_cuda/src/cfitsio/include
```

And for the cfitsio library:
```
/<dir>/healpix_cuda/src/cfitsio/lib
```