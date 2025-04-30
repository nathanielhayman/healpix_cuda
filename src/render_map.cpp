#include "healpix_map.h"
#include "healpix_base.h"
#include "healpix_map_fitsio.h"
#include "alloc_utils.h"
#include "arr.h"

#include <string>
#include <fstream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("error: filename not supplied\n");
        return -1;
    }

    const std::string filename = argv[1];

    std::vector<int> arr;
    std::ifstream inputFile(filename, std::ios::binary);

    int order;
    long size;

    if (inputFile.is_open()) {
        inputFile.read(reinterpret_cast<char*>(&order), sizeof(order)); // read order from first line
        inputFile.read(reinterpret_cast<char*>(&size), sizeof(size)); // get data array size
        arr.resize(size);
        inputFile.read(reinterpret_cast<char*>(arr.data()), size * sizeof(int)); // Read array elements
        inputFile.close();
    } else {
        std::cerr << "Error opening file" << std::endl;
    }

    Healpix_Map<int>* map = new Healpix_Map<int>(order, RING);

    // copy data to the map
    for (int i = 0; i < size; i++) {
        (*map)[i] = arr[i];
    }

    write_Healpix_map_to_fits(
        "output.fits", *map, PLANCK_FLOAT64
    );

    printf("Wrote map to file!\n");

    free(map);
    
    return 0;
}