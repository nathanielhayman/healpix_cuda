void interpolate_spherical(float *image, rangeset<double> &hp_map, Camera cam, int width, int height) {
    rangeset<double> corners(4);

    // register pixel corner destinations
    mapToHealpix(image, hp_map, cam, width, height);

    vector<double> pixset;

    // identify overlapping spherical pixels
    query_polygon(vertex, pixset);

    // perform bilinear interpolation over spherical pixels given neighboring 2D pixels
    for (int i = 0; i < pixset.size(); i++) {
        
    }
}

int main(int argc, char** argv) {

    return 0;
}