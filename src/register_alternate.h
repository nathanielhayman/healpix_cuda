struct pixeldata {
    double theta;
    double phi;
    int value;
};

struct rgb_t {
    char r;
    char g;
    char b;
};

struct angle_t {
    double theta;
    double phi;
};

struct camera_t {
    angle_t fov;
    angle_t off;
    int n_channels;
};

struct hprange_t {
    long start;
    long end;
    long ring;
};

struct hpbound_t {
    int rstart;
    int rend;
    int data_size;
    hprange_t* pixels;
};