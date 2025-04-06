#include <stdio.h>
#include <math.h>

#define NPIX  1000000   // Number of HEALPix points
#define NSIDE 512       // HEALPix NSIDE parameter
#define PI 3.14159265358979323846
#define INV_PI (1.0 / PI)

// Camera Parameters
struct Camera {
    float elevation;  // Rotation in elevation (degrees)
    float azimuth;    // Rotation in azimuth (degrees)
    float x, y, z;    // Absolute position in space
    float fov;        // Field of View (degrees)
};

// Converts degrees to radians
__device__ float degToRad(float degrees) {
    return degrees * (PI / 180.0);
}

// CUDA Kernel: Maps 2D Image to HEALPix
__global__ void mapToHealpix(float *image, int *healpix_output, Camera cam, int width, int height) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height) return;

    // Normalize coordinates (-1 to 1)
    float u = (2.0f * x / width - 1.0f);
    float v = (2.0f * y / height - 1.0f);

    // Convert to spherical coordinates
    float theta = cam.fov * v;  // Vertical FOV scaling
    float phi   = cam.azimuth + cam.fov * u; // Azimuthal shift

    // Convert to Cartesian 3D space
    float px = cos(degToRad(theta)) * cos(degToRad(phi));
    float py = sin(degToRad(theta));
    float pz = cos(degToRad(theta)) * sin(degToRad(phi));

    // Apply camera position offset
    px += cam.x;
    py += cam.y;
    pz += cam.z;

    // Convert to HEALPix index (simplified)
    int healpix_idx = (int)((theta / PI) * height) * width + (int)((phi / (2 * PI)) * width);
    healpix_idx = healpix_idx % (width * height);

    // Store result
    healpix_output[y * width + x] = healpix_idx;
}

// CUDA Kernel: Parallel ang2pix
__global__ void ang2pix(float *theta, float *phi, int *pix, int nside, int npix) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= npix) return;

    float t = theta[i];
    float p = phi[i];

    // HEALPix resolution
    int npface = 4 * nside * nside;  
    float z = cos(t);  
    float za = fabs(z);
    
    int face_num, ix, iy;
    if (za <= 2.0 / 3.0) {  // Equatorial region
        float temp = nside * (0.5 + p * INV_PI);
        int j = (int)(nside * (0.75 + 0.5 * z));
        face_num = (int)(temp / nside);
        ix = (int)(temp - face_num * nside);
        iy = (int)(j);
    } else {  // Polar caps
        float tt = nside * sqrt(3 * (1 - za));
        int jp = (int)(tt);
        int jm = (int)(2 * tt - jp);
        face_num = (z > 0) ? 2 : 3;
        ix = (int)(p * INV_PI * jm);
        iy = (int)(jp);
    }

    // Compute HEALPix index
    pix[i] = face_num * npface + ix * nside + iy;
}

__global__ void pix2ang(float *theta, float *phi, int *pix, int nside, int npix) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= npix) return;

    int face_num = pix[i] / (4 * nside * nside);
    int ix = (pix[i] % (4 * nside * nside)) / nside;
    int iy = (pix[i] % nside);

    // Convert HEALPix index to angle
    float za = 2.0 / 3.0;
    if (face_num < 2) {  // Equatorial region
        float j = iy / (float)nside - 0.75;
        float temp = ix / (float)nside;
        theta[i] = acos(j);
        phi[i] = temp * 2 * PI;
    } else {  // Polar caps
        float tt = sqrt(3 * iy / (float)nside);
        theta[i] = (face_num == 2) ? acos(1 - tt) : acos(tt - 1);
        phi[i] = (ix / (float)nside) * 2 * PI;
    }
}

// Host function
int main() {

    /* MAIN IDEA FOR NOW :    
        CONVERT MAP FROM IMAGE -> HEALPix Map
        Then HEALPix Map -> Angle !    
    */
    
    int npix = NPIX;

    // Allocate host memory
    float *h_theta = (float*)malloc(npix * sizeof(float));
    float *h_phi   = (float*)malloc(npix * sizeof(float));
    int *h_pix     = (int*)malloc(npix * sizeof(int));

    // Initialize random angles
    for (int i = 0; i < npix; i++) {
        h_theta[i] = degToRad((rand() % 180));  // 0 to 180 degrees
        h_phi[i]   = degToRad((rand() % 360));  // 0 to 360 degrees
    }

    // Allocate device memory
    float *d_theta, *d_phi;
    int *d_pix;
    cudaMalloc(&d_theta, npix * sizeof(float));
    cudaMalloc(&d_phi, npix * sizeof(float));
    cudaMalloc(&d_pix, npix * sizeof(int));

    // Copy data to device
    cudaMemcpy(d_theta, h_theta, npix * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_phi, h_phi, npix * sizeof(float), cudaMemcpyHostToDevice);

    // Launch kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (npix + threadsPerBlock - 1) / threadsPerBlock;
    ang2pix<<<blocksPerGrid, threadsPerBlock>>>(d_theta, d_phi, d_pix, NSIDE, npix);
    
    // pix2ang<<<blocksPerGrid, threadsPerBlock>>>(d_theta, d_phi, d_pix, NSIDE, npix);

    /*
    float *h_image = (float*)malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(float));
    int *h_healpix = (int*)malloc(IMG_WIDTH * IMG_HEIGHT * sizeof(int));

    // Initialize image with some values
    for (int i = 0; i < IMG_WIDTH * IMG_HEIGHT; i++) {
        h_image[i] = (float)(i % 256) / 255.0;  // Example grayscale intensity
    }

    // Allocate device memory
    float *d_image;
    int *d_healpix;
    cudaMalloc(&d_image, IMG_WIDTH * IMG_HEIGHT * sizeof(float));
    cudaMalloc(&d_healpix, IMG_WIDTH * IMG_HEIGHT * sizeof(int));

    // Copy data to device
    cudaMemcpy(d_image, h_image, IMG_WIDTH * IMG_HEIGHT * sizeof(float), cudaMemcpyHostToDevice);

    // Set up Camera
    Camera h_cam = {45.0f, 30.0f, 0.0f, 0.0f, 0.0f, 90.0f};
    Camera *d_cam;
    cudaMalloc(&d_cam, sizeof(Camera));
    cudaMemcpy(d_cam, &h_cam, sizeof(Camera), cudaMemcpyHostToDevice);

    // Launch kernel
    dim3 blockSize(16, 16);
    dim3 gridSize((IMG_WIDTH + blockSize.x - 1) / blockSize.x, (IMG_HEIGHT + blockSize.y - 1) / blockSize.y);
    mapToHealpix<<<gridSize, blockSize>>>(d_image, d_healpix, *d_cam, IMG_WIDTH, IMG_HEIGHT);

    // Copy result back
    cudaMemcpy(h_healpix, d_healpix, IMG_WIDTH * IMG_HEIGHT * sizeof(int), cudaMemcpyDeviceToHost);

    // Print some sample results
    for (int i = 0; i < 10; i++) {
        printf("Pixel %d -> HEALPix index %d\n", i, h_healpix[i]);
    }
 
    */


    // Copy result back
    cudaMemcpy(h_pix, d_pix, npix * sizeof(int), cudaMemcpyDeviceToHost);

    // Print some results
    for (int i = 0; i < 10; i++) {
        printf("Theta: %.2f, Phi: %.2f -> HEALPix Index: %d\n", h_theta[i], h_phi[i], h_pix[i]);
    }

    // Cleanup
    free(h_theta);
    free(h_phi);
    free(h_pix);
    cudaFree(d_theta);
    cudaFree(d_phi);
    cudaFree(d_pix);

    return 0;
}
