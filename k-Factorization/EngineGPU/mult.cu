#include <cuda_runtime.h>
#include <iostream>

__global__ void logical_mul_kernel(const long long* __restrict__ A,
                                   const long long* __restrict__ B,
                                   long long* __restrict__ C,
                                   size_t N) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        // Logical multiplication = bitwise AND
        C[idx] = A[idx] & B[idx];
    }
}

__declspec(dllexport) void logical_mul(const long long* h_A, const long long* h_B, long long* h_C, int N, int nc) {
    long long* d_A, * d_B, * d_C;

    size_t bytes = N * sizeof(long long);
    cudaMalloc(&d_A, bytes);
    cudaMalloc(&d_B, bytes);
    cudaMalloc(&d_C, bytes);

    cudaMemcpy(d_A, h_A, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, bytes, cudaMemcpyHostToDevice);

    // Configure grid
    int blockSize = 32;//  256;
    int gridSize = (N + blockSize - 1) / blockSize;

    printf("Starting CUDA gpu  gridSize = %d  blockSize = %d\n", gridSize, blockSize);
    auto tm = clock();
    for (int i = 0; i < nc; i++) {
        logical_mul_kernel << <gridSize, blockSize >> > (d_A, d_B, d_C, N);
        cudaDeviceSynchronize();
    }
    printf("CUDA gpu x %d iterations  lenArray = %d: %ld ms\n", nc, N, clock() - tm);

    cudaMemcpy(h_C, d_C, bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}