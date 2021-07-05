#include <cuda.h>
#include <cmath>
#include <cstdlib>
#include <cub/cub.cuh>

#include "helper_cuda.h"
#include "core.h"
#include "parameters.h"

__host__ __device__ uint div_ceil(uint a, uint b)
{
    return (a + b - 1) / b;
}

__device__ double minimum_image(double cordi, const double cell_length)
{
    if (cordi <= -0.5 * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5 * cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}

__global__ void forces_cuda_1(const double* rxyz, double* fxyz, const double L,
                                const int i, double* epot_array, double* pres_vir_array) //, double* flop)
{
    uint gtid = blockIdx.x * blockDim.x + threadIdx.x;

    double epot_gtid = 0;
    double pres_vir_gtid = 0;

    if ((i < N) && (gtid < N) && (gtid != i)) {
    // if ((i < N-1) && (gtid < N) && (i < gtid)) {
        double xi = rxyz[i];
        double yi = rxyz[i + N];
        double zi = rxyz[i + 2*N];

        double xj = rxyz[gtid];
        double yj = rxyz[gtid + N];
        double zj = rxyz[gtid + 2*N];

        double rx = xi - xj;
        rx = minimum_image(rx, L);
        double ry = yi - yj;
        ry = minimum_image(ry, L);
        double rz = zi - zj;
        rz = minimum_image(rz, L);

        double rij2 = rx * rx + ry * ry + rz * rz;

        // *flop += 14;

        double rcut2 = RCUT * RCUT;

        if (rij2 <= rcut2) {
            double r2inv = 1.0 / rij2;
            double r6inv = r2inv * r2inv * r2inv;
            double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

            fxyz[gtid] -= fr * rx;
            fxyz[gtid + N] -= fr * ry;
            fxyz[gtid + 2*N] -= fr * rz;

            if (gtid>i){
                epot_gtid = 4.0 * r6inv * (r6inv - 1.0) - ECUT;
                pres_vir_gtid = fr * rij2;
            }
            // epot_gtid = 4.0 * r6inv * (r6inv - 1.0) - ECUT;
            // pres_vir_gtid = fr * rij2;
            // *flop += 27; // 24 operacioens ?
        }
    }

    typedef cub::BlockReduce<double, BLOCK_SIZE> ReduceDobles;
    __shared__ typename ReduceDobles::TempStorage  temp_storage;

    double epot_blocks = ReduceDobles(temp_storage).Sum(epot_gtid);
    __syncthreads();
    double pres_vir_blocks = ReduceDobles(temp_storage).Sum(pres_vir_gtid);

    if (threadIdx.x == 0) {
        uint bid = blockIdx.x;
        epot_array[bid] = epot_blocks;
        pres_vir_array[bid] = pres_vir_blocks;
    }

}

void launch_forces_cuda_1(const double* rxyz, double* fxyz, const double L,
                            const int i, double* epot_array, double* pres_vir_array) // float elapsed_k, double flop)
{
    // cudaEvent_t start_k, end_k;
    // checkCudaCall(cudaEventCreate(&start_k));
    // checkCudaCall(cudaEventCreate(&end_k));

    dim3 block(BLOCK_SIZE);
    dim3 grid(div_ceil(N, BLOCK_SIZE));

    // checkCudaCall(cudaEventRecord(start_k));
    forces_cuda_1<<<grid, block>>>(rxyz, fxyz, L, i, epot_array, pres_vir_array); // &flop);
    checkCudaCall(cudaGetLastError());
    // checkCudaCall(cudaEventRecord(end_k));

    checkCudaCall(cudaDeviceSynchronize());

    // checkCudaCall(cudaEventElapsedTime(&elapsed_k, start_k, end_k));

    // checkCudaCall(cudaEventDestroy(start_k));
    // checkCudaCall(cudaEventDestroy(end_k));

}

__global__ void forces_cuda_2(double* epot_array, double* pres_vir_array,
                                double* epot_qdah, double* pres_vir_qdah) //, double* flop)
{
    // for (int idx = 0; idx< div_ceil(N, BLOCK_SIZE); idx++){
    //     *epot_qdah += epot_array[idx];
    //     *pres_vir_qdah += pres_vir_array[idx];
    // }

    uint gtid = blockIdx.x * blockDim.x + threadIdx.x;

    double epot_gtid = 0;
    double pres_vir_gtid = 0;

    if (gtid < div_ceil(N, BLOCK_SIZE)) {
        epot_gtid = epot_array[gtid];
        pres_vir_gtid = pres_vir_array[gtid];
    }
    // hasta acá viene bien
    typedef cub::BlockReduce<double, BLOCK_SIZE> ReduceDobles;
    __shared__ typename ReduceDobles::TempStorage  temp_storage;

    double epot_aux = ReduceDobles(temp_storage).Sum(epot_gtid);
    __syncthreads();
    double pres_vir_aux = ReduceDobles(temp_storage).Sum(pres_vir_gtid);

    if (threadIdx.x == 0) {
        *epot_qdah = epot_aux;
        *pres_vir_qdah = pres_vir_aux;
    }
}

void launch_forces_cuda_2(double* epot_array, double* pres_vir_array,
                                double* epot_qdah, double* pres_vir_qdah) //, double* flop)
{
    // cudaEvent_t start_k, end_k;
    // checkCudaCall(cudaEventCreate(&start_k));
    // checkCudaCall(cudaEventCreate(&end_k));

    dim3 block(BLOCK_SIZE);
    dim3 grid(div_ceil(div_ceil(N, BLOCK_SIZE), BLOCK_SIZE));

    // checkCudaCall(cudaEventRecord(start_k));
    forces_cuda_2<<<grid, block>>>(epot_array, pres_vir_array, epot_qdah, pres_vir_qdah); //, flop)
    checkCudaCall(cudaGetLastError());
    // checkCudaCall(cudaEventRecord(end_k));

    checkCudaCall(cudaDeviceSynchronize());

    // checkCudaCall(cudaEventElapsedTime(&elapsed_k, start_k, end_k));

    // checkCudaCall(cudaEventDestroy(start_k));
    // checkCudaCall(cudaEventDestroy(end_k));

}
