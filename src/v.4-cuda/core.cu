#include <cuda.h>
#include <cmath>
#include <cstdlib>
#include <cub/cub.cuh>

#include "helper_cuda.h"
#include "core.h"
#include "kernel.h"
#include "parameters.h"
// #include "wtime.h"

void init_pos(double* rxyz, const double rho)
{
    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
                // Eje x
                rxyz[idx + 0] = i * a;
                rxyz[idx + 1] = (i + 0.5) * a;
                rxyz[idx + 2] = (i + 0.5) * a;
                rxyz[idx + 3] = i * a;

                // Eje y
                rxyz[idx + 0 + N] = j * a;
                rxyz[idx + 1 + N] = (j + 0.5) * a;
                rxyz[idx + 2 + N] = j * a;
                rxyz[idx + 3 + N] = (j + 0.5) * a;

                // Eje z
                rxyz[idx + 0 + 2*N] = k * a;
                rxyz[idx + 1 + 2*N] = k * a;
                rxyz[idx + 2 + 2*N] = (k + 0.5) * a;
                rxyz[idx + 3 + 2*N] = (k + 0.5) * a;

                idx += 4;
            }
        }
    }
}

void init_vel(double* vxyz, double* temp, double* ekin)
{
    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    for (int i = 0; i < 3 * N; i += 3) {
        vxyz[i] = rand() / (double)RAND_MAX - 0.5;
        vxyz[i + N] = rand() / (double)RAND_MAX - 0.5;
        vxyz[i + 2*N] = rand() / (double)RAND_MAX - 0.5;

        sumvx += vxyz[i];
        sumvy += vxyz[i + N];
        sumvz += vxyz[i + 2*N];
        sumv2 += vxyz[i] * vxyz[i] + vxyz[i + N] * vxyz[i + N]
                + vxyz[i + 2*N] * vxyz[i + 2*N];
    }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < 3 * N; i += 3) {
        vxyz[i] = (vxyz[i] - sumvx) * sf;
        vxyz[i + N] = (vxyz[i + N] - sumvy) * sf;
        vxyz[i + 2*N] = (vxyz[i + 2*N] - sumvz) * sf;
    }
}

void forces(const double* rxyz, double* fxyz, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)//, double *f_force)
{
    // double start_w = wtime();
    // double flop = 0;
    // float elapsed_k = 0;

    double pres_vir = 0;
    *epot = 0;

    for (int i = 0; i < 3 * N; i++) {
        fxyz[i] = 0.0;
    }

    // for (int i = 0; i < N; i++) {
    for (int i = 0; i < N-1; i++) {

        size_t total_blocks_size = div_ceil(N, BLOCK_SIZE) * sizeof(double);
        double* epot_array = nullptr;
        double* pres_vir_array = nullptr;

        double* epot_qdah = nullptr;
        double* pres_vir_qdah = nullptr;

        checkCudaCall(cudaMallocManaged(&epot_array, total_blocks_size));
        checkCudaCall(cudaMallocManaged(&pres_vir_array, total_blocks_size));

        checkCudaCall(cudaMallocManaged(&epot_qdah, sizeof(double)));
        checkCudaCall(cudaMallocManaged(&pres_vir_qdah, sizeof(double)));

        launch_forces_cuda_1(rxyz, fxyz, L, i, epot_array, pres_vir_array); // elapsed_k, flop);

        launch_forces_cuda_2(epot_array, pres_vir_array, epot_qdah, pres_vir_qdah);

        *epot += *epot_qdah;
        pres_vir += *pres_vir_qdah;

        checkCudaCall(cudaFree(epot_array));
        checkCudaCall(cudaFree(pres_vir_array));

        checkCudaCall(cudaFree(epot_qdah));
        checkCudaCall(cudaFree(pres_vir_qdah));
    }

    pres_vir /= (V * 3.0);
    *pres = *temp * rho + pres_vir;

    // double elapsed_w = wtime()-start_w;

    // *f_force += ((flop + 4) / elapsed_w)*0.000000001;
}

static double pbc(double cordi, const double cell_length)
{
    if (cordi <= 0) {
        cordi += cell_length;
    } else if (cordi > cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}

void velocity_verlet(double* rxyz, double* vxyz, double* fxyz, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L) //, double *f_force)
{
    for (int i = 0; i < N; i++) {
        rxyz[i] += vxyz[i] * DT + 0.5 * fxyz[i] * DT * DT;
        rxyz[i + N] += vxyz[i + N] * DT + 0.5 * fxyz[i + N] * DT * DT;
        rxyz[i + 2*N] += vxyz[i + 2*N] * DT + 0.5 * fxyz[i + 2*N] * DT * DT;

        rxyz[i] = pbc(rxyz[i], L);
        rxyz[i + N] = pbc(rxyz[i + N], L);
        rxyz[i + 2*N] = pbc(rxyz[i + 2*N], L);

        vxyz[i] += 0.5 * fxyz[i] * DT;
        vxyz[i + N] += 0.5 * fxyz[i + N] * DT;
        vxyz[i + 2*N] += 0.5 * fxyz[i + 2*N] * DT;
    }

    forces(rxyz, fxyz, epot, pres, temp, rho, V, L);//, f_force);

    double sumv2 = 0.0;

    for (int i = 0; i < N; i++) {
        vxyz[i] += 0.5 * fxyz[i] * DT;
        vxyz[i + N] += 0.5 * fxyz[i + N] * DT;
        vxyz[i + 2*N] += 0.5 * fxyz[i + 2*N] * DT;

        sumv2 += vxyz[i] * vxyz[i] + vxyz[i + N] * vxyz[i + N]
                + vxyz[i + 2*N] * vxyz[i + 2*N];
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
