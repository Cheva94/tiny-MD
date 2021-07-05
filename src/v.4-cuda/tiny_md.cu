#include <cuda.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "helper_cuda.h"
#include "core.h"
#include "parameters.h"
#include "kernel.h"
#include "wtime.h"

int main()
{
    // FILE *file_metrica;
    // file_metrica = fopen("metrica.temp","w");
    // fprintf(file_metrica,"GFLOPS\n");

    double Ekin, Epot, Temp, Pres;
    double Rho, cell_V, cell_L, tail, Etail, Ptail;
    size_t array_size = 3 * N * sizeof(double);
    double* rxyz = nullptr;
    double* vxyz = nullptr;
    double* fxyz = nullptr;

    checkCudaCall(cudaMallocManaged(&rxyz, array_size));
    checkCudaCall(cudaMallocManaged(&vxyz, array_size));
    checkCudaCall(cudaMallocManaged(&fxyz, array_size));

    printf("# Tamaño del bloque:         %d\n", BLOCK_SIZE);
    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", TEQ);
    printf("# Pasos de medición:         %d\n", TRUN - TEQ);
    printf("# (mediciones cada %d pasos)\n", TMES);
    printf("# densidad, volumen, energía potencial media, presión media\n");

    srand(SEED);
    double t = 0.0, sf;
    double Rhob;
    Rho = RHOI;

    init_pos(rxyz, Rho);

    double start = wtime();

    // double f_force = 0;

    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = RHOI - 0.1 * (double)m;
        cell_V = (double)N / Rho;
        cell_L = cbrt(cell_V);
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9)
                - pow(RCUT, -3)) / 3.0;
        Etail = tail * (double)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrt(Rhob / Rho);

        for (int k = 0; k < 3 * N; k++) {
            rxyz[k] *= sf;
        }

        init_vel(vxyz, &Temp, &Ekin);

        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L); //, &f_force);

        for (i = 1; i < TEQ; i++) {
            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho,
                            cell_V, cell_L); // &f_force);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) {
                vxyz[k] *= sf;
            }
        }

        int mes = 0;
        double epotm = 0.0, presm = 0.0;

        for (i = TEQ; i < TRUN; i++) {
            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho,
                            cell_V, cell_L); // &f_force);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) {
                vxyz[k] *= sf;
            }

            if (i % TMES == 0) {
                Epot += Etail;
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;
            }
            t += DT;
        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (double)mes,
                presm / (double)mes);
    }

    double elapsed = wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);

    // double gf_force = f_force/(9+TRUN); // Contador de flops de forces

    // fprintf(file_metrica,"%f\n",gf_force);

    checkCudaCall(cudaFree(rxyz));
    checkCudaCall(cudaFree(vxyz));
    checkCudaCall(cudaFree(fxyz));

    return 0;
}
