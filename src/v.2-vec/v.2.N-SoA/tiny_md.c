#include "parameters.h"
#include "core.h"
#include "wtime.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
    FILE *file_metrica;
    file_metrica = fopen("metrica.temp","w");
    fprintf(file_metrica,"GFLOPS\n");

    double Ekin, Epot, Temp, Pres;
    double Rho, cell_V, cell_L, tail, Etail, Ptail;

    struct SoA * rxyz = malloc(sizeof(struct SoA));
    struct SoA * fxyz = malloc(sizeof(struct SoA));
    struct SoA * vxyz = malloc(sizeof(struct SoA));

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

    // NO AUTOVECTORIZADO: loops anidados.
    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = RHOI - 0.1 * (double)m;
        cell_V = (double)N / Rho;
        cell_L = cbrt(cell_V);
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0;
        Etail = tail * (double)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrt(Rhob / Rho);

        // AUTOVECTORIZADO
        for (int k = 0; k < N; k++) { // reescaleo posiciones a nueva densidad
            rxyz-> x[k] *= sf;
            rxyz-> y[k] *= sf;
            rxyz-> z[k] *= sf;
        }

        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        for (i = 1; i < TEQ; i++) { // loop de equilibracion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);

            // AUTOVECTORIZADO
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                vxyz-> x[k] *= sf;
                vxyz-> y[k] *= sf;
                vxyz-> z[k] *= sf;
            }
        }

        int mes = 0;
        double epotm = 0.0;
        double presm = 0.0;

        // NO AUTOVECTORIZADO: control flow
        for (i = TEQ; i < TRUN; i++) { // loop de medicion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);

            // AUTOVECTORIZADO
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                vxyz-> x[k] *= sf;
                vxyz-> y[k] *= sf;
                vxyz-> z[k] *= sf;
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

        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (double)mes, presm / (double)mes);

    }

    double elapsed = wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);

    double gf_force = f_force/(9+TRUN);    // flop/s

    fprintf(file_metrica,"%f\n",gf_force);

    return 0;
}

