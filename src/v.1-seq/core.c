#include "core.h"
#include "parameters.h"
#include "wtime.h"

#include <math.h>
#include <stdlib.h> // rand()

void init_pos(double* rxyz, const double rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    double a = cbrt(4.0 / rho); // cbrt div -> 2
    int nucells = ceil(cbrt((double)N / 4.0)); // cbrt div -> 2
    int idx = 0;

    for (int i = 0; i < nucells; i++) {         // nucells iteraciones
        for (int j = 0; j < nucells; j++) {     // nucells iteraciones
            for (int k = 0; k < nucells; k++) { // nucells iteraciones
                rxyz[idx + 0] = i * a; // x mult
                rxyz[idx + 1] = j * a; // y mult
                rxyz[idx + 2] = k * a; // z mult
                    // del mismo átomo
                rxyz[idx + 3] = (i + 0.5) * a; //sum mult -> 2
                rxyz[idx + 4] = (j + 0.5) * a; //sum mult -> 2
                rxyz[idx + 5] = k * a; // mult

                rxyz[idx + 6] = (i + 0.5) * a; // sum mult -> 2
                rxyz[idx + 7] = j * a; // mult
                rxyz[idx + 8] = (k + 0.5) * a; // sum mult -> 2

                rxyz[idx + 9] = i * a; // mult
                rxyz[idx + 10] = (j + 0.5) * a; // mult sum -> 2
                rxyz[idx + 11] = (k + 0.5) * a; // mult sum -> 2

                idx += 12; // sum
            }
        }
    }
}

// N*15 + 8 + N*6 = N*21 + 8 Flop
void init_vel(double* vxyz, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    for (int i = 0; i < 3 * N; i += 3) { // N iteraciones
        vxyz[i + 0] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2
        vxyz[i + 1] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2
        vxyz[i + 2] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2

        sumvx += vxyz[i + 0]; // sum
        sumvy += vxyz[i + 1]; // sum
        sumvz += vxyz[i + 2]; // sum
        sumv2 += vxyz[i + 0] * vxyz[i + 0] + vxyz[i + 1] * vxyz[i + 1] // 3mult 3sum -> 6
            + vxyz[i + 2] * vxyz[i + 2];
    } // N * 15

    sumvx /= (double)N; // div
    sumvy /= (double)N; // div
    sumvz /= (double)N; // div
    *temp = sumv2 / (3.0 * N); // div mult -> 2
    *ekin = 0.5 * sumv2; // mult
    sf = sqrt(T0 / *temp); // sqrt div -> 2

    // N iteraciones
    for (int i = 0; i < 3 * N; i += 3) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        vxyz[i + 0] = (vxyz[i + 0] - sumvx) * sf; // rest mult -> 2
        vxyz[i + 1] = (vxyz[i + 1] - sumvy) * sf; // rest mult -> 2
        vxyz[i + 2] = (vxyz[i + 2] - sumvz) * sf; // rest mult -> 2
    } // -> N*6
}

// 2 Flop
static double minimum_image(double cordi, const double cell_length)
{
    // imagen más cercana

    if (cordi <= -0.5 * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5 * cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}

double f_force = 0;

// N(N-1)/2 * (41) + 5 flop
void forces(const double* rxyz, double* fxyz, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)
{
    // calcula las fuerzas LJ (12-6)
    double start = wtime();
    double flop = 0;
    for (int i = 0; i < 3 * N; i++) {
        fxyz[i] = 0.0;
    }
    double pres_vir = 0.0;
    double rcut2 = RCUT * RCUT; // mult
    *epot = 0.0;

    for (int i = 0; i < 3 * (N - 1); i += 3) { //(N-1) iteraciones

        double xi = rxyz[i + 0];
        double yi = rxyz[i + 1];
        double zi = rxyz[i + 2];

        for (int j = i + 3; j < 3 * N; j += 3) { // (N-i/3-1) iteraciones

            double xj = rxyz[j + 0];
            double yj = rxyz[j + 1];
            double zj = rxyz[j + 2];

            // distancia mínima entre r_i y r_j
            // (rest flop m_i)*3 -> 3*3 = 9
            double rx = xi - xj;
            rx = minimum_image(rx, L);
            double ry = yi - yj;
            ry = minimum_image(ry, L);
            double rz = zi - zj;
            rz = minimum_image(rz, L);

            double rij2 = rx * rx + ry * ry + rz * rz; // mult sum mult sum mult -> 5

            flop += 14;

            if (rij2 <= rcut2) {
                double r2inv = 1.0 / rij2; // div
                double r6inv = r2inv * r2inv * r2inv; // mult mult -> 2

                double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0); // 4mult rest -> 5

                fxyz[i + 0] += fr * rx; // sum mult -> 2
                fxyz[i + 1] += fr * ry; // sum mult -> 2
                fxyz[i + 2] += fr * rz; // sum mult -> 2

                fxyz[j + 0] -= fr * rx; // rest mult -> 2
                fxyz[j + 1] -= fr * ry; // rest mult -> 2
                fxyz[j + 2] -= fr * rz; // rest mult -> 2

                *epot += 4.0 * r6inv * (r6inv - 1.0) - ECUT; //sum 2mult 2rest -> 5
                pres_vir += fr * rij2; // sum mult -> 2

                flop += 27;
            }
        }
    }
    pres_vir /= (V * 3.0); // div mult -> 2
    *pres = *temp * rho + pres_vir; // mult sum -> 2

    f_force += ((flop + 4) / (wtime()-start))*0.000000001;
}


static double pbc(double cordi, const double cell_length)
{
    // condiciones periodicas de contorno coordenadas entre [0,L)
    if (cordi <= 0) {
        cordi += cell_length;
    } else if (cordi > cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}


void velocity_verlet(double* rxyz, double* vxyz, double* fxyz, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L)
{

    // N iteraciones
    for (int i = 0; i < 3 * N; i += 3) { // actualizo posiciones
        rxyz[i + 0] += vxyz[i + 0] * DT + 0.5 * fxyz[i + 0] * DT * DT; // 2sum 4mult -> 6
        rxyz[i + 1] += vxyz[i + 1] * DT + 0.5 * fxyz[i + 1] * DT * DT; // 2sum 4mult -> 6
        rxyz[i + 2] += vxyz[i + 2] * DT + 0.5 * fxyz[i + 2] * DT * DT; // 2sum 4mult -> 6

        rxyz[i + 0] = pbc(rxyz[i + 0], L); // sum o rest -> 1
        rxyz[i + 1] = pbc(rxyz[i + 1], L); // sum o rest -> 1
        rxyz[i + 2] = pbc(rxyz[i + 2], L); // sum o rest -> 1

        vxyz[i + 0] += 0.5 * fxyz[i + 0] * DT; // sum mult -> 2
        vxyz[i + 1] += 0.5 * fxyz[i + 1] * DT; // sum mult -> 2
        vxyz[i + 2] += 0.5 * fxyz[i + 2] * DT; // sum mult -> 2
    } // N*27 flop

    // N(N-1)/2 * (41) + 5 flop
    forces(rxyz, fxyz, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < 3 * N; i += 3) { // actualizo velocidades
        vxyz[i + 0] += 0.5 * fxyz[i + 0] * DT; // sum 2mult -> 3
        vxyz[i + 1] += 0.5 * fxyz[i + 1] * DT; // sum 2mult -> 3
        vxyz[i + 2] += 0.5 * fxyz[i + 2] * DT; // sum 2mult -> 3

        sumv2 += vxyz[i + 0] * vxyz[i + 0] + vxyz[i + 1] * vxyz[i + 1]
            + vxyz[i + 2] * vxyz[i + 2]; // 3sum 3mult -> 6
    } // N*15

    *ekin = 0.5 * sumv2; // mult
    *temp = sumv2 / (3.0 * N); // div mult -> 2
} // N*27 + N*15 + N(N-1)/2 * (41) + 5 + 3
