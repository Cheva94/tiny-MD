#include "core.h"
#include "parameters.h"
#include "wtime.h"
#include "vec.h"

#include <math.h>
#include <stdlib.h> // rand()
#include <stdio.h>

void init_pos(double* rxyz, const double rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    double a = cbrt(4.0 / rho); // cbrt div -> 2
    int nucells = ceil(cbrt((double)N / 4.0)); // cbrt div -> 2
    int idx = 0;

    // Loop NO autovectorizado: control de flujo (if) y los costos pueden no valer la pena
    for (int i = 0; i < nucells; i++) {         // nucells iteraciones
        for (int j = 0; j < nucells; j++) {     // nucells iteraciones
        // optimized: basic block part vectorized using 16 byte vectors
            for (int k = 0; k < nucells; k++) { // nucells iteraciones
                rxyz[idx + 0] = i * a; // x mult
                rxyz[idx + 1] = (i + 0.5) * a; //sum mult -> 2
                rxyz[idx + 2] = (i + 0.5) * a; // sum mult -> 2
                rxyz[idx + 3] = i * a; // mult
                
                rxyz[idx + 0 + N] = j * a; // y mult
                rxyz[idx + 1 + N] = (j + 0.5) * a; //sum mult -> 2
                rxyz[idx + 2 + N] = j * a; // mult
                rxyz[idx + 3 + N] = (j + 0.5) * a; // mult sum -> 2
                
                rxyz[idx + 0 + 2*N] = k * a; // z mult
                rxyz[idx + 1 + 2*N] = k * a; // mult
                rxyz[idx + 2 + 2*N] = (k + 0.5) * a; // sum mult -> 2
                rxyz[idx + 3 + 2*N] = (k + 0.5) * a; // mult sum -> 2

                idx += 4; // sum
            }
        }
    }
}

// N*15 + 8 + N*6 = N*21 + 8 Flop
void init_vel(double* vxyz, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    // Loop NO autovectorizado: rand statement clobbers memory.
    for (int i = 0; i < N; i++) { // N iteraciones
        vxyz[i] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2
        vxyz[i + N] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2
        vxyz[i + 2*N] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2

        sumvx += vxyz[i]; // sum
        sumvy += vxyz[i + N]; // sum
        sumvz += vxyz[i + 2*N]; // sum
        sumv2 += vxyz[i] * vxyz[i] + vxyz[i + N] * vxyz[i + N] // 3mult 3sum -> 6
            + vxyz[i + 2*N] * vxyz[i + 2*N];
    } // N * 15

    sumvx /= (double)N; // div
    sumvy /= (double)N; // div
    sumvz /= (double)N; // div
    *temp = sumv2 / (3.0 * N); // div mult -> 2
    *ekin = 0.5 * sumv2; // mult
    sf = sqrt(T0 / *temp); // sqrt div -> 2

    // N iteraciones
    // Loop autovectorizado usando vectores de 32 bytes
    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        vxyz[i] = (vxyz[i] - sumvx) * sf; // rest mult -> 2
        vxyz[i + N] = (vxyz[i + N] - sumvy) * sf; // rest mult -> 2
        vxyz[i + 2*N] = (vxyz[i + 2*N] - sumvz) * sf; // rest mult -> 2
    } // -> N*6
}

// 2 Flop

double f_force = 0;

// NO autovectorizado: memset statement clobbers memory.
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

    // Loop NO autovectorizado: contiene loops anidados.
    for (int i = 0; i < (N - 1); i++) { //(N-1) iteraciones

        // Loop NO autovectorizado:
        flop += dist(rxyz,fxyz,epot,&pres_vir,L,i,rcut2);

    }
    
    //printf("PRESION = %f\n",pres_vir);
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
    // Loop autovectorizado usando vectores de 32 bytes (reversionado por posible aliasing)
    for (int i = 0; i < N; i++) { // actualizo posiciones
        rxyz[i] += vxyz[i] * DT + 0.5 * fxyz[i] * DT * DT; // 2sum 4mult -> 6
        rxyz[i + N] += vxyz[i + N] * DT + 0.5 * fxyz[i + N] * DT * DT; // 2sum 4mult -> 6
        rxyz[i + 2*N] += vxyz[i + 2*N] * DT + 0.5 * fxyz[i + 2*N] * DT * DT; // 2sum 4mult -> 6

        rxyz[i] = pbc(rxyz[i], L); // sum o rest -> 1
        rxyz[i + N] = pbc(rxyz[i + N], L); // sum o rest -> 1
        rxyz[i + 2*N] = pbc(rxyz[i + 2*N], L); // sum o rest -> 1

        vxyz[i] += 0.5 * fxyz[i] * DT; // sum mult -> 2
        vxyz[i + N] += 0.5 * fxyz[i + N] * DT; // sum mult -> 2
        vxyz[i + 2*N] += 0.5 * fxyz[i + 2*N] * DT; // sum mult -> 2
    } // N*27 flop

    // NO autovectorizado: forces statement clobbers memory.
    // N(N-1)/2 * (41) + 5 flop
    forces(rxyz, fxyz, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;

    // Loop autovectorizado usando vectores de 32 bytes (reversionado por posible aliasing)
    for (int i = 0; i < N; i++) { // actualizo velocidades
        vxyz[i] += 0.5 * fxyz[i] * DT; // sum 2mult -> 3
        vxyz[i + N] += 0.5 * fxyz[i + N] * DT; // sum 2mult -> 3
        vxyz[i + 2*N] += 0.5 * fxyz[i + 2*N] * DT; // sum 2mult -> 3

        sumv2 += vxyz[i] * vxyz[i] + vxyz[i + N] * vxyz[i + N]
            + vxyz[i + 2*N] * vxyz[i + 2*N]; // 3sum 3mult -> 6
    } // N*15

    *ekin = 0.5 * sumv2; // mult
    *temp = sumv2 / (3.0 * N); // div mult -> 2
} // N*27 + N*15 + N(N-1)/2 * (41) + 5 + 3

