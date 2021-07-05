#include "parameters.h"
#include "wtime.h"
#include "forces.h"
#include <math.h>

#include <stdlib.h>

void init_pos(struct SoA * rxyz,
              const double rho)
{
  double a = cbrt(4.0 / rho); // cbrt div -> 2
  int idx = 0, nucells = ceil(cbrt((double)N / 4.0)); // cbrt div -> 2

  for (int i = 0; i < nucells; i++) {
    for (int j = 0; j < nucells; j++) {
      for (int k = 0; k < nucells; k++) {
        rxyz-> x[idx] = a * i;
        rxyz-> x[idx + 1] = a * (i + 0.5);
        rxyz-> x[idx + 2] = a * (i + 0.5);
        rxyz-> x[idx + 3] = a * i;

        rxyz-> y[idx] = a * j;
        rxyz-> y[idx + 1] = a * (j + 0.5);
        rxyz-> y[idx + 2] = a * j;
        rxyz-> y[idx + 3] = a * (j + 0.5);

        rxyz-> z[idx] = a * k;
        rxyz-> z[idx + 1] = a * k;
        rxyz-> z[idx + 2] = a * (k + 0.5);
        rxyz-> z[idx + 3] = a * (k + 0.5);

        idx += 4;
      }
    }
  }
}

// N*15 + 8 + N*6 = N*21 + 8 Flop
void init_vel(struct SoA * vxyz,
              double* temp,
              double* ekin)
{
  double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

  for (int i = 0; i < N; i++) { // N iteraciones
    vxyz-> x[i] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2
    vxyz-> y[i] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2
    vxyz-> z[i] = rand() / (double)RAND_MAX - 0.5; // rest div -> 2

    sumvx += vxyz-> x[i]; // sum
    sumvy += vxyz-> y[i]; // sum
    sumvz += vxyz-> z[i]; // sum
    sumv2 += vxyz-> x[i] * vxyz-> x[i] + vxyz-> y[i] * vxyz-> y[i] // 3mult 3sum -> 6
      + vxyz-> z[i] * vxyz-> z[i];
  } // N * 15

  sumvx /= (double)N; // div
  sumvy /= (double)N; // div
  sumvz /= (double)N; // div
  *temp = sumv2 / (3.0 * N); // div mult -> 2
  *ekin = 0.5 * sumv2; // mult
  sf = sqrt(T0 / *temp); // sqrt div -> 2

  // N iteraciones
  for (int i = 0; i < N; i++) {
    vxyz-> x[i] = (vxyz-> x[i] - sumvx) * sf; // rest mult -> 2
    vxyz-> y[i] = (vxyz-> y[i] - sumvy) * sf; // rest mult -> 2
    vxyz-> z[i] = (vxyz-> z[i] - sumvz) * sf; // rest mult -> 2
  } // -> N*6
}

double f_force = 0;

// N(N-1)/2 * (41) + 5 flop
void forces(const struct SoA * rxyz,
            struct SoA * fxyz,
            struct SoA * dxyz,
            double* epot,
            double* pres,
            const double* temp,
            const double rho,
            const double V,
            const double L)
{
  double start = wtime();
  double flop = 0;

  to_zero(fxyz);

  double pres_vir = 0.0;
  double rcut2 = RCUT * RCUT; // mult
  *epot = 0.0;

  for (int i = 0; i < N-1; i++) { //(N-1) iteraciones
    forces_contrib(rxyz, fxyz, dxyz, epot, L, i);
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

void velocity_verlet(struct SoA * rxyz,
                      struct SoA * vxyz,
                      struct SoA * fxyz,
                      struct SoA * dxyz,
                      double* epot,
                      double* ekin,
                      double* pres,
                      double* temp,
                      const double rho,
                      const double V,
                      const double L)
{

  // N iteraciones
  for (int i = 0; i < N; i++) { // actualizo posiciones
    rxyz-> x[i] += vxyz-> x[i] * DT + 0.5 * fxyz-> x[i] * DT * DT; // 2sum 4mult -> 6
    rxyz-> y[i] += vxyz-> y[i] * DT + 0.5 * fxyz-> y[i] * DT * DT; // 2sum 4mult -> 6
    rxyz-> z[i] += vxyz-> z[i] * DT + 0.5 * fxyz-> z[i] * DT * DT; // 2sum 4mult -> 6

    rxyz-> x[i] = pbc(rxyz-> x[i], L); // sum o rest -> 1
    rxyz-> y[i] = pbc(rxyz-> y[i], L); // sum o rest -> 1
    rxyz-> z[i] = pbc(rxyz-> z[i], L); // sum o rest -> 1

    vxyz-> x[i] += 0.5 * fxyz-> x[i] * DT; // sum mult -> 2
    vxyz-> y[i] += 0.5 * fxyz-> y[i] * DT; // sum mult -> 2
    vxyz-> z[i] += 0.5 * fxyz-> z[i] * DT; // sum mult -> 2
  } // N*27 flop

  // N(N-1)/2 * (41) + 5 flop
  forces(rxyz, fxyz, dxyz, epot, pres, temp, rho, V, L); // actualizo fuerzas

  double sumv2 = 0.0;

  for (int i = 0; i < N; i++) { // actualizo velocidades
    vxyz-> x[i] += 0.5 * fxyz-> x[i] * DT; // sum 2mult -> 3
    vxyz-> y[i] += 0.5 * fxyz-> y[i] * DT; // sum 2mult -> 3
    vxyz-> z[i] += 0.5 * fxyz-> z[i] * DT; // sum 2mult -> 3

    sumv2 += vxyz-> x[i] * vxyz-> x[i] + vxyz-> y[i] * vxyz-> y[i]
      + vxyz-> z[i] * vxyz-> z[i]; // 3sum 3mult -> 6
  } // N*15

  *ekin = 0.5 * sumv2; // mult
  *temp = sumv2 / (3.0 * N); // div mult -> 2
} // N*27 + N*15 + N(N-1)/2 * (41) + 5 + 3
