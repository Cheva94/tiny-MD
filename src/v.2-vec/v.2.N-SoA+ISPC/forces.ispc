#include "parameters.h"
#include "core.h"

double flop = 0;
double rcut2 = RCUT * RCUT;
double pres_vir = 0.0;

export void to_zero(uniform struct SoA xyz[])
{
  foreach (i = 0 ... N){
    xyz-> x[i] = 0.0;
    xyz-> y[i] = 0.0;
    xyz-> z[i] = 0.0;
  }
}

double minimum_image(double cordi, const double cell_length)
{
  if (cordi <= -0.5 * cell_length) {
    cordi += cell_length;
  } else if (cordi > 0.5 * cell_length) {
    cordi -= cell_length;
  }
  return cordi;
}

export void forces_contrib(const uniform struct SoA rxyz[],
                            uniform struct SoA fxyz[],
                            uniform struct SoA dxyz[],
                            uniform double epot[],
                            const uniform double L,
                            uniform int i)
{
  double parc_epot = 0.0, parc_pres_vir = 0.0;
  // to_zero(sxyz);
  double sumx = 0, sumy = 0, sumz = 0;

  // foreach (i = 0 ... N-1, j = 1 ... N) {
  foreach (j = 1 ... N) {
    // cif (i<j) {
    if (i<j) {
      dxyz-> x[j] = rxyz-> x[i] - rxyz-> x[j];
      dxyz-> x[j] = minimum_image(dxyz-> x[j], L);

      dxyz-> y[j] = rxyz-> y[i] - rxyz-> y[j];
      dxyz-> y[j] = minimum_image(dxyz-> y[j], L);

      dxyz-> z[j] = rxyz-> z[i] - rxyz-> z[j];
      dxyz-> z[j] = minimum_image(dxyz-> z[j], L);

      double rij2 = dxyz-> x[j] * dxyz-> x[j] + dxyz-> y[j] * dxyz-> y[j] + dxyz-> z[j] * dxyz-> z[j]; // mult sum mult sum mult -> 5

      flop += 14;

      // cif (rij2 <= rcut2) {
      if (rij2 <= rcut2) {
        double r2inv = 1.0 / rij2; // div
        double r6inv = r2inv * r2inv * r2inv; // mult mult -> 2
        double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0); // 4mult rest -> 5

        parc_epot += 4.0 * r6inv * (r6inv - 1.0) - ECUT; //sum 2mult 2rest -> 5
        parc_pres_vir += fr * rij2; // sum mult -> 2

        sumx += fr * dxyz-> x[j];
        sumy += fr * dxyz-> y[j];
        sumz += fr * dxyz-> z[j];

        fxyz-> x[j] -= fr * dxyz-> x[j]; // rest mult -> 2
        fxyz-> y[j] -= fr * dxyz-> y[j]; // rest mult -> 2
        fxyz-> z[j] -= fr * dxyz-> z[j]; // rest mult -> 2

        flop += 27;
      } else {
      continue;
      }
    }
  }

  fxyz-> x[i] += reduce_add(sumx); // sum mult -> 2
  fxyz-> y[i] += reduce_add(sumy); // sum mult -> 2
  fxyz-> z[i] += reduce_add(sumz); // sum mult -> 2
  epot[programIndex] += reduce_add(parc_epot); //sum 2mult 2rest -> 5
  pres_vir += reduce_add(parc_pres_vir); // sum mult -> 2
  // return 0 ;
}
