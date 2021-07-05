#pragma once

struct SoA {
  double x[N];
  double y[N];
  double z[N];
};

void init_pos(struct SoA * rxyz, const double rho);
void init_vel(struct SoA * vxyz, double* temp, double* ekin);
void forces(const struct SoA * rxyz, struct SoA * fxyz, double* epot,
  						double* pres, const double* temp, const double rho,
              const double V, const double L);
void velocity_verlet(struct SoA * rxyz, struct SoA * vxyz, struct SoA * fxyz,
  										double* epot, double* ekin, double* pres, double* temp,
  										const double rho, const double V, const double L);
double f_force;
