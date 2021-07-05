#pragma once

void launch_forces_cuda_1(const double* rxyz, double* fxyz, const double L,
                            const int i, double* epot_array, double* pres_vir_array); // float elapsed_k, double flop)

void launch_forces_cuda_2(double* epot_array, double* pres_vir_array,
                            double* epot, double* pres_vir); //, double* flop)

uint div_ceil(uint a, uint b);
