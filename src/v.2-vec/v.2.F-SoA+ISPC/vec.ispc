#include "parameters.h"

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

export uniform double dist(uniform const double rxyz[],uniform double fxyz[],
                uniform double* uniform epot, uniform double* uniform pres_vir,uniform const double L,
                uniform int i, uniform double rcut2){

    double flop = 0;

    double epot_v = 0.0;
    double pres_vir_v = 0.0;
    double sumx = 0.0;
    double sumy = 0.0;
    double sumz = 0.0;

    uniform double xi = rxyz[i];
    uniform double yi = rxyz[i + N];
    uniform double zi = rxyz[i + 2*N];

    foreach(j = i + 1 ... N) { // (N-i-1) iteraciones

        double xj = rxyz[j];
        double yj = rxyz[j + N];
        double zj = rxyz[j + 2*N];
        
        // distancia mínima entre r_i y r_j
        // (rest + flop m_i)*3 -> 3*3 = 9
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
            
            epot_v += 4.0 * r6inv * (r6inv - 1.0) - ECUT; //sum 2mult 2rest -> 5
            pres_vir_v += fr * rij2; // sum mult -> 2
            
                        
            sumx += fr * rx; // sum mult -> 2
            sumy += fr * ry; // sum mult -> 2
            sumz += fr * rz; // sum mult -> 2
            
            fxyz[j] -= fr * rx; // rest mult -> 2
            fxyz[j + N] -= fr * ry; // rest mult -> 2
            fxyz[j + 2*N] -= fr * rz; // rest mult -> 2

            flop += 27;
        }
    }
    
    *epot += reduce_add(epot_v);
    fxyz[i] += reduce_add(sumx);
    fxyz[i+N] += reduce_add(sumy);
    fxyz[i+2*N] += reduce_add(sumz);
    *pres_vir += reduce_add(pres_vir_v);

    return reduce_add(flop);
}

