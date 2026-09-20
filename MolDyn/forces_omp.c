/* forces_omp.c - OpenMP PARALLEL version
 *
 * The outer loop over particles i is distributed with
 * #pragma omp parallel for.
 *
 *  private   : j, xi,yi,zi, fxi,fyi,fzi, xx,yy,zz, rd, rrd..rrd7, r148,
 *              forcex,forcey,forcez  - scratch values for ONE pair, each
 *              thread needs its own copy.
 *  reduction : epot and vir (the 2 reduction variables) - every pair adds
 *              to them, so each thread accumulates a private partial sum
 *              that OpenMP adds together at the end of the loop.
 *  critical  : the force array f[].  Thread A working on row i writes
 *              f[j] for every j > i, and thread B may own row j (or write
 *              the same f[j] from its own row).  Two threads can therefore
 *              update the same f[] element at the same time, so the
 *              updates are wrapped in #pragma omp critical.
 *
 * SCHEDULE: the amount of work per i is NOT constant (row i only visits
 * j > i, so early rows do far more work than late rows).  A plain
 * static schedule gives thread 0 the heavy rows; schedule(static,n)
 * with a small chunk n interleaves rows between threads and balances the
 * load.  CHUNK can be set at compile time: -DCHUNK=n  (0 = plain static)
 */
#include <omp.h>
#include "moldyn.h"

#ifndef CHUNK
#define CHUNK 0
#endif

void forces(int npart_, double x[], double f[], double side, double rcoff) {
  int    i, j;
  double sideh, rcoffs;
  double xi, yi, zi, fxi, fyi, fzi, xx, yy, zz;
  double rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148;
  double forcex, forcey, forcez;

  vir    = 0.0;
  epot   = 0.0;
  sideh  = 0.5 * side;
  rcoffs = rcoff * rcoff;

#if CHUNK > 0
  #pragma omp parallel for schedule(static,CHUNK) default(none) \
          shared(npart_, x, f, side, sideh, rcoffs) \
          private(j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, \
                  rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, \
                  forcex, forcey, forcez) \
          reduction(+:epot, vir)
#else
  #pragma omp parallel for schedule(static) default(none) \
          shared(npart_, x, f, side, sideh, rcoffs) \
          private(j, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, \
                  rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, \
                  forcex, forcey, forcez) \
          reduction(+:epot, vir)
#endif
  for (i = 0; i < npart_*3; i += 3) {
    xi = x[i]; yi = x[i+1]; zi = x[i+2];
    fxi = 0.0; fyi = 0.0; fzi = 0.0;

    for (j = i+3; j < npart_*3; j += 3) {
      xx = xi - x[j];
      yy = yi - x[j+1];
      zz = zi - x[j+2];
      if (xx < -sideh) xx += side;
      if (xx >  sideh) xx -= side;
      if (yy < -sideh) yy += side;
      if (yy >  sideh) yy -= side;
      if (zz < -sideh) zz += side;
      if (zz >  sideh) zz -= side;
      rd = xx*xx + yy*yy + zz*zz;

      if (rd <= rcoffs) {
        rrd  = 1.0 / rd;
        rrd2 = rrd  * rrd;
        rrd3 = rrd2 * rrd;
        rrd4 = rrd2 * rrd2;
        rrd6 = rrd2 * rrd4;
        rrd7 = rrd6 * rrd;
        epot += (rrd6 - rrd3);           /* reduction variable 1 */
        r148  = rrd7 - 0.5 * rrd4;
        vir  -= rd * r148;               /* reduction variable 2 */
        forcex = xx * r148; fxi += forcex;
        forcey = yy * r148; fyi += forcey;
        forcez = zz * r148; fzi += forcez;

        /* f[j] is shared between threads -> must be updated atomically */
        #pragma omp critical
        {
          f[j]   -= forcex;
          f[j+1] -= forcey;
          f[j+2] -= forcez;
        }
      }
    }

    /* f[i] can also be written by another thread as its f[j] */
    #pragma omp critical
    {
      f[i] += fxi; f[i+1] += fyi; f[i+2] += fzi;
    }
  }
}
