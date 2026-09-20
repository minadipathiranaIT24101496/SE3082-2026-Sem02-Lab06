/* forces.c - SERIAL version
 * Compute forces on all particles and accumulate the virial (vir) and
 * the potential energy (epot).  Newton 3rd law is used: each pair (i,j)
 * is visited once, the force is added to i and subtracted from j.
 */
#include "moldyn.h"

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

  for (i = 0; i < npart_*3; i += 3) {
    xi = x[i]; yi = x[i+1]; zi = x[i+2];
    fxi = 0.0; fyi = 0.0; fzi = 0.0;

    for (j = i+3; j < npart_*3; j += 3) {
      xx = xi - x[j];
      yy = yi - x[j+1];
      zz = zi - x[j+2];
      /* minimum image convention (periodic boundaries) */
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
        epot += (rrd6 - rrd3);
        r148  = rrd7 - 0.5 * rrd4;
        vir  -= rd * r148;
        forcex = xx * r148; fxi += forcex; f[j]   -= forcex;
        forcey = yy * r148; fyi += forcey; f[j+1] -= forcey;
        forcez = zz * r148; fzi += forcez; f[j+2] -= forcez;
      }
    }
    f[i] += fxi; f[i+1] += fyi; f[i+2] += fzi;
  }
}
