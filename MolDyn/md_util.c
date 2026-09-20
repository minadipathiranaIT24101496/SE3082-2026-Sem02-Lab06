/* md_util.c - helper routines: lattice generation, initial velocities,
 * particle moves, kinetic energy, velocity average and output.
 * (fcc.c, mxwell.c, domove.c, mkekin.c, velavg.c, prnout.c, dfill.c and
 *  dscal.c of the original EPCC code merged into one file.)
 */
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "moldyn.h"

/* ---- portable drand48() replacement (same 48-bit LCG as POSIX) ----- */
static uint64_t rng_state;
static void   my_srand48(long seed) { rng_state = ((uint64_t)seed << 16) | 0x330E; }
static double my_drand48(void) {
  rng_state = (rng_state * 0x5DEECE66DULL + 0xBULL) & ((1ULL << 48) - 1);
  return (double)rng_state / (double)(1ULL << 48);
}

/* ---- generate fcc lattice for atoms inside the box ------------------ */
void fcc(double x[], int npart_, int mm_, double a) {
  int ijk = 0, i, j, k, lg;
  (void)npart_;
  for (lg = 0; lg < 2; lg++)
    for (i = 0; i < mm_; i++)
      for (j = 0; j < mm_; j++)
        for (k = 0; k < mm_; k++) {
          x[ijk]   = i*a + lg*a*0.5;
          x[ijk+1] = j*a + lg*a*0.5;
          x[ijk+2] = k*a;
          ijk += 3;
        }
  for (lg = 1; lg < 3; lg++)
    for (i = 0; i < mm_; i++)
      for (j = 0; j < mm_; j++)
        for (k = 0; k < mm_; k++) {
          x[ijk]   = i*a + (2-lg)*a*0.5;
          x[ijk+1] = j*a + (lg-1)*a*0.5;
          x[ijk+2] = k*a + a*0.5;
          ijk += 3;
        }
}

/* ---- sample Maxwell distribution at temperature tref ---------------- */
void mxwell(double vh[], int n3, double h, double tref) {
  int i, npart_ = n3/3;
  double r, tscale, v1, v2, s, ekin = 0.0, sp = 0.0, sc;

  my_srand48(4711);
  tscale = 16.0 / ((double)npart_ - 1.0);
  for (i = 0; i < n3; i += 2) {
    s = 2.0;
    while (s >= 1.0) {
      v1 = 2.0*my_drand48() - 1.0;
      v2 = 2.0*my_drand48() - 1.0;
      s  = v1*v1 + v2*v2;
    }
    r = sqrt(-2.0*log(s)/s);
    vh[i]   = v1*r;
    vh[i+1] = v2*r;
  }
  /* remove centre-of-mass motion in each direction, sum k.e. */
  for (int d = 0; d < 3; d++) {
    sp = 0.0;
    for (i = d; i < n3; i += 3) sp += vh[i];
    sp /= (double)npart_;
    for (i = d; i < n3; i += 3) { vh[i] -= sp; ekin += vh[i]*vh[i]; }
  }
  sc = h * sqrt(tref / (tscale*ekin));
  for (i = 0; i < n3; i++) vh[i] *= sc;
}

void dfill(int n, double val, double a[], int ia) {
  for (int i = 0; i < n*ia; i += ia) a[i] = val;
}

void dscal(int n, double sa, double sx[], int incx) {
  for (int i = 0; i < n*incx; i += incx) sx[i] *= sa;
}

/* ---- move particles and partially update velocities ----------------- */
void domove(int n3, double x[], double vh[], double f[], double side) {
  for (int i = 0; i < n3; i++) {
    x[i] += vh[i] + f[i];
    /* periodic boundary conditions */
    if (x[i] < 0.0)  x[i] += side;
    if (x[i] > side) x[i] -= side;
    /* partial velocity update */
    vh[i] += f[i];
    /* initialise forces for the next iteration */
    f[i] = 0.0;
  }
}

/* ---- scale forces, update velocities and compute k.e. --------------- */
double mkekin(int npart_, double f[], double vh[], double hsq2, double hsq) {
  double sum = 0.0;
  for (int i = 0; i < 3*npart_; i++) {
    f[i]  *= hsq2;
    vh[i] += f[i];
    sum   += vh[i]*vh[i];
  }
  return sum / hsq;
}

/* ---- compute average velocity ---------------------------------------- */
double velavg(int npart_, double vh[], double vaver, double h) {
  double vaverh = vaver*h, vel = 0.0, sq;
  count = 0.0;
  for (int i = 0; i < npart_*3; i += 3) {
    sq = sqrt(vh[i]*vh[i] + vh[i+1]*vh[i+1] + vh[i+2]*vh[i+2]);
    if (sq > vaverh) count++;
    vel += sq;
  }
  return vel / h;
}

/* ---- print out energies etc. ---------------------------------------- */
void prnout(int move, double ekin, double epot_, double tscale, double vir_,
            double vel, double count_, int npart_, double den) {
  double ek, etot, temp, pres, rp;
  ek    = 24.0*ekin;
  epot_ *= 4.0;
  etot  = ek + epot_;
  temp  = tscale*ekin;
  pres  = den*16.0*(ekin - vir_)/npart_;
  vel  /= npart_;
  rp    = (count_/(double)npart_)*100.0;
  printf(" %6d%12.4f%12.4f%12.4f%10.4f%10.4f%10.4f%6.1f\n",
         move, ek, epot_, etot, temp, pres, vel, rp);
}
