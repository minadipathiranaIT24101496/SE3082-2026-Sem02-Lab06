/* Exercise 7 - Area of the Mandelbrot Set (PARALLEL version)
 * Compile: gcc -O2 -fopenmp area_parallel.c -o area_parallel
 * Run:     OMP_NUM_THREADS=4 ./area_parallel [mode]
 *
 *   mode 0 (default) : BLOCK mapping   - done "by hand" with the OpenMP
 *                      library routines omp_get_num_threads() and
 *                      omp_get_thread_num().  Thread t gets rows
 *                      [t*NPOINTS/nthreads , (t+1)*NPOINTS/nthreads).
 *                      Every thread gets an EQUAL number of points as
 *                      the exercise asks.
 *   mode 1           : CYCLIC mapping  - also by hand: thread t gets
 *                      rows t, t+nthreads, t+2*nthreads, ...
 *   mode 2           : omp for schedule(dynamic,10) - the runtime hands
 *                      out chunks of 10 rows on demand.
 *
 * Modes 1 and 2 are the "extra exercise" (different ways of mapping
 * iterations to threads).  They matter because the work per row is
 * very uneven: points INSIDE the set always run the full MAXITER
 * iterations while points far outside escape after a handful.
 *
 * Data sharing:
 *   shared    : nothing is written except numoutside (reduction)
 *   private   : i, j, iter, ztemp, z, c  (declared inside the region,
 *               so they are automatically private to each thread)
 *   reduction : numoutside (+)  - each thread counts its own escapes,
 *               OpenMP sums them when the region ends
 */
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define NPOINTS 2000
#define MAXITER 2000

struct complex { double real; double imag; };

/* Mandelbrot test for one row of the grid; returns number of points
 * in that row that are OUTSIDE the set */
static int do_row(int i) {
  int j, iter, out = 0;
  double ztemp;
  struct complex z, c;

  for (j = 0; j < NPOINTS; j++) {
    c.real = -2.0 + 2.5 * (double)(i) / (double)(NPOINTS) + 1.0e-7;
    c.imag =        1.125 * (double)(j) / (double)(NPOINTS) + 1.0e-7;
    z = c;
    for (iter = 0; iter < MAXITER; iter++) {
      ztemp  = (z.real * z.real) - (z.imag * z.imag) + c.real;
      z.imag = z.real * z.imag * 2 + c.imag;
      z.real = ztemp;
      if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
        out++;
        break;
      }
    }
  }
  return out;
}

int main(int argc, char **argv) {
  int mode = (argc > 1) ? atoi(argv[1]) : 0;
  int numoutside = 0;
  double area, error;
  double tstart, tstop, tcalc;
  const char *name[] = { "block (manual)", "cyclic (manual)", "dynamic,10" };

  tstart = omp_get_wtime();

  /* 1. Parallel region started BEFORE the main loop */
  #pragma omp parallel default(none) shared(mode) reduction(+:numoutside)
  {
    int nthreads = omp_get_num_threads();   /* OpenMP library routines */
    int tid      = omp_get_thread_num();
    int i;

    if (mode == 0) {
      /* 2. BLOCK distribution: equal contiguous share of rows per thread */
      int lo = (tid     * NPOINTS) / nthreads;
      int hi = ((tid+1) * NPOINTS) / nthreads;
      for (i = lo; i < hi; i++)
        numoutside += do_row(i);
    }
    else if (mode == 1) {
      /* CYCLIC distribution: rows dealt out like a pack of cards */
      for (i = tid; i < NPOINTS; i += nthreads)
        numoutside += do_row(i);
    }
    else {
      /* Let the runtime balance the load */
      #pragma omp for schedule(dynamic,10)
      for (i = 0; i < NPOINTS; i++)
        numoutside += do_row(i);
    }
  }

  tstop = omp_get_wtime();
  tcalc = tstop - tstart;

  area  = 2.0 * 2.5 * 1.125 * (double)(NPOINTS * NPOINTS - numoutside)
        / (double)(NPOINTS * NPOINTS);
  error = area / (double)NPOINTS;

  printf("Parallel : threads=%d  mapping=%-16s Area = %12.8f +/- %12.8f  time = %.3f s\n",
         omp_get_max_threads(), name[mode], area, error, tcalc);
  return 0;
}
