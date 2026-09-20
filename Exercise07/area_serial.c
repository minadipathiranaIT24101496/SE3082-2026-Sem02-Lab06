/* Exercise 7 - Area of the Mandelbrot Set (SERIAL version)
 * Compile: gcc -O2 -fopenmp area_serial.c -o area_serial
 *
 * A grid of NPOINTS x NPOINTS points c is generated in the box that
 * contains the upper half of the Mandelbrot set.  Each point is
 * iterated z = z*z + c up to MAXITER times; if |z| > 2 the point is
 * outside the set.  The area is estimated from the fraction of points
 * that never escape.
 */
#include <omp.h>
#include <stdio.h>

#define NPOINTS 2000
#define MAXITER 2000

struct complex { double real; double imag; };

int main(void) {
  int i, j, iter, numoutside = 0;
  double area, error, ztemp;
  double tstart, tstop, tcalc;
  struct complex z, c;

  tstart = omp_get_wtime();

  /* Outer loops iterate over the grid points, inner loop is the
   * Mandelbrot iteration for one point */
  for (i = 0; i < NPOINTS; i++) {
    for (j = 0; j < NPOINTS; j++) {
      c.real = -2.0 + 2.5 * (double)(i) / (double)(NPOINTS) + 1.0e-7;
      c.imag =        1.125 * (double)(j) / (double)(NPOINTS) + 1.0e-7;
      z = c;
      for (iter = 0; iter < MAXITER; iter++) {
        ztemp  = (z.real * z.real) - (z.imag * z.imag) + c.real;
        z.imag = z.real * z.imag * 2 + c.imag;
        z.real = ztemp;
        if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
          numoutside++;
          break;
        }
      }
    }
  }

  tstop = omp_get_wtime();
  tcalc = tstop - tstart;

  /* Box is 2.5 x 1.125, doubled because only the upper half was done */
  area  = 2.0 * 2.5 * 1.125 * (double)(NPOINTS * NPOINTS - numoutside)
        / (double)(NPOINTS * NPOINTS);
  error = area / (double)NPOINTS;

  printf("Serial   : Area of Mandelbrot set = %12.8f +/- %12.8f  time = %.3f s\n",
         area, error, tcalc);
  printf("Correct answer should be around 1.510659\n");
  return 0;
}
