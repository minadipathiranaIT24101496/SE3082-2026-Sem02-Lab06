/* Exercise 5 - Serial Fibonacci (reference version)
 * Compile: gcc -O2 -fopenmp fib_serial.c -o fib_serial
 * Run:     ./fib_serial 35
 */
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int fib(int n) {
  int i, j;
  if (n<2)
    return n;
  else {
    i=fib(n-1);
    j=fib(n-2);
    return i+j;
  }
}

int main(int argc, char **argv) {
  int n = (argc > 1) ? atoi(argv[1]) : 35;
  double tstart, tstop, tcalc;

  tstart = omp_get_wtime();
  int result = fib(n);
  tstop  = omp_get_wtime();
  tcalc  = tstop - tstart;      /* seconds */

  printf("Serial   fib(%d) = %d   time = %.4f s\n", n, result, tcalc);
  return 0;
}
