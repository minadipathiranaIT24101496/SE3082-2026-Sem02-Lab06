/* Exercise 5 - Parallel Fibonacci using OpenMP task parallelism
 * Compile: gcc -O2 -fopenmp fib_parallel.c -o fib_parallel
 * Run:     OMP_NUM_THREADS=4 ./fib_parallel 35
 *
 * Each recursive call fib(n-1) and fib(n-2) is spawned as an
 * independent OpenMP task. `shared(i,j)` lets the child tasks write
 * their result back to the parent's variables, and `taskwait` makes
 * the parent block until both children are finished before it adds
 * the two results together.
 *
 * Creating a task for every single call would drown the program in
 * scheduling overhead (fib(1), fib(0) ... are trivial), so below a
 * cut-off value the code falls back to the plain serial recursion.
 */
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define CUTOFF 20   /* below this n we stop creating tasks */

int fib_serial(int n) {
  if (n<2) return n;
  return fib_serial(n-1) + fib_serial(n-2);
}

int fib(int n) {
  int i, j;
  if (n<2)
    return n;
  if (n < CUTOFF)               /* granularity control */
    return fib_serial(n);
  else {
    #pragma omp task shared(i)
    i=fib(n-1);

    #pragma omp task shared(j)
    j=fib(n-2);

    #pragma omp taskwait          /* wait for both child tasks */
    return i+j;
  }
}

int main(int argc, char **argv) {
  int n = (argc > 1) ? atoi(argv[1]) : 35;
  int result;
  double tstart, tstop, tcalc;

  tstart = omp_get_wtime();

  #pragma omp parallel            /* create the team of threads   */
  {
    #pragma omp single            /* only ONE thread starts the   */
    result = fib(n);              /* recursion; the tasks it      */
  }                               /* creates are picked up by all */

  tstop = omp_get_wtime();
  tcalc = tstop - tstart;         /* seconds */

  printf("Parallel fib(%d) = %d   threads = %d   time = %.4f s\n",
         n, result, omp_get_max_threads(), tcalc);
  return 0;
}
