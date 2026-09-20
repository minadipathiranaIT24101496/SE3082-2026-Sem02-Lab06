/* Exercise 6 - Serial element-wise multiplication  C[i] = A[i] * B[i]
 * Compile: gcc -O2 -fopenmp strip_serial.c -o strip_serial
 *
 * A single pass over 1,000,000 floats takes well under a millisecond,
 * which is below the resolution of omp_get_wtime() on some systems, so
 * the kernel is repeated REPS times and the total time is reported.
 */
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N     1000000
#define REPS  200

int main(void) {
  float *A = malloc(N * sizeof(float));
  float *B = malloc(N * sizeof(float));
  float *C = malloc(N * sizeof(float));
  double tstart, tstop, tcalc;

  for (int i = 0; i < N; i++) { A[i] = i * 0.5f; B[i] = 2.0f; }

  tstart = omp_get_wtime();
  for (int rep = 0; rep < REPS; rep++)
    for (int i = 0; i < N; i++)
      C[i] = A[i] * B[i];
  tstop = omp_get_wtime();
  tcalc = tstop - tstart;

  /* checksum so the compiler cannot optimise the loop away */
  double sum = 0.0;
  for (int i = 0; i < N; i++) sum += C[i];

  printf("Serial   : C[10]=%.1f C[N-1]=%.1f checksum=%.1f  time(%d reps) = %.6f s\n",
         C[10], C[N-1], sum, REPS, tcalc);

  free(A); free(B); free(C);
  return 0;
}
