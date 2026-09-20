/* Exercise 6 - Parallel element-wise multiplication using OpenMP and
 *              STRIP MINING.
 * Compile: gcc -O2 -fopenmp strip_parallel.c -o strip_parallel
 * Run:     OMP_NUM_THREADS=4 ./strip_parallel
 *
 * Strip mining: the single loop 0..N-1 is split into two loops.
 *   - The OUTER loop walks over the array one strip at a time
 *     (i = 0, STRIP, 2*STRIP, ...).  This is the loop OpenMP
 *     distributes across the threads.
 *   - The INNER loop processes the elements inside one strip.  Its
 *     trip count is a fixed, small, power of two, so the compiler can
 *     turn it into SIMD instructions without any remainder handling.
 *
 * Strip size:  AVX2 registers are 256 bits = 8 floats, AVX-512 is
 * 512 bits = 16 floats.  STRIP = 1024 floats (4 KB) is a multiple of
 * both 8 and 16, so every strip maps onto whole vector registers, and
 * it is small enough to stay in L1 cache.  N = 1,000,000 is not a
 * multiple of 1024 so the last strip is clamped to N.
 *
 * As in the serial version the kernel is repeated REPS times so that
 * the timing is measurable.
 */
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N      1000000
#define STRIP  1024        /* multiple of the SIMD width (8/16 floats) */
#define REPS   200

int main(void) {
  float *A = malloc(N * sizeof(float));
  float *B = malloc(N * sizeof(float));
  float *C = malloc(N * sizeof(float));
  double tstart, tstop, tcalc;

  for (int i = 0; i < N; i++) { A[i] = i * 0.5f; B[i] = 2.0f; }

  tstart = omp_get_wtime();

  for (int rep = 0; rep < REPS; rep++) {
    /* outer loop over strips is shared between the threads */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i += STRIP) {
      int end = (i + STRIP < N) ? i + STRIP : N;   /* clamp last strip */

      /* inner loop over one strip is vectorised */
      #pragma omp simd
      for (int j = i; j < end; j++)
        C[j] = A[j] * B[j];
    }
  }

  tstop = omp_get_wtime();
  tcalc = tstop - tstart;

  double sum = 0.0;
  for (int i = 0; i < N; i++) sum += C[i];

  printf("Parallel : C[10]=%.1f C[N-1]=%.1f checksum=%.1f  threads=%d strip=%d  time(%d reps) = %.6f s\n",
         C[10], C[N-1], sum, omp_get_max_threads(), STRIP, REPS, tcalc);

  free(A); free(B); free(C);
  return 0;
}
