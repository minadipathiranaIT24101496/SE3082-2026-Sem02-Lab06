/* Molecular Dynamics exercise - shared declarations */
#ifndef MOLDYN_H
#define MOLDYN_H

#define mm     15
#define npart  (4*mm*mm*mm)      /* 13500 argon atoms */

/* globals accumulated inside forces() */
extern double epot, vir, count;

void   fcc    (double x[], int npart_, int mm_, double a);
void   mxwell (double vh[], int n3, double h, double tref);
void   dfill  (int n, double val, double a[], int ia);
void   dscal  (int n, double sa, double sx[], int incx);
void   domove (int n3, double x[], double vh[], double f[], double side);
void   forces (int npart_, double x[], double f[], double side, double rcoff);
double mkekin (int npart_, double f[], double vh[], double hsq2, double hsq);
double velavg (int npart_, double vh[], double vaver, double h);
void   prnout (int move, double ekin, double epot_, double tscale, double vir_,
               double vel, double count_, int npart_, double den);

#endif
