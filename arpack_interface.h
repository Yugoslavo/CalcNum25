#pragma once
#include <stddef.h>

int arpackSym(
    int n,
    const int *ia, const int *ja, const double *a,
    int nev, int ncv, double tol, int maxit, const char *which,
    double *evals, double *evecs, int *nconv_out);

double arpackResidu(
    int n, const int *ia, const int *ja, const double *a,
    const double *v, double lambda);

static inline int arpack_solve_smallest_sym_simple(
    int n, const int *ia, const int *ja, const double *a,
    int nev, double tol, int maxit,
    double *evals, double *evecs)
{
    return arpackSym(
        n, ia, ja, a,
        nev, /*ncv*/0, tol, maxit, "SA",
        evals, evecs, /*nconv_out*/NULL);
}
