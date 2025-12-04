#ifndef JADAMILU_INTERFACE_H
#define JADAMILU_INTERFACE_H

int jadamilu_solve(
    int n,
    const int *ia,
    const int *ja,
    const double *a,
    int nev,
    double *evals,
    double *evecs
);

#endif
