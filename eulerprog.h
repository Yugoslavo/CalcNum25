#ifndef EULERPROG_H
#define EULERPROG_H

void matvec_csr(int n, const int *ia, const int *ja, const double *a,
                const double *x, double *y);

void euler_prog(int n, const int *ia, const int *ja, const double *a,
                double bheta2, int it, double h,
                double *phi_init, double *phi_prog);
                
#endif
