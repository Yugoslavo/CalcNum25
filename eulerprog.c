#include <stdlib.h>
#include <string.h> /*memcpy*/
#include <math.h>
/*Calcul matriciel*/
static void matvec_csr(int n, const int *ia, const int *ja, const double *a,
                const double *x, double *y) { /*y = Av*/
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int p = ia[i]; p < ia[i+1]; ++p) s += a[p] * x[ ja[p] ];
    y[i] = s;
  }

}
/*Euler Progressif*/
void euler_prog(int n, const int *ia, const int *ja, const double *a,
    double bheta2, int it, double dt,
    double *phi_init, double *phi_prog){

    int tau = 300;
    double c = tau * dt;

   double *y = (double*)malloc((size_t)n*sizeof(double));
   double *phi = (double*)malloc((size_t)n*sizeof(double));
    /* J'utilise memcpy pour copier la mémoire sans reasigner les pointeurs*/
   memcpy(phi, phi_init, (size_t)n*sizeof(double));
   
    for (int t = 0; t < it; ++t){
      matvec_csr(n,ia,ja,a,phi,y); /*A*phi*/
      for (int i = 0; i < n ; ++i){
         double r = y[i] - bheta2*phi[i];
         phi[i] -= c*r;
      }
      } /*y est egale a h*tau(Aphi - Bheta2phi)*/
    
    memcpy(phi_prog, phi, (size_t)n*sizeof(double));
    free(y);
    free(phi);
}
