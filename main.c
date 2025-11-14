#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "time.h"
#include "interface_primme.h"
#include "csr_io.h"
#include "plotflux.h"
#include <math.h>
#include "eulerprog.h"



/*Fonction pour produit matrice vecteur Av*/
/*J'ai préferé ne pas utiliser matvec_primme pour  ne pas devoir initialiser les var*/
static void matvec_csr(int n, const int *ia, const int *ja, const double *a,
                const double *x, double *y) {
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int p = ia[i]; p < ia[i+1]; ++p) s += a[p] * x[ ja[p] ];
    y[i] = s;
  }
}

/* Norm et produit scalaire,  finalement pas utilisé*/
double dot(int n, const double *x, const double *y) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += x[i] * y[i];
  return s;
}

double norm2(int n, const double *x) {
  return sqrt(dot(n, x, x));
}
/* Calcul du residu*/

double residual_ratio(int n, const int *ia, const int *ja, const double *a, const double *v, double bheta){
  double *Av = (double*)malloc((size_t)n * sizeof(double));
  if (!Av) return -1.0; // OOM, error
  matvec_csr(n, ia, ja, a, v, Av);
  double num = 0.0, den = 0.0;
  for (int i = 0; i < n; ++i) {
    double r = Av[i] - bheta * v[i];
    num += r*r;
    den += v[i]*v[i];
  }
  free(Av);
  return sqrt(num) / sqrt(den);
}

int main(int argc, char *argv[])
{
  /* déclarer les variables */
  double beta2 = 4;
  double lam = 2.949094; /*Jai enlevé evals[0] pour linstant*/
  double alpha = sqrt(lam/beta2) ; /*sqrt(lam/beta2)*/

  int m = 90, nev = 1;
  int n, *ia, *ja; 
  double *a;
  double *evals, *evecs;
  double tc1, tc2, tw1, tw2;

  /* générer le problème */
  if (prob(m, alpha ,&n, &ia, &ja, &a))
     return 1;

   /* Copier la matrice CSR*/
  if  (!write_csr_arrays("Mat_CSR",n,ia,ja,a)) return 1;
  printf("\nPROBLÈME: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
  
  /* allouer la mémoire pour vecteurs & valeurs propres */
  evals = calloc(nev, sizeof(double));
  evecs = calloc(nev * n, sizeof(double));

  if (evals == NULL || evecs == NULL) {
      printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
      return 1;
  }

  /* primme - résolution */
  tc1 = mytimer_cpu(); tw1 = mytimer_wall();
  if(primme(n, ia, ja, a, nev, evals, evecs))
     return 1;
  tc2 = mytimer_cpu(); tw2 = mytimer_wall();

  /* temps de solution */
  printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1);
  printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1);
  printf("\nValeur propre minimale calculée: %5.1f\n",evals[0]);
  
  /* Affichage de la val propre et vecteur propre*/
  
  /* Affichage de la plus petite propre des valeurs propres */
  for (int k = 0; k < nev; ++k) { 
      printf("eval[%d] = %.15g\n", k, evals[k]);
  }
  /* Affichage des 24 premieres valeurs du vecteur propre */
  
  int head = (n <= 24) ? n : 24;
  printf("taille du vecteur propre : %d \n",n);
  for (int i = 0; i < head; ++i) {
      printf("evec0[%d] = %.6e\n", i, evecs[i]);  // evecs contient le vecteur 0 en tête
  }  

  /*Initialiser les vecteur et val prop*/
  const double *v0 = evecs; /*Vecteur propre*/
  double bheta0 = evals[0];
  double *y = (double*)malloc((size_t)n * sizeof(double));
  matvec_csr(n, ia, ja, a, v0, y); /*Output y = Av0*/
  double res = residual_ratio(n, ia, ja, a, v0, bheta0);
  printf("Résidu relatif  ||Av0 - bheta*v0|| / ||v0|| = %.3e\n", res);
  free(y);

  /*Calcul des dimentions pour cas stationnaire*/
  /*Alpha étant l'homothétie a realiser pour ajuster la valeur propre a 4*/

  printf("Les dimentions fois alpha = %.3e\n",alpha); 
  /*La valeur de h doit varier, en fixant m, h' = alpha*h */
  double h = 3.5/(m - 1);
  printf("Nouveau pad : h' = %.3e\n Ancien pad h = %3.e \n", h*alpha, h);

  /*Plot grace a gnuplot*/
  const double Lx = 3.5, Ly = 3.5;

  /*Normalisation pour optimiser l'affichage*/
  
  //double maxabs = 0.0;
  //for(int k = 0 ; k < n; ++k) if (fabs(evecs[k])>maxabs) maxabs = fabs(evecs[k]); /*fabs retourne la valeur absolue*/
  //if (maxabs > 0 ) for (int k = 0;k<n;++k) evecs[k] /= maxabs; /*Tout diviser par le max*/
  //printf("maxabs : %d \n",maxabs);

  double maxabs = 0.0;
  int pmax = 0;
  for (int p = 0; p < n; ++p) {
      double ap = fabs(evecs[p]);
      if (ap > maxabs) { maxabs = ap; pmax = p; }
  }
  if (maxabs > 0.0) {
    for (int p = 0; p < n; ++p) evecs[p] /= maxabs; // scale to [-1,1]
  }
  // flip so the peak is positive
  if (evecs[pmax] < 0.0) for (int p = 0; p < n; ++p) evecs[p] = -evecs[p];
  printf("maxabs : %.6e\n", maxabs);

  printf("taille du vecteur propre : %d \n",n);
  for (int i = 0; i < head; ++i) {
      printf("evec0[%d] = %.6e\n", i, evecs[i]);  // evecs contient le vecteur 0 en tête
  }

/*Approximation avec euler prog*/
double lam2 = 8.0/(alpha*h*alpha*h);
double dt = 2/(300*(lam2-beta2));
int it = 3000;
double *phi_init = (double*)malloc(n*sizeof(double));
for(int i = 0; i < n ; ++i)phi_init[i] = 1.0;w

double *phi_prog = (double*)malloc(n*sizeof(double));

euler_prog(n, ia, ja, a, beta2,  it,  dt,
     phi_init,  phi_prog);

double *Z = malloc(sizeof(double)*m*m);
/*Decompresser le vecteur phi en Z de taille mm*/
decompress(phi_prog, m, alpha, 3.5, Z);
/*Ecrire sur format TXT pour etre lut par gp*/  
write_phi("phiB.txt", Z, m, alpha,3.5,3.5);

plot_phi("phiB.txt", m, alpha,3.5,3.5, "Mode fondamental", 0);


  /* libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs); free(Z);
  return 0;
  
}