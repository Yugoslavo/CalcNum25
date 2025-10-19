#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "time.h"
#include "interface_primme.h"
#include "csr_io.h"
#include <math.h>


/*Fonction pour produit matrice vecteur Av*/
/*J'ai préferé ne pas utiliser matvec_primme pour pour ne pas devoir initialiser les var*/
void matvec_csr(int n, const int *ia, const int *ja, const double *a,
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

  int m = 8, nev = 1;
  int n, *ia, *ja; 
  double *a;
  double *evals, *evecs;
  double tc1, tc2, tw1, tw2;

  /* générer le problème */
  if (prob(m, &n, &ia, &ja, &a))
     return 1;

   /* Copier la matrice CSR*/
   if  (!write_csr_arrays("Mat_CSR",n,ia,ja,a)) return 1;
  printf("\nPROBLÈME: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
  



  
  /* allouer la mémoire pour vecteurs & valeurs propres */
  evals = malloc(nev * sizeof(double));
  evecs = malloc(nev * n * sizeof(double));

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
  const double *v0 = evecs;
  double bheta0 = evals[0];
  double *y = (double*)malloc((size_t)n * sizeof(double));
  matvec_csr(n, ia, ja, a, v0, y); /*Output y = Av0*/
  double res = residual_ratio(n, ia, ja, a, v0, bheta0);
  printf("Résidu relatif  ||Av0 - bheta*v0|| / ||v0|| = %.3e\n", res);
  free(y);


  /* libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs);
  return 0;
  
}

