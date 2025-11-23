#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "time.h"
#include "interface_primme.h"
#include "csr_io.h"
#include "plotflux.h"
#include <math.h>
#include "eulerprog.h"

/*Calcul du résidus*/
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
  double beta2 = 4.0;
  int m = 90, nev = 1;
  double h = 3.5/(m - 1);
  int n, *ia, *ja; 
  double *a;
  double *evals = calloc(nev, sizeof(double));
  double *evecs = calloc((size_t)nev *  (size_t)/*n après prob*/ 1, sizeof(double)); // pas oublier de réalluer aprés prob

  /*Sollution pour alpha = 1*/
  double alpha = 1.0;
  if (prob(m, alpha, &n, &ia, &ja, &a)) return 1;
  evecs = realloc(evecs, (size_t)nev * (size_t)n * sizeof(double));
  if (primme(n, ia, ja, a, nev, evals, evecs)) return 1;
  double lam = evals[0];

  /*Nouvel alpha critique*/
  double alpha_crit = sqrt(lam / beta2); /*Mettre a 1 pour le test 1 et m = 8*/
  printf("lambda1(alpha = 1) = %.12g => alpha_crit %.12g \n", lam, alpha_crit);

  free(ia); free(ja); free(a); 
  double tc1, tc2, tw1, tw2;

  /* générer le problème */
  if (prob(m, alpha_crit ,&n, &ia, &ja, &a)) return 1;
  evecs = realloc(evecs, (size_t)nev * (size_t)n * sizeof(double));
  if (!evecs) { perror("realloc evecs"); return 1; } /*vérifier ca*/

   /* Copier la matrice CSR*/
  if  (!write_csr_arrays("Mat_CSR",n,ia,ja,a)) return 1;
  printf("\nPROBLÈME: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
  /* allouer la mémoire pour vecteurs & valeurs propres */
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
  printf("Check: lambda1(alpha_crit) = %.12g (cible = %.12g)\n", evals[0], beta2); /*Comparaison avec lambda*/

  
  /* Affichage de la plus petite propre des valeurs propres */
  for (int k = 0; k < nev; ++k) { 
      printf("eval[%d] = %.15g\n", k, evals[k]);
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

  printf("Les dimentions fois alpha = %.3e\n",alpha_crit); 
  /*La valeur de h doit varier, en fixant m, h' = alpha*h */
  printf("Nouveau pad : h' = %.3e\n Ancien pad h = %3.e \n", h*alpha_crit, h);

  /*Normalisation pour optimiser l'affichage*/

  double maxabs = 0.0;
  int pmax = 0;
  for (int p = 0; p < n; ++p) {
      double ap = fabs(evecs[p]);
      if (ap > maxabs) { maxabs = ap; pmax = p; }
  }
  if (maxabs > 0.0) {
    for (int p = 0; p < n; ++p) evecs[p] /= maxabs; // Normalisation de evecs
  }
  if (evecs[pmax] < 0.0) for (int p = 0; p < n; ++p) evecs[p] = -evecs[p];
  printf("maxabs : %.6e\n", maxabs);


/*Boucle avec les scénarios possible en variant alpha_crit*/

double scales[] = {0.99,1.01};
for (int si = 0; si < 2; ++si){
  double alpha_run = scales[si] * alpha_crit;

  if (ia) { free(ia); ia = NULL; }
  if (ja) { free(ja); ja = NULL; }
  if (a)  { free(a);  a  = NULL; }
  if (prob(m, alpha_run, &n, &ia, &ja, &a)) return 1;

  double h_progg = alpha_run*h;
  double lam_max = 8.0/(h_progg*h_progg);
  double tau = 300.0;
  double dt = 0.8 * 2.0 / (tau * fmax(lam_max - beta2, 1e-12)); /*fmax pour ne pas renvoyer 0 et j'ai pris 80% de la borne */


  // Nombre d'itération : Quand le flux est doublé
  double lam1_run = lam / (alpha_run * alpha_run);
  double Delta = fabs(beta2 - lam1_run);
  double R = 2.0; // 2 car flux double
  double Ttarget = log(R) / (300.0 * fmax(Delta, 1e-12));
  int it = (int)ceil(Ttarget / dt);


  double *phi_init = (double*)malloc(n*sizeof(double));
  for(int i = 0; i < n ; ++i)phi_init[i] = 1.0;
  double *phi_prog = (double*)malloc(n*sizeof(double));

  euler_prog(n, ia, ja, a, beta2,  it,  dt,
      phi_init,  phi_prog);

  double *Z = (double*)malloc((size_t)m*(size_t)m*sizeof(double));
  decompress(phi_prog, m, alpha_run, 3.5, Z);
  const char* ftxt = (si==0) ? "phi_099.txt" : "phi_101.txt";
  write_phi(ftxt, Z, m, alpha_run, 3.5, 3.5);
  plot_phi(ftxt, m, alpha_run, 3.5, 3.5,
             (si==0) ? "eulerprogg alpha=0.99" : "eulerprogg alpha=1.01",
             0);
  free(Z);
  free(phi_init); free(phi_prog);
} 

  double *Z = (double*)malloc((size_t)m*(size_t)m*sizeof(double));
  decompress(evecs, m, alpha_crit, 3.5, Z);
  write_phi("evecs.txt", Z, m, alpha_crit, 3.5, 3.5);
  plot_phi("evecs.txt", m, alpha_crit, 3.5, 3.5,"evecs alpha=1",
             0);


  /* libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs); free(Z);
  return 0;
  
}