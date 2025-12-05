#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "prob.h"
#include "time.h"
#include "interface_primme.h"
#include "csr_io.h"
#include "plotflux.h"
#include "eulerprog.h"
#include "arpack_interface.h"
#include <arpack/arpack.h>
#include "jadamilu_interface.h"

/* =========================================
 * Question 2 : calcul du résidu relatif
 *  ||A v - beta v||_2 / ||v||_2
 * =========================================
 */
double residual_ratio(int n,
                      const int    *ia,
                      const int    *ja,
                      const double *a,
                      const double *v,
                      double beta)
{
    double *Av = (double*)malloc((size_t)n * sizeof(double));
    if (!Av) {
        fprintf(stderr, "ERREUR: malloc(Av) dans residual_ratio\n");
        return -1.0; /* OOM */
    }

    matvec_csr(n, ia, ja, a, v, Av);

    double num = 0.0;
    double den = 0.0;

    for (int i = 0; i < n; ++i) {
        double r = Av[i] - beta * v[i];
        num += r * r;
        den += v[i] * v[i];
    }

    free(Av);

    if (den == 0.0) {
        fprintf(stderr, "ATTENTION: ||v|| = 0 dans residual_ratio\n");
        return -1.0;
    }

    return sqrt(num) / sqrt(den);
}

/* =========================================
 *                  main
 * =========================================
 */
int main(int argc, char *argv[])
{
    /* paramètres physiques / géométriques (adapter si besoin) */
    const double beta2 = 4.0;
    const double L     = 3.5;  /* taille du côté */
    const double tau   = 300.0;

    int m   = 80;   /* A regler pour valeur voulue */
    int nev = 1;    /* nombre de valeurs propres */

    double h = L / (m - 1);

    /* matrices et vecteurs */
    int    n   = 0;
    int   *ia  = NULL;
    int   *ja  = NULL;
    double *a  = NULL;

    double *evals  = calloc((size_t)nev, sizeof(double));
    double *evals2 = calloc((size_t)nev, sizeof(double));
    double *evecs  = NULL;  /* alloués après prob() */
    double *evecs2 = NULL;

    if (!evals || !evals2) {
        fprintf(stderr, "ERREUR: pas assez de mémoire pour evals/evals2\n");
        return 1;
    }

    double tc1, tc2, tw1, tw2;

    /* =========================================
     * Question 3 : alpha = 1, calcul de lambda1,
     * puis alpha_crit = sqrt(lambda1 / beta2)
     * =========================================
     */
    double alpha = 1.0;

    if (prob(m, alpha, &n, &ia, &ja, &a)) {
        fprintf(stderr, "ERREUR: prob(m, alpha=1)\n");
        return 1;
    }

    evecs  = (double*)malloc((size_t)nev * (size_t)n * sizeof(double));
    evecs2 = (double*)malloc((size_t)nev * (size_t)n * sizeof(double));
    if (!evecs || !evecs2) {
        fprintf(stderr, "ERREUR: malloc(evecs/evecs2)\n");
        return 1;
    }

    /* PRIMME (alpha = 1) */
    if (primme(n, ia, ja, a, nev, evals2, evecs2)) {
        fprintf(stderr, "ERREUR: primme(alpha=1)\n");
        return 1;
    }

    printf("Valeur propre PRIMME (alpha=1) = %.12g\n", evals2[0]);

    double lam = evals2[0];
    double alpha_crit = sqrt(lam / beta2);   /* alpha critique */

    printf("lambda1(alpha = 1) = %.12g => alpha_crit = %.12g\n",
           lam, alpha_crit);

    /* on repart avec la géométrie homothétique */
    free(ia); ia = NULL;
    free(ja); ja = NULL;
    free(a);  a  = NULL;

    /* =========================================
     * Question 6 :
     *   - générer le problème pour alpha_crit
     *   - PRIMME, ARPACK, JADAMILU
     *   - écrire la matrice CSR
     * 
     * INITIALEMENT j'avais utilisé le solveur ARPACK pour repondre a cette question 
     * cependant je décidé d'utiliser JADAMILU, par concequent j'ai comparé avec JADAMILU 
     * Mais le code pour ARPACK est totalement fontionnel, je l'ai pas effacé par peine
     * =========================================
     */ 

    if (prob(m, alpha_crit, &n, &ia, &ja, &a)) {
        fprintf(stderr, "ERREUR: prob(m, alpha_crit)\n");
        return 1;
    }

    evecs = realloc(evecs, (size_t)nev * (size_t)n * sizeof(double));
    evecs2 = realloc(evecs2, (size_t)nev * (size_t)n * sizeof(double));
    if (!evecs || !evecs2) {
        perror("realloc(evecs/evecs2)");
        return 1;
    }

/* debuger PROB*/
    if (!write_csr_arrays("Mat_CSR", n, ia, ja, a)) {
        fprintf(stderr, "ERREUR: write_csr_arrays\n");
        return 1;
    }

    printf("\nPROBLÈME (alpha_crit): m = %5d   n = %8d   nnz = %9d\n",
           m, n, ia[n]);

    if (!evals) {
        fprintf(stderr, "ERREUR: evals est NULL\n");
        return 1;
    }

    /* PRIMME - résolution */
    tc1 = mytimer_cpu();
    tw1 = mytimer_wall();

    if (primme(n, ia, ja, a, nev, evals2, evecs2)) {
        fprintf(stderr, "ERREUR: primme(alpha_crit)\n");
        return 1;
    }
   
    /* ARPACK */
    int nconv = 0;
    int err = arpackSym(
        n, ia, ja, a,
        /* nev  */ nev,
        /* ncv  */ 0,
        /* tol  */ 1e-10,
        /* it   */ 1000,
        /* which*/ "SA",
        evals, evecs,
        &nconv
    );

    /* JADAMILU */
    int nev_jada = 1;
    double evals_jada[1];
    double *evecs_jada = (double*)malloc((size_t)nev_jada * (size_t)n * sizeof(double));
    if (!evecs_jada) {
        fprintf(stderr, "ERREUR evec_jada\n");
        return 1;
    }

    if (jadamilu_solve(n, ia, ja, a,
                       nev_jada,
                       evals_jada,
                       evecs_jada)) {
        fprintf(stderr, "JADAMILU a échoué\n");
        free(evecs_jada);
        return 1;
    }

    if (err) {
        fprintf(stderr, "ARPACK failed (%d)\n", err);
        free(evecs_jada);
        return 1;
    }

    tc2 = mytimer_cpu();
    tw2 = mytimer_wall();

    printf("\nTemps de solution (CPU)    : %5.1f sec\n", tc2 - tc1);
    printf("Temps de solution (horloge): %5.1f sec\n", tw2 - tw1);
    
    printf("Valeur propre JADAMILU = %.15e\n", evals_jada[0]);
    printf("Valeur propre ARPACK : %15e\n", evals[0]);
    printf("Valeur propre PRIMME : %15e\n", evals2[0]);

    /*COMPARAISON*/
    double diff_abs = fabs(evals2[0] - evals_jada[0]);
    double diff_rel = diff_abs / fabs(evals_jada[0]);

    printf("Différence absolue PRIMME - JADAMILU = %.3e\n", diff_abs);
    printf("Différence relative = %.3e\n", diff_rel);
    
    printf("Check: lambda1(alpha_crit) = %.12g (bheta² = %.12g)\n",
           evals[0], beta2);

    /* affichage de la plus petite VP */
    for (int k = 0; k < nev; ++k) {
        printf("eval[%d] = %.15g\n", k, evals[k]);
    }

    /* =========================================
     * Question 2 : résidu relatif pour la VP
     * =========================================
     */
    const double *v0 = evecs2;
    double beta0 = evals2[0];

    double *y = (double*)malloc((size_t)n * sizeof(double));
    if (!y) {
        fprintf(stderr, "ERREUR: malloc(y)\n");
        free(evecs_jada);
        return 1;
    }
    matvec_csr(n, ia, ja, a, v0, y);

    double res = residual_ratio(n, ia, ja, a, v0, beta0);
    printf("Résidu relatif PRIMME ||A v0 - lambda v0|| / ||v0|| = %.3e\n", res);
    free(y);

    /* =========================================
     * Question 3 : dimensions pour cas stationnaire
     * =========================================
     */
    printf(" dimensions critique alpha_crit = %.3e\n", alpha_crit);
    printf("Nouveau pas : h' = %.3e\nAncien pas : h  = %.3e\n",
           h * alpha_crit, h);

    /* Normalisation pour affichage */
    double fluxmax = 0.0;
    int pmax = 0;
    for (int p = 0; p < n; ++p) {
        double ap = fabs(evecs[p]);
        if (ap > fluxmax) {
            fluxmax = ap;
            pmax = p;
        }
    }
    if (fluxmax > 0.0) {
        for (int p = 0; p < n; ++p)
            evecs[p] /= fluxmax;
    }
    if (evecs[pmax] < 0.0) {
        for (int p = 0; p < n; ++p)
            evecs[p] = -evecs[p];
    }
    printf("fluxmax : %.6e\n", fluxmax);

    /* =========================================
     * Question 5 : Euler progressif,
     *              alpha = 0.99*alpha_crit et 1.01*alpha_crit
     * =========================================
     */
    double crit[] = {0.99, 1.01};

    for (int si = 0; si < 2; ++si) {
        double alpha_run = crit[si] * alpha_crit;

        /* régénérer le problème pour alpha_run */
        if (ia) { free(ia); ia = NULL; }
        if (ja) { free(ja); ja = NULL; }
        if (a)  { free(a);  a  = NULL; }

        if (prob(m, alpha_run, &n, &ia, &ja, &a)) {
            fprintf(stderr, "ERREUR: prob(m, alpha_run)\n");
            free(evecs_jada);
            return 1;
        }

        double h_progg = alpha_run * h;
        double lam_max = 8.0 / (h_progg * h_progg);
        double dt = 0.8 * 2.0 / (tau * fmax(lam_max - beta2, 1e-12));

        /* temps pour que le flux double */
        double lam1_run = lam / (alpha_run * alpha_run);
        double Delta = fabs(beta2 - lam1_run);
        double R = 2.0;  /* flux doublé */
        double Ttarget = log(R) / (tau * fmax(Delta, 1e-12));
        int it = (int)ceil(Ttarget / dt);

        double *phi_init = (double*)malloc((size_t)n * sizeof(double));
        double *phi_prog = (double*)malloc((size_t)n * sizeof(double));
        if (!phi_init || !phi_prog) {
            fprintf(stderr, "ERREUR: malloc(phi_init/phi_prog)\n");
            free(phi_init);
            free(phi_prog);
            free(evecs_jada);
            return 1;
        }

        for (int i = 0; i < n; ++i)
            phi_init[i] = 1.0;

        euler_prog(n, ia, ja, a, beta2, it, dt,
                   phi_init, phi_prog);

        double *Z = (double*)malloc((size_t)m * (size_t)m * sizeof(double));
        if (!Z) {
            fprintf(stderr, "ERREUR: malloc(Z)\n");
            free(phi_init);
            free(phi_prog);
            free(evecs_jada);
            return 1;
        }

        decompress(phi_prog, m, alpha_run, L, Z);
        const char *ftxt  = (si == 0) ? "phi_099.txt" : "phi_101.txt";
        const char *title = (si == 0)
                            ? "eulerprog alpha=0.99"
                            : "eulerprog alpha=1.01";
        write_phi(ftxt, Z, m, alpha_run, L, L);
        plot_phi(ftxt, m, alpha_run, L, L, title, 0);

        free(Z);
        free(phi_init);
        free(phi_prog);
    }

    /* =========================================
     * Question 4 : visualisation du cas stationaire
     * =========================================
     */
    {
        double *Z = (double*)malloc((size_t)m * (size_t)m * sizeof(double));
        if (!Z) {
            fprintf(stderr, "ERREUR: malloc(Z) pour evecs\n");
            free(evecs_jada);
            return 1;
        }

        decompress(evecs, m, alpha_crit, L, Z);
        write_phi("evecs.txt", Z, m, alpha_crit, L, L);
        plot_phi("evecs.txt", m, alpha_crit, L, L,
                 "evecs alpha=alpha_crit", 0);

        free(Z);
    }

    /* libérer la mémoire */
    free(ia);
    free(ja);
    free(a);
    free(evals);
    free(evecs);
    free(evals2);
    free(evecs2);
    free(evecs_jada);

    return 0;
}
