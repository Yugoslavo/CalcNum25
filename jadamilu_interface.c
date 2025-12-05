#include "jadamilu_interface.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* D’après DPJD_prototype.f */
extern void dpjd_(
    int    *N,
    double *A,
    int    *JA,
    int    *IA,
    double *EIGS,
    double *RES,
    double *X,
    int    *LX,
    int    *NEIG,
    double *SIGMA,
    int    *ISEARCH,
    int    *NINIT,
    int    *MADSPACE,
    int    *ITER,
    double *TOL,
    double *SHIFT,
    double *DROPTOL,
    double *MEM,
    int    *ICNTL,
    int    *IPRINT,
    int    *INFO,
    double *GAP
);

/*CSR end 1-forme fortran*/
/* Le code doit envoyer ladresse des tables CSR a cause de fortran
*/

int buildUpper_CSR(
    int n,
    const int *ia, const int *ja, const double *a,
    int **iaF_out, int **jaF_out, double **aF_out,
    int *nnzF_out
){
    int i, k;

    /*Compter les nnz dans la partie triangulaire sup */
    int nnzF = 0;
    for (i = 0; i < n; ++i) {
        for (k = ia[i]; k < ia[i+1]; ++k) {
            int j = ja[k];
            if (j >= i) nnzF++;
        }
    }

    /* Allouer*/
    int *iaF = (int*)malloc((size_t)(n+1) * sizeof(int));
    int *jaF = (int*)malloc((size_t)nnzF * sizeof(int));
    double *aF = (double*)malloc((size_t)nnzF * sizeof(double));
    if (!iaF || !jaF || !aF) {
        free(iaF); free(jaF); free(aF);
        return 1;
    }

    /*Remplir en 1-based pour Fortran pas oublier que en fortran ca commence en 1 */
    int pos = 0;
    iaF[0] = 1;  /* ligne 1 commence à l’entrée 1 */

    for (i = 0; i < n; ++i) {
        for (k = ia[i]; k < ia[i+1]; ++k) {
            int j = ja[k];
            if (j >= i) {
                jaF[pos] = j + 1;   /*  1-based */
                aF[pos]  = a[k];
                pos++;
            }
        }
        iaF[i+1] = pos + 1; /* début de la ligne suivante (1-based) */
    }

    *iaF_out  = iaF;
    *jaF_out  = jaF;
    *aF_out   = aF;
    *nnzF_out = nnzF;
    return 0;
}


int jadamilu_solve(
    int n, const int *ia, const int *ja, const double *a,
    int nev,
    double *evals,  /* taille >= nev */
    double *evecs   /* taille >= nev * n (stockage colonnes successives) */
){
    int *iaF = NULL, *jaF = NULL;
    double *aF = NULL;
    int nnzF;
    int info;

    /* CSR 1-based */
    if (buildUpper_CSR(n, ia, ja, a, &iaF, &jaF, &aF, &nnzF)) {
        fprintf(stderr, "[JADAMILU] erreur alloc / conversion CSR\n");
        return 1;
    }

    /* PARAMETRES */

    int N = n;
    int NVP = nev;      /* nombre VP à demander */

    /* tableau de résidus (normes des résidus) */
    double *RES = (double*)calloc((size_t)NVP, sizeof(double));
    if (!RES) {
        fprintf(stderr, "[JADAMILU] OOM RES\n");
        free(iaF); free(jaF); free(aF);
        return 1;
    }

    /* Paramètres de contrôle */
    /*Quasi tout a 0*/

    double SIGMA   = 0.0;      /* shift  */
    int    ISEARCH = 0;        /* 0, + petite VP */
    int    NINIT   = 0;        /* JADAMILU choisis les vecteurs */

    int MADSPACE   = 20;       /* dimension espace de krylov */
    int ITER       = 1000;     /* nb max de multiplications matrice-vecteur */
    double TOL     = 1e-8;     /* norme du résidu que je veut, 10^-8 ca me parait resonable */

    double SHIFT   = 0.0;      /* voir doc, on met 0 pour démarrer */
    double DROPTOL = 1e-3;     /* drop tolerance ILU */
    double MEM     = 10.0;     /*factorisation de 10 parait pas trop lourd*/

    int ICNTL[5];
    for (int i = 0; i < 5; ++i) ICNTL[i] = 0; /*Tout a 0, mit a default*/

    int IPRINT = 0;            /* Cas derreur importante */
    int INFO   = 0;
    double GAP = 0.0;

    /* Taille LX pour X : on prend la borne recommandée MADSPACE>=3 : */
   /*Calcul de la taille du Workspace X*/

    int ms2 = MADSPACE * MADSPACE;
    int LX  = N * (3* MADSPACE + NVP + 1)
              + 3 * ms2
              + (ms2 > NVP ? ms2 : NVP);

    double *X = (double*)malloc((size_t)LX * sizeof(double));
    if (!X) {
        fprintf(stderr, "[JADAMILU] OOM X (LX=%d)\n", LX);
        free(RES);
        free(iaF); free(jaF); free(aF);
        return 1;
    }

    /* EIGS = evals (JADAMILU écrira dedans) */
    /* X contiendra les vecteurs propres, 
    Pour les scalaires je dois déferencer car fortran accepte que des adresse
    */
    dpjd_(
        &N,
        aF,
        jaF,
        iaF,
        evals,
        RES,
        X,
        &LX,
        &NVP,
        &SIGMA,
        &ISEARCH,
        &NINIT,
        &MADSPACE,
        &ITER,
        &TOL,
        &SHIFT,
        &DROPTOL,
        &MEM,
        ICNTL,
        &IPRINT,
        &INFO,
        &GAP
    );

    info = INFO;

    /* Copier les vecteurs propres de X vers evecs (1..NEIG) */
    if (info == 0 && NVP >= 1) {
        for (int k = 0; k < NVP; ++k) {
            /* en Fortran : x(i) = X(1+N*(k):N*(k+1)) ; en C on a X[0..LX-1] */
            for (int i = 0; i < N; ++i) {
                evecs[k*(size_t)N + i] = X[k*(size_t)N + i];
            }
        }
    }

    free(X);
    free(RES);
    free(iaF); free(jaF); free(aF);

    if (info != 0) {
        fprintf(stderr,
                "[JADAMILU] DPJD a retourné INFO=%d (NEIG=%d, GAP=%.3e)\n",
                info, NVP, GAP);
        return 1;
    }

    return 0;
}
