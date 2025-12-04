#include "jadamilu_interface.h"
#include <stdlib.h>
#include <stdio.h>

/* Tu n'as pas besoin de <string.h> ici pour l'instant */

/* === 1) Prototype Fortran JADAMILU ===================== */
/* À adapter si DPJD_prototype.f a plus/d'autres arguments */
extern void dpjd_(
    int *N,
    int *IA,
    int *JA,
    double *A,
    int *NEV,
    double *EV,
    double *X,
    double *WORK,
    int *LWORK,
    int *IPARAM,
    double *RPARAM,
    int *IFLAG
);

/* === 2) Conversion CSR vers upper-tri Fortran 1-based ==== */

int buildUpper_CSR(
    int n,
    const int *ia, const int *ja, const double *a,
    int **iaF_out, int **jaF_out, double **aF_out,
    int *nnzF_out
){
    int i, k;

    /* 1) Compter les éléments non nuls dans la partie triangulaire sup */
    int nnzF = 0;
    for (i = 0; i < n; ++i) {
        for (k = ia[i]; k < ia[i+1]; ++k) {
            int j = ja[k];
            if (j >= i) nnzF++;
        }
    }

    /* 2) Allouer */
    int *iaF = (int*)malloc((size_t)(n+1) * sizeof(int));
    int *jaF = (int*)malloc((size_t)nnzF * sizeof(int));
    double *aF = (double*)malloc((size_t)nnzF * sizeof(double));
    if (!iaF || !jaF || !aF) {
        free(iaF); free(jaF); free(aF);
        return 1;
    }

    /* 3) Remplir en 1-based pour Fortran */
    int pos = 0;
    iaF[0] = 1;  /* ligne 1 commence à l’entrée 1 */

    for (i = 0; i < n; ++i) {
        for (k = ia[i]; k < ia[i+1]; ++k) {
            int j = ja[k];
            if (j >= i) {
                jaF[pos] = j + 1;   /* 0-based -> 1-based */
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

/* === 3) Wrapper C pour JADAMILU ========================= */

int jadamilu_solve(
    int n, const int *ia, const int *ja, const double *a,
    int nev,
    double *evals,  /* taille >= nev */
    double *evecs   /* taille >= nev * n (colonnes) */
){
    int *iaF = NULL, *jaF = NULL;
    double *aF = NULL;
    int nnzF;          /* <-- manquait chez toi */
    int info;

    if (buildUpper_CSR(n, ia, ja, a, &iaF, &jaF, &aF, &nnzF)) {
        fprintf(stderr, "[JADAMILU] erreur alloc / conversion CSR\n");
        return 1;
    }

    /* -------- paramètres JADAMILU -------- */

    int N   = n;
    int NEV = nev;

    int IFLAG = 0;
    int    IPARAM[40];   /* vérifier la taille exacte dans la doc */
    double RPARAM[40];

    for (int i = 0; i < 40; ++i) {
        IPARAM[i] = 0;
        RPARAM[i] = 0.0;
    }

    /* exemple : tu pourras fixer des options ici :
       IPARAM[0] = 0; etc. selon le manuel JADAMILU */

    int LWORK = 10 * N * NEV;  /* valeur prudente; voir doc pour mieux */
    double *WORK = (double*)malloc((size_t)LWORK * sizeof(double));
    if (!WORK) {
        fprintf(stderr, "[JADAMILU] erreur alloc WORK\n");
        free(iaF); free(jaF); free(aF);
        return 1;
    }

    /* appel Fortran */
    dpjd_(
        &N,
        iaF, jaF, aF,
        &NEV,
        evals,
        evecs,
        WORK, &LWORK,
        IPARAM, RPARAM, &IFLAG
    );

    info = IFLAG;

    free(WORK);
    free(iaF); free(jaF); free(aF);

    if (info != 0) {
        fprintf(stderr, "[JADAMILU] DPJD a retourné IFLAG=%d\n", info);
        return 1;
    }

    return 0;
}
