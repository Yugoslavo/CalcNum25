#include "arpack_interface.h"
#include <arpack/arpack.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* CSR SpMV (0-based CSR) */
static void matvec_csr(int n, const int *ia, const int *ja, const double *a,
                             const double *x, double *y)
{
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int p = ia[i]; p < ia[i+1]; ++p) s += a[p] * x[ ja[p] ];
        y[i] = s;
    }
}

double arpackResidu(
    int n, const int *ia, const int *ja, const double *a,
    const double *v, double lambda)
{
    double *Av = (double*)malloc((size_t)n * sizeof(double));
    if (!Av) return -1.0;
    matvec_csr(n, ia, ja, a, v, Av);
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
        double r = Av[i] - lambda * v[i];
        num += r*r; den += v[i]*v[i];
    }
    free(Av);
    return (den > 0.0) ? sqrt(num)/sqrt(den) : -1.0;
}

int arpackSym(
    int n,
    const int *ia, const int *ja, const double *a,
    int nev, /*nbr de val prop*/
    int ncv, /*Quel vecteur ritz extraire*/
    double tol, int maxit,
    const char *which,               
    double *evals, double *evecs, int *nconv_out)
{
    if(n<=0 || nev <= 0 || !ia || !ja || !a || !evals || !evecs){
        fprintf(stderr, "Problem d'argument\n");
        return 1;
    }

    /* ARPACK int type _ai, 32-bit*/
    a_int n_ai   = (a_int)n;
    a_int nev_ai = (a_int)nev;

    /* taille Arnoldi space ((source doc ARPACK))*/
    if (ncv <= 0) ncv = nev + 20; /*régle pratique, stackoverflow*/
    if (ncv < nev + 2) ncv = nev + 2;
    if (ncv > n) ncv = n;
    a_int ncv_ai = (a_int)ncv;

    /*Format _ia de ARPACK*/
    a_int ldv_ai     = n_ai;
    a_int ldz_ai     = n_ai;
    a_int lworkl_ai  = ncv_ai * (ncv_ai + 8); /*taille minimale demandé par ARPACK du workspace*/

    /* Workspace */
    double *V     = (double*)calloc((size_t)n_ai * (size_t)ncv_ai, sizeof(double)); 
    double *workd = (double*)calloc(3*(size_t)n_ai, sizeof(double));
    double *workl = (double*)calloc((size_t)lworkl_ai, sizeof(double));
    double *resid = (double*)calloc((size_t)n_ai, sizeof(double));
    a_int  *select= (a_int*)  calloc((size_t)ncv_ai, sizeof(a_int));
    if (!V || !workd || !workl || !resid || !select) {
        fprintf(stderr, "Problem avec workspace\n");
        free(V); free(workd); free(workl); free(resid); free(select);
        return 2;
    }

    /*"tableau de bord pour l'utilisation du reverse communication"*/
    a_int iparam[11] = {0};
    a_int ipntr[14]  = {0};     /*Pointeur 1-based au workd */
    iparam[0] = 1;              /* ISHIFT = 1 (valeur exact) */
    iparam[2] = (a_int)maxit;   /* max iterations */
    iparam[6] = 1;              /*standard Av = λv */
    /*le reste des valeurs sont calculée de ARPACK*/

    a_int ido = 0;
    a_int info = 0; /*Status code proposé par ARPACK*/
    const char bmat[] = "I";    /*problem : Ax = λx pour calculer les val propres*/

    /* 
    Arnoldi Loop
    Il utilise le Lanczos algorithm, avec la "reverse communication"
    "Reverse communication": a chaque fois qu'il a besoin dun produit matrice vecteur 
    le process pause
    */
    
    for (;;) {
        dsaupd_c(&ido, /*reverse communication flag*/
                 bmat,
                n_ai, 
                which, 
                nev_ai, 
                tol,
                resid,
                ncv_ai,
                V,
                ldv_ai,
                iparam,
                ipntr, 
                workd, 
                workl,
                lworkl_ai,
                &info
            );

        if (ido == -1 || ido == 1) {
            /* y = A x  */
            double *x = &workd[ipntr[0]-1];
            double *y = &workd[ipntr[1]-1];
            matvec_csr(n, ia, ja, a, x, y);
            continue;
        }
        break;
    }

    if (info != 0) {
        fprintf(stderr,"Problem dsaupd \n", (int)info);
        free(V); free(workd); free(workl); free(resid); free(select);
        return 3;
    }

    /* extract eigenvalues */
    a_int rvec = 1;                /* return Ritz vectors */
    const char howmny[] = "A";     
    double sigma = 0.0;            /* unused in mode 1 */
    a_int ierr = 0;

    dseupd_c(
        rvec,   /*retourne vecprop*/
        howmny, /*que ceux qui converge*/
        select, /*ignoré car howmny = A*/ 
        evals,
        evecs,
        ldz_ai, /*dimention the evec (n)*/
        sigma,  /*pas utilisé*/
        bmat,   /*I car mode 1*/
        n_ai,   /*dimention matrice*/
        which,
        nev_ai,
        tol,
        resid,  /*workspace crée par dsaupd*/
        ncv_ai, /*taille subspace*/
        V,      /*base*/
        ldv_ai, /*dimention base (n)*/
        iparam,
        ipntr,
        workd,  /*workspace*/
        workl,  /*taille workspace*/
        lworkl_ai,
        &ierr
    );

    if (ierr != 0) {
        fprintf(stderr,"[ARPACK] dseupd failed (ierr=%d)\n", (int)ierr);
        free(V); free(workd); free(workl); free(resid); free(select);
        return 4;
    }

    if (nconv_out) *nconv_out = (int)iparam[4]; /* #converged Ritz values */

    free(V); free(workd); free(workl); free(resid); free(select);
    return 0;
}