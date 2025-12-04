#include <stdlib.h>
#include <stdio.h>


/* =========================================
                Question 1
   =========================================
*/
/* test “point dans le trou” (coordonnées physiques) */
static inline int in_hole(double x, double y, double alpha) {
    return (x >= 1.5*alpha && x <= 3.5*alpha &&
            y >= 2.0*alpha && y <= 3.5*alpha);
}
int prob(int m, double alpha ,int *n, int **ia, int **ja, double **a)

/*
  Arguments
  =========
  m (input)   - nombre de points par direction dans la grille 
  n  (output) - pointeur vers le nombre d'inconnues dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
 
*/
{
    int  nnz, ix, iy, nx, ind = 0;
    double invh2;
    nx = m - 2;

    /* Domaine aprés homothétie */
        double L = alpha*3.5;
        double h = L / (m - 1);
        invh2 = 1.0 / (h*h);

        /*vecteur Mapping de taille nx*nx avec les points numeroté et les trou en -1 bord exclus*/
        /*Mapping une facon dexprimer les coord 2D en un vecteur 1D bord exclus*/
        int Ngrid = nx * nx;
        int *map = (int*)malloc((size_t)Ngrid * sizeof(int));
        if (!map) return 1;

        int N = 0;
        for (iy = 0; iy < nx; ++iy) {
            for (ix = 0; ix < nx; ++ix) {
                ind = ix + nx * iy;
                /* check si c dans le trou */
                double x = (ix + 1) * h;
                double y = (iy + 1) * h;
                map[ind] = in_hole(x, y, alpha) ? -1 : N++;
            }
        }
        /*Finit mapping*/
        
        *n  = N;                     /* nombre réel d'inconnues (hors trou) */
        nnz = 5 * nx * nx - 4 * nx; /* nombre d'éléments non nuls */

        /* allocation des tableaux */

        *ia  = malloc((*n + 1) * sizeof(int));
        *ja  = malloc(nnz * sizeof(int));
        *a   = malloc(nnz * sizeof(double));
            
        /* allocation réussite? */

        if (*ia == NULL || *ja == NULL || *a == NULL ) {
            printf("\n ERREUR : pas assez de mémoire pour générer la matrice\n\n");
            return 1;
        }

        /* partie principale : remplissage de la matrice */
        
        nnz = 0;
        for (iy = 0; iy < nx; iy++) {
            for (ix = 0; ix < nx; ix++) {
                /* numéro de l'équation */
                ind = ix + nx * iy;

                /* ligne compacte r ; si -1 noeud suprimé */
                int r = map[ind];
                if (r < 0) continue;
                /* marquer le début de la ligne r dans 'ia' */
                (*ia)[r] = nnz;

                /* voisin sud */
                if (iy > 0)  {
                    int c = map[ind - nx];
                    if (c >= 0) {
                        (*a)[nnz]  = -invh2; /* pour D=1 */
                        (*ja)[nnz] = c;
                        nnz++;
                    }
                }

                /* voisin ouest */
                if (ix > 0)  { 
                    int c = map[ind - 1];
                    if (c >= 0) {
                        (*a)[nnz]  = -invh2; /* pour D=1 */
                        (*ja)[nnz] = c;
                        nnz++;
                    }
                }

                /*  diagonal */
                (*a)[nnz]  = 4.0*invh2; /* pour D=1 */
                (*ja)[nnz] = r;
                nnz++;

                /* voisin est */
                if (ix < nx - 1) {
                    int c = map[ind + 1];
                    if (c >= 0) {
                        (*a)[nnz]  = -invh2; /* pour D=1 */
                        (*ja)[nnz] = c;
                        nnz++;
                    }
                }

                /*voisin nord */
                if (iy < nx - 1) {
                    int c = map[ind + nx];
                    if (c >= 0) {
                        (*a)[nnz]  = -invh2; /* pour D=1 */
                        (*ja)[nnz] = c;
                        nnz++;
                    }
                }
            }
        }

        /* dernier élément du tableau 'ia' = nnz réel écrit */
        (*ia)[*n] = nnz;

        free(map);

        /* retour habituel de fonction */
        return 0;
    }

