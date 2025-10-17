#include <stdlib.h>
#include <stdio.h>

/* test “point dans le trou” (coordonnées physiques) */
static inline int in_hole(double x, double y) {
    return (x >= 1.5 && x <= 3.5 && y >= 2.0 && y <= 3.5);
}

int prob(int m, int *n, int **ia, int **ja, double **a)
/*
   But
   ===
   Génère la matrice n x n qui correspond à la discrétisation sur une grille 
   cartésienne régulière m x m de l'opérateur de Laplace à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )        sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec la fonction u qui satisfait les conditions aux limites de Dirichlet
         
         u = 0  sur (0,y), (1,y), (x,0) et (x,1), avec 0 <= x,y <= 1 .
  
  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice est retournée dans le format CRS
  qui est défini par le scalaire 'n' et les trois tableaux 'ia, 'ja' et 'a'.

  Arguments
  =========
  m (input)   - nombre de points par direction dans la grille 
  n  (output) - pointeur vers le nombre d'inconnues dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
 
  Sortie
  ======
  0 - exécution avec succès
  1 - erreurs
*/
{
    int  nnz, ix, iy, nx, ind = 0;
    double invh2;

    nx = m - 2; /* nœuds de Dirichlet ne sont pas pris en compte */

    /* Domaine (0,3.5)x(0,3.5) pour exprimer le trou en mètres */
        double L = 3.5;
        double h = L / (m - 1);
        invh2 = 1.0 / (h*h);

        /* ---- Pass 1 : renumérotation compacte (suppression des DOF du trou) ---- */
        int Ngrid = nx * nx;
        int *map = (int*)malloc((size_t)Ngrid * sizeof(int));
        if (!map) return 1;

        int rows = 0;
        for (iy = 0; iy < nx; ++iy) {
            for (ix = 0; ix < nx; ++ix) {
                ind = ix + nx * iy;
                /* coords physiques des nœuds intérieurs : ((ix+1)h, (iy+1)h) */
                double x = (ix + 1) * h;
                double y = (iy + 1) * h;
                map[ind] = in_hole(x, y) ? -1 : rows++;
            }
        }

        *n  = rows;                     /* nombre réel d'inconnues (hors trou) */
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

                /* ligne compacte r ; si -1 => nœud supprimé (dans le trou) */
                int r = map[ind];
                if (r < 0) continue;

                /* marquer le début de la ligne r dans 'ia' */
                (*ia)[r] = nnz;
        
                /* remplissage de la ligne : voisin sud */
                if (iy > 0)  {
                    int c = map[ind - nx];
                    if (c >= 0) {
                        (*a)[nnz]  = -invh2; /* pour D=1 */
                        (*ja)[nnz] = c;
                        nnz++;
                    }
                }

                /* remplissage de la ligne : voisin ouest */
                if (ix > 0)  {
                    int c = map[ind - 1];
                    if (c >= 0) {
                        (*a)[nnz]  = -invh2; /* pour D=1 */
                        (*ja)[nnz] = c;
                        nnz++;
                    }
                }

                /* remplissage de la ligne : élément diagonal */
                (*a)[nnz]  = 4.0*invh2; /* pour D=1 */
                (*ja)[nnz] = r;
                nnz++;

                /* remplissage de la ligne : voisin est */
                if (ix < nx - 1) {
                    int c = map[ind + 1];
                    if (c >= 0) {
                        (*a)[nnz]  = -invh2; /* pour D=1 */
                        (*ja)[nnz] = c;
                        nnz++;
                    }
                }

                /* remplissage de la ligne : voisin nord */
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

