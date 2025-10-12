#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* inside-rectangle test for the hole */
static inline int in_hole(double x, double y) {
    return (x > 1.5 && x < 3.5 && y > 2.0 && y < 3.5);
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

    nx = m - 2;                 /* nœuds de Dirichlet ne sont pas pris en compte */

    /* geometry: (0,3.5)x(0,3.5) */
    double L = 3.5;
    double h = L / (m - 1);
    double invh2 = 1.0 / (h*h);

    //invh2 = (m-1)*(m-1);    /* h^-2 pour L=3.5 */
    
    *n  = nx * nx;              /* nombre d'inconnues */
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

    nnz = 0; /*element non nul*/ 
    for (iy = 0; iy < nx; iy++) {
        for (ix = 0; ix < nx; ix++) {
            /* numéro de l'équation */
            ind = ix + nx * iy;

            /* marquer le début de la ligne suivante dans le tableau 'ia' */
            (*ia)[ind] = nnz;

            /* physical coords of this interior node */
            double x = (ix + 1) * h;
            double y = (iy + 1) * h;

            /* If the node is inside the hole: make an identity row and skip neighbors */
//            if (in_hole(x, y)) {
 //               (*a)[nnz]  = 1.0;
  //              (*ja)[nnz] = ind;   /* self-column */
   //             nnz++;
    //            continue;
     //       }
   
            /* SOUTH neighbor (iy-1) if valid and not in hole */
            if (iy > 0) {
                double xs = x;
                double ys = (iy    ) * h;      /* iy-1 -> (iy)*h on interior indexing */
                if (!in_hole(xs, ys)) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = ind - nx;
                    nnz++;
                }
            }

            /* WEST neighbor (ix-1) */
            if (ix > 0) {
                double xw = (ix    ) * h;
                double yw = y;
                if (!in_hole(xw, yw)) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = ind - 1;
                    nnz++;
                }
            }

            /* remplissage de la ligne : élément diagonal */
            (*a)[nnz] = 4.0*invh2; /* pour D=1 */
            (*ja)[nnz] = ind;
            nnz++;

            /* EAST neighbor (ix+1) */
            if (ix < nx - 1) {
                double xe = (ix + 2) * h;
                double ye = y;
                if (!in_hole(xe, ye)) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = ind + 1;
                    nnz++;
                }
            }

            /* NORTH neighbor (iy+1) */
            if (iy < nx - 1) {
                double xn = x;
                double yn = (iy + 2) * h;
                if (!in_hole(xn, yn)) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = ind + nx;
                    nnz++;
                }
            }
        }
    }

    /* dernier élément du tableau 'ia' */
    (*ia)[ind + 1] = nnz;

    /* retour habituel de fonction */
    return 0;
}
