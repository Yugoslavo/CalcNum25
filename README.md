#include <stdlib.h>
#include <stdio.h>

static inline int in_hole(double x, double y) {
    return (x >= 1.5 && x <= 3.5 && y >= 2.0 && y <= 3.5);
}

int prob(int m, int *n, int **ia, int **ja, double **a) {
    if (m < 3) return 1;

    int ix, iy, nx = m - 2;
    int Ngrid = nx * nx;

    double L = 3.5;                    // adapte selon ton domaine
    double h = L / (m - 1);
    double invh2 = 1.0 / (h*h);

    // --- Pass 1: mapping des DOF valides
    int *map = (int*)malloc((size_t)Ngrid * sizeof(int));
    if (!map) return 1;

    int rows = 0;
    for (iy = 0; iy < nx; ++iy) {
        for (ix = 0; ix < nx; ++ix) {
            int ind = ix + nx*iy;
            double x = (ix + 1) * h;
            double y = (iy + 1) * h;
            if (in_hole(x, y)) map[ind] = -1;   // DOF supprimé
            else               map[ind] = rows++;
        }
    }
    *n = rows;
    if (*n == 0) { free(map); return 1; }

    // --- Compter nnz exactement
    int nnz = 0;
    for (iy = 0; iy < nx; ++iy) {
        for (ix = 0; ix < nx; ++ix) {
            int ind = ix + nx*iy;
            int r = map[ind];
            if (r < 0) continue;       // nœud supprimé

            nnz += 1;                  // diagonale
            if (iy > 0   && map[ind-nx] >= 0) ++nnz;  // S
            if (ix > 0   && map[ind-1]  >= 0) ++nnz;  // W
            if (ix < nx-1&& map[ind+1]  >= 0) ++nnz;  // E
            if (iy < nx-1&& map[ind+nx] >= 0) ++nnz;  // N
        }
    }

    // --- Allocation CSR
    *ia = (int*)   malloc((size_t)(*n + 1) * sizeof(int));
    *ja = (int*)   malloc((size_t)nnz       * sizeof(int));
    *a  = (double*)malloc((size_t)nnz       * sizeof(double));
    if (!*ia || !*ja || !*a) {
        free(map); free(*ia); free(*ja); free(*a);
        *ia = NULL; *ja = NULL; *a = NULL;
        return 1;
    }

    // --- Pass 2: remplissage CSR (ordre S, W, diag, E, N)
    int k = 0;
    for (iy = 0; iy < nx; ++iy) {
        for (ix = 0; ix < nx; ++ix) {
            int ind = ix + nx*iy;
            int r = map[ind];
            if (r < 0) continue;

            (*ia)[r] = k;

            if (iy > 0)   { int c = map[ind - nx]; if (c >= 0) { (*a)[k] = -invh2; (*ja)[k] = c; ++k; } } // S
            if (ix > 0)   { int c = map[ind - 1 ]; if (c >= 0) { (*a)[k] = -invh2; (*ja)[k] = c; ++k; } } // W

            (*a)[k] = 4.0*invh2; (*ja)[k] = r; ++k;                                                     // diag

            if (ix < nx-1){ int c = map[ind + 1 ]; if (c >= 0) { (*a)[k] = -invh2; (*ja)[k] = c; ++k; } } // E
            if (iy < nx-1){ int c = map[ind + nx]; if (c >= 0) { (*a)[k] = -invh2; (*ja)[k] = c; ++k; } } // N
        }
    }
    (*ia)[*n] = k;

    free(map);
    return 0;
}
