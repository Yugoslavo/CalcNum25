// plotflux.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "plotflux.h"

int in_hole(double x, double y, double alpha) {
    return (x >= 1.5*alpha && x <= 3.5*alpha &&
            y >= 2.0*alpha && y <= 3.5*alpha);
}
/*Systeme d'indexation pour creer le phi qui serra lut par gnuplot*/
/*Meme indexisation que utilisé en prob*/
int build_map(int m, double alpha, double Lx, int *map){ 
    const int nx = m - 2;
    const double L = alpha * Lx;
    const double h = L / (m - 1);
    int N = 0; /*Valeurs non nulles*/
    for (int j = 1; j <= nx; ++j) {
        double y = j*h;
        for (int i = 1; i <= nx; ++i) {
            double x = i*h;
            map[(j-1)*nx + (i-1)] = in_hole(x,y,alpha) ? -1 : N++;
        }
    }
    return N;
}

void decompress(const double *evec, int m, double alpha, double Lx,
                                 double *Z)
{   
    /*Tout mettre a 0*/
    const int nx = m - 2;
    for (int K = 0; K < m*m; ++K) Z[K] = 0.0;
    
    /*Generer le Map*/
    int *map = (int*)malloc((size_t)nx*nx*sizeof(int));
    (void)build_map(m, alpha, Lx, map);

    
    for (int j = 1; j <= nx; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = map[(j-1)*nx + (i-1)]; /*Bord exclus*/
            if (k >= 0) Z[j*m + i] = evec[k];
            else        Z[j*m + i] = NAN; /*Trou visuel, rien a la place de 0*/
        }
    }
    free(map);

    // signe global positif au cas ou
    double s = 0.0;
    for (int K = 0; K < m*m; ++K) if (!isnan(Z[K])) s += Z[K];
    if (s < 0.0) for (int K = 0; K < m*m; ++K) if (!isnan(Z[K])) Z[K] = -Z[K];
}

/*Vecteur Z, taille m*m lut par gnuplot*/
int write_phi(const char *path, const double *Z, int m,
                          double alpha, double Lx, double Ly)
{
    FILE *f = fopen(path, "w");
    if (!f) return -1;

    const double Lx_eff = alpha*Lx, Ly_eff = alpha*Ly;
    const double hx = Lx_eff / (m - 1);
    const double hy = Ly_eff / (m - 1);

    for (int J = 0; J < m; ++J) {
        double y = J*hy;
        for (int I = 0; I < m; ++I) {
            double x = I*hx;
            double z = Z[J*m + I];
            if (isnan(z)) fprintf(f, "% .16g % .16g NaN\n", x, y); /*Sinon gnuplot ne lit pas les NAN*/
            else          fprintf(f, "% .16g % .16g % .16g\n", x, y, z);
        }
        fputc('\n', f); // séparateur de blocs pour pm3d
    }
    fclose(f);
    return 0;
}
/*Option plot Heatmap et 3D*/
int plot_phi(const char *path, int m, double alpha, double Lx, double Ly,
              const char *title, int as3d)
{
    FILE *gp = popen("gnuplot -persistent", "w");
    if (!gp) return -1;

    double Lx_eff = alpha * Lx, Ly_eff = alpha * Ly;
    double hx = Lx_eff / (m - 1), hy = Ly_eff / (m - 1);

    fprintf(gp, "set term qt noraise\n");
    fprintf(gp, "set mouse on\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "set datafile missing 'NaN'\n");
    fprintf(gp, "set pm3d missingcolor rgb 'white'\n");
    fprintf(gp, "set pm3d corners2color median\n");
    fprintf(gp, "set pm3d interpolate 1,1\n");
    fprintf(gp, "set palette rgb 33,13,10\n");
    fprintf(gp, "set xlabel 'x (%s)'\n", "m");
    fprintf(gp, "set ylabel 'y (%s)'\n", "m");
    fprintf(gp, "set zlabel 'phi (%s)'\n", "a.u.");
    fprintf(gp, "set xrange [0:%.12g]\n", Lx_eff);
    fprintf(gp, "set yrange [0:%.12g]\n", Ly_eff);

    /*delimiter contour*/
    fprintf(gp, "set object 1 rect from 0,0 to %.12g,%.12g front fs empty border lc rgb 'black' lw 1.5\n",
            Lx_eff, Ly_eff);
    /*delimiter contour trou*/
    fprintf(gp, "set object 2 rect from %.12g,%.12g to %.12g,%.12g front fs empty border lc rgb 'black' lw 1.5\n",
            1.5*alpha, 2.0*alpha, 3.5*alpha, 3.5*alpha);
    /*Titre*/
    fprintf(gp, "set title '%s (m=%d, hx=%.6g %s, hy=%.6g %s, %s)'\n",
            title ? title : "Flux neutronique", m, hx, "m", hy, "m", as3d ? "3D" : "heatmap"); /*Si as3d 3D alors heatmap*/
    
    /*Lire le .TXT*/
    if (as3d) {
        fprintf(gp, "set view 60,35\n");
        fprintf(gp, "set pm3d explicit\n");
        fprintf(gp, "set xyplane 0\n");
        fprintf(gp, "splot '%s' using 1:2:3 with pm3d\n", path);
    } else {
        fprintf(gp, "set view map\n");
        fprintf(gp, "set pm3d map explicit\n");
        fprintf(gp, "splot '%s' using 1:2:3 with pm3d\n", path);
    }
    fprintf(gp, "pause mouse close\n");
    fflush(gp); pclose(gp);
    return 0;
}
