// plotflux.c
#include <stdio.h>
#include <math.h>

static inline int in_hole(double x, double y) {
    return (x >= 1.5 && x <= 3.5 && y >= 2.0 && y <= 3.5);
}
/*Systeme d'indexation pour creer le phi qui serra lut par gnuplot*/
static inline int idx_border(int i, int j, int nx){
    return j*nx + i;
}
static inline int idx_interior(int i, int j, int m){
    return (j-1)*m + (i-1);    
}
/*Constructio de phi qui serra lut par gnuplot*/
void build_borders(const double *phi, int m, double *phiB){
    int nx = m + 2, ny = m + 2;
    for(int J = 0; J < ny; ++J)
        for (int I = 0; I < nx ; ++I)
            phiB[idx_border(I,J,nx)] = 0.0;
    
    for(int j=1; j<=m; ++j)
        for(int i=1; i <= m; ++i){
            int k = idx_interior(i,j,m);
            int KB = idx_border(i,j,nx);
            phiB[KB] = phi[k];
        }
}

void decompress_to_with_borders(
    const double *evec_compact, int m, double Lx, double Ly,
    double *phiB_out
){
    const int nx = m+2, ny = m+2;
    const double hx = Lx/(m+1), hy = Ly/(m+1);  

    for (int J=0; J<ny; ++J) /*mise a 0*/
        for (int I=0; I<nx; ++I)
            phiB_out[idx_border(I,J,nx)] = 0.0;
    
    int k = 0;
    for (int j=1; j<=m; ++j) {
        double y = j*hy;
        for (int i=1; i<=m; ++i) {
            double x = i*hx;
            int KB = idx_border(i,j,nx);
            if (in_hole(x,y)) {
                phiB_out[KB] = NAN;      // pas d’inconnue, trou visuel
            } else {
                phiB_out[KB] = evec_compact[k++];
            }
        }
    }
}

void mask_hole_inplace(double *phiB, int m, double Lx, double Ly){
    double hx = Lx/(m+1), hy = Ly/(m+1);
    for (int j = 0; j <= m+1; ++j) {
        double y = j*hy;
        for (int i = 0; i <= m+1; ++i) {
            double x = i*hx;
            if (in_hole(x,y)) {
                phiB[idx_border(i,j,m+2)] = NAN;
            }
        }
    }
}

int plot_phi(const double *phiB, int m, double Lx, double Ly,
             const char *title, const char *unit_label)
{
    FILE *gp = popen("gnuplot -persistent", "w");
    if (!gp) return -1;

    int nx = m + 2, ny = m + 2;
    double hx = Lx/(m+1), hy = Ly/(m+1);
    const char *u = (unit_label && unit_label[0]) ? unit_label : "";

    fprintf(gp, "set term qt noraise\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set view map\n");
    fprintf(gp, "set size ratio -1\n");
    fprintf(gp, "set pm3d map explicit interpolate 2,2\n");
    fprintf(gp, "set pm3d missingcolor rgb 'white'\n"); // trous (NaN) en blanc
    fprintf(gp, "set palette rgb 33,13,10\n");
    fprintf(gp, "set xlabel 'x (%s)'\n", u);
    fprintf(gp, "set ylabel 'y (%s)'\n", u);
    fprintf(gp, "set cblabel 'phi (a.d)'\n"); /*adimentionnal*/
    fprintf(gp, "set object 1 rect from 0,0 to %.12g,%.12g front fs empty border lc rgb 'black' lw 1.5\n", Lx, Ly);
    fprintf(gp, "set object 2 rect from 1.5,2.0 to 3.5,3.5 front fs empty border lc rgb 'black' lw 1.5\n");
    fprintf(gp, "set title '%s (m=%d, hx=%.6g %s, hy=%.6g %s)'\n",
            title ? title : "Flux neutronique", m, hx, u, hy, u);

    // Data block: x y z (NaN for hole)
    fprintf(gp, "$D << EOD\n");
    for (int j = 0; j < ny; ++j) {
        double y = j*hy;
        for (int i = 0; i < nx; ++i) {
            double x = i*hx;
            double z = phiB[idx_border(i, j, nx)];
            if (isnan(z)) fprintf(gp, "%.12g %.12g NaN\n", x, y);
            else          fprintf(gp, "%.12g %.12g %.12g\n", x, y, z);
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "EOD\n");

    // IMPORTANT: use splot (not plot) with pm3d
    fprintf(gp, "splot $D using 1:2:3 with pm3d\n");

    fflush(gp);
    pclose(gp);
    return 0;
}

/*
In the bash 
gcc plot_inline.c -lm -o plot_inline
./plot_inline
*/