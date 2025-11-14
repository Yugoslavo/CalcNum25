#ifndef PLOTFLUX_H
#define PLOTFLUX_H
int   in_hole(double x, double y, double alpha);
int   build_map(int m, double alpha, double Lx, int *map);
void  decompress(const double *evec,
                                  int m, double alpha, double Lx,
                                  double *Z);
int   write_phi(const char *path, const double *Z, int m,
                            double alpha, double Lx, double Ly);
int   plot_phi(const char *path, int m, double alpha, double Lx, double Ly,
                const char *title, int as3d);

#endif