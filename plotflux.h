void decompress_to_with_borders(
    const double *evec_compact, int m, double Lx, double Ly,
    double *phiB_out
);
static inline int idx_border(int i, int j, int nx); 
static inline int idx_interior(int i, int j, int m); 
static inline int in_hole(double x, double y); 
int plot_phi(const double *phiB, int m, double Lx, double Ly, const char *title, const char *unit_label);
