void read_input(char filename[100], int *NUM_POINTS, double X[], int *GSI_MODEL);

void defzeros(int i, char **zeros);

double testparticle(char filename[100], double Umag, double theta, double phi, double Ts, double Ta, double epsilon, double alpha, double alphan, double sigmat, double X[], int GSI_MODEL, int iTOT);

int read_num_lines(char filename[100]);

void facet_properties(char filename[100], double Ux, double Uy, double Uz, int nfacets, struct facet_struct *pfacet, double Umag);

void domain_boundary(double *XMIN, double *YMIN, double *ZMIN, double *XMAX, double *YMAX, double *ZMAX, int nfacets, struct facet_struct *pfacet);

void max_facet_dimension(double XMIN, double YMIN, double ZMIN, double XMAX, double YMAX, double ZMAX, struct facet_struct *pfacet, int nfacets, int nc[3]);

void create_grid(double XMIN, double YMIN, double ZMIN, double XMAX, double YMAX, double ZMAX, struct cell_struct *pcell, struct facet_struct *pfacet, int nfacets, int nc[3]);

void sort_facets(struct facet_struct *pfacet, struct cell_struct *pcell, int nfacets, int nc[3]);

double projected_area(double XMIN, double YMIN, double ZMIN, double XMAX, double YMAX, double ZMAX, double Ux, double Uy, double Uz, int nfacets, struct facet_struct *pfacet);

void pposition(double XMIN, double YMIN, double ZMIN, double XMAX, double YMAX, double ZMAX, struct particle_struct *pparticle, int ipart, int particle_surf);

void pvelocity(double s[NSURF][NSPECIES], double Usurf[], int direction[], double Cmp[], int particle_surf, struct particle_struct *pparticle, int ipart, int species);

void cell_track(struct particle_struct *pparticle, int ipart, struct cell_struct *pcell, int nc[3], int particle_cell_track[MAXCELLS]);

void compile_facet_list(struct facet_struct *pfacet, int nfacets, struct cell_struct *pcell, int particle_cell_track[MAXCELLS], int *total_facet_list, int *fcount, int ipart);

void pf_intercept(int *min_fc, struct facet_struct *pfacet, struct particle_struct *pparticle, struct samp_struct *psample, int ipart, int particle_surf, int nfacets, struct cell_struct *pcell, int particle_cell_track[MAXCELLS], int *total_facet_list, int fcount);

void gsi(struct particle_struct *pparticle, int ipart, struct facet_struct *pfacet, int min_fc, double Vw, double pveln[], int GSI_MODEL, double epsilon, double alpha, double alphan, double sigmat);

void maxwell(double pvelf[], double pvelr[], double Vw, double epsilon);

void dria(double pvelf[], double pvelr[], double Vw, int GSI_MODEL, double alpha, double alphan, double sigmat);

void cll(double pvelf[], double pvelr[], double Vw, int GSI_MODEL, double alphan, double sigmat, double alpha);

double compute_Cd(double Umag, double theta, double phi, double Fx, double Fy, double Fz, double proj_area, double n, double m_avg);

double dot(double V[], double W[]);

void cross(double V[], double W[], double VEC[]);

double ranf0(gsl_rng *rr);
