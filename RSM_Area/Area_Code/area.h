void defzeros(int i, char **zeros);

double compute_area(char *meshfilename, double Ux, double Uy, double Uz);

int read_num_lines(char *meshfilename);

void facet_properties(char *meshfilename, double Ux, double Uy, double Uz, int nfacets, 
		      struct facet_struct *pfacet, double Umag);

void domain_boundary(double *XMIN, double *YMIN, double *ZMIN, double *XMAX, double *YMAX, double *ZMAX, int nfacets, struct facet_struct *pfacet);

double projected_area(double XMIN, double YMIN, double ZMIN, double XMAX, double YMAX, double ZMAX, double Ux, double Uy, double Uz, int nfacets, struct facet_struct *pfacet);

double dot(double V[], double W[]);

void cross(double V[], double W[], double VEC[]);
