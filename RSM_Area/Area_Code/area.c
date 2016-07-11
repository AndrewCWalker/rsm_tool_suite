/* Calculates the drag coefficient of a given object (from a .stl mesh file) using */
/* the test particle method */

/* INPUTS: */
/* Ux, Uy, Uz = satellite speed relative to atmosphere [m/s] */
/* f = directory path for the mesh file to be used e.g "./Mesh Files/sphere_res05_ascii.stl" */

/* OUTPUT: */
/* A = Projected Area [m^2] */

#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PARALLEL  1          /* 1 = PARALLEL COMPUTATION; 0 = SERIAL */

struct facet_struct {
  double normal[3];
  double vertex1[3];
  double vertex2[3];
  double vertex3[3];
  double area;
};


#include "area.h"
#if PARALLEL
#include "mpi.h"
#endif /* PARALLEL */

#define M_PI 3.141593 /* PI */

int nProcs;

#if PARALLEL
int rank;
MPI_Comm io_comm;
#endif /* PARALLEL */
  
int main(int argc, char *argv[])

{

#if PARALLEL
  /***************************** MPI Initialization *******************************/
  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  /*Create communicator for the I/O functions */
  MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
  /* Find my ID number */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /* Find the total number of procs */
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs); 

#endif /* PARALLEL */

  double A;
  char outfilename[100];
  char meshfilename[100];
  int NUM_POINTS = 0;

  char *temp;
  char *data;
 
#if PARALLEL 
  /************* Set up number of 0 to use before processor number ****************/
  char *zeros="p";
  if (nProcs>9) {
    if (rank<10) zeros="p0";
  }
  if (nProcs>99) {
    if (rank<10) {
      zeros="p00";
    }
    else if (rank<100) {
      zeros="p0";
    }
  }
#endif /*PARALLEL*/

  
  char line[1024];
#if PARALLEL
  sprintf(outfilename, "Aout_%s%d.dat", zeros, rank);
#else /* PARALLEL */
  sprintf(outfilename, "Aout.dat");
#endif /* PARALLEL */

  FILE *fout = fopen(outfilename, "w");
#if PARALLEL
  FILE *ftot = fopen("../Aout.dat", "w");
#endif /* PARALLEL */
  
  int i;

#if PARALLEL
  double ppN;          /* Particles per processor */
  double pc;
#endif /* PARALLEL */

  /* GENERATE ENSEMBLE ARRAYS */
  double *pUx, *pUy, *pUz;

  /* Define new data variables */
  double *LHS_Umag, *LHS_theta, *LHS_phi;

  /* Open ensemble file */
  FILE *lhs_file = fopen("./area.ensemble", "r");

  /* Read the top two ensemble file headers */
  /* Number of Ensemble Members */
  fgets(line, 1024, lhs_file);
  temp = strtok(line, "#");
  data = strtok(NULL, "#");
  sscanf(data, "%d\n", &NUM_POINTS);

  /* Allocate Memory */
  pUx = (double *) calloc(NUM_POINTS, sizeof(double));
  pUy = (double *) calloc(NUM_POINTS, sizeof(double));
  pUz = (double *) calloc(NUM_POINTS, sizeof(double));
  LHS_Umag = (double *) calloc(NUM_POINTS, sizeof(double));
  LHS_theta = (double *) calloc(NUM_POINTS, sizeof(double));
  LHS_phi = (double *) calloc(NUM_POINTS, sizeof(double));

  /* Mesh Filename */
  fgets(line, 1024, lhs_file);
  temp = strtok(line, "#");
  data = strtok(NULL, "#");
  sscanf(data, "%s\n", meshfilename);

  /* Read the third header line */
  fgets(line, 1024, lhs_file);

  /* Read the rest of the data */
  for(i=0; i<NUM_POINTS; i++) {
    fscanf(lhs_file, "%lf %lf %lf\n", &LHS_Umag[i], &LHS_theta[i], &LHS_phi[i]);
  }

  /* Overwrite the ensemble arrays with the LHS data */
  for(i=0; i<NUM_POINTS; i++) {
    pUx[i] = LHS_Umag[i]*cos(LHS_theta[i])*cos(LHS_phi[i]);
    pUy[i] = LHS_Umag[i]*sin(LHS_theta[i])*cos(LHS_phi[i]);
    pUz[i] = LHS_Umag[i]*sin(LHS_phi[i]);
  }

#if PARALLEL
  ppN = (double)NUM_POINTS/(double)nProcs;
#endif /* PARALLEL */


#if PARALLEL
  if(rank==0) {
    printf("Starting simulation\n");
  }

  int start = rank*ppN;
  int end = (rank+1)*ppN;

  for(i=start; i<end; i++) {
    if(rank==0) {
      pc = 100.0*(double)i / ppN;
      printf("%2.1lf percent complete\r", pc);
    }
    A = compute_area(meshfilename, pUx[i], pUy[i], pUz[i]);
    fprintf(fout, "%e %e %e\n", LHS_theta[i], LHS_phi[i], A);
  }
#else /* PARALLEL */

  for(i=0; i<NUM_POINTS; i++) {

    A = compute_area(meshfilename, pUx[i], pUy[i], pUz[i]);
    fprintf(fout, "%e %e %e\n", LHS_theta[i], LHS_phi[i], A);

  }
#endif /* PARALLEL */

  fclose(fout);

#if PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);

  /* Merge Files */
  if(rank==0) {
    fprintf(ftot, "Theta in X-Y plane [radians]     Phi from X-Y plane [radians]      Area [m^2]\n");
    for(i=0; i<nProcs; i++) {
      defzeros(i, &zeros);
      sprintf(outfilename, "Aout_%s%d.dat", zeros, i);
      fout = fopen(outfilename, "r");
      while(!feof(fout)) {
	if(fgets(line, 1024, fout)) {
	  fprintf(ftot, line);
	}
      }
      fclose(fout);
    }
    fclose(ftot);
  }
#endif /* PARALLEL */
   
  free(pUx);
  free(pUy);
  free(pUz);

  free(LHS_Umag);
  free(LHS_theta);
  free(LHS_phi);

  fclose(lhs_file);


#if PARALLEL
  MPI_Finalize();
#endif /* PARALLEL */

#if PARALLEL
  if(rank==0) {
    printf("Simulation complete\n");
  }
#endif /* PARALLEL */

  return(0);

}


/***************************** DEFINE ZEROS FORMATTING ***************************/ 
 void defzeros(int i, char **zeros)
 {

   if (nProcs>9) {
     if (i<10) {
       *zeros="p0";
     } else {
       *zeros = "p";
     }
   }
   if (nProcs>99) {
     if (i<10) {
       *zeros="p00";
     }
     else if (i<100) {
       *zeros="p0";
     }
   }
   
 }


/***************************** TEST PARTICLE MONTE CARLO ***************************/ 
double compute_area(char *meshfilename, double Ux, double Uy, double Uz)

{

  struct facet_struct *pfacet=NULL;

  int nfacets = 0;
  double proj_area = 0.0;

/* INITIALIZE MIN AND MAX POSITIONS FOR EACH DIMENSION */
  double XMIN = 0.0;
  double YMIN = 0.0;
  double ZMIN = 0.0;

  double XMAX = 0.0;
  double YMAX = 0.0;
  double ZMAX = 0.0;

  double Umag = sqrt(Ux*Ux + Uy*Uy + Uz*Uz);

/* DETERMINE NUMBER OF FACETS FROM MESH FILE */
  nfacets = read_num_lines(meshfilename);
  pfacet = (struct facet_struct *) calloc(nfacets, sizeof(struct facet_struct));

  /* READ IN FACET PROPERTIES FROM MESH FILE */
  facet_properties(meshfilename, Ux, Uy, Uz, nfacets, pfacet, Umag);

  /* DYNAMIC COMPUTATION OF DOMAIN BOUNDARIES */
  domain_boundary(&XMIN, &YMIN, &ZMIN, &XMAX, &YMAX, &ZMAX, nfacets, pfacet);

  /* COMPUTE PROJECTED AREA OF MESH */
  proj_area = projected_area(XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX, Ux, Uy, Uz, nfacets, pfacet);


  free(pfacet);

  return(proj_area);

}

  
 /******************************** MESH READ ***************************/
int read_num_lines(char *meshfilename)
{

  /* Read in the STL mesh file and determine number of lines*/
  /* Input: STL mesh filename */
  /* Output: Number of Facets */

  int nfacets = 0;
  int num_lines = 0;
  int ch;
  FILE *f = fopen(meshfilename, "r");
  
  /*Check that file exists*/
  if(!f) {
    printf("Mesh File does not exist\n");
    exit(1);
  }
  
  while (EOF != (ch=fgetc(f))) 
    if (ch=='\n')
      num_lines++;

  fclose(f);

  /* DETERMINE NUMBER OF FACETS */
  nfacets = (num_lines-2)/7;

  return(nfacets);

}


/******************************** FACET PROPERTIES ***************************/
void facet_properties(char *meshfilename, double Ux, double Uy, double Uz, int nfacets, 
		      struct facet_struct *pfacet, double Umag)
{

  struct facet_struct *facet;

  /* Input: STL mesh filename */
  /* Output: Facet Properties Stucture Containing: */
  /*         Facet Normal [x, y, z]  */
  /*         Vertex1 [x, y, z] */
  /*         Vertex2 [x, y, z] */
  /*         Vertex3 [z, y, z] */

  FILE *f = fopen(meshfilename, "r");

  /* READ IN STL FILE HEADER */
  char header[1024];
  char line1[1024], line2[1024], line3[1024], line4[1024], line5[1024], line6[1024], line7[1024];
  char *vert1x, *vert1y, *vert1z;
  char *vert2x, *vert2y, *vert2z;
  char *vert3x, *vert3y, *vert3z;
  char *normx, *normy, *normz;
  char *temp;
  int HeaderSize = 1;
  int i, ifacet;

  double v[3];
  double dist1[3], dist2[3];
  double tc[3];

  for(i=0; i<3; i++) {
    dist1[i] = 0.0;
    dist2[i] = 0.0;
    tc[i] = 0.0;
  }

  for(i=0; i<HeaderSize; i++) {
    fgets(header, 1024, f);
  }

  for(ifacet=0; ifacet<nfacets; ifacet++) {
    /* Read the first line and assign normal components */
    fgets(line1, 1024, f);
    temp = strtok(line1, " ");
    temp = strtok(NULL, " ");
    normx = strtok(NULL, " ");
    normy = strtok(NULL, " ");
    normz = strtok(NULL, " ");
    /* Throw away second line */
    fgets(line2, 1024, f);
    /* Read the third line and assign vertex #1 position */
    fgets(line3, 1024, f);
    temp = strtok(line3, " ");
    vert1x = strtok(NULL, " ");
    vert1y = strtok(NULL, " ");
    vert1z = strtok(NULL, " ");
    /* Read the fourth line and assign vertex #2 position */
    fgets(line4, 1024, f);
    temp = strtok(line4, " ");
    vert2x = strtok(NULL, " ");
    vert2y = strtok(NULL, " ");
    vert2z = strtok(NULL, " ");
    /* Read the fifth line and assign vertex #3 position */
    fgets(line5, 1024, f);
    temp = strtok(line5, " ");
    vert3x = strtok(NULL, " ");
    vert3y = strtok(NULL, " ");
    vert3z = strtok(NULL, " ");
    /* Throw away the sixth and seventh lines */
    fgets(line6, 1024, f);
    fgets(line7, 1024, f);

    facet = pfacet + ifacet;

    /* Assign facet normal vector to "facet" structure */
    facet->normal[0] = atof(normx);
    facet->normal[1] = atof(normy);
    facet->normal[2] = atof(normz);

    /* Assign vertex #1 position to "facet" structure */
    facet->vertex1[0] = atof(vert1x);
    facet->vertex1[1] = atof(vert1y);
    facet->vertex1[2] = atof(vert1z);

    /* Assign vertex #2 position to "facet" structure */
    facet->vertex2[0] = atof(vert2x);
    facet->vertex2[1] = atof(vert2y);
    facet->vertex2[2] = atof(vert2z);

    /* Assign vertex #3 position to "facet" structure */
    facet->vertex3[0] = atof(vert3x);
    facet->vertex3[1] = atof(vert3y);
    facet->vertex3[2] = atof(vert3z);

    /* COMPUTE UNIT VELOCITY VECTOR */
    v[0] = Ux/Umag;
    v[1] = Uy/Umag;
    v[2] = Uz/Umag;

    /* Find length of 1st side of the triangle */
    dist1[0] = facet->vertex2[0] - facet->vertex1[0];
    dist1[1] = facet->vertex2[1] - facet->vertex1[1];
    dist1[2] = facet->vertex2[2] - facet->vertex1[2];
    
    /* Find length of 2st side of the triangle */
    dist2[0] = facet->vertex1[0] - facet->vertex3[0];
    dist2[1] = facet->vertex1[1] - facet->vertex3[1];
    dist2[2] = facet->vertex1[2] - facet->vertex3[2];

    /* Compute the cross product between the two sides of the triangle */
    cross(dist1, dist2, tc);

    facet->area = 0.5*sqrt(dot(tc, tc));
     
  }

  fclose(f);
}


/***************************** PROJECTED AREA **********************************/
double projected_area(double XMIN, double YMIN, double ZMIN, double XMAX, double YMAX, double ZMAX, double Ux, double Uy, double Uz, int nfacets, struct facet_struct *pfacet)
{

  /* Calculate the projected area of the mesh */
  /* Translated from F. Lumpkin's Fortran Area routine for DAC */
  /* Inputs: facet_normal = Facet Normal [x, y, z] */
  /*         vertex1 = x, y, z positions of 1st vertex */
  /*         vertex2 = x, y, z positions of 2nd vertex */
  /*         vertex3 = x, y, z positions of 3rd vertex */
            
  /* Outputs: facet_area = Facet Area [m^2] */
  /*          proj_area = Projected Area of the Mesh [m^2] */ 

  struct facet_struct *facet;
  double proj_area = 0.0;
  double v[3];
  double Umag = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);
  double theta, phi;
  double **fv1, **fv2, **fv3;
  double transform[3][3];
  double sp, cp, st, ct;

  int ifacet, i, j, k, n;
  
  int idir0 = 0;
  int msort = 64;
  int mgrid = 4096;
  int melms = 2000000;
  int melm2 = 2*melms;
  int msort2 = msort*msort;
  int icount;
  int j1, j2, j3, k1, k2, k3;
  int maxj, minj, maxk, mink;
  int jj;
  int iarea;
  int ix[3], iy[3];
  int *isort, *j1sort, *j2sort;

  double a = 2.0e-7;
  double eps;
  double xmax[3] = {0.0, 0.0, 0.0};
  double xmin[3] = {0.0, 0.0, 0.0};
  double fact;
  double dx[3], dy[3];
  double dxi[3], dyi[3];
  double **xp, **yp;
  double *x1, *x2, *x3, *y1, *y2, *y3;
  double *dx21, *dx32, *dx13, *dy21, *dy32, *dy13;
  double cx1, cx2, cx3;
  double area, diffarea;
  double ddiv, ddx, ddy, ddxy;

  /* ADD SMALL VALUE TO BOUNDARIES */
  /* NOTE - SHOULD REMOVE 1.1. IN FAVOR OF A CONSTANT VARIABLE DEFINED IN MAIN */
  /* CURRENTLY 1.1. COMES FROM THE VALUE SET IN DOMAIN_BOUNDARY */
  double fxmin = 1.0e6;
  double fxmax = -1.0e6;
  double fymin = 1.0e6;
  double fymax = -1.0e6;
  double fzmin = 1.0e6;
  double fzmax = -1.0e6;

  /* ALLOCATE MEMORY */
  xp = (double **) calloc(3, sizeof(double));
  yp = (double **) calloc(3, sizeof(double));

  for(i=0; i<3; i++) {
    xp[i] = (double *) calloc(mgrid, sizeof(double));
    yp[i] = (double *) calloc(mgrid, sizeof(double));
  }

  isort = (int *) calloc(msort2, sizeof(int));
  j1sort = (int *) calloc(melm2, sizeof(int));
  j2sort = (int *) calloc(melm2, sizeof(int));

  x1 = (double *) calloc(nfacets, sizeof(double));
  x2 = (double *) calloc(nfacets, sizeof(double));
  x3 = (double *) calloc(nfacets, sizeof(double));
  y1 = (double *) calloc(nfacets, sizeof(double));
  y2 = (double *) calloc(nfacets, sizeof(double));
  y3 = (double *) calloc(nfacets, sizeof(double));

  dx21 = (double *) calloc(nfacets, sizeof(double));
  dx32 = (double *) calloc(nfacets, sizeof(double));
  dx13 = (double *) calloc(nfacets, sizeof(double));
  dy21 = (double *) calloc(nfacets, sizeof(double));
  dy32 = (double *) calloc(nfacets, sizeof(double));
  dy13 = (double *) calloc(nfacets, sizeof(double));

  fv1 = (double **) calloc(nfacets, sizeof(double));
  fv2 = (double **) calloc(nfacets, sizeof(double));
  fv3 = (double **) calloc(nfacets, sizeof(double));
  for(ifacet=0; ifacet<nfacets; ifacet++) {
    fv1[ifacet] = (double *) calloc(3, sizeof(double));
    fv2[ifacet] = (double *) calloc(3, sizeof(double));
    fv3[ifacet] = (double *) calloc(3, sizeof(double));
  }

  /* Define orientation angles in radians*/
  if(Ux==0.0 && Uy==0.0) {
    theta = 0.0;
  } else {
    theta = atan2(Uy, Ux);
  }
  phi = asin(Uz/sqrt(Ux*Ux+Uy*Uy+Uz*Uz));

  /* Transform vertices based on theta and phi */
  if(theta != 0.0 || phi != 0.0) {
    sp = sin(phi);
    cp = cos(phi);
    st = sin(theta);
    ct = cos(theta);
    transform[0][0] = ct*cp;
    transform[0][1] = st*cp;
    transform[0][2] = -sp;
    transform[1][0] = -st;
    transform[1][1] = ct;
    transform[1][2] = 0.0;
    transform[2][0] = ct*sp;
    transform[2][1] = st*sp;
    transform[2][2] = cp;
    for(ifacet=0; ifacet<nfacets; ifacet++) {
      facet = pfacet + ifacet;
      for(i=0; i<3; i++) {
  	fv1[ifacet][i] = 0.0;
  	fv2[ifacet][i] = 0.0;
  	fv3[ifacet][i] = 0.0;
  	for(j=0; j<3; j++) {
  	  fv1[ifacet][i] += transform[j][i]*facet->vertex1[j];
  	  fv2[ifacet][i] += transform[j][i]*facet->vertex2[j];
  	  fv3[ifacet][i] += transform[j][i]*facet->vertex3[j];
  	}
      }
    }
  } else {
    for(ifacet=0; ifacet<nfacets; ifacet++) {
      facet = pfacet + ifacet;
      for(i=0; i<3; i++) {
  	fv1[ifacet][i] = facet->vertex1[i];
  	fv2[ifacet][i] = facet->vertex2[i];
  	fv3[ifacet][i] = facet->vertex3[i];
      }
    }
  }

  for(ifacet=0; ifacet<nfacets; ifacet++) {

    /* Determine X mininimum and maximum for this facet */
    if(fv1[ifacet][0] < fxmin) fxmin = fv1[ifacet][0];
    if(fv2[ifacet][0] < fxmin) fxmin = fv2[ifacet][0];
    if(fv3[ifacet][0] < fxmin) fxmin = fv3[ifacet][0];

    if(fv1[ifacet][0] > fxmax) fxmax = fv1[ifacet][0];
    if(fv2[ifacet][0] > fxmax) fxmax = fv2[ifacet][0];
    if(fv3[ifacet][0] > fxmax) fxmax = fv3[ifacet][0];

    /* Determine Y mininimum and maximum for this facet */
    if(fv1[ifacet][1] < fymin) fymin = fv1[ifacet][1];
    if(fv2[ifacet][1] < fymin) fymin = fv2[ifacet][1];
    if(fv3[ifacet][1] < fymin) fymin = fv3[ifacet][1];

    if(fv1[ifacet][1] > fymax) fymax = fv1[ifacet][1];
    if(fv2[ifacet][1] > fymax) fymax = fv2[ifacet][1];
    if(fv3[ifacet][1] > fymax) fymax = fv3[ifacet][1];

    /* Determine Z mininimum and maximum for this facet */
    if(fv1[ifacet][2] < fzmin) fzmin = fv1[ifacet][2];
    if(fv2[ifacet][2] < fzmin) fzmin = fv2[ifacet][2];
    if(fv3[ifacet][2] < fzmin) fzmin = fv3[ifacet][2];

    if(fv1[ifacet][2] > fzmax) fzmax = fv1[ifacet][2];
    if(fv2[ifacet][2] > fzmax) fzmax = fv2[ifacet][2];
    if(fv3[ifacet][2] > fzmax) fzmax = fv3[ifacet][2];

  }

  eps = a*((xmax[0] + xmax[1] + xmax[2]) - (xmin[0] + xmin[1] + xmin[2]));

  xmin[0] = fxmin - eps;
  xmax[0] = fxmax + eps;
  xmin[1] = fymin - eps;
  xmax[1] = fymax + eps;
  xmin[2] = fzmin - eps;
  xmax[2] = fzmax + eps;

  /* CALCULATE DOMAIN SIZE AND INVERSE DOMAIN GRID SIZE */
  for(i=0; i<3; i++) {
    ix[i] = (int)((i+1) % 3);
    iy[i] = (int)(((i+2) % 3));
    dx[i] = xmax[ix[i]] - xmin[ix[i]];
    dy[i] = xmax[iy[i]] - xmin[iy[i]];
    dxi[i] = (double)msort/dx[i];
    dyi[i] = (double)msort/dy[i];
  }

  /* CALCULATE GRID POINT CENTERS */
  fact = 1.0/(double)mgrid;
  for(i=0; i<3; i++) {
    for(n=0; n<mgrid; n++) {
      xp[i][n] = fact*dx[i]*(double)(n+0.5) + xmin[ix[i]];
      yp[i][n] = fact*dy[i]*(double)(n+0.5) + xmin[iy[i]];
    }
  }

  /* INITIALIZE SORTING ARRAYS */
  icount = 0;
  for(i=0; i<msort2; i++) {
    isort[i] = 0;
  }
  for(i=0; i<melm2; i++) {
    j1sort[i] = 0;
    j2sort[i] = 0;
  }
  
  for(ifacet=0; ifacet<nfacets; ifacet++) { 

    facet = pfacet + ifacet;

    x1[ifacet] = fv1[ifacet][ix[idir0]];
    x2[ifacet] = fv2[ifacet][ix[idir0]];
    x3[ifacet] = fv3[ifacet][ix[idir0]];
    y1[ifacet] = fv1[ifacet][iy[idir0]];
    y2[ifacet] = fv2[ifacet][iy[idir0]];
    y3[ifacet] = fv3[ifacet][iy[idir0]];

    dx21[ifacet] = x2[ifacet] - x1[ifacet];
    dx32[ifacet] = x3[ifacet] - x2[ifacet];
    dx13[ifacet] = x1[ifacet] - x3[ifacet];
    dy21[ifacet] = y2[ifacet] - y1[ifacet];
    dy32[ifacet] = y3[ifacet] - y2[ifacet];
    dy13[ifacet] = y1[ifacet] - y3[ifacet];
    
    j1 = (int)((x1[ifacet] - xmin[ix[idir0]])*dxi[idir0]);
    j2 = (int)((x2[ifacet] - xmin[ix[idir0]])*dxi[idir0]);
    j3 = (int)((x3[ifacet] - xmin[ix[idir0]])*dxi[idir0]);

    k1 = (int)((y1[ifacet] - xmin[iy[idir0]])*dyi[idir0]);
    k2 = (int)((y2[ifacet] - xmin[iy[idir0]])*dyi[idir0]);
    k3 = (int)((y3[ifacet] - xmin[iy[idir0]])*dyi[idir0]);

    
  
    maxj = (int)fmax((double)j1, (double)j2);
    maxj = (int)fmax((double)j3, (double)maxj);

    minj = (int)fmin((double)j1, (double)j2);
    minj = (int)fmin((double)j3, (double)minj);

    maxk = (int)fmax((double)k1, (double)k2);
    maxk = (int)fmax((double)k3, (double)maxk);

    mink = (int)fmin((double)k1, (double)k2);
    mink = (int)fmin((double)k3, (double)mink);

    //TEMPORARY FIX
    if(maxj > 63) {
      maxj = 63;
    }
    if(maxk > 63) {
      maxk = 63;
    }
    if(minj < 0) {
      minj = 0;
    }
    if(mink < 0) {
      mink = 0;
    }
   
    for(j=minj; j<=maxj; j++) {
      for(k=mink; k<=maxk; k++) {
	icount = icount + 1;
	jj = msort*(k) + j;
	j2sort[icount] = ifacet;
	j1sort[icount] = isort[jj];
	isort[jj] = icount;
      }
    }

  }
   
  area = 0.0;
  for(i=0; i<mgrid; i++) {
    iarea = 0;
    for(j=0; j<mgrid; j++) {
      jj = msort*(int)((yp[idir0][j]-xmin[iy[idir0]])*dyi[idir0]) + (int)((xp[idir0][i]-xmin[ix[idir0]])*dxi[idir0]);
      k = isort[jj];
      while(k!=0) {
	ifacet = j2sort[k];
	cx1 = dx21[ifacet]*(yp[idir0][j]-y1[ifacet]) - dy21[ifacet]*(xp[idir0][i]-x1[ifacet]);
	cx2 = dx32[ifacet]*(yp[idir0][j]-y2[ifacet]) - dy32[ifacet]*(xp[idir0][i]-x2[ifacet]);
	if(cx1*cx2 > 0.0) {
	  cx3 = dx13[ifacet]*(yp[idir0][j]-y3[ifacet]) - dy13[ifacet]*(xp[idir0][i]-x3[ifacet]);
	  if(cx1*cx3 > 0.0) {
	    iarea++;
	    break;
	  }
	}
	k = j1sort[k];
      }
    }
    diffarea = (double)iarea;
    area += diffarea;
  }
  ddiv = (double)mgrid;
  ddx = (double)(dx[idir0])/ddiv;
  ddy = (double)(dy[idir0])/ddiv;
  ddxy = ddx*ddy;
  proj_area = ddxy*area;

  for(i=0; i<3; i++) {
    free(xp[i]);
    free(yp[i]);
  }

  free(xp);
  free(yp);
  
  free(isort);
  free(j1sort);
  free(j2sort);

  free(x1);
  free(x2);
  free(x3);
  free(y1);
  free(y2);
  free(y3);

  free(dx21);
  free(dx32);
  free(dx13);
  free(dy21);
  free(dy32);
  free(dy13);

  for(ifacet=0; ifacet<nfacets; ifacet++) {
    free(fv1[ifacet]);
    free(fv2[ifacet]);
    free(fv3[ifacet]);
  }

  free(fv1);
  free(fv2);
  free(fv3);

  return(proj_area);

}


/***************************** DOMAIN BOUNDARY **********************************/
void domain_boundary(double *XMIN, double *YMIN, double *ZMIN, double *XMAX, double *YMAX, double *ZMAX, int nfacets, struct facet_struct *pfacet)
{

  /* Calculate the domain boundaries based on mesh */
  /* Inputs: nfacets = Number of Facets */
  /*         pfacet = Facet structure pointer          */
            
  /* Outputs: = Xmin = X-direction domain minimum boundary [m] */
  /*            Xmax = X-direction domain maximum boundary [m] */
  /*            Ymin = Y-direction domain minimum boundary [m] */
  /*            Ymax = Y-direction domain maximum boundary [m] */
  /*            Zmin = Z-direction domain minimum boundary [m] */
  /*            Zmax = Z-direction domain maximum boundary [m] */

  struct facet_struct *facet;
  int ifacet;
  double fxmin = 1.0e6;
  double fxmax = -1.0e6;
  double fymin = 1.0e6;
  double fymax = -1.0e6;
  double fzmin = 1.0e6;
  double fzmax = -1.0e6;

  double A = 1.1; /* Safety factor on domain bounds */

  for(ifacet=0; ifacet<nfacets; ifacet++) {
    
    facet = pfacet + ifacet;

    /* Determine X mininimum and maximum for this facet */
    if(facet->vertex1[0] < fxmin) fxmin = facet->vertex1[0];
    if(facet->vertex2[0] < fxmin) fxmin = facet->vertex2[0];
    if(facet->vertex3[0] < fxmin) fxmin = facet->vertex3[0];

    if(facet->vertex1[0] > fxmax) fxmax = facet->vertex1[0];
    if(facet->vertex2[0] > fxmax) fxmax = facet->vertex2[0];
    if(facet->vertex3[0] > fxmax) fxmax = facet->vertex3[0];

    /* Determine Y mininimum and maximum for this facet */
    if(facet->vertex1[1] < fymin) fymin = facet->vertex1[1];
    if(facet->vertex2[1] < fymin) fymin = facet->vertex2[1];
    if(facet->vertex3[1] < fymin) fymin = facet->vertex3[1];

    if(facet->vertex1[1] > fymax) fymax = facet->vertex1[1];
    if(facet->vertex2[1] > fymax) fymax = facet->vertex2[1];
    if(facet->vertex3[1] > fymax) fymax = facet->vertex3[1];

    /* Determine Z mininimum and maximum for this facet */
    if(facet->vertex1[2] < fzmin) fzmin = facet->vertex1[2];
    if(facet->vertex2[2] < fzmin) fzmin = facet->vertex2[2];
    if(facet->vertex3[2] < fzmin) fzmin = facet->vertex3[2];

    if(facet->vertex1[2] > fzmax) fzmax = facet->vertex1[2];
    if(facet->vertex2[2] > fzmax) fzmax = facet->vertex2[2];
    if(facet->vertex3[2] > fzmax) fzmax = facet->vertex3[2];

  }

  *XMIN = A*fxmin;
  *XMAX = A*fxmax;
  *YMIN = A*fymin;
  *YMAX = A*fymax;
  *ZMIN = A*fzmin;
  *ZMAX = A*fzmax;

}


/*********************************** DOT PRODUCT *************************************/
double dot(double V[], double W[]) {
   
  /* Computes dot product between vectors V and W */
  /* Inputs: V = Vector #1 [Vx, Vy, Vz] */
  /*         W = Vector #2 [Wx, Wy, Wz] */
  /* Output: a = Dot Product of V and W */

  double a = 0.0;

  a = V[0]*W[0] + V[1]*W[1] + V[2]*W[2];
  return(a);

}


/*********************************** DOT PRODUCT *************************************/
void cross(double V[], double W[], double VEC[]) {
  /* Computes cross product between vectors V and W */

  /* Inputs: V = Vector #1 [Vx, Vy, Vz] */
  /*         W = Vector #2 [Wx, Wy, Wz] */

  /* Output: b = Cross Product of V and W */

  int i;

  VEC[0] = V[1]*W[2] - V[2]*W[1];
  VEC[1] = V[2]*W[0] - V[0]*W[2];
  VEC[2] = V[0]*W[1] - V[1]*W[0];

}


 
