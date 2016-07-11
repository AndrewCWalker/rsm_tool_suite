#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <time.h>

#define NSPECIES 6
#define M_PI 3.141593
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

struct initCd_struct {
  double X[NSPECIES];
  double n;
  double U;
  double Ta;
  double Ts;
  int ads_model;
  double pitch;
  double yaw;
  double surf_mass;
  double sat_mass;
  char rsmfilename[256];
  char areafilename[256];
};

struct RSM_struct {
  double beta[6][7];
  double xmin[6][7];
  double xrange[6][7];
  double mean[6];
  double sd[6];
  double lamz[6];
  double lamws[6];
  double x[6][1000][7];
  double w[6][1000];
  double **area;
};

struct area_struct {
  int NPITCH;
  int NYAW;
  double MinPitch;
  double MaxPitch;
  double MinYaw;
  double MaxYaw;
  double **area;
};

#include "RSM.h"

int main(int argc, char *argv[]) {

  int i;

  /* Initialize variables */
  double RSMCd = 0.0;     /* Satellite Drag Coefficient from RSM */
  double RSMA = 0.0;      /* Satellite Area from Area Look-up Table */
  double RSMBC = 0.0;     /* Ballistic Coefficient from RSM & Look-up Table */

  char *adsstr;           /* String for Adsorption Model Name */

  /* Initiliaze and allocate Drag Coefficient structure */
  struct initCd_struct *initCd = NULL;
  initCd = (struct initCd_struct *) calloc(1, sizeof(struct initCd_struct));

  /* Initialize and allocated Response Surface Model structure */
  struct RSM_struct *RSMdata = NULL;
  RSMdata = (struct RSM_struct *) calloc(1, sizeof(struct RSM_struct));

  /* Initialize and allocated look-up table area structure */
  struct area_struct *AreaData = NULL;
  AreaData = (struct area_struct *) calloc(1, sizeof(struct area_struct));

  /* Read the input file */
  read_input_file(initCd);

  /* Read the RSM file */
  read_RSM(initCd, RSMdata);

  /* Read the area look-up table */
  read_area(initCd, AreaData);

  /* Calculate the drag coefficient from the RSM */
  RSMCd = CdSetup(initCd, RSMdata);

  /* Calculate the area from the look-up table */
  RSMA = AreaLookUp(AreaData, initCd);

  /* Calculate the ballistic coefficient */
  RSMBC = RSMCd*(RSMA/initCd->sat_mass);

  /* Write the outputs to the screen */
  WriteOutput(initCd, RSMCd, RSMA, RSMBC, adsstr);

  /* Free memory */
  for(i=0; i<AreaData->NPITCH; i++) {
    free(AreaData->area[i]);
  }
  free(AreaData->area);

  free(initCd);
  free(RSMdata);
  free(AreaData);

}


double CdSetup(struct initCd_struct *initCd, struct RSM_struct *RSMdata) {


  /* Setup all the input properties to calculate the drag coefficient */

  /****** INPUTS: ******/
  /* initCd structure containing all independent Cd variables */
  /* initCd->X[6] is a list of the species mole fractions in the order [H, He, N, N2, O, O2] */
  /* initCd->n is the total number density in m^-3 */
  /* initCd->U is the satellite velocity relative to the atmosphere in km/s */
  /* initCd->Ta is the atmospheric translational temperature in K */
  /* initCd->Ts is the satellite surface temperature in K */
  /* initCd->ads_model is the flag to determine the adsorption model (0 = Langmuir; 1 = Freundlich) */
  /* initCd->pitch is the pitch angle (phi measured from the x-y plane) in radians */
  /* initCd->yaw is the yaw angle (theta measured in the x-y plane from the x-axis) in radians */
  /* initCd->surf_mass is the mass of a particle that composes the surface lattice in kg */

  /****** OUTPUTS ******/
  /* Cd_TOTAL is the satellite drag coefficient based on the atmospheric and satellite properties */

  int i;
  double yaw; /* Yaw angle - Defined due to RSM symmetry */

  double ads_input[7], srf_input[7];   /* Input parameters to the RSM emulator [U, Ts, Ta, alphan, sigmat, yaw, pitch] */
  double m_avg = 0.0;                  /* Average particle mass of the atmosphere [kg] */

  double m_surf = initCd->surf_mass;   /* Define particle mass of satellite surface [kg] */

  double sigmat = 1.0;  /* Tangential momentum accomodation of 1.0 based on Comsa (1980)/Porodnov (1974) [unitless] */
  double kB = 1.38e-23; /* Boltzmann's constant [J/K] */

  double kL_CLL = 2.89e6;       /* CLL - Define Langmuir Adsorpate Constant [unitless] */
  double aF_CLL = 0.089;        /* CLL - Define Freundlich Alpha Constant [unitless] */
  double kF_CLL = 2.35;         /* CLL - Define Freundlich K constant [unitless] */

  /* Alpha for atomic oxygen covered surface is unity */
  double alpha_ads = 1.0;       /* Define energy accommodation coefficient (DRIA) */
  double alphan_ads = 1.0;      /* Define the normal energy accommodation coefficient (CLL) */

  double mu;                    /* Ratio of atmospheric particle mass to satellite surface particle mass [unitless] */
  double fsc = 0.0;             /* Fraction of Surface Covered by atomic oxygen [unitless]*/
  double alpha_sub, alphan_sub; /* Energy accommodation coefficients for clear surface [unitless] */
  double Cd_ads, Cd_srf;        /* Drag coefficients for adsorbate and clean satellite surfaces [unitless] */
  double Po;                    /* Partial pressure of atomic oxygen [Pa] */
  double Cd_TOTAL;              /* Total Satellite Drag Coefficient [unitless] */

  /* Define atomic/molecular masses of each species */
  double mass[NSPECIES];               /* The masses of the six species [kg] */
  mass[0] = 1.674e-27;  /* hydrogen */
  mass[1] = 6.646e-27;  /* helium */
  mass[2] = 2.326e-26;  /* atomic nitrogen */
  mass[3] = 4.652e-26;  /* diatomic nitrogen */
  mass[4] = 2.657e-26;  /* atomic oxygen */
  mass[5] = 5.314e-26;  /* diatomic oxygen */
  
  /* Compute the average mass */
  for(i=0; i<NSPECIES; i++) {
    m_avg += initCd->X[i]*mass[i];
  }
    
  /* alpha for clean surface is based on Goodman's (1966) empirical formula */
  mu = m_avg/initCd->surf_mass;
  alpha_sub = 2.4*mu/pow(1+mu, 2.0);
  alphan_sub = 2.0*alpha_sub-1.0;

  /* Can't allow alphan below 0.0 */
  if(alphan_sub < 0.0) {
    alphan_sub = 0.0;
  }

  /* Define response surface emulator inputs */
  ads_input[0] = srf_input[0] = initCd->U*1.0e3; /* Convert speed to m/s */
  ads_input[1] = srf_input[1] = initCd->Ts;
  ads_input[2] = srf_input[2] = initCd->Ta;
  ads_input[4] = srf_input[4] = sigmat;
  ads_input[5] = srf_input[5] = fabs(initCd->yaw);
  ads_input[6] = srf_input[6] = initCd->pitch;
  
  /* Specific inputs for adsorbate and clean surfaces for accommodation coefficients */
  ads_input[3] = alphan_ads;
  srf_input[3] = alphan_sub;
  
  /* Compute surface-specific drag coefficients from RSM */
  emu(RSMdata, ads_input, initCd->X, &Cd_ads);
  emu(RSMdata, srf_input, initCd->X, &Cd_srf);

  /* Compute the partial pressure of atomic oxygen */
  Po = initCd->n*initCd->X[4]*kB*initCd->Ta;

  /* Compute fsc: fractional coverage of surface by atomic oxygen */
  if(initCd->ads_model) {                                   /* Freundlich */
    fsc = kF_CLL*pow(Po, aF_CLL);
    if(fsc > 1.0) {
      fsc = 1.0;
    }
  } else {                                          /* Langmuir */
    fsc = (kL_CLL*Po)/(1+kL_CLL*Po);
  }

  /* Compute satellite drag coefficient based on weighted sum of surface-specific Cd */
  Cd_TOTAL = fsc*Cd_ads + (1.0-fsc)*Cd_srf;

  return(Cd_TOTAL);

}



// Species order for emu will always be H, He, N, N2, O, O2

// Number of simulations, number of inputs, number of chemical species.
static int m=1000, p=7, nspec=6;
// Kriging basis computed by emuInit, sizes should be nspec and m
static double KrigBasis[7][1000];


// Initialization function that computes the Kriging basis
void emuInit(struct RSM_struct *RSMdata) {
    int i,j,k,l;
    double cov;
    gsl_matrix *SigmaSim = gsl_matrix_alloc(m,m);
    gsl_vector *b = gsl_vector_alloc(m);
    
    // Do these one principal component at a time
    for(i=0; i<nspec; i++) {
		
        // Fill in the covariance matrix for the principals components
        // Also make a gsl_vector with the weights.
        for(j=0; j<m; j++) {
            // Diagonal
            gsl_matrix_set(SigmaSim, j, j, (1.0/RSMdata->lamz[i]) + (1.0/RSMdata->lamws[i]));
            // Off-diagonals
            for(k=0; k<j; k++) {
                // Compute the covariance
                cov = 0.0;
                for(l=0; l<p; l++) {
                    cov -= RSMdata->beta[i][l]*pow(RSMdata->x[i][j][l]-RSMdata->x[i][k][l], 2.0);
                }
                cov = exp(cov)/RSMdata->lamz[i];
                gsl_matrix_set(SigmaSim, j, k, cov);
                gsl_matrix_set(SigmaSim, k, j, cov);
            } // for(k=0; k<j; k++)
            gsl_vector_set(b, j, RSMdata->w[i][j]);
        } // for(j=0; j<m; j++)
        
        // Cholesky and solve
        gsl_linalg_cholesky_decomp(SigmaSim);
        gsl_linalg_cholesky_svx(SigmaSim, b);
        
        // Copy into the Kriging Basis
        for(j=0; j<m; j++) {
            KrigBasis[i][j] = gsl_vector_get(b, j);
        }
    } // for(i=0; i<nspec; i++)
    gsl_matrix_free(SigmaSim);
    gsl_vector_free(b);
    
}

// The actual emulation
// input parameters, species weights, output
void emu(struct RSM_struct *RSMdata, double *xstar, double *specw, double *ystar) {
    static int inited=0;
    int i, j, k;
    double xstarstd[nspec][p], wstar[nspec], Sigmastar[nspec][m], logc;

    double mass[nspec];
    double mavg = 0.0;
    /* Define atomic/molecular masses of each species */
    mass[0] = 1.674e-27;  /* hydrogen */
    mass[1] = 6.646e-27;  /* helium */
    mass[2] = 2.326e-26;  /* atomic nitrogen */
    mass[3] = 4.652e-26;  /* diatomic nitrogen */
    mass[4] = 2.657e-26;  /* atomic oxygen */
    mass[5] = 5.314e-26;  /* diatomic oxygen */

    /* Compute average mass */
    /* Compute the average mass */
    for(i=0; i<nspec; i++) {
      mavg += specw[i]*mass[i];
    }
    		
    // Iinitialize if necessary
    if(inited==0) {
        emuInit(RSMdata);
        inited=1;
    }
    
    // Standardize the inputs
    // The different species require slightly different standardizations.
    for(i=0; i<nspec; i++) {
        for(j=0; j<p; j++) {
            xstarstd[i][j] = (xstar[j] - RSMdata->xmin[i][j]) / RSMdata->xrange[i][j];
        }
    }
    
    // Compute the covariances between the new input and sims for all PCs
    for(i=0; i<nspec; i++) {
        for(j=0; j<m; j++) {
            logc = 0.0;
            for(k=0; k<p; k++) {
                logc -= RSMdata->beta[i][k]*pow(RSMdata->x[i][j][k]-xstarstd[i][k], 2.0);
            }
            Sigmastar[i][j] = exp(logc)/RSMdata->lamz[i];
        }
    }
    
    // Compute wstar, the predicted species results for the new input
    for(i=0; i<nspec; i++) {
        wstar[i]=0.0;
        for(j=0; j<m; j++) {
            wstar[i] += Sigmastar[i][j] * KrigBasis[i][j];
        }
		wstar[i] = wstar[i]*RSMdata->sd[i] + RSMdata->mean[i];
    }
    
    // Compute ystar, the new output
    *ystar = 0.0;
    for(j=0; j<nspec; j++) {
      *ystar += specw[j]*wstar[j]*mass[j]/mavg;
    }

}


/******************************** READ INPUT FILE ********************************/
double AreaLookUp(struct area_struct *AreaData, struct initCd_struct *initCd) {

  /* Find nearest pitch and yaw and linearly interpolate projected area */
  /* Input: RSM areas */
  /* Output: Projected area */

  double proj_area;
  int i, j;
  int ipitch, iyaw;
  double yaw;

  /* Define orientation bounds */
  double MinPitch = 0.0;
  double MaxPitch = M_PI/2.0;
  double MinYaw = 0.0; /* Note by symmetry, this extends to -5.0 */
  double MaxYaw = 5.0*M_PI/180.0;

  /* Define dummy pitch and yaw arrays */
  double *pitch_array = calloc(AreaData->NPITCH, sizeof(double));
  double *yaw_array = calloc(AreaData->NYAW, sizeof(double));

  /* Define area array for bilinear interpolation */
  double area_array[4];

  /* Populate dummy arrays for pitch and yaw */
  for(i=0; i<AreaData->NPITCH; i++) {
    pitch_array[i] = (AreaData->MaxPitch - AreaData->MinPitch)*i/(AreaData->NPITCH-1);
  }
  
  for(i=0; i<AreaData->NYAW; i++) {
    yaw_array[i] = (AreaData->MaxYaw - AreaData->MinYaw)*i/(AreaData->NYAW-1);
  }

  /* Take absolute value of yaw (by symmetry) */
  yaw = fabs(initCd->yaw);

  /* Find narest pitch and yaw indices */
  ipitch = (int) (((initCd->pitch - AreaData->MinPitch)/(AreaData->MaxPitch - AreaData->MinPitch))
		  *((double)AreaData->NPITCH-1.0));
  iyaw = (int) (((yaw - AreaData->MinYaw)/(AreaData->MaxYaw - AreaData->MinYaw))*((double)AreaData->NYAW-1.0));

  area_array[0] = AreaData->area[ipitch][iyaw];
  area_array[1] = AreaData->area[ipitch][iyaw+1];
  area_array[2] = AreaData->area[ipitch+1][iyaw];
  area_array[3] = AreaData->area[ipitch+1][iyaw+1];

  /* BiLinear interpolation */
  proj_area = BilinearInterp(pitch_array, yaw_array, initCd->pitch, yaw, area_array, ipitch, iyaw);

  free(pitch_array);
  free(yaw_array);
  
  return(proj_area);
      
}


/******************************** READ INPUT FILE ********************************/
void read_input_file(struct initCd_struct *initCd) {

  /* Read in the RSM input file and assign variables */
  /* Input: RSM.inp */
  /* Output: Drag coefficient input structure, *initCd */

  char *filepath = "./RSM.inp";

  FILE *f = fopen(filepath, "r");
  if (f == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }

  char line[1024];
  char *temp;
  char *data;
  int i;

  double Xsum = 0.0;

  /* READ IN HEADER AND EMPTY LINES */
  for(i=0; i<4; i++) {
    fgets(line, 1024, f);
  }

  /* READ IN DRAG COEFFICIENT CONTROLS */
  for(i=0; i<12; i++) {
    fgets(line, 1024, f);
    temp = strtok(line, "#");
    data = strtok(NULL, "#");

    if(i==0) /* SPECIES MOLE FRACTIONS */
      sscanf(data, "%lf %lf %lf %lf %lf %lf\n", &initCd->X[0], &initCd->X[1], &initCd->X[2], 
	                                        &initCd->X[3], &initCd->X[4], &initCd->X[5]);
    if(i==1) /* SATELLITE GEOMETRY */
      sscanf(data, "%lf\n", &initCd->n);
    if(i==2) /* SATELLITE VELOCITY RELATIVE TO ATMOSPHERE */
      sscanf(data, "%lf\n", &initCd->U);
    if(i==3) /* ATMOSPHERIC TRANSLATIONAL TEMPERATURE */
      sscanf(data, "%lf\n", &initCd->Ta);
    if(i==4) /* SATELLITE SURFACE TEMPERATURE */
      sscanf(data, "%lf\n", &initCd->Ts);
    if(i==5) /* ADSORPTION MODEL */
      sscanf(data, "%d\n", &initCd->ads_model);
    if(i==6) /* PITCH ANGLE */
      sscanf(data, "%lf\n", &initCd->pitch);
    if(i==7) /* YAW ANGLE */
      sscanf(data, "%lf\n", &initCd->yaw);
    if(i==8) /* SURFACE MATERIAL PARTICLE MASS */
      sscanf(data, "%lf\n", &initCd->surf_mass);
    if(i==9) /* SATELLITE MASS */
      sscanf(data, "%lf\n", &initCd->sat_mass);
    if(i==10) /* SATELLITE RSM FILENAME */
      sscanf(data, "%s\n", &initCd->rsmfilename);
    if(i==11) /* SATELLITE AREA FILENAME */
      sscanf(data, "%s\n", &initCd->areafilename);

  }

  /* CHECK MOLE FRACTION BOUNDS */
  for(i=0; i<NSPECIES; i++) {
    Xsum += initCd->X[i];
    if(initCd->X[i] < 0.0 || initCd->X[i] > 1.0) {
      printf("Species Mole Fraction is outside maximal bounds for species %d!\n", i);
      printf("Min X = 0.0\n");
      printf("Max X = 1.0\n");
      exit(1);
    }
  }

  /* CHECK MOLE FRACTION SUM */
  if(Xsum != 1.0) {
    printf("Species Mole Fractions don't sum to 1.0!\n");
    printf("Xsum = %e\n", Xsum);
    exit(1);
  }

  /* CHECK NUMBER DENSITY BOUNDS */
  if(initCd->n < 0.0 || initCd->n > 1.0e16) {
    printf("Total Number Density is outside maximal bounds!\n");
    printf("Min n = 0 m^-3\n");
    printf("Max n = 1.0e16 m^-3\n");
    exit(1);
  }

  /* CHECK VELOCITY BOUNDS */
  if(initCd->U < 5.5 || initCd->U > 9.5) {
    printf("Velocity is outside maximal bounds!\n");
    printf("Min U = 5.5 km/s\n");
    printf("Max U = 9.5 km/s\n");
    exit(1);
  }

  /* CHECK ATMOSPHERIC TEMPERATURE BOUNDS */
  if(initCd->Ta < 200.0 || initCd->Ta > 2000.0) {
    printf("Atmospheric Translational Temperature is outside maximal bounds!\n");
    printf("Min Ta = 200 K\n");
    printf("Max Ta = 2000 K\n");
    exit(1);
  }

  /* CHECK SATELLITE SURFACE TEMPERATURE BOUNDS */
  if(initCd->Ts < 100.0 || initCd->Ts > 500.0) {
    printf("Satellite Surface Temperature is outside maximal bounds!\n");
    printf("Min Ts = 100 K\n");
    printf("Max Ts = 500 K\n");
    exit(1);
  }

  /* CHECK ADSORPTION MODEL BOUNDS */
  if(initCd->ads_model < 0 || initCd->ads_model > 1) {
    printf("Adsorption Model is outside maximal bounds!\n");
    printf("Allowed values = 0 or 1\n");
    printf("0 = Langmuir; 1 = Fruendlich\n");
    exit(1);
  }

    /* Check orientation bounds */
  if(initCd->pitch < 0.0 || initCd->pitch > 90.0) {
    printf("Pitch is outside maximal bounds!\n");
    printf("Min Pitch = 0 degrees\n");
    printf("Max Pitch = +90 degrees\n");
    exit(1);
  }

  if(initCd->yaw < -5.0 || initCd->yaw > 5.0) {
    printf("Yaw is outside maximal bounds!\n");
    printf("Min Yaw = -5 degrees\n");
    printf("Max Yaw = +5 degrees\n");
    exit(1);
  }

  /* CONVERT PITCH AND YAW TO RADIANS */
  initCd->pitch *= DEG2RAD;
  initCd->yaw *= DEG2RAD;
      
  fclose(f);

}



/******************************** READ RSM FILE ********************************/
void read_RSM(struct initCd_struct *initCd, struct RSM_struct *RSMdata) {

  /* Read the Response Surface Model data */
  /* Input: Area File (e.g. CYGNSS_area.dat) */
  /* Output: beta[6][7]    */
  /*         xmin[6][7]    */
  /*         xrange[6][7]  */
  /*         mean[6]       */
  /*         sd[6]         */
  /*         lamz[6]       */
  /*         lamws[6]      */
  /*         x[6][1000][7] */
  /*         w[6][1000]    */

  char filepath[256];
  char filename[256];

  strcpy(filepath, "../data/");
  strcpy(filename, initCd->rsmfilename);
  strcat(filepath, filename); 

  FILE *f = fopen(filepath, "r");
  if (f == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }

  char line[1024];
  int i, j;

  double wt[1000][6];

  /* READ IN HEADER AND EMPTY LINES */
  for(i=0; i<3; i++) {
    fgets(line, 1024, f);
  }

  /* READ IN BETA */
  for(i=0; i<6; i++) {
    fscanf(f, "%lf %lf %lf %lf %lf %lf %lf\n", 
	   &RSMdata->beta[i][0], &RSMdata->beta[i][1], &RSMdata->beta[i][2], &RSMdata->beta[i][3], 
	   &RSMdata->beta[i][4], &RSMdata->beta[i][5], &RSMdata->beta[i][6]);
  }

  fgets(line, 1024, f);

  /* READ IN XMIN */
  for(i=0; i<6; i++) {
    fscanf(f, "%lf %lf %lf %lf %lf %lf %lf\n", 
	   &RSMdata->xmin[i][0], &RSMdata->xmin[i][1], &RSMdata->xmin[i][2], &RSMdata->xmin[i][3], 
	   &RSMdata->xmin[i][4], &RSMdata->xmin[i][5], &RSMdata->xmin[i][6]);
  }

  fgets(line, 1024, f);

  /* READ IN XRANGE */
  for(i=0; i<6; i++) {
    fscanf(f, "%lf %lf %lf %lf %lf %lf %lf\n", 
	   &RSMdata->xrange[i][0], &RSMdata->xrange[i][1], &RSMdata->xrange[i][2], &RSMdata->xrange[i][3], 
	   &RSMdata->xrange[i][4], &RSMdata->xrange[i][5], &RSMdata->xrange[i][6]);
  }

  fgets(line, 1024, f);

  /* READ IN RSM MEAN */
  fscanf(f, "%lf %lf %lf %lf %lf %lf\n", 
	 &RSMdata->mean[0], &RSMdata->mean[1], &RSMdata->mean[2], 
	 &RSMdata->mean[3], &RSMdata->mean[4], &RSMdata->mean[5]);

  fgets(line, 1024, f);

  /* READ IN SD */
  fscanf(f, "%lf %lf %lf %lf %lf %lf\n", 
	 &RSMdata->sd[0], &RSMdata->sd[1], &RSMdata->sd[2], 
	 &RSMdata->sd[3], &RSMdata->sd[4], &RSMdata->sd[5]);

  fgets(line, 1024, f);

  /* READ IN LAMZ */
  fscanf(f, "%lf %lf %lf %lf %lf %lf\n", 
	 &RSMdata->lamz[0], &RSMdata->lamz[1], &RSMdata->lamz[2], 
	 &RSMdata->lamz[3], &RSMdata->lamz[4], &RSMdata->lamz[5]);

  fgets(line, 1024, f);

  /* READ IN LAMWS */
  fscanf(f, "%lf %lf %lf %lf %lf %lf\n", 
	 &RSMdata->lamws[0], &RSMdata->lamws[1], &RSMdata->lamws[2], 
	 &RSMdata->lamws[3], &RSMdata->lamws[4], &RSMdata->lamws[5]);

  /* READ IN HEADER AND EMPTY LINES */
  for(i=0; i<2; i++) {
    fgets(line, 1024, f);
  }

  /* READ IN X */
  for(i=0; i<6; i++) {
    for(j=0; j<1000; j++) {
      fscanf(f, "%lf %lf %lf %lf %lf %lf %lf\n", 
	   &RSMdata->x[i][j][0], &RSMdata->x[i][j][1], &RSMdata->x[i][j][2], &RSMdata->x[i][j][3], 
	   &RSMdata->x[i][j][4], &RSMdata->x[i][j][5], &RSMdata->x[i][j][6]);
    }
    fgets(line, 1024, f);
  }

  /* READ IN W TRANSPOSE*/
  for(i=0; i<1000; i++) {
    fscanf(f, "%lf %lf %lf %lf %lf %lf\n", 
	   &wt[i][0], &wt[i][1], &wt[i][2], 
	   &wt[i][3], &wt[i][4], &wt[i][5]);
  }

  /* TRANSPOSE W */
  for(i=0; i<1000; i++) {
    for(j=0; j<6; j++) {
      RSMdata->w[j][i] = wt[i][j];
    }
  }

  fclose(f);

}


/******************************** READ NUMBER OF LINES ***************************/
int read_num_lines(char filename[256]) {

  /* Read in the EGM96 coefficient file and determine number of lines*/
  /* Input:  */
  /* Output: Number of Facets */

  int num_lines = 0;
  int ch;

  char filepath[256];

  strcpy(filepath, "../data/");
  strcat(filepath, filename);
  FILE *f = fopen(filepath, "r");
  
  while (EOF != (ch=fgetc(f))) 
    if (ch=='\n')
      num_lines++;

  fclose(f);

  /* Subtract one for header */
  return(num_lines-1);

}


/******************************** READ RSM FILE ********************************/
void read_area(struct initCd_struct *initCd, struct area_struct *AreaData) {

  /* Read the Response Surface Model data */
  /* Input: CYGNSS_area.dat */
  /* Output: area[NPITCH][NYAW] */

  char filepath[256];
  char filename[256];

  int num_lines = 0;

  char line[1024];
  double junk;
  int i, j;

  /* Define the file path */
  strcpy(filepath, "../data/");
  strcpy(filename, initCd->areafilename);
  strcat(filepath, filename); 

  /* Calculate the total number of lines in the file */
  num_lines = read_num_lines(filepath);
  double *dyaw = calloc(num_lines, sizeof(double));
  double *dpitch = calloc(num_lines, sizeof(double));

  /* Open the file. Print an error if the file can't be opened */
  FILE *f1 = fopen(filepath, "r");
  if (f1 == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }

  /* Read in the header */
  fgets(line, 1024, f1);

  /* Loop over the file and assign a dummy angle array */
  for(i=0; i<num_lines; i++) {
    fscanf(f1, "%lf %lf %lf\n", &dyaw[i], &dpitch[i], &junk);
  }

  /* Loop over the dummy angle array to find npitch and nyaw */
  for(i=0; i<num_lines; i++) {
    if((dyaw[i] == dyaw[0]) && i>0) {
      AreaData->NYAW = i;
      AreaData->NPITCH = num_lines/i;
      AreaData->MinYaw = dyaw[0];
      AreaData->MaxYaw = dyaw[i-1];
      break;
    }
  }

  /* Define Min and Max Pitch Angles */
  AreaData->MinPitch = dpitch[0];
  AreaData->MaxPitch = dpitch[num_lines-1];

  /* Define Memory for **area */
  AreaData->area = calloc(AreaData->NPITCH, sizeof(double));
  for(i=0; i<AreaData->NPITCH; i++) {
    AreaData->area[i] = calloc(AreaData->NYAW, sizeof(double));
  }

  fclose(f1);

  /* Reopen the file */
  FILE *f2 = fopen(filepath, "r");
  if (f2 == NULL) {
    printf("%s: %s\n", strerror(errno), filepath);
    exit(1);
  }

  /* Read in header */
  fgets(line, 1024, f2);

  /* Read in areas from look-up table */
  for(i=0; i<AreaData->NPITCH; i++) {
    for(j=0; j<AreaData->NYAW; j++) {
      fscanf(f2, "%lf %lf %lf\n",  &junk, &junk, &AreaData->area[i][j]);
    }
  }

  fclose(f2);

  free(dyaw);
  free(dpitch);

}


double BilinearInterp(double *pitch_array, double *yaw_array, double pitch, double yaw, double param[], int ipitch, int iyaw) {

  double x1, x2, y1, y2;
  double xp, yp;
  double f11, f12, f21, f22;
  double vi;

  x1 = yaw_array[iyaw];
  x2 = yaw_array[iyaw+1];
  y1 = pitch_array[ipitch];
  y2 = pitch_array[ipitch+1];
        
  xp = yaw;
  yp = pitch;

  f11 = param[0];
  f12 = param[1];
  f21 = param[2];
  f22 = param[3];
        
  vi = (1.0/((x2-x1)*(y2-y1)))*(f11*(x2-xp)*(y2-yp)+f21*(xp-x1)*(y2-yp)+f12*(x2-xp)*(yp-y1)+f22*(xp-x1)*(yp-y1));

  return(vi);

}


void WriteOutput(struct initCd_struct *initCd, double RSMCd, double RSMA, double RSMBC, char *adsstr) {

  /* Write the inputs and outputs to the screen */
  printf("\n");
  printf("************** Inputs ****************\n");
  printf("RSM filename                              = %s\n", initCd->rsmfilename);
  printf("Area filename                             = %s\n", initCd->areafilename);
  printf("H  Mole Fraction                          = %e\n", initCd->X[0]);
  printf("He Mole Fraction                          = %e\n",  initCd->X[1]);
  printf("N  Mole Fraction                          = %e\n",  initCd->X[2]);
  printf("N2 Mole Fraction                          = %e\n", initCd->X[3]);
  printf("O  Mole Fraction                          = %e\n", initCd->X[4]);
  printf("O2 Mole Fractoin                          = %e\n", initCd->X[5]);
  printf("Total Number Density, n                   = %e m^-3\n", initCd->n);
  printf("Satellite Velocity, U                     = %e m/s\n", initCd->U);
  printf("Atmospheric Translational Temperature, Ta = %e K\n", initCd->Ta);
  printf("Satellite Surface Temperature, Ts         = %e K\n", initCd->Ts);
  if(initCd->ads_model==0) {
    adsstr = "Langmuir";
  } else {
    adsstr = "Freundlich";
  }
  printf("Adsorption Model                          = %s\n", adsstr);
  printf("Satellite Pitch Angle, Phi                = %e degrees\n", initCd->pitch*RAD2DEG);
  printf("Satellite Yaw Angle, Theta                = %e degrees\n", initCd->yaw*RAD2DEG);
  printf("Surface Material Particle Mass, ms        = %e kg\n", initCd->surf_mass);
  printf("\n");
  printf("************** Outputs ***************\n");
  printf("Drag Coefficient, Cd      = %e\n", RSMCd);
  printf("Projected Area, A         = %e m^2\n", RSMA);
  printf("Satellite Mass, m         = %e kg\n", initCd->sat_mass);
  printf("Ballistic Coefficient, BC = %e m^2/kg\n", RSMBC);
  printf("\n");

}



