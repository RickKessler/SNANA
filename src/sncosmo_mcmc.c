//Maximum number of SN in input file
#define MAXSN 100000
//Maximum number of mcmc samples
#define MAXCHAIN 100000
//Cosmological parameters for SN=4
#define NCOSPAR 4
#define NGITER 10
#define NCUTS 4
#define TRUE 1
#define FALSE 0
//Number of bins for 1D plots
#define BIN1D 200
//Number of bins in each variable for 2D plots
#define BIN2D 100

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sntools.h"
#include "sntools_output.h"

//sncosmo_mcmc functions
double cmb_point(double params[NCOSPAR],double *redshift);
double confidVal(int ix, double limfrac, int upper) ;
void deriv(double p[NCOSPAR],int ipr,double del[NCOSPAR]);
int domcmc(double parini[NCOSPAR],double errini[NCOSPAR][NCOSPAR],int chain_length,int chainmult[MAXCHAIN],
double chainlike[MAXCHAIN],double chainpar[4][MAXCHAIN]);
void getCovMatrix(int numsamp,double params[NCOSPAR],double covmatrix[NCOSPAR][NCOSPAR],char fname[]);
double like(double p[NCOSPAR]);
int readinputfile(FILE* funit);
void sortChainData();
void writeSqMatrix(char fname[],double *mat,int n);

int read_fitres(char filnam[ ]) ;

//cosmological distance functions
double cosmodl(double z);
double distance(double p[NCOSPAR],double z);
double inc(double z);

//Generic mathematical functions
void gengauss(double r[2]);
void invert(double matrix[],int n);
void jacobi(double A[],int N,double D[],double V[],int *NROT);
void ludcmp(double* a, const int n, const int ndim, int* indx, 
	    double* d, int* icon);
void lubksb(const double* a, const int n, const int ndim, 
       const int* indx, double* b);
void qSort(int gindex[],double value[], int last);
double rombint(double f(double z),double a,double b,double tol);




/*
! Supernovae cosmology parameter analysis code [i.e., "fitter"]
 
 Program to calculate the luminosity distance to supernovae
 and determine cosmological parameters and uncertianties
 npar gives number of cosmological parameters to be varied
 Any combination of Omega_DE and w0, wa, and Omega_k can be varied

! for SNANA light curve simulations based on
! Monte Carlo Markov Chains (mcmc).
!
!    USAGE
!  =================
!
!   sncosmo_mcmc  sncosmo_mcmc.input
!   where sncosmo_mcmc.input is a file of fit parameters including
!   input and output file names and locations   
!
March 5, 2012
Completely rewritten in C
Added gradient search feature at end of fit
version 2.00 

General information of usage
The program fits for 4 parametes: omega_l, w0, wa, omega_k
The input values of the parameters are taken to be the prior values
There is no provision for using more than 1 prior with inconsistent prior values
The best way to include priors is with the "Planck" matrix, which is a general Fischer matrix on the variables.
If the priors are not well described by a Fischer matrix, the results may be incaccurate.
The BAO is the Eisenstein A parameter.  The prior value is given by the starting values of the parameters.
The CMB reference is unknown, but the constraint also is calculated from the starting values of the parameters.

The starting values of 
!
! version 1.35 
! by David Parkinson and Pia Mukherjee
! 30th January 2008
!
!
!    SNANA HISTORY
!  ==================
!
! Aug 7, 2008 R.Kessler: all files combined into one file for snana.
!
! Aug 13, 2008 R.R.R.Reis: 
!   - All information in the input file is dumped to stdout.
!
!   - Parameters names and order for the
!     Planck Fisher matrix specified in output.
!                         
!   - two additional entries in the input file:
!       * User must choose if the program will use the light curve fit
!         result for distance modulus or the actual value from 
!         simulation (keyword 'USE_TRUE_MU: T' or F)
! 
!       * User must set the scale factor (keyword 'SCALE_ERR: 1.0') 
!         if  USE_TRUE_MU: T.
!
!   - still working with old style fitres format
!     if the fitres file is in the new format, convert it to the old
!     using temporary files. All tmp files are deleted after 
!     read all data
!
!
! Aug 14, 2008 R.R.R.Reis: 
!   - User now must specify in the input file the directory where
!     the plot data files should be created. It can be any path
!     up to 56 characters.
!
! Aug 19, 2008 R.R.R.Reis: 
!   - Included the subroutines rdfitres, parse_words and madabort
!     from R. Kessler, to read the new style fitres format, instead
!     using system calls. The program still accepts the old style
!
! Aug 21, 2008 R.R.R.Reis: 
!  - Bug found. The code was creating the output directory but
!    was not write the files on it. Bug fixed.
!  - The user must provide a valid output directory ('PLOT_DIR:' keyword).
!    if the directory does not exist the program will abort.
!  - 'PLOT_DIR: ' now accepts up to 200 characters.

! Aug 25, 2008 J. Marriner 
!  - Reworked read_input routine.  All variables have default 
!    values set in read_input.  Coding was simplified.  New
!    parameters added.
!    Added code to domcmc that steps in n gaussian dimensions
!    including arbitrary coefficients
!    Code added to domcmc automatically compute optimum step size
!    Some vestigal code that was commented out and no longer seemed
!    useful was removed.
!
! Aug 26, 2008 J. Marriner 
!  - Changed parameter selection so that any combination 
!    of Omega_DE, w0, wa, and Omega_k can be fit
!    NPAR: keyword is obsolete.  NPAR is calculated
!    from the "USE_parameter" key words
!
! Sep 22, 2008 R.R.R.Reis
!  - Fixed a bug when reading the fitres files produced by the
!    latest version of snlc_fit.exe
!
! March 2012: f90 -> C (J.Marriner)
!
! March 30, 2012
! Fixed a misplaced {} around the if (USE_PLANCK) sequence.
! Fixed a diagnostic printout when using cmb + planck prior
! Fixed problem with CMB prior (missing initialization)
! Fixed problem with BAO prior (missing scale factors)
!
! April 4, 2012
! Fixed a problem with parsing plot directory names
! 1.  Strip off leading blanks
! 2.  Strip off trailing blanks, /'s, CR and LF
! version 2.1.0
!
! Feb 21 2013 RK - udpate to compile with c++.
!                  Mainly replaec try* with Try* inside confidVal().
!
! Oct 27 2014 RK - switch to refactored table-read functions, SNTABLE_xxx
*/

//Global variables
FILE *fout;
//Speed of light (km/sec)
double cvel=2.99792458e5;
char  data_fname[256]={""};
char planck_fname[256]={""};
int use_lowz_sn;
int use_cmb,use_bao, use_planck;
int make_1d_plots, make_2d_plots;
int debug;
double cmb_redshift,cmb_prior, cmb_error;
double nscalar, Abao0, dAbao, zbao, ns_factor;
char plot_data_dir[256];
int cheat;
int sucheat;
int ipar[NCOSPAR];
int npar;
int ninit, ntrain, nsamps;
double in_params[NCOSPAR], in_parvar[NCOSPAR], in_parmin[NCOSPAR], in_parmax[NCOSPAR];
double zlim_mn, zlim_mx;
double scerr, sigint;
double ptol[NCOSPAR], delta[NCOSPAR];
double params[NCOSPAR];
int npt;
double zdata[MAXSN];
double mudata[MAXSN];
double sigma[MAXSN];
double delmu[MAXSN];
double z_val[MAXSN];
double z_sig[MAXSN];
double mu_val[MAXSN];
double mu_sig[MAXSN];
int simok;
double mu_sim[MAXSN];
int nchain;
int inch[MAXCHAIN];
int chainmult[MAXCHAIN];
double chainlike[MAXCHAIN];
double chainpar[4][MAXCHAIN];
int npoints;
double planck_cov_matrix[NCOSPAR][NCOSPAR];

double RNORM;
double omega_k, omega_l, omega_m, wde, wa;
double H0;
int nz = 2000;
double dtab[2000];
  double cmb_chisq;
  double sn_chisq;
  double bao_chi;
  double pl_chisq;

int main(int argc,char* argv[])
{
  FILE *finput;
  FILE *fout1;
  FILE *fout2;

  int num_samps;
  int np;
  int i, j, k;
  double coeff;
  int nlzsn;
  double ddl_lzsn;
  double z_lzsn;
  double mu_lzsn;
  double dl_lzsn;
  int npt_dat, npt_sim;
  double mu_sig_z;
  double z, dl;
  double chi;
  int idum, jdum;
  double Planck_fish[8][8];
  char planck_par_names[8][12];
  double sigmap, corr[8];
  double grand[2];
  double errini[NCOSPAR][NCOSPAR];
  char filename[128];
  int numsamp;
  int n;
  double covmatrix[NCOSPAR][NCOSPAR];
  double meanlike;
  int nbin[NCOSPAR];
  double bin[NCOSPAR], ledge[NCOSPAR];
  double limmin[NCOSPAR], limmax[NCOSPAR];
  int num_contours;
  double contours[2];
  int ie, je;
  int ix1, ix2;
  int bincounts[BIN1D];
  double binlikes[BIN1D];
  double binmaxlikes[BIN1D];
  double bins2D[BIN2D][BIN2D], bin2Dlikes[BIN2D][BIN2D], bin2Dmax[BIN2D][BIN2D];
  double binV;
  double limfrac;
  double cont_lines[NCOSPAR][2][2];
  double maxbin;
  int l;
  double spar[NCOSPAR];
  int ind;
  double maxlike;
  char rootname[256];
  double chinew, dchi, chitol;
  double diff;
  int nconv;
  int nread;
  double resultb;

  printf("sncosmo_mcmc version 2.1.0 \n");
  printf("omega_m at start=%f\n",omega_m);
  debug = 0;
  RNORM = RAND_MAX;

  printf("Start sncosmo with %s\n",argv[0]);
  if (argc<2)
    {
      printf("You must specify an input file.\n");
      exit(1);
    }
  finput = fopen(argv[1],"r");
  if (!finput)
    {
      printf("Error in opening input file=%s\nABORTING.\n",argv[1]);
      exit(2);
    }

  readinputfile(finput);
  printf("omega_m after input=%f \n",omega_m);
  
  ptol[0] = 1.e-6;
  ptol[1] = 1.e-4;
  ptol[2] = 1.e-4;
  ptol[3] = 1.e-6;

  delta[0] = 0.001;
  delta[1] = 0.01;
  delta[2] = 0.03;
  delta[3] = 0.001;

  num_samps = nsamps;

  np = npar;
  for (i=0;i<NCOSPAR;++i) params[i] = in_params[i];

  printf("No. variable parameters = %i\n", npar);  
  if (ipar[0]) 
    {
      printf("Fitting for Omega_DE\n");
      printf("Minimum=%f  Maximum=%f' Estimated variance=%f\n",
	     in_parmin[0],in_parmax[0],in_parvar[0]);
    }
  if (ipar[1]) 
    {
      printf("Fitting for w0\n");
      printf("Minimum=%f  Maximum=%f' Estimated variance=%f\n",
	     in_parmin[1],in_parmax[1],in_parvar[1]);
    }
  if (ipar[2]) 
    {
      printf("Fitting for wa\n");
      printf("Minimum=%f  Maximum=%f' Estimated variance=%f\n",
	     in_parmin[2],in_parmax[2],in_parvar[2]);
    }
  if (ipar[3]) 
    {
      printf("Fitting for Omega_k.\n");
      printf("Minimum=%f  Maximum=%f' Estimated variance=%f\n",
	     in_parmin[3],in_parmax[3],in_parvar[3]);
    }
  printf("Initial values:\n");
  printf("Omega_DE = %10.4f\n", params[0]);
  printf("w0       = %10.4f\n", params[1]);
  printf("wa       = %10.4f\n", params[2]);
  printf("Omega_K  = %10.4f\n", params[3]);
  printf("H0       = %10.4f\n", H0);
  printf ("Input fitres file: %s\n",data_fname);

  if (use_bao) printf("Using BAO prior? YES\n");
  else printf("Using BAO prior? NO\n");
 
  if (use_cmb) printf("Using CMB prior? YES\n");
  else printf("Using CMB prior? NO\n");
  
  if (use_lowz_sn) printf("Using low z SN prior? YES\n");
  else printf("Using low z SN prior? NO\n");
  
  if (use_planck) printf("Using Planck prior? YES\n");
  else printf("Using Planck prior? NO\n");
  
  printf("Planck Fisher matrix file: %s\n", planck_fname);
  printf("Initial sample size: %i\n", ninit);
  printf("Training sample size: %i\n", ntrain);
  printf("Number of samples: %i\n", num_samps);
  printf("DEBUG LEVEL = %i\n", debug);

  if (make_1d_plots) printf("Generating 1D plots? YES\n");
  else printf ("Generating 1D plots? NO\n");

  if (make_2d_plots) printf("Generating 2D plots? YES\n");
  else  printf("Generating 2D plots? NO\n");
 
  strcpy(filename,data_fname);
  strcpy(rootname,strtok(filename,"."));

  if (make_1d_plots || make_2d_plots) 
    {
      printf("The output plot data files will be created with the root file name of %s\n",rootname);
    }     

  if (cheat) 
    {
      printf("Using actual distance modulus in place of lc fit? YES\n");
      printf("Scale factor for error = %f\n", scerr);
      if (sucheat) printf("Computing distance modulus from input redshift.\n");
    }
  else
    {
      sucheat = FALSE;
      printf("Using actual distance modulus in place of lc fit? NO\n");
      printf("Ignoring scale factor.");
    }
  
  // zmin and zmax
  printf("Data redshift range limits %.4f<z<%.4f\n", zlim_mn,zlim_mx);
  printf("Intrinsic error=%f\n",sigint);


  // read in data from file ...
  coeff = 5.0/log(10.0);

  // Initialize point count
  npt = 0;
  // Add low_z SN point just like other points
  //  This is supposed to account for systematic errors.
  if(use_lowz_sn) 
    {
      nlzsn = 100;
      z_lzsn = 0.055;
      ddl_lzsn = 0.15/sqrt((double) nlzsn);
      dl_lzsn = distance(params,z_lzsn);
      mu_lzsn = 5.0*log10(dl_lzsn)+25.0;
      zdata[0] = 0.055;
      mudata[0] = mu_lzsn; 
      sigma[0] = ddl_lzsn;
      mu_sim[0] = mudata[0];
      if(debug>=3) printf("%f %f \n",z_lzsn, mu_lzsn);
      npt = 1;
    }

  npt_dat = read_fitres(data_fname);
  npt_sim = 0;
  
  for (i=0;i<npt_dat;++i)
    {
      if((z_val[i]>zlim_mn) && (z_val[i]<zlim_mx)) 
	{
	  if (npt>MAXSN) 
	    {
	      printf("Maximum number of input points (%i) exceeded.\nAborting.\n",MAXSN);
	      exit(3);
	    }
	  zdata[npt]=z_val[i];
	  mudata[npt]=mu_val[i];
	  mu_sig_z = coeff * ((1.0+z_val[i])/(z_val[i]*(1.0+z_val[i]/2.0))) * z_sig[i];
	  sigma[npt]=sqrt(mu_sig[i]*mu_sig[i]+sigint*sigint*mu_sig_z*mu_sig_z);
	  if (cheat) 
	    {
	      if (sucheat)
		{
		  z = z_val[i];
		  dl = distance(params,z);
		  mu_sim[i] = 5.0*log10(dl) + 25.0;
		}
	      gengauss(grand);
	      mudata[npt] = mu_sim[i] + grand[0]*scerr*sigma[npt];
	    }
	  if (debug>=3) printf("Data %5i %8.4f %8.4f %7.3f %7.3f %7.3f \n",
			       i,zdata[i],z_sig[i],mudata[i],sigma[i],mu_sim[i]);
	  ++npt;
	}
    }
  
  npoints = npt;
  if(debug>=1) printf("Read in %6i data points %6.3f < z < %6.3f \n",npoints,zlim_mn,zlim_mx);
  if (debug>=1 && use_lowz_sn) printf("First data point is lowz pseudo-data.\n");
  
  // initialise cmb distance prior ...
  if(use_cmb)
    {
      cmb_prior = cmb_point(params,&cmb_redshift);
      cmb_error = 0.007*cmb_prior;
      printf("CMB redshift=%.2f distance=%.4e dist_err=%.4e \n",
	      cmb_redshift,cmb_prior,cmb_error);
    }
  
  // initialise BAO prior ...
  if(use_bao) 
    {
      // Abao0 from Eisenstein et al. 2005:
      nscalar=0.95; //wmap 3 year data
      Abao0=0.469;    //(nscalar/0.98)**0.35
      dAbao=0.017;
      zbao=0.35;
      ns_factor=1.0/pow((nscalar/0.98),0.35);
      resultb = H0*distance(params,zbao)/(cvel*(1.0+zbao));
      Abao0=sqrt(omega_m)*pow(sqrt(inc(zbao))*resultb/zbao,2.0/3.0);
      printf("BAO A parameter=%.4f error=%.4f \n",Abao0,dAbao);
    }
  // read in Planck file and prepare planck cov matrix ...
				      
  if(use_planck) 
    {
      finput = fopen(planck_fname,"r");
      if (!finput)
	{
	  printf("Error in opening planck file=%s\nABORTING",planck_fname);
	  exit(2);
	}

      for(i=0;i<8;++i)
	{
	  for(j=0;j<8;++j)
	    {
              nread = fscanf(finput,"%i %i %lf ",&idum,&jdum,&Planck_fish[i][j]);
	      if (nread!=3)
		{
		  printf("Error in scanning planck file=%s.\nAborting\n",planck_fname);
		  exit(2);
		}
	    }
        }
      fclose(finput);
    

      writeSqMatrix("Planck_fisher",&Planck_fish[0][0],8);
      invert(&Planck_fish[0][0],8);
      
      strcpy(&planck_par_names[0][0],"w0");
      strcpy(&planck_par_names[1][0],"wa");
      strcpy(&planck_par_names[2][0],"Omega_DE");
      strcpy(&planck_par_names[3][0],"Omega_K");
      strcpy(&planck_par_names[4][0],"Omega_m*h^2");
      strcpy(&planck_par_names[5][0],"Omega_b*h^2");
      strcpy(&planck_par_names[6][0],"n_s");
      strcpy(&planck_par_names[7][0],"ln P");
      
      printf("Planck Fisher matrix\n");
      for (i=0;i<8;++i)
	{
	  sigmap = sqrt(Planck_fish[i][i]);
	  printf(" %3i  %11s %7.4f ",i,&planck_par_names[i][0],sigmap);
	  for (j=0;j<8;++j)
	    {
	      corr[j] = Planck_fish[j][i]/sqrt(Planck_fish[j][j])/sigmap;
	      printf(" %6.3f",corr[j]);
	    }
	  printf("\n");
	}
      
      //  Build NCOSPARxNCOSPAR Planck Prior Matrix regardless of number of parameters fit
      planck_cov_matrix[0][0] = Planck_fish[2][2]; // Omega_DE^2
      planck_cov_matrix[0][1] = Planck_fish[2][0]; // Omega_DE.w0
      planck_cov_matrix[0][2] = Planck_fish[2][1]; // Omega_DE.wa
      planck_cov_matrix[0][3] = Planck_fish[2][3]; // Omega_DE.Omega_K
      
      planck_cov_matrix[1][0] = planck_cov_matrix[0][1];
      planck_cov_matrix[1][1] = Planck_fish[0][0]; // w0^2
      planck_cov_matrix[1][2] = Planck_fish[0][1]; // w0.wa
      planck_cov_matrix[1][3] = Planck_fish[0][3]; // w0.Omega_K
      
      planck_cov_matrix[2][0] = planck_cov_matrix[0][2];
      planck_cov_matrix[2][1] = planck_cov_matrix[1][2];
      planck_cov_matrix[2][2] = Planck_fish[1][1]; // wa^2
      planck_cov_matrix[2][3] = Planck_fish[1][3]; // wa.Omega_K
      
      planck_cov_matrix[3][0] = planck_cov_matrix[0][3];
      planck_cov_matrix[3][1] = planck_cov_matrix[1][3];
      planck_cov_matrix[3][2] = planck_cov_matrix[2][3];
      planck_cov_matrix[3][3] = Planck_fish[3][3]; // Omega_K.Omega_K
      
      writeSqMatrix("Planck_matrix",&planck_cov_matrix[0][0],NCOSPAR);
      invert(&planck_cov_matrix[0][0],NCOSPAR);
    }

  // ... calculate chi-squared for input cosmology ...
  chi = like(params);
  
  if(debug>=1) printf("Initial chi-squared =%12.3e\n", chi);
  
		   
  // do mcmc runs ...

  //Initial error matrix from input values
  for (i=0;i<NCOSPAR;++i)
    {
      for (j=0;j<NCOSPAR;++j) errini[i][j] = 0.0;
      errini[i][i] = in_parvar[i]*in_parvar[i];
    }

  printf("\nBeginning initial error estimate with %i samples.\n",ninit);
  nchain = domcmc(params,errini,ninit,chainmult,chainlike,chainpar);
  numsamp = 0;
  for (n=0;n<nchain;++n) numsamp += chainmult[n];

  strcpy(filename,rootname);
  strcat(filename,"_0.covmat");
  getCovMatrix(numsamp,params,covmatrix,filename);

  printf( "\nBeginning training stage with %i samples.\n",ntrain);
  nchain = domcmc(params,covmatrix,ntrain,chainmult,chainlike,chainpar);
  numsamp = 0;
  for (n=0;n<nchain;++n) numsamp += chainmult[n];
  //  printf("nchain=%i numsamp=%i\n",nchain,numsamp);

  strcpy(filename,rootname);
  strcat(filename,"_1.covmat");
  getCovMatrix(numsamp,params,covmatrix,filename);

  printf("\nBeginning final stage with %i samples.\n",num_samps);
  nchain = domcmc(params,covmatrix,num_samps,chainmult,chainlike,chainpar);
  //Sort chain data
  for (i=0;i<nchain;++i) inch[i] = i;
  qSort(inch,chainlike,nchain-1);
  numsamp = 0;
  for (n=0;n<nchain;++n) numsamp += chainmult[n];

  //   printf("nchain=%i numsamp=%i\n",nchain,numsamp);
  strcpy(filename,rootname);
  strcat(filename,".covmat");
  getCovMatrix(numsamp,params,covmatrix,filename);
	       
  // analyse chains ...

  // Find best fit, and mean likelihood     
  // Best fit is first row since we have sorted the lines
  ind = inch[0];
  maxlike = chainlike[ind];
  printf("\nBest fit -Ln(like) = %e \n",maxlike);

  meanlike = 0.0;  
  ind = inch[nchain-1];
  if (chainlike[ind] - maxlike < 30) 
    {
      for (n=0;n<nchain;++n)
	{
	  //mean 1/like
	  meanlike += chainmult[n]*exp(chainlike[n]-maxlike);
	}
      meanlike = log(meanlike/numsamp) + maxlike;
      printf(" -Ln(Mean(1/like)) = %e\n",meanlike);
    }
  
  for (n=0;n<nchain;++n)
    {
      meanlike += chainmult[n]*chainlike[n];
    }
  meanlike /= numsamp;
  printf("   Mean(-Ln(like)) = %e\n",meanlike);
  
  for (n=0;n<nchain;++n)
    {
      //mean like
      meanlike += chainmult[n]*exp(maxlike-chainlike[n]);
    }
  meanlike = -log(meanlike/numsamp) + maxlike;
  printf("   -Ln(mean like)  = %e\n\n",meanlike);

  //   Just enter binning as fixed for now. 
  //   ledge=lower edge, bin=bin size, nbin=no. of bins
  //   2D binning is the same variables, but redefined to coarser bins (2x) below
  ledge[0] = 0.68;
  bin[0] = 0.001;
  nbin[0] = 100;
  ledge[1] = -1.6;
  bin[1] = 0.01;
  nbin[1] = 120;
  ledge[2] = -2.0;
  bin[2] = .05;
  nbin[2] = 80;
  ledge[3] = -0.03;
  bin[3] = 0.0005;
  nbin[3] = 120;
     
  // Output files for 1D plots
  // Do 1D bins

  // ... find limits ...
  for (i=0;i<NCOSPAR;++i)
    {
      if (ipar[i]<=0) continue;
      limmin[i] = chainpar[i][0];
      limmax[i] = chainpar[i][0];
      for (n=1;n<nchain;++n)
	{
	  if (chainpar[i][n]<limmin[i]) limmin[i]=chainpar[i][n];
	  if (chainpar[i][n]>limmax[i]) limmax[i]=chainpar[i][n]; 
	}
    }
  
  num_contours=2;
  contours[0]=0.68;
  contours[1]=0.95;

  //!  allocate(cont_lines(n,2,num_contours))
  for (j=0;j<NCOSPAR;++j)
    {
      je = ipar[j]-1;
      if (ipar[j]<=0) continue;
      for (ix1=0;ix1<120;++ix1)
	{
	  bincounts[ix1] = 0;
	  binlikes[ix1] = 0;
	  binmaxlikes[ix1] = 0;
	}
      //     !   binsraw = 0
      //!** Add debug print of contour limits
      printf("Contour limits\n");
      for (ix1=0;ix1<num_contours;++ix1)
	{
	  limfrac = (1.0-contours[ix1])/2.0;
	  cont_lines[j][1][ix1] = confidVal(j,limfrac,1) ;
	  cont_lines[j][0][ix1] = confidVal(j,limfrac,0) ;
	  // Print contour limits.
	  printf("Param=%i level=%4.2f %8.3f %8.3f \n",j,contours[ix1],
		 cont_lines[j][0][ix1],cont_lines[j][1][ix1]);
	}

      for (n=0;n<nchain;++n)
	{
	  ix2 =  (int) ((chainpar[j][n]-ledge[j])/bin[j]+1.0 ) ;
	  if (ix2>=0 && ix2<nbin[j]) 
	  {
	    bincounts[ix2] = bincounts[ix2] + chainmult[n];
	    binlikes[ix2] = binlikes[ix2] + chainmult[n]*exp(meanlike-chainlike[n]);
	    binV = exp(maxlike-chainlike[n]);
	    if (binmaxlikes[ix2]<binV) binmaxlikes[ix2]=binV;
	  }
	}

      for(i=0;i<nbin[j];++i)
	{
	  if (bincounts[i]>0) binlikes[i] = binlikes[i]/bincounts[i];
	}

     if (make_1d_plots) 
       {
	 if(debug>=1) printf("Constructing 1-d prob. distribution for param %i\n",j);
	 sprintf(filename,"%s/%s_p%i.txt",plot_data_dir,rootname,j);
	 printf("Writing=%s\n",filename);
	 fout1= fopen(filename,"w");
	 if (!fout1)
	   {
	     printf("Failed to open file=%s.\nABORTING.\n",filename);
	     exit(21);
	   }
	 //is loop counter i ok here
	 maxbin = bincounts[0];
	 for (i=1;i<nbin[j];++i) if (maxbin<bincounts[i]) maxbin=bincounts[i];

	 if (maxbin==0) 
	   {
	     printf("no samples in bin, param:%i\n",j);
	     exit(66);
	   }
	 for (i=0;i<nbin[j];++i)
	   {
	     fprintf(fout1,"%16.7e %16.7e %16.7e \n",
		     ledge[j] + (i-0.5)*bin[j], bincounts[i]/maxbin, binmaxlikes[i]);
	   }
	 fclose(fout1);
       }
     
    }

  // construct two dimensional probability distributions ...

  
    //  Change definition to be the same as 1D plots 
    //  (bin2Dmax=max liklihood, not min chisq
  for (i=0;i<BIN2D;++i)
    {
      for (j=0;j<BIN2D;++j) bin2Dmax[i][j] = 0.0;
    }

  bin[0] = 2.0*bin[0];
  bin[1] = 2.0*bin[1];
  bin[2] = 2.0*bin[2];
  bin[3] = 2.0*bin[3];
  nbin[0] = nbin[0]/2;
  nbin[1] = nbin[1]/2;
  nbin[2] = nbin[2]/2;
  nbin[3] = nbin[3]/2;

  for(i=0;i<NCOSPAR-1;++i)
    {
      ie = ipar[i]-1;
      if (ipar[i]<=0) continue;
      for (j=i+1;j<NCOSPAR;++j)
	{
	  je = ipar[j]-1;
	  if (ipar[j]<=0) continue;
	  if(debug>=1) printf("Constructing 2-d prob. distribution for params %i and %i \n",i,j);
	  for (k=0;k<BIN2D;++k)
	    {
	      for (l=0;l<BIN2D;++l)
		{
		  bins2D[k][l] = 0.0;
		  bin2Dlikes[k][l] = 0.0;
		  bin2Dmax[k][l] = -999.99;
		}
	    }
	  for (n=0;n<nchain;++n)
	    {
	      ix1= (int) ( (chainpar[i][n]-ledge[i])/bin[i]+1.0 ) ;
	      ix2= (int) ( (chainpar[j][n]-ledge[j])/bin[j]+1.0 ) ;
	      if (ix1>=0 && ix1<nbin[i] && ix2>=0 && ix2<nbin[j]) 
		{
		  bins2D[ix1][ix2] = bins2D[ix1][ix2] + chainmult[n];
		  bin2Dlikes[ix1][ix2] = bin2Dlikes[ix1][ix2] + chainmult[n]*exp(meanlike-chainlike[n]);
		  binV = exp(maxlike-chainlike[n]);
		  if (bin2Dmax[ix1][ix2]<binV)  bin2Dmax[ix1][ix2]=binV;
		}
	    }

	  for (ix1=0;ix1<nbin[i];++ix1)
	    {
	      for (ix2=0;ix2<nbin[j];++ix2)
		{
		  if (bins2D[ix1][ix2]>0.0) bin2Dlikes[ix1][ix2] = bin2Dlikes[ix1][ix2]/bins2D[ix1][ix2];
		}
	    }

	  maxbin = bins2D[0][0];
	  for (ix1=0;ix1<nbin[i];++ix1)
	    {
	      for (ix2=0;ix2<nbin[j];++ix2) if(maxbin<bins2D[ix1][ix2]) maxbin=bins2D[ix1][ix2];
	    } 
	  if (make_2d_plots) 
	    {
	      sprintf(filename,"%s/%s_2D_p%i_p%i.txt",plot_data_dir,rootname,i,j);
	      printf("Writing file=%s\n",filename);
	      fout2 = fopen(filename,"w");
	      if (!fout2)
		{
		  printf("Failed to open file=%s.\nABORTING.\n",filename);
		  exit(21);
		}
	    //!** Attempt to calculate DETF FoM from 2D plots.
	    //!        area1 = 0.0
	    //!        area2 = 0.0
	    //!        iprob = 0.0
	    //!        totprob = 0.0
	    //!** end
	      for(ix1=0;ix1<nbin[i];++ix1)
		{
		  for(ix2=0;ix2<nbin[j];++ix2)
		    {
		    // Output bin2Dmax for 2D plots (same as 1D plots)
		      fprintf(fout2,"%16.7e %16.7e %16.7e %16.7e \n",ledge[i]+(ix1-0.5)*bin[i], 
			      ledge[j]+(ix2-0.5)*bin[j],bins2D[ix1][ix2]/maxbin, bin2Dmax[ix1][ix2]);
		      //!** Begin 1 sigma and 2 sigma contour calculations
		    //!                 if (bin2Dmax(ix1,ix2)>0.607) then
		    //!                    area1 = area1 + 1.0
		    //!                    iprob = iprob + bin2Dmax(ix1,ix2)
		    //!                    endif
		    //!                 if (bin2Dmax(ix1,ix2)>0.135) area2 = area2 + 1.0
		    //!                 totprob = totprob + bin2Dmax(ix1,ix2)
		    //!** end area calcs
		    }
		}
	     /*
!** Print out area calculation results--still in development
!           iprob=iprob/totprob
!           write (*,'(a,i2,i2,2f10.4,2f10.4)') &
!                'Areas', i, j, area1, area2, iprob, totprob
!           area1 = area1*width(i)*width(j)
!           area2 = area2*width(i)*width(j)
!           write (*,'(a,i2,i2,2f10.4)') 'Areas', i, j, area1, area2
!** end
	     */
	      fclose(fout2);
	    }
	}
    }

  fout = fopen("resid.dat","w");
  if (!fout)
    {
      printf("Failed to open file=resid.dat.\nABORTING.\n");
      exit(11);
    }
  for (i=0;i<npoints;++i) fprintf(fout,"sn data %i %f %f %f %f \n",i,zdata[i],mudata[i],delmu[i],sigma[i]);
  fclose(fout);
  //Do gradient search starting from mean parameters
				      
  printf("Starting gradient search\n");
  for (n=0;n<NGITER;++n)
    {
      for (i=0;i<NCOSPAR;++i) spar[i] = params[i];
      deriv(params,0,delta);
	     for (i=0;i<NCUTS;++i)
	       {
		 chinew = like(params);
		 // Check for improved chi-squared
		 dchi = chi - chinew;
		 chitol = 1.e-4;
		 if (dchi>-chitol) break;
	 
		 //No. Cut step size
		 for  (j=0;j<NCOSPAR;++j) params[j] = 0.5*params[j] + 0.5*spar[j];
	       }
	     if (i==NCUTS)  
	       {
		 printf("Gradient search failed.  Too many cut steps\n");
		 n = NGITER;
		 break;
	       }
	    nconv = 0;
	    if (fabs(dchi)>1.e-6) nconv=1; 
	    for (i=0;i<NCOSPAR;++i)
	      {
		diff = params[i] - spar[i];
		if (fabs(diff)>ptol[i]) nconv=nconv+1;
	      }
	    chi = chinew;
	    if (nconv==0) break;
    }
  if (n<NGITER) printf("Gradient search converged.\n");
  else printf("Gradient search failed.  Too many steps or cut steps.\n");

  //** final call to print parameters around solution
  deriv(params,1,delta);
  //chi = like(params,1);
  
  //! Ad-hoc Calculation
  // use_bao = FALSE;
  //use_cmb = FALSE;
  //use_planck = FALSE;
  //H0 = 65.0;
  //params[0] = 0.73;
  //params[1] = -1.0;
  //params[2] = 0.0;
  //params[3] = 0.0;
  //delta[0] = 0.001;
  //delta[1] = 0.01;
  //delta[2] = 0.03;
  //delta[3] = 0.001;
  
  //delta[0] = 0.00001;
  //delta[1] = 0.00002;
  //delta[2] = 0.0002;
  //delta[3] = 0.00001;
  
  //deriv(params,2,delta);
  //!End Ad-hoc calculation
   return(0);
}


double cmb_point(double params[NCOSPAR], double *redshift)
{
  double ombh2, omdmh2, dist, h, z;
  printf("cmb_point omega_m=%f \n",omega_m);
  h = (65.0/100.0);

  ombh2 = 0.02229;
  omdmh2 = params[0]*h*h - ombh2;

  //! initialise last scaterring surface measurement
  //Redshift copied from the original but looks wrong?  Reference?
  z = (1.0+ (0.0783*pow(ombh2,-0.238)/(1.0+39.5*pow(ombh2,0.763)))*(pow(omdmh2+ombh2,0.560/(1.+21.1*pow(ombh2,1.81)))));
  z *= 1048.0*(1.0+0.00124*pow(ombh2,-0.738));
		
  dist = distance(params,z);
  *redshift = z;
  return(dist);
}


double confidVal ( int ix, double limfrac, int upper ) {

  // Feb 2013: try -> Try so that it compiles with c++.
  //   (RK)     Lower-case 'try' seems to conflict with something ??
  //

  double Try, lasttry, Try_m, Try_b, Try_t ;
  double samps ;
  int n ;


  Try_b = chainpar[ix][0] ;
  Try_t = chainpar[ix][0] ;
  samps = chainmult[0];

  for (n=1;n<nchain;++n)
    {
      if (Try_b>chainpar[ix][n]) Try_b=chainpar[ix][n];
      if (Try_t<chainpar[ix][n]) Try_t=chainpar[ix][n];
      samps += chainmult[n];
    }

  Try = -1.0;
  if (upper) 
    {
      do 
	{
	  lasttry = Try;
	  Try_m = (Try_b+Try_t)/2.0;
	  Try = 0.0;
	  for (n=0;n<nchain;++n) if (chainpar[ix][n]>Try_m) Try+=chainmult[n];
	  if (Try > samps*limfrac) Try_b = Try_m;
	  else Try_t = Try_m;
	} while (Try!=lasttry);
    }
  else
    {
      do 
	{
	  lasttry = Try;
	  Try_m = (Try_b+Try_t)/2.0;
	  Try = 0.0;
	  for (n=0;n<nchain;++n) if (chainpar[ix][n]<Try_m) Try+=chainmult[n];
	  if (Try > samps*limfrac) Try_t = Try_m;
	  else Try_b = Try_m;
	} while (Try!=lasttry);
    }
  return(Try_t);
} 


void deriv(double p[NCOSPAR],int ipr,double del[NCOSPAR])
{
  double delta[NCOSPAR], pp[NCOSPAR];
  double l0, lp, lm, lpp, lpm, lmp, lmm;
  double  dl1[NCOSPAR], dl2[NCOSPAR][NCOSPAR], err[NCOSPAR*NCOSPAR];
  double  corr[NCOSPAR], derr[NCOSPAR];
  int i, j, k;
  int ie, je;

  if (ipr>0) 
    {
      printf(" Computing derivatives around current solution.\n");
      printf("Initial: Omega_DE=%15.5e w0=%15.5e wa=%15.5e Omega_K=%15.5e\n",p[0],p[1],p[2],p[3]);
    }
  
  l0 = like(p);
  
  for (i=0;i<NCOSPAR;++i)
    {
      for (j=0;j<NCOSPAR;++j) pp[j] = p[j];
      pp[i] = p[i] + del[i];
      lp = like(pp);
      pp[i] = p[i] - del[i];
      lm = like(pp);
      dl1[i] = (lp-lm)/(2.0*del[i]);
      dl2[i][i] = (lp+lm-2.0*l0)/(del[i]*del[i]);
    }
  
  for (i=0;i<NCOSPAR;++i)
    {
      for (j=0;j<NCOSPAR;++j)
	{
	  if (i>=j) continue;
	  for (k=0;k<NCOSPAR;++k) pp[k] = p[k];
	  pp[i] = p[i] + del[i];
	  pp[j] = p[j] + del[j];
	  lpp = like(pp);
	  pp[j] = p[j] - del[j];
	  lmp = like(pp);
	  pp[i] = p[i] - del[i];
	  lmm = like(pp);
	  pp[j] = p[j] + del[j];
	  lpm = like(pp);
	  dl2[i][j] = (lpp-lpm-lmp+lmm)/(4.0*del[i]*del[j]);
	  dl2[j][i] = dl2[i][j];
	}
    }

  
  for (i=0;i<NCOSPAR;++i)
    {
      ie = ipar[i] - 1;
      if (ie<0) continue;
      for (j=0;j<NCOSPAR;++j) 
	{
	  je = ipar[j] - 1;
	  if (je<0) continue;
	  err[ie*npar+je] = 0.5*dl2[i][j];
	}
    }

  invert(err,npar);
  
  for (i=0;i<NCOSPAR;++i)
    {
      delta[i] = 0.0;
      derr[i] = 0.0;
      ie = ipar[i] - 1;
      if (ie<0) continue;
      for (j=0;j<NCOSPAR;++j) 
	{
	  je = ipar[j] - 1;
	  if (je<0) continue;
	  delta[i] += 0.5*err[ie*npar+je]*dl1[j];
	}
      derr[i] = sqrt(err[ie*npar+ie]);
    }
  
  for (i=0;i<NCOSPAR;++i)
    {
      pp[i] = p[i];
      p[i] = p[i] - delta[i];
    }
  
  if (ipr>0) 
    {
      printf("Derivs:  %8.5f %8.5f %8.5f %8.5f \n",dl1[0],dl1[1],dl1[2],dl1[3]);
      printf("Deltas:  %8.5f %8.5f %8.5f %8.5f \n",delta[0],delta[1],delta[2],delta[3]);
      printf("Updated: %8.5f %8.5f %8.5f %8.5f \n ",p[0],p[1],p[2],p[3]);
      printf("Log(likelihood)=%14.4e\n",l0);
      for (i=0;i<NCOSPAR;++i)
	{
	  ie = ipar[i] - 1;
	  for (j=0;j<i;++j)
	    {
	      corr[j] = 0;
	      je = ipar[j] - 1;
	      if (ie<0 || je<0) continue;
	      corr[j] = err[ie*npar+je]/(derr[i]*derr[j]);
	    }
	  if (i==0) printf("Omega_DE %12.4e +/- %-12.4e\n",pp[0],derr[0]);
	  if (i==1) printf("   w0    %12.4e +/- %-12.4e %7.3f \n",pp[1],derr[1],corr[0]);
	  if (i==2) printf("   wa    %12.4e +/- %-12.4e %7.3f %7.3f \n",pp[2],derr[2],corr[0],corr[1]);
	  if (i==3) printf("Omega_K  %12.4e +/- %-12.4e %7.3f %7.3f %7.3f \n",pp[3],derr[3],corr[0],corr[1],corr[2]);
	}
    }
  
  if (ipr>=2) 
    {
      printf("Error Matrix\n");
      for (i=0;i<npar;++i)
	{
	  for (j=0;j<npar;++j) printf("%15.5e ",err[i*npar+j]);
	  printf("\n");
	}
    }
  
  if (ipr>=1) 
    {
      printf("Fisher Matrix\n");
      for (i=0;i<NCOSPAR;++i) 
	{
	  if (ipar[i]<=0) continue;
	  for (j=0;j<NCOSPAR;++j) 
	    {
	      if (ipar[j]>0) printf("%15.5e ",0.5*dl2[i][j]);
	    } 
	  printf("\n");
	}
    }
  
  return;
}    

int domcmc(double parini[NCOSPAR],double errini[NCOSPAR][NCOSPAR],int chain_length,int chainmult[MAXCHAIN],double chainlike[MAXCHAIN],
double chainpar[4][MAXCHAIN])
{
  double eigval[NCOSPAR];
  double parc[NCOSPAR], g3[NCOSPAR], parn[NCOSPAR], maxpar[NCOSPAR];
  double parvar[NCOSPAR*NCOSPAR], eigvec[NCOSPAR*NCOSPAR];
  
  double logzero,maxlike,curlike,loglike,u,oldlike;
  int numaccept,acceptnum,numpropose,mult ;
  int i, j, ie, je;
  int nrot;
  int nprint;
  nprint = 0;
  int gotnew;
  for(i=0;i<NCOSPAR;++i)
    {
      ie = ipar[i]-1;
      if (ie<0) continue;
      for (j=0;j<NCOSPAR;++j)
	{
	  je = ipar[j]-1;
	  if (je<0) continue;
	  parvar[ie*npar+je] = errini[i][j];
	}
    }
  //check if jacobi works for nxn matrix
 jacobi(parvar,npar,eigval,eigvec,&nrot);

  for (i=0;i<npar;++i) eigval[i] = sqrt(eigval[i]);
  
  logzero = 1.0e30;
  maxlike=logzero;
  curlike=logzero;
  acceptnum=chain_length;
  numaccept=0;
  mult=1;
  numpropose=0;

  for (i=0;i<NCOSPAR;++i) parc[i] = parini[i];

  loglike= 0.5*like(parc);
  curlike = loglike;
  //!   write(*,*) 'typical loglike', loglike
  //!   pause

   while (numaccept < acceptnum)
     {
       ++numpropose;
      
       if (debug==1 && numpropose%5000==0) printf("numaccepted=%i, numproposed=%i \n", numaccept,numpropose);
       
       // get new proposed point and likelihood
       gotnew=0;
       for (i=0;i<NCOSPAR;++i) parn[i]=parc[i];
       while(!gotnew)
	 {
	   gotnew=1;
	   gengauss(g3);
	   gengauss(&g3[2]);
	   for (i=0;i<NCOSPAR;++i)
	     {
	       parn[i] = parc[i];
	       ie = ipar[i]-1;
	       if (ie<0) continue;
	       for (j=0;j<NCOSPAR;++j)
		 {
		   je = ipar[j]-1;
		   if (je<0) continue;
		   //check ie, je order.  Does it matter?
		   parn[i]=parn[i]+eigvec[ie*npar+je]*eigval[je]*g3[je];
		 }
	     }
	   
	   for (i=0;i<NCOSPAR;++i) if (parn[i]<in_parmin[i] || parn[i]>in_parmax[i]) gotnew=0;
	 }
       loglike = 0.5*like(parn);
      //!     write(*,*) 'got new parameter point', parn, 'likelihood', loglike

      //! implement metropolis hastings

       u = rand()/RNORM;
        if ((loglike != logzero) && (curlike>loglike || u<exp(-loglike+curlike))) 
	 {
	   chainmult[numaccept] = mult;
	   chainlike[numaccept] = curlike;
	   for (i=0;i<NCOSPAR;++i) chainpar[i][numaccept] = parc[i];
	   ++numaccept;

	   oldlike=curlike;
	   curlike=loglike;
	   if(curlike<maxlike) 
	     {
	       maxlike=curlike;
	       for (i=0;i<NCOSPAR;++i) maxpar[i]=parc[i];
	     }
	   if(debug>=2) 
	     {
	       //!            write(*,*) 'point accepted'
	       //!            write(*,*) 'new point, likelihood', parn, curlike
	       //!            write(*,*) 'moving from: point, likelihood', parc, oldlike
	       //!            write(*,*) 'max like so far: point, likelihood', maxpar, maxlike
	       //temporary disable
	       //    printf("numaccepted=%i numproposed=%i \n", numaccept,numpropose);
	     }
	   for (i=0;i<NCOSPAR;++i) parc[i]=parn[i];
	   mult=1;
	 }
       else
	 {
	   mult=mult+1;
//!         write(*,*) 'staying on', curlike, parc
//!         write(*,*) 'rejected point', loglike, parn
//!         write(*,*) 'rejection ratio, proposed like, current like ', &
//!              1.-real(numaccept)/real(numpropose),loglike,curlike
	 }
    //   ### End to Do */
     }
   return(numaccept);
}
void getCovMatrix(int numsamp,double params[NCOSPAR],double covmatrix[NCOSPAR][NCOSPAR],char fname[])
{
  int i, j, n;
  int ie, je;
  int ncor;
  double sigres[NCOSPAR], corr[NCOSPAR];
  //Is this needed???
  double corrmatrix[NCOSPAR][NCOSPAR];
  double a_p, z_p, sigwp, DETF_FoM;
  double mat[8][8];
  double scale;
  double chisq;

  for (i=0;i<NCOSPAR;++i)
    {
      params[i] = in_params[i];
      for (j=0;j<NCOSPAR;++j) covmatrix[i][j] = 0.0;
      if (ipar[i]<=0) continue;
      params[i] = 0;
      for (n=0;n<nchain;++n) params[i] += (chainmult[n]*chainpar[i][n]);
      params[i] /= numsamp;
    }

  for (i=0;i<NCOSPAR;++i)
    {
      covmatrix[i][i] = 1.e-20;
      ie = ipar[i]-1;
      if (ie<0) continue;
      for (j=i;j<NCOSPAR;++j)
	{
	  je = ipar[j] - 1;
	  if (je<0) continue;
	  covmatrix[i][j] = 0.0;
	  for (n=0;n<nchain;++n) 
	    {
	      covmatrix[i][j] += (chainmult[n]*(chainpar[i][n]-params[i])*(chainpar[j][n]-params[j]));
	    }
	  
	  covmatrix[i][j] /= numsamp;
	  covmatrix[j][i] = covmatrix[i][j];
	}
    }

  chisq = like(params);
  printf("\nChi-squared at current solution=%12.3e\n",chisq);
  printf("Contribution from SN=%12.3e\n",sn_chisq);
  if (use_cmb) printf("Contribution from CMB prior=%12.3e\n",cmb_chisq);
  if (use_bao) printf("Contribution from BAO prior=%12.3e\n",bao_chi);
  if (use_planck) printf("Contribution from Planck prior=%12.3e\n",pl_chisq);

    printf("\nWriting covariance matrix for %i parameters\n\n",npar);
    for (i=0;i<NCOSPAR;++i)
      {
	sigres[i] = sqrt(covmatrix[i][i]);
	if (ipar[i]<=0) continue;
	ncor = 0;
	for (j=0;j<i;++j)
	  {
	    if (ipar[j]<=0) continue;
	    corr[ncor] = covmatrix[i][j]/(sigres[i]*sigres[j]);
	    ++ncor;
	  }
	if (i==0) printf("Omega_DE %9.4f +/-%-9.4f", params[0],sigres[0]);
	if (i==1) printf("w_0      %9.3f +/-%-9.3f", params[1],sigres[1]);
	if (i==2) printf("w_a      %9.2f +/-%-9.2f", params[2],sigres[2]);
	if (i==3) printf("Omega_K  %9.4f +/-%-9.4f", params[3],sigres[3]);
	if (ncor>0) 
	  {
	    printf(" corrs");
	    for (j=0;j<ncor;++j) printf("%6.2f ",corr[j]);
	  }
	printf("\n");
      }
    
    if (ipar[1]>0 && ipar[2]>0) 
    {
      a_p = 1.0 + covmatrix[1][2]/covmatrix[2][2];
      z_p = 1.0/a_p - 1.0;
      if (z_p<0.0) z_p = 0.0;
      if (z_p>9.999) z_p = 9.999;
      sigwp = (1.0-a_p)*sigres[2];
      sigwp = sqrt(sigres[1]*sigres[1] - sigwp*sigwp);
      DETF_FoM = 1.0/(sigwp*sigres[2]);
      printf("a_p=%8.3f   z_p=%5.3f\n",a_p,z_p);
      printf("sigma(w_p)  %8.3f\n",sigwp);
      printf("DETF FoM    %8.3f\n",DETF_FoM);
    }

    for (i=0;i<NCOSPAR;++i)
      {
	ie = ipar[i] - 1;
	if (ie<0) continue;
	for (j=0;j<NCOSPAR;++j)
	  {
	    je = ipar[j]-1;
	    if (je<0) continue;
	    mat[ie][je] = covmatrix[i][j];
	  }
      }

    writeSqMatrix(fname,&mat[0][0],npar);  
    if(debug>=3) printf("Covariance matrix written to %s\n",fname);

    for (i=0;i<NCOSPAR;++i)
      {
        scale = sqrt(covmatrix[i][i]);
        for (j=0;j<NCOSPAR;++j) corrmatrix[i][j] = covmatrix[i][j]/(sigres[i]*sigres[j]);
      }
 
    for (i = 0;i<NCOSPAR;++i)
      {
	ie = ipar[i]-1;
	if (ie<0) continue;
	for (j=0;j<NCOSPAR;++j)
	  {
	    je = ipar[j]-1;
	    if (je<0) continue;
	    mat[ie][je] = corrmatrix[i][j];
	  }
      }
    writeSqMatrix(fname,&mat[0][0],npar);
    if(debug>=3) printf("Correlation matrix written to %s\n",fname);
    return;
}
double like(double p[NCOSPAR])
{
  double chisq;
// Update cosmology parameters for current step
  double mu0;
  int i, j;
  int sn_marge;
  double aprima, bprima, cprima;
  double del, z, delz ;
  int iz;
  double dl, dz;
  double sig2;
  double a1;
  double resultb;
  double cmb_red, cmb_dist;
  double dev;
  double di, dj;
  double Abao;

  omega_l = p[0];
  wde = p[1];
  wa = p[2];
  omega_k = p[3];
  omega_m = 1.0 - omega_l - omega_k;
  if(fabs(omega_k)<1.0e-8) omega_k = 0.0;
			       
  //** Build distance table.  The idea here is to speed up the program
  //** distances are calculated only once for a given cosmology.  
  //** Individual SN found via lookup table dtab
			       
  dtab[0] = 0.0;
  z = 0.0;
  delz = 2.0/nz;
  del = inc(0.0);
  for (i=1;i<=nz;++i)
    {
      dtab[i] = del + dtab[i-1];
      dtab[i-1] = 0.5*delz*dtab[i-1];
      if(omega_k==0.0) dtab[i-1] = (1.0+z)*(cvel/H0)*dtab[i-1];
      else if(omega_k<0.0) dtab[i-1] = (1.0+z)*(cvel/H0)*sin(sqrt(-omega_k)*dtab[i-1])/sqrt(-omega_k);
      else dtab[i-1] = (1.0+z)*(cvel/H0)*sinh(sqrt(omega_k)*dtab[i-1])/sqrt(omega_k);      
     z = z + delz;
     del = inc(z);
     dtab[i] = dtab[i] + del;
    }
 
  dtab[nz] = 0.5*delz*dtab[nz];
  if(omega_k==0.0) dtab[nz] = (1.0+z)*(cvel/H0)*dtab[nz];
  else if(omega_k<0.0) dtab[nz] = (1.0+z)*(cvel/H0)*sin(sqrt(-omega_k)*dtab[nz])/sqrt(-omega_k);
  else dtab[nz] = (1.0+z)*(cvel/H0)*sinh(sqrt(omega_k)*dtab[nz])/sqrt(omega_k);
    
  sn_marge = TRUE;
					       
  // SN data - may include lowz "super" point
  for (i=0;i<npoints;++i)
    {
      z=zdata[i];
      // get interpolated distance
      iz = (int)(z/delz);
      dz = z/delz-iz;
      dl = dtab[iz]*(1.-dz) + dtab[iz+1]*dz;
      mu0 = 5.0*log10(dl) + 25.0;
      delmu[i] = mu0-mudata[i];
    }
  
  aprima = 0.0;
  bprima = 0.0;
  cprima = 0.0;
  for (i=0;i<npoints;++i)
    {
      sig2 = sigma[i]*sigma[i];
      aprima += delmu[i]*delmu[i]/sig2;
      bprima += delmu[i]/sig2;
      cprima += 1.0/sig2;
    }
  if(sn_marge) 
    {
      // to calculate the abslolute chisquare
      a1=log(cprima/(2.0*M_PI));
      sn_chisq=a1+aprima-((bprima*bprima)/(cprima));
    }
  else
    {
      a1 = 0.0;
      sn_chisq=aprima;
    }
  

  //  printf("sn chisq %13.4e  minimum %13.4e \n",sn_chisq, a1);
  chisq = sn_chisq;

  if(use_bao) 
    {
      resultb = H0*distance(p,zbao)/(cvel*(1.0+zbao));
      Abao=sqrt(omega_m)*pow(sqrt(inc(zbao))*resultb/zbao,2.0/3.0);
      dev=(Abao-Abao0)/dAbao;
      bao_chi = dev*dev;
      chisq += bao_chi;
    }


  if(use_cmb) 
    {
      // cmb point
      cmb_dist = cmb_point(p,&cmb_red);
      dev = (cmb_dist-cmb_prior)/cmb_error;
      cmb_chisq = dev*dev;
      //      if (pflag>0) printf("cmb chisq %15.6e", cmb_chisq);
      chisq += cmb_chisq;
    }

  if(use_planck) 
    {
      pl_chisq = 0.0;
      for (i=0;i<NCOSPAR;++i)
	{
	  di = p[i] - in_params[i];
	  for (j=0;j<NCOSPAR;++j)
	    {
	      dj = p[j] - in_params[j];
	      pl_chisq  += di*planck_cov_matrix[i][j]*dj;
	    }
	}
      //      if (pflag>0) printf("planck chisq %15.6e \n", pl_chisq);
      chisq += pl_chisq;
    }

  //  if (pflag>0) printf("total chisq %15.6e \n", chisq);
  
return(chisq);
}

int readinputfile(FILE* finput)
{
  int i;
  int nlen, use_w0, use_wa ;
  char instring[128];
  double in_H0 ;
  int use_omega_de, use_omega_k;
  //! Preset defaults.  Overwrite as set in parameter file 
  //! No default for data file
  char logical;
  strcpy(planck_fname,"planck_wayne_stage_II.dat");
  make_1d_plots = TRUE;
  make_2d_plots = TRUE;
  strcpy(plot_data_dir,"plots");
  // Debug level
  debug = 0;
  // Number of samples
  nsamps = 50000;
  ninit = 1000;
  ntrain = 5000;
  // intrinsic error
  sigint = 0.08;
  // cheat=.true. means use simulated value (not light curve fit value)
  cheat = FALSE;
  sucheat = FALSE;
  // If cheat=.true., the simulated values are randomly smeared by scerr
  // standard deviations 
  scerr = 0.0;
  // Redshift limits on input data
  zlim_mn = 0.02;
  zlim_mx = 1.40;
  // Default priors
  use_bao = FALSE;
  use_cmb = FALSE;
  use_lowz_sn = FALSE;
  use_planck = TRUE;
  // Default parameters = fit all
  use_omega_de = TRUE;
  use_w0 = TRUE;
  use_wa = TRUE;
  use_omega_k = TRUE;
  // Hubble constant fixed value
  in_H0 = 65.0;
  // Starting values of parameters
  in_params[0] = 0.73;
  in_params[1] = -1.0;
  in_params[2] = 0.0;
  in_params[3] = 0.0;
  // Parameter error estimates
  in_parvar[0] = 0.002;
  in_parvar[1] = 0.02;
  in_parvar[2] = 0.1;
  in_parvar[3] = 0.002;
  // Parameter minimum values
  in_parmin[0] = 0.05;
  in_parmin[1] = -5.0;
  in_parmin[2] = -5.0;
  in_parmin[3] = -1.0;
  // Parameter maximum values
  in_parmax[0] = 1.0;
  in_parmax[1] = 5.0;
  in_parmax[2] = 5.0;
  in_parmax[3] = 5.0;

  // read in input file
  while (fgets(instring,111,finput))
    {
      instring[112] = 0;
      //      printf("read string=%s\n",instring);
      if(!strncmp(instring,"FITRES_FILE:",12)) 
	{
	  sscanf(&instring[12],"%s",data_fname);
	  printf("data=%s\n",data_fname);
	}

      if(!strncmp(instring,"PLANCK_FILE:",12)) 
	{
	  sscanf(&instring[12],"%s",planck_fname);
	  printf("planck=%s\n",planck_fname);
	}

      if(!strncmp(instring,"BAO_PRIOR:",10)) 
	{
	  logical=*strtok(&instring[10]," ");
	  if(logical == 'T') use_bao = TRUE;
	  else if (logical == 'F') use_bao = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}

      if(!strncmp(instring,"CMB_PRIOR:",10)) 
	{
	  logical=*strtok(&instring[10]," ");
	  if(logical == 'T') use_cmb = TRUE;
	  else if (logical == 'F') use_cmb = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}

      if(!strncmp(instring,"LZ_SN_PRIOR:",12)) 
	{
	  logical=*strtok(&instring[12]," ");
	  if(logical == 'T') use_lowz_sn = TRUE;
	  else if (logical == 'F') use_lowz_sn = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}

      if(!strncmp(instring,"PLANCK_PRIOR:",13)) 
	{
	  logical=*strtok(&instring[13]," ");
	  if(logical == 'T') use_planck = TRUE;
	  else if (logical == 'F') use_planck = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}


      if (!strncmp(instring,"DEBUG:",6)) sscanf(&instring[6],"%i",&debug);

      if (!strncmp(instring,"NUM_SAMPLES:",12)) sscanf(&instring[12],"%i",&nsamps);
      if (!strncmp(instring,"NUM_INITIAL:",12)) sscanf(&instring[12],"%i",&ninit);
      if (!strncmp(instring,"NUM_TRAIN:",10)) sscanf(&instring[10],"%i",&ntrain);


      if (!strncmp(instring,"H0:",3)) sscanf(&instring[4],"%lf",&in_H0);
	

      if(!strncmp(instring,"OMEGA_DE:",9)) sscanf(&instring[9],"%lf",&in_params[0]);
      if(!strncmp(instring,"W0:",3))  sscanf(&instring[3],"%le",&in_params[1]);
      if(!strncmp(instring,"WA:",3))  sscanf(&instring[3],"%le",&in_params[2]);
      if(!strncmp(instring,"OMEGA_K:",9))  sscanf(&instring[9],"%le",&in_params[3]);
      if(!strncmp(instring,"USE_OMEGA_DE:",13)) 
	{
	  logical=*strtok(&instring[13]," ");
	  if(logical == 'T') use_omega_de = TRUE;
	  else if (logical == 'F') use_omega_de = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	   }
	}
	
      if(!strncmp(instring,"USE_W0:",7)) 
	{
	  logical=*strtok(&instring[7]," ");
	  if(logical == 'T') use_w0 = TRUE;
	  else if (logical == 'F') use_w0 = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}
	  
      if(!strncmp(instring,"USE_WA:",7)) 
	{
	  logical=*strtok(&instring[7]," ");
	  if(logical == 'T') use_wa = TRUE;
	  else if (logical == 'F') use_wa = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}

      if(!strncmp(instring,"USE_OMEGA_K:",12)) 
	{
	  logical=*strtok(&instring[12]," ");
	  if(logical == 'T') use_omega_k = TRUE;
	  else if (logical == 'F') use_omega_k = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}
      
      if (!strncmp(instring,"OMEGA_DE_RANGE:",15)) sscanf(&instring[15],"%lf %lf %lf",&in_parvar[0],&in_parmin[0],&in_parmax[0]);
      if (!strncmp(instring,"W0_RANGE:",9)) sscanf(&instring[9],"%lf %lf %lf",&in_parvar[1],&in_parmin[1],&in_parmax[1]);
      if (!strncmp(instring,"WA_RANGE:",9)) sscanf(&instring[9],"%lf %lf %lf",&in_parvar[2],&in_parmin[2],&in_parmax[2]);
      if (!strncmp(instring,"OMEGA_K_RANGE:",14)) sscanf(&instring[14],"%lf %lf %lf",&in_parvar[3],&in_parmin[3],&in_parmax[3]);
 
					      
      if (!strncmp(instring,"ZLIM_MN:",8)) sscanf(&instring[8],"%lf",&zlim_mn);
      if (!strncmp(instring,"ZLIM_MX:",8)) sscanf(&instring[8],"%lf",&zlim_mx);


      if(!strncmp(instring,"1D_PLOTS:",9)) 
	{
	  logical=*strtok(&instring[9]," ");
	  if(logical == 'T') make_1d_plots = TRUE;
	  else if (logical == 'F') make_1d_plots = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}
      if(!strncmp(instring,"2D_PLOTS:",9)) 
	{
	  logical=*strtok(&instring[9]," ");
	  if(logical == 'T') make_2d_plots = TRUE;
	  else if (logical == 'F') make_2d_plots = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}

      if (!strncmp(instring,"PLOT_DIR:",9)) 
	{
	  //skip leading whitespec
	  for (i=9;i<200;++i) if (instring[i]!=' ') break;
	  strncpy(plot_data_dir,&instring[i],200-i);
	  //eliminate trailing whitespece, /, CR, LF
	  nlen = strlen(plot_data_dir);
	  for (i=nlen-1;i>=0;--i)
	    {
	      if (plot_data_dir[i]==' '  ||
		  plot_data_dir[i]=='/'  ||
		  plot_data_dir[i]=='\n' ||
		  plot_data_dir[i]=='\r') 
		  plot_data_dir[i]=0;
	      else break;
	    }
	  //eliminate leading whitespace
	  nlen = strlen(plot_data_dir);
	}

      if(!strncmp(instring,"USE_TRUE_MU:",12)) 
	{
	  logical=*strtok(&instring[12]," ");
	  if(logical == 'T') cheat = TRUE;
	  else if (logical == 'F') cheat = FALSE;
	  else
	    {
	      printf("non-logical input for %s\n",instring);
	      exit(1);
	    }
	}

      if(!strncmp(instring,"USE_EXACT_MU:",13)) 
	{
	  logical=*strtok(&instring[13]," ");
	  if(logical == 'T') sucheat = TRUE;
	  else if (logical == 'F') cheat = FALSE;
	  else
	    {
	      printf("non logical input for %s",instring);
	      exit(1);
	    }
	}

      if (!strncmp(instring,"SCALE_ERR:",10)) sscanf(&instring[10],"%lf",&scerr);
      if (!strncmp(instring,"INTRINSIC_ERR:",14)) sscanf(&instring[14],"%lf",&sigint);

    }

  H0 = in_H0;

  npar = 0;
  for (i=0;i<NCOSPAR;++i) ipar[i] = 0;
  if (use_omega_de) 
    {
      ++npar;
      ipar[0] = npar;
    }
  if (use_w0)
    {
      ++npar;
      ipar[1] = npar;
    }
  
  if (use_wa)
    {
      ++npar;
      ipar[2] = npar;
    }
  
  if (use_omega_k) 
    {
      ++npar;
      ipar[3] = npar;
    }
  
  if(npar==0 || npar>NCOSPAR)
    {
      printf("Number of parameters to be varied out of range. npar=%i\n", npar);
      exit(2);
    }
  
  if(!strlen(data_fname)) 
    {
      printf("No datafile specified.\n");
      exit(3);
    }
  
  if(use_planck && use_cmb) 
    {
      printf("Cannot use both cmb distance prior and planck Fisher matrix.\n");
      exit(1);
    }
  
  
  if(use_planck  && !strlen(planck_fname)) 
    {
      printf("No Planck file specified.\n");
      exit(3);
    }
  

  return(0);
}



void writeSqMatrix(char fname[],double *mat,int n)
{
  FILE *fout;
  int i, j, in;

  fout = fopen(fname,"w");
  if (!fout)
    {
      printf("Failed to open file=%s.\nABORTING.\n",fname);
      exit(11);
    }


  for (i=0;i<n;++i)
    {
      for (j=0;j<n;++j) 
	{
	  in = 8*i + j;
	  fprintf(fout," %15.6e",mat[in]);
	}
      fprintf(fout,"\n");
    }  
  fclose(fout);
  return;
}

double cosmodl(double z)
{
  const double tol=1.e-6;
  double dflat;
  double dist;

  if(fabs(omega_k)<tol) omega_k = 0.0;
  omega_m = 1.0 - omega_l - omega_k;
  //comoving distance to redshift
  dflat = rombint(inc,0.0,z,tol);
  //  printf("z %f dflat %f \n",z,dflat);

  if(omega_k==0.0) dist = cvel*(1.0/H0)*dflat;
  else if(omega_k<0.0) 
    dist = cvel*(1.0/H0)*(1.0/sqrt(-omega_k))*sin(sqrt(-omega_k)*dflat);
  else dist = cvel*(1.0/H0)*(1.0/sqrt(omega_k))*sinh(sqrt(omega_k)*dflat);
  return((1.0+z)*dist);
}
double distance(double p[NCOSPAR],double z)
{
  omega_l = p[0];
  wde = p[1];
  wa = p[2];
  omega_k = p[3];
  omega_m = 1.0 - omega_l - omega_k;
  return(cosmodl(z));
}
double inc(double z)
{
  double hubble, rhode;
  double zz;

  zz = 1.0 + z;
  rhode = omega_l*pow(zz,3.0*(1.0+wde+wa));
  rhode = rhode*exp(-3.0*(wa*z/zz));
  hubble = sqrt((omega_m*(zz*zz*zz))+rhode+(omega_k*(zz*zz)));
  return(1.0/hubble);					  
}

int read_fitres(char filnam[ ])
{
  int nsn;
  char vartmp[20];

  //Open input file
  printf("Reading fitres file=%s.",filnam);

  TABLEFILE_INIT();
  int IFILETYPE = TABLEFILE_OPEN(filnam,"read");
  SNTABLE_READPREP(IFILETYPE,"FITRES");
  //  fitresFile_init(filnam);


  //Get redshift: key = Z, z, or REDSHIFT
  sprintf(vartmp,"Z z REDSHIFT");
  SNTABLE_READPREP_VARDEF(vartmp,z_val,MAXSN,3);  // 3-> abort if missing


  //Get redshift error: key=ZERR
  sprintf(vartmp,"ZERR");
  SNTABLE_READPREP_VARDEF(vartmp,z_sig,MAXSN,3);

  //Get distance modulus: key=DLMAG or MU
  sprintf(vartmp,"DLMAG MU");
  SNTABLE_READPREP_VARDEF(vartmp,mu_val,MAXSN,3);


  //Get distance modulus: key=DLMAGERR or MUERR
  sprintf(vartmp,"DLMAGERR MUERR");
  SNTABLE_READPREP_VARDEF(vartmp,mu_sig,MAXSN,3); 

  
  //Get simulation mu: key=SIMDLMAG or SIMMU
  simok = 0 ;
  sprintf(vartmp,"SIMDLMAG SIMMU");
  if ( SNTABLE_READPREP_VARDEF(vartmp,mu_sim,MAXSN,1) >0 )
    { simok = 1; }

  
  nsn = SNTABLE_READ_EXEC();
  
  if (nsn > MAXSN) 
    {
      printf("Number of SN=%i greater than maximum=%i.\n", nsn,MAXSN);
      exit(11);
    }

return(nsn);
}



void invert(double matrix[],int n)
{
  //Inverts a nxn matrix[n*n]
  //n must be less than 10 (but no check)
  int i, j;
  int indx[10];
  double matinv[100];
  double d;
  int icon;

  ludcmp(matrix,n,n,indx,&d,&icon);
  for (i=0;i<n;++i)
    {
      for (j=0;j<n;++j) matinv[n*i+j]=0.0;
      matinv[n*i+i] = 1.0;
      lubksb(matrix,n,n,indx,&matinv[n*i]);
    }


  //Replace input matrix with inverse
  for (i=0;i<n*n;++i) matrix[i] = matinv[i];
  return;
}

void jacobi(double A[],int N,double D[],double V[],int* NROT)
{
 /*
!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real symmetric square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
Converted from FORTRAN by J. Marriner
 */
  int i, j, ip, iq;
  double c, g, h, s, sm, t, tau, theta, tresh;
  double B[100], Z[100];

  //!initialize V to identity matrix
  for (ip=0; ip<N; ++ip)    
    {
      for (iq=0;iq<N;++iq) V[ip*N+iq]=0.0; 
      V[ip*N+ip]=1.0;
    }
    
  for (ip=0; ip<N; ++ip)
    {
      B[ip]=A[N*ip+ip];
      D[ip]=B[ip];
      Z[ip]=0.0;    
    }
  *NROT=0;

  for (i=0;i<50;++i)
    {
      sm=0.0;
      //  !sum off-diagonal elements
      for (ip=0;ip<N-1;++ip)
	{
	  for (iq=ip+1;iq<N;++iq) sm += fabs(A[ip*N+iq]);
	}
      if(sm==0.0) return;  //!normal return
      if(i<3) tresh=0.2*sm*sm;
      else tresh=0.0;
    
      for (ip=0;ip<N-1;++ip)
	{
	  for (iq=ip+1;iq<N;++iq)
	    {
	      g=100.0*fabs(A[ip*N+iq]);
	      //! after 4 sweeps, skip the rotation if the off-diagonal element is small
	      if((i>3) && (fabs(D[ip])+g==fabs(D[ip])) && (fabs(D[iq])+g==fabs(D[iq]))) A[ip*N+iq]=0.0;
	      else if(fabs(A[ip*N+iq])>tresh) 
		{
		  h=D[iq]-D[ip];
		  if(fabs(h)+g==fabs(h)) t=A[ip*N+iq]/h;
		  else
		    {
		      theta=0.5*h/A[ip*N+iq];  
		      t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
		      if(theta<0.0) t=-t;
		    }
		  c=1.0/sqrt(1.0+t*t);
		  s=t*c;
		  tau=s/(1.0+c);
		  h=t*A[ip*N+iq];
		  Z[ip]=Z[ip]-h;
		  Z[iq]=Z[iq]+h;
		  D[ip]=D[ip]-h;
		  D[iq]=D[iq]+h;
		  A[ip*N+iq]=0.0;
		  for(j=0;j<=ip-1;++j)
		    {
		      g=A[j*N+ip];
		      h=A[j*N+iq];
		      A[j*N+ip]=g-s*(h+g*tau);
		      A[j*N+iq]=h+s*(g-h*tau);
		    }
		  for (j=ip+1;j<=iq-1;++j)
		    {
		      g=A[ip*N+j];
		      h=A[j*N+iq];
		      A[ip*N+j]=g-s*(h+g*tau);
		      A[j*N+iq]=h+s*(g-h*tau);
		    }		      
		  for (j=iq+1;j<N;++j)
		    {
		      g=A[ip*N+j];
		      h=A[iq*N+j];
		      A[ip*N+j]=g-s*(h+g*tau);
		      A[iq*N+j]=h+s*(g-h*tau);
		    }		  
		  for (j=0;j<N;++j)
		    {
		      g=V[j*N+ip];
		      h=V[j*N+iq];
		      V[j*N+ip]=g-s*(h+g*tau);
		      V[j*N+iq]=h+s*(g-h*tau);
		    }		  
		  *NROT=*NROT+1;
		  //end of if if ((i>4)
		}
	      //end iq loop
	    }
	  // end ip loop
	}
      for (ip=0;ip<N;++ip)
	{
	  B[ip]=B[ip]+Z[ip];
	  D[ip]=B[ip];
	  Z[ip]=0.0;
	}
      // end do !main i loop
    }
  printf("Jacobi 50 iterations !\n");
  exit(66);
}

void gengauss(double r[2])
{
  // Generates two random numbers with a Gaussian distribution
  double radius, phi;
  radius = sqrt(-2.0*log(rand()/RNORM));
  phi = TWOPI*rand()/RNORM;
  //  printf("radius %f phi %f \n",radius,phi);
  r[0] = radius*cos(phi);
  r[1] = radius*sin(phi);
  return;
}
void ludcmp(double* a, const int n, const int ndim, int* indx, 
       double* d, int* icon)
{
  /* System generated locals */
  int i1=0, i2=0, i3;

  /* Local variables */
  int imax = 0;
  int  i, j, k;
  double aamax, dum, sum;
  double vv[100] = { 0.0 };
    
  /* Function Body */
  *icon = 1;
  *d = 1.0;


  for (i = 0; i < n; ++i) 
  {
    aamax = 0.0;
  
    for (j = 0; j < n; ++j) 
    {
      if (fabs(a[i + j * ndim]) > aamax) 
      {
	aamax = fabs(a[i + j * ndim]);
      }
    }

    /* MH       IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.' */
    if (aamax == 0.0) 
    {
      printf("LU decomposition (ludcmp.c) : Singular matrix, but continue. %d %d %d %d %d\n",i,j,i1,i2,ndim);
      *icon = -10;
      return;
      
    }

    vv[i] = 1.0 / aamax;
  }

  /* OK to HERE */
  i1 = n;
  for (j = 0; j < n; ++j) 
  {
    if (j > 0) 
    {
      for (i = 0; i < j; ++i) 
      {
	sum = a[i + j * ndim];
	if (i > 0) 
	{
	  i3 = i - 1;
	  for (k = 0; k < i; ++k) 
	  {
	    sum -= a[i + k * ndim] * a[k + j * ndim];
	  }
	  a[i + j * ndim] = sum;
	}
      }
    }

    aamax = 0.0;

    for (i = j; i < n; ++i) 
    {
      sum = a[i + j * ndim];
      if (j > 0) 
      {
     
	for (k = 0; k < j; ++k) 
	{
	  sum -= a[i + k * ndim] * a[k + j * ndim];
	}
	a[i + j * ndim] = sum;
      }

      dum = vv[i] * fabs(sum);
      if (dum >= aamax) 
      {
	imax = i;
	aamax = dum;
      }
    }

    if (j != imax) 
    {
    
      for (k = 0; k < n; ++k) 
      {
	dum = a[imax + k * ndim];
	a[imax + k * ndim] = a[j + k * ndim];
	a[j + k * ndim] = dum;
      }
      *d = - *d;
      vv[imax] = vv[j];
    }

    indx[j] = imax;
    if (j != (n-1)) 
    {
      /*          IF(A(J,J).EQ.0.)A(J,J)=TINY */
      if (fabs(a[j + j * ndim]) <= 1.e-20) 
      {
	if (a[j + j * ndim] > 0.0) 
	{
	  a[j + j * ndim] = 1.0e-20;
	} 
	else 
	{
	  a[j + j * ndim] = -1.0e-20;
	}
      }

      dum = 1.0 / a[j + j * ndim];
     
      for (i = j + 1; i < n; ++i) 
      {
	a[i + j * ndim] *= dum;
      }
    }
  }

  /*      IF(A(N,N).EQ.0.)A(N,N)=TINY */
  if (fabs(a[n-1 + (n-1) * ndim]) <= 1.0e-20) 
  {
    if (a[n-1 + (n-1) * ndim] > 0.0) 
    {
      a[n-1 + (n-1) * ndim] = 1.0e-20;
    }
    else 
    {
      a[n-1 + (n-1) * ndim] = -1.0e-20;
    }
  }
  return;
} /* ludcmp */

void lubksb(const double* a, const int n, const int ndim, 
	    const int* indx, double* b)
{
  int i1;
  int i, j, ii, ll;
  double sum;

  
  /* Function Body */
  ii = -1;
  i1 = n;
  for (i = 0; i < i1; ++i) 
  {
    ll = indx[i];
    sum = b[ll];
    b[ll] = b[i];
    if (ii != -1) 
    {
      for (j = ii; j < i; ++j) 
      {
	sum -= a[i + j * ndim] * b[j];
      }
    } 
    else 
    {
      if (sum != 0.0) 
      {
	ii = i;
      }
    }
    b[i] = sum;
  }

  for (i = n-1; i >= 0; --i) 
  {
    sum = b[i];
    if (i < n-1) 
    {
   
      for (j = i + 1; j < n; ++j) 
      {
	sum -= a[i + j * ndim] * b[j];
      }
    }
    b[i] = sum / a[i + i * ndim];
  }
  return;
} /* lubksb */


#define INSERTION_SORT_BOUND 16 /* boundary point to use insertion sort */
 void qSort(int gindex[],double value[], int last)
{
  //Quick sort routine
  int stack_pointer = 0;
  int first_stack[32];
  int last_stack[32];
  int ifirst, ilast, imed, idown, iup;
  int first=0;
  for (;;)
  {
    if (last - first <= INSERTION_SORT_BOUND)
    {
      /* for small sort, use insertion sort */
      int indx;
      int prev_val = gindex[first];
      int cur_val;

      for (indx = first + 1; indx <= last; ++indx)
      {
        cur_val = gindex[indx];
        if (value[prev_val]>value[cur_val])
        {
          /* out of order: array[indx-1] > array[indx] */
          int indx2;
          gindex[indx] = prev_val; /* move up the larger item first */

          /* find the insertion point for the smaller item */
          for (indx2 = indx - 1; indx2 > first; )
          {
            int temp_val = gindex[indx2 - 1];
            if (value[temp_val]>value[cur_val])
            {
              gindex[indx2--] = temp_val;
              /* still out of order, move up 1 slot to make room */
            }
            else
              break;
          }
          gindex[indx2] = cur_val; /* insert the smaller item right here */
        }
        else
        {
          /* in order, advance to next element */
          prev_val = cur_val;
        }
      }
    }
    else
    {
      int pivot;
 
      /* try quick sort */
      {
        int temp;
        int med = (first + last) >> 1;
        /* Choose pivot from first, last, and median position. */
        /* Sort the three elements. */
        temp = gindex[first];
	ilast = gindex[last];
        if (value[temp]>value[ilast])
	  {
	    gindex[first] = gindex[last]; gindex[last] = temp;
	  }
        temp = gindex[med];
	ifirst = gindex[first];
        if (value[ifirst]>value[temp])
        {
          gindex[med] = gindex[first]; gindex[first] = temp;
        }
        temp = gindex[last];
	    imed = gindex[med];
        if (value[imed]>value[temp])
        {
          gindex[last] = gindex[med]; gindex[med] = temp;
        }
        pivot = gindex[med];
      }
      {
        int up;
        {
	  int down;
          /* First and last element will be loop stopper. */
	  /* Split array into two partitions. */
	  down = first;
	  up = last;
	  for (;;)
	  {
	    do
	    {
	      ++down;
	      idown = gindex[down];
	    } while (value[pivot]>value[idown]); 
	    do
	    {
	      --up;
	      iup = gindex[up];
	    } while (value[iup]>value[pivot]);
 
	    if (up > down)
	    {
	      int temp;
	      /* interchange L[down] and L[up] */
	      temp = gindex[down]; gindex[down]= gindex[up]; gindex[up] = temp;
	    }
	    else
	      break;
	  }
	}
	{
	  int len1; /* length of first segment */
	  int len2; /* length of second segment */
	  len1 = up - first + 1;
	  len2 = last - up;
	  /* stack the partition that is larger */
	  if (len1 >= len2)
	  {
	    first_stack[stack_pointer] = first;
	    last_stack[stack_pointer++] = up;
 
	    first = up + 1;
	    /*  tail recursion elimination of
	     *  Qsort(gindex,fun_ptr,up + 1,last)
	     */
	  }
	  else
	  {
	    first_stack[stack_pointer] = up + 1;
	    last_stack[stack_pointer++] = last;

	    last = up;
	    /* tail recursion elimination of
	     * Qsort(gindex,fun_ptr,first,up)
	     */
	  }
	}
        continue;
      }
      /* end of quick sort */
    }
    if (stack_pointer > 0)
    {
      /* Sort segment from stack. */
      first = first_stack[--stack_pointer];
      last = last_stack[stack_pointer];
    }
    else
      break;
  } /* end for */
}
double rombint(double f(double z),double a,double b,double tol)
{
  /*
!  Rombint returns the integral from a to b of using Romberg integration.
                                                                               
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be double precision and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
  */
#define MAXITER 23
#define MAXJ 5
  double g[MAXJ+2];
  double rombint=0.0, h,gmax,error,g0=0, fourj,g1;
  int nint,i,k,j,jmax;
   h = 0.5*(b-a);
   gmax = h*(f(a)+f(b));
   g[0] = gmax;
   nint = 1;
   error = 1.0e20;
   for (i=1;i<=MAXITER;++i)
     {
       if (fabs(error)<tol) break;
   //Calculate next trapezoidal rule approximation to integral.
       g0=0.0;
       for (k=0;k<nint;++k) g0+=f(a+(2*k+1)*h);
       g0=0.5*g[0]+h*g0;
       h=0.5*h;
       nint=nint+nint;
       if (i<MAXJ) jmax=i;
       else jmax=MAXJ;
       fourj=1.0;
    for (j=0;j<jmax;++j)
      {
	//  ! Use Richardson extrapolation.
       fourj=4.0*fourj;
       g1=g0+(g0-g[j])/(fourj-1.0);
       g[j]=g0;
       g0=g1;
      }
       if (fabs(g0)>tol) error=1.0-gmax/g0;
       else error=gmax;
 
       gmax=g0;
       g[jmax]=g0;
     }

   if (i==MAXITER && fabs(error)>tol) 
     printf("Rombint failed to converge; integral, error=%e %e \n",rombint,error);
   return(g0);
}


