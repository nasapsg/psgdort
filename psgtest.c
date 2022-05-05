#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

// Physical constants
#define CS        299792458.0  // Speed of light [m/s]
#define HP    6.626070040e-34  // Planck's constant [W s2] or [J s]
#define KB     1.38064852e-23  // Boltzmann' constant [J/K]
#define C2         1.43877736  // Second radiation constant [K/cm-1]
#define VERBOSE             1  // Verbosely print results

// Types of phases for getphase
#define ISOTROPIC            1
#define RAYLEIGH             2
#define HENYEY_GREENSTEIN    3
#define HAZE_GARCIA_SIEWERT  4
#define CLOUD_GARCIA_SIEWERT 5

// ---------------------------------------------
// Main function - Read PSG input file and call scatter
// ---------------------------------------------
// Memory allocation functions
double *array1D(long m1, int clear) {
  double *ptr; long i;
  if (clear==1) {
    ptr = calloc(m1, sizeof(double));
  } else {
    ptr = malloc(m1*sizeof(double));
    if (clear==2) for (i=0;i<m1;i++) ptr[i] = 1.0;
  }
  return ptr;
}
double **array2D(long m1, long m2, int clear) {
  double **ptr = malloc(m1*sizeof(double*)); long i;
  for (i=0;i<m1;i++) ptr[i] = array1D(m2, clear);
  return ptr;
}
double ***array3D(long m1, long m2, long m3, int clear) {
  double ***ptr = malloc(m1*sizeof(double**)); long i;
  for (i=0;i<m1;i++) ptr[i] = array2D(m2, m3, clear);
  return ptr;
}
double ****array4D(long m1, long m2, long m3, long m4, int clear) {
  double ****ptr = malloc(m1*sizeof(double***)); long i;
  for (i=0;i<m1;i++) ptr[i] = array3D(m2, m3, m4, clear);
  return ptr;
}
void free2D(double **ptr, long m1) {
  long i; for (i=0;i<m1;i++) free(ptr[i]);
  free(ptr);
}
void free3D(double ***ptr, long m1, long m2) {
  long i; for (i=0;i<m1;i++) free2D(ptr[i], m2);
  free(ptr);
}
void free4D(double ****ptr, long m1, long m2, long m3) {
  long i; for (i=0;i<m1;i++) free3D(ptr[i], m2, m3);
  free(ptr);
}

// Include main routines
#include "psgdort.c"

// Prototypes
void psgtest01(void);
void psgtest02(void);
void psgtest03(void);
void psgtest04(void);
void psgtest05(void);
void psgtest06(void);
void psgtest07(void);
void psgtest08(void);
void psgtest09(void);
void psgtest10(void);
void psgtest11(void);
void psgtest12(void);
void psgtest14(void);
void getphase(int iphas, double gg, double albedo, int nmom, int ilyr, double ****phfunc);
void testprint(int nphi, int ntau, int numu, double *phi, double *utau, double *umu, double *rfldir, double *rfldn, double *flup, double *dfdt, double ****uu, double *grfldir, double *grfldn, double *gflup, double *gdfdt, double ****guu);

// Main function
int main() {

  psgtest01();
  psgtest02();
  psgtest03();
  psgtest04();
  psgtest05();
  psgtest06();
  psgtest07();
  psgtest08();
  psgtest09();
  psgtest10();
  psgtest11();
  psgtest12();
  psgtest14();
  return 0;
}

// ----------------------------------------------------------------
// Test Problem 1:  Isotropic Scattering                       ----
// (Compare to Ref. VH1, Table 12)                             ----
// ----------------------------------------------------------------
void psgtest01(void) {

  int nstr   = 16;
  int nlyr   = 1;
  int ntau   = 2;
  int numu   = 6;
  int nphi   = 1;
  int nmom   = nstr;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0=0.1;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=NULL;
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1); //  Corrected intensity field

  umu[0] = -1.0;
  umu[1] = -0.5;
  umu[2] = -0.1;
  umu[3] =  0.1;
  umu[4] =  0.5;
  umu[5] =  1.0;

  utau[0] = 0.0;

  tsrf[0] = 0.0;
  esrf[0] = 1.0 - tsrf[0];

  for (icas=1; icas<=6; icas++) {
    switch(icas) {
      case 1:
        utau[1] = 0.03125;
        getphase(ISOTROPIC, 0.0, 0.2, nstr, 0, phfunc);

        fstar[0] = 1.0/umu0;

        // Correct answers
        grfldir[0] = 3.14159;       grfldir[1] = 2.29844;
        grfldn[0]  = 0.0;           grfldn[1]  = 7.94108E-02;
        gflup[0]   = 7.99451E-02;   gflup[1]   = 0.0;
        gdfdt[0]   = 2.54067E+01;   gdfdt[1]   = 1.86531E+01;
        guu[0][0][0][0] = 0.0;         guu[0][1][0][0]=0.0;         guu[0][2][0][0]=0.0;         guu[0][3][0][0]=1.17771E-01; guu[0][4][0][0]=2.64170E-02; guu[0][5][0][0]=1.34041E-02;
        guu[0][0][1][0] = 1.33826E-02; guu[0][1][1][0]=2.63324E-02; guu[0][2][1][0]=1.15898E-01; guu[0][3][1][0]=0.0;         guu[0][4][1][0]=0.0;         guu[0][5][1][0]=0.0;
      break;
      case 2:
        utau[1]  =  .03125;
        getphase(ISOTROPIC, 0.0, 1.0, nstr, 0, phfunc);

        fstar[0] = 1.0/umu0;
        fisot = 0.0;

        // Correct answers
        grfldir[0]=3.14159,    grfldir[1]=2.29844;
        grfldn[0] =0.,         grfldn[1] =4.20233E-01;
        gflup[0]  =4.22922E-01,gflup[1]  =0.;
        gdfdt[0]  =0.,         gdfdt[1]  =0.;
        guu[0][0][0][0]=0.;          guu[0][1][0][0]=0.;          guu[0][2][0][0]=0.;          guu[0][3][0][0]=6.22884E-01; guu[0][4][0][0]=1.39763E-01; guu[0][5][0][0]=7.09192E-02;
        guu[0][0][1][0]=7.08109E-02; guu[0][1][1][0]=1.39337E-01; guu[0][2][1][0]=6.13458E-01; guu[0][3][1][0]=0.;          guu[0][4][1][0]=0.;          guu[0][5][1][0]=0.;
      break;
      case 3:
        utau[1]  = .03125;
        getphase(ISOTROPIC, 0.0, 0.99, nstr, 0, phfunc);

        fstar[0] = 0.0;
        fisot = 1.;

        // Correct answers
        grfldir[0]=0.,         grfldir[1]=0.;
        grfldn[0] =3.14159,    grfldn[1] =3.04897;
        gflup[0]  =9.06556E-02,gflup[1]  =0.;
        gdfdt[0]  =6.66870E-02,gdfdt[1]  =5.88936E-02;
        guu[0][0][0][0]=1.;          guu[0][1][0][0]=1.;          guu[0][2][0][0]=1.;          guu[0][3][0][0]=1.33177E-01; guu[0][4][0][0]=2.99879E-02; guu[0][5][0][0]=1.52233E-02;
        guu[0][0][1][0]=9.84447E-01; guu[0][1][1][0]=9.69363E-01; guu[0][2][1][0]=8.63946E-01; guu[0][3][1][0]=0.;          guu[0][4][1][0]=0.;          guu[0][5][1][0]=0.;
      break;
      case 4:
        utau[1]  = 32.;
        getphase(ISOTROPIC, 0.0, 0.2, nstr, 0, phfunc);

        fstar[0] = 1.0/umu0;
        fisot = 0.0;

        // Correct answers
        grfldir[0]=3.14159,    grfldir[1]=0.;
        grfldn[0] =0.,         grfldn[1] =0.;
        gflup[0]  =2.59686E-01,gflup[1]  =0.;
        gdfdt[0]  =2.57766E+01,gdfdt[1]  =0.;
        guu[0][0][0][0]=0.;          guu[0][1][0][0]=0.;          guu[0][2][0][0]=0.;          guu[0][3][0][0]=2.62972E-01; guu[0][4][0][0]=9.06967E-02; guu[0][5][0][0]=5.02853E-02;
        guu[0][0][1][0]=1.22980E-15; guu[0][1][1][0]=1.30698E-17; guu[0][2][1][0]=6.88840E-18; guu[0][3][1][0]=0.;          guu[0][4][1][0]=0.;          guu[0][5][1][0]=0.;
      break;
      case 5:
        utau[1]  = 32.;
        getphase(ISOTROPIC, 0.0, 1.0, nstr, 0, phfunc);

        fstar[0] = 1.0/umu0;
        fisot = 0.0;

        // Correct answers
        grfldir[0]=3.14159,grfldir[1]=0.;
        grfldn[0] =0.,     grfldn[1] =6.76954E-02;
        gflup[0]  =3.07390,gflup[1]  =0.;
        gdfdt[0]  =0.,     gdfdt[1]  =0.;
        guu[0][0][0][0]=0.;          guu[0][1][0][0]=0.;          guu[0][2][0][0]=0.;          guu[0][3][0][0]=1.93321E+00; guu[0][4][0][0]=1.02732E+00; guu[0][5][0][0]=7.97199E-01;
        guu[0][0][1][0]=2.71316E-02; guu[0][1][1][0]=1.87805E-02; guu[0][2][1][0]=1.16385E-02; guu[0][3][1][0]=0.;          guu[0][4][1][0]=0.;          guu[0][5][1][0]=0.;
      break;
      case 6:
        utau[1]  = 32.;
        getphase(ISOTROPIC, 0.0, 0.99, nstr, 0, phfunc);

        fstar[0] = 0.;
        fisot = 1.;

        // Correct answers
        grfldir[0]=0.,         grfldir[1]=0.;
        grfldn[0] =3.14159,    grfldn[1] =4.60048E-03;
        gflup[0]  =2.49618,    gflup[1]  =0.;
        gdfdt[0]  =1.14239E-01,gdfdt[1]  =7.93633E-05;
        guu[0][0][0][0]=1.;          guu[0][1][0][0]=1.;          guu[0][2][0][0]=1.;          guu[0][3][0][0]=8.77510E-01; guu[0][4][0][0]=8.15136E-01; guu[0][5][0][0]=7.52715E-01;
        guu[0][0][1][0]=1.86840E-03; guu[0][1][1][0]=1.26492E-03; guu[0][2][1][0]=7.79280E-04; guu[0][3][1][0]=0.;          guu[0][4][1][0]=0.;          guu[0][5][1][0]=0.;
      break;

    } dtau[0][0][0] = utau[1];

    if (VERBOSE) printf("\n\nTest Case No. 1-%02d:  Isotropic Scattering, Ref. VH1, Table 12:  b =%9.5f, a =%5.2f\n",icas,utau[2],phfunc[0][0][0][0]);

    // Call the scattering module
    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
            &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    // Compare results
    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
  return;
}


//========================== psgtest02() ==============================
//---------------------------------------------------------------------
//---  Test Problem 2:  Rayleigh Scattering, Beam Source           ----
//---  (Compare To Ref. SW, Table 1)                               ----
//---------------------------------------------------------------------
void psgtest02(void) {

  int nstr   = 16;
  int nlyr   = 1;
  int ntau   = 2;
  int numu   = 6;
  int nphi   = 1;
  int nmom   = nstr;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0=0.1;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, iod, iss;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=NULL;
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field

  utau[0]   = 0.0;

  umu[0]    = -0.981986;
  umu[1]    = -0.538263;
  umu[2]    = -0.018014;
  umu[3]    = 0.018014;
  umu[4]    = 0.538263;
  umu[5]    = 0.981986;

  phi[0]    = 0.0;
  fstar[0]  = 1.0;
  umu0      = 0.080442;
  fisot     = 0.0;
  tsrf[0]   = 0.0;
  esrf[0]   = 1.0-tsrf[0];

  for (iod = 1, icas=0; iod <= 2; iod++) {
    if (iod == 1) {
      utau[1] = 0.2;
    }
    else {
      utau[1] = 5.0;
    }
    dtau[0][0][0] = utau[1];
    for (iss = 1; iss <= 2; iss++) {
      if(iss == 1) {
        getphase(RAYLEIGH, 0.0, 0.5, nstr, 0, phfunc);
      } else {
        getphase(RAYLEIGH, 0.0, 1.0, nstr, 0, phfunc);
      }

      icas++;

      switch(icas) {
        case 1:
          // Correct answers
          grfldir[0]=2.52716E-01,grfldir[1]=2.10311E-02;
          grfldn[0] =0.,         grfldn[1] =4.41791E-02;
          gflup[0]  =5.35063E-02,gflup[1]  =0.0;
          gdfdt[0]  =1.66570E+00,gdfdt[1]  =1.89848E-01;
          guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=1.61796E-01,guu[0][4][0][0]=2.11501E-02,guu[0][5][0][0]=7.86713E-03,
          guu[0][0][1][0]=7.71897E-03,guu[0][1][1][0]=2.00778E-02,guu[0][2][1][0]=2.57685E-02,guu[0][3][1][0]=0.,         guu[0][4][1][0]=0.,         guu[0][5][1][0]=0.0;
        break;
        case 2:
          // Correct answers
          grfldir[0]=2.52716E-01,grfldir[1]=2.10311E-02;
          grfldn[0] =0.,         grfldn[1] =1.06123E-01;
          gflup[0]  =1.25561E-01,gflup[1]  =0.0;
          gdfdt[0]  =0.,         gdfdt[1]  =0.0;
          guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=3.47678E-01,guu[0][4][0][0]=4.87120E-02,guu[0][5][0][0]=1.89387E-02,
          guu[0][0][1][0]=1.86027E-02,guu[0][1][1][0]=4.64061E-02,guu[0][2][1][0]=6.77603E-02,guu[0][3][1][0]=0.,         guu[0][4][1][0]=0.,         guu[0][5][1][0]=0.0;
        break;
        case 3:
          // Correct answers
          grfldir[0]=2.52716E-01,grfldir[1]=2.56077E-28;
          grfldn[0] =0.,         grfldn[1] =2.51683E-04;
          gflup[0]  =6.24730E-02,gflup[1]  =0.0;
          gdfdt[0]  =1.67462E+00,gdfdt[1]  =1.75464E-04;
          guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=1.62566E-01,guu[0][4][0][0]=2.45786E-02,guu[0][5][0][0]=1.01498E-02,
          guu[0][0][1][0]=1.70004E-04,guu[0][1][1][0]=3.97168E-05,guu[0][2][1][0]=1.32472E-05,guu[0][3][1][0]=0.,         guu[0][4][1][0]=0.,         guu[0][5][1][0]=0.0;
        break;
        case 4:
          // Correct answers
          grfldir[0]=2.52716E-01,grfldir[1]=0.0;
          grfldn[0] =0.,         grfldn[1] =2.68008E-02;
          gflup[0]  =2.25915E-01,gflup[1]  =0.0;
          gdfdt[0]  =0.,         gdfdt[1]  =0.0;
          guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=3.64010E-01,guu[0][4][0][0]=8.26993E-02,guu[0][5][0][0]=4.92370E-02,
          guu[0][0][1][0]=1.05950E-02,guu[0][1][1][0]=7.69337E-03,guu[0][2][1][0]=3.79276E-03,guu[0][3][1][0]=0.,         guu[0][4][1][0]=0.,         guu[0][5][1][0]=0.0;
        break;
      }

      if (VERBOSE) printf("\n\nTest Case No. 2-%02d, Rayleigh Scattering, Ref. SW, Table 1:  tau =%5.2f, mu0 =%9.6f, ss-albedo =%4.2f",
                        icas,utau[1],umu0,phfunc[0][0][0][0]);

      psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
            &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

      testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
    }
  }
  return;
}

//========================== psgtest03() ==============================
//---------------------------------------------------------------------
//---  Test Problem 3:  Henyey-Greenstein Scattering               ----
//---  (Compare To Ref. VH2, Table 37)                             ----
//---------------------------------------------------------------------
void psgtest03(void) {

  int nstr   = 16;
  int nlyr   = 1;
  int ntau   = 2;
  int numu   = 6;
  int nphi   = 1;
  int nmom   = 32;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0=0.1;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, iod, iss;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=NULL;
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field

  double gg = .75;
  getphase(HENYEY_GREENSTEIN, gg, 1.0, nmom, 0, phfunc);

  utau[0]   = 0.0;

  umu[0]    = -1.0;
  umu[1]    = -0.5;
  umu[2]    = -0.1;
  umu[3]    = 0.1;
  umu[4]    = 0.5;
  umu[5]    = 1.;

  phi[0]    = 0.0;

  umu0      = 1.;
  fstar[0]  = 1.0/umu0;
  fisot     = 0.0;
  tsrf[0]   = 0.0;
  esrf[0]   = 1.0 - tsrf[0];

  for (icas = 1; icas <= 2; icas++) {
    switch(icas) {
      case 1:
        utau[1] = 1.0;

        // Correct answers
        grfldir[0]=3.14159,    grfldir[1]=1.15573;
        grfldn[0] =0.,         grfldn[1] =1.73849;
        gflup[0]  =2.47374E-01,gflup[1]  =0.0;
        gdfdt[0]  =0.,         gdfdt[1]  =0.0;
        guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=1.51159E-01,guu[0][4][0][0]=1.01103E-01,guu[0][5][0][0]=3.95460E-02,
        guu[0][0][1][0]=3.05855E+00,guu[0][1][1][0]=2.66648E-01,guu[0][2][1][0]=2.13750E-01,guu[0][3][1][0]=0.,         guu[0][4][1][0]=0.,         guu[0][5][1][0]=0.0;
      break;
      case 2:
        utau[1] = 8.0;

        // Correct answers
        grfldir[0]=3.14159,    grfldir[1]=1.05389E-03;
        grfldn[0] =0.,         grfldn[1] =1.54958;
        gflup[0]  =1.59096E+00,gflup[1]  =0.0;
        gdfdt[0]  =0.,         gdfdt[1]  =0.0;
        guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=3.79740E-01,guu[0][4][0][0]=5.19598E-01,guu[0][5][0][0]=4.93302E-01,
        guu[0][0][1][0]=6.69581E-01,guu[0][1][1][0]=4.22350E-01,guu[0][2][1][0]=2.36362E-01,guu[0][3][1][0]=0.,         guu[0][4][1][0]=0.,         guu[0][5][1][0]=0.0;
      break;
    }
    dtau[0][0][0] = utau[1];

    if (VERBOSE) printf("\n\nTest Case No. 3-%02d,  Henyey-Greenstein Scattering, Ref. VH2, Table 37, g = %3.2f, b =%9.5f, a =%5.2f\n\n",
                      icas,gg,utau[1],phfunc[0][0][0][0]);

    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
          &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
  return;
}

//========================== psgtest04() ==============================
//---------------------------------------------------------------------
//---  Test Problem 4:  Haze-L Scattering, Beam Source             ----
//---  (Compare to Ref. GS, Tables 12-16)                          ----
//---------------------------------------------------------------------
void psgtest04(void) {

  int nstr   = 32;
  int nlyr   = 1;
  int ntau   = 3;
  int numu   = 6;
  int nphi   = 3;
  int nmom   = 32;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0=0.1;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, iod, iss;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu = array4D(nlam,numu,ntau,nphi,1);

  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field
  char title[100];

  fstar[0]  = 1.0;
  tsrf[0]   = 0.0;
  tsrf[0]   = 0.0;
  esrf[0]   = 1.0 - tsrf[0];

  for (icas = 1; icas <= 3; icas++) {
    sprintf(title,"Test Case No. 4-%02d, Haze-L Scattering, Ref. GS, Table ",icas);
    switch(icas) {
      case 1:
        nphi = 1;

        getphase(HAZE_GARCIA_SIEWERT, 0.0, 1.0, nmom, 0, phfunc);

        dtau[0][0][0]  = 1.;

        utau[0]   = 0.0;
        utau[1]   = 0.5;
        utau[2]   = 1.;

        umu[0]    = -1.;
        umu[1]    = -0.5;
        umu[2]    = -0.1;
        umu[3]    = 0.1;
        umu[4]    = 0.5;
        umu[5]    = 1.;

        phi[0]   = 0.0;

        umu0  = 1.0;
        if (VERBOSE) printf("\n\n%s 12\n",title);

        // Correct answers
        grfldir[0]=3.14159,    grfldir[1]=1.90547,    grfldir[2]=1.15573;
        grfldn[0] =0.,         grfldn[1] =1.17401,    grfldn[2] =1.81264;
        gflup[0]  =1.73223E-01,gflup[1]  =1.11113E-01,gflup[2]  =0.0;
        gdfdt[0]  =0.,         gdfdt[1]  =0.,         gdfdt[2]  =0.0;
        guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=9.26837E-02,guu[0][4][0][0]=6.59569E-02,guu[0][5][0][0]=3.64755E-02,
        guu[0][0][1][0]=2.51608E+00,guu[0][1][1][0]=1.19287E-01,guu[0][2][1][0]=1.34962E-01,guu[0][3][1][0]=1.23887E-01,guu[0][4][1][0]=4.02058E-02,guu[0][5][1][0]=1.77746E-02,
        guu[0][0][2][0]=3.37302E+00,guu[0][1][2][0]=2.19835E-01,guu[0][2][2][0]=1.56893E-01,guu[0][3][2][0]=0.,         guu[0][4][2][0]=0.,         guu[0][5][2][0]=0.0;
      break;
      case 2:
        nphi = 1;

        getphase(HAZE_GARCIA_SIEWERT, 0.0, 0.9, nmom, 0, phfunc);

        dtau[0][0][0] = 1.;

        utau[0]  = 0.0;
        utau[1]  = 0.5;
        utau[2]  = 1.;

        umu[0]   = -1.;
        umu[1]   = -0.5;
        umu[2]   = -0.1;
        umu[3]   = 0.1;
        umu[4]   = 0.5;
        umu[5]   = 1.;

        phi[0]   = 0.0;

        umu0  = 1.0;
        if (VERBOSE) printf("\n\n%s 13\n",title);

        // Correct answers
        grfldir[0]=3.14159,    grfldir[1]=1.90547,    grfldir[2]=1.15573;
        grfldn[0] =0.,         grfldn[1] =1.01517,    grfldn[2] =1.51554;
        gflup[0]  =1.23665E-01,gflup[1]  =7.88690E-02,gflup[2]  =0.0;
        gdfdt[0]  =3.43724E-01,gdfdt[1]  =3.52390E-01,gdfdt[2]  =3.19450E-01;
        guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=6.53056E-02,guu[0][4][0][0]=4.55144E-02,guu[0][5][0][0]=2.82693E-02,
        guu[0][0][1][0]=2.24258E+00,guu[0][1][1][0]=9.66049E-02,guu[0][2][1][0]=9.61335E-02,guu[0][3][1][0]=8.43278E-02,guu[0][4][1][0]=2.79473E-02,guu[0][5][1][0]=1.38835E-02,
        guu[0][0][2][0]=2.97057E+00,guu[0][1][2][0]=1.67698E-01,guu[0][2][2][0]=1.08115E-01,guu[0][3][2][0]=0.,         guu[0][4][2][0]=0.,         guu[0][5][2][0]=0.0;
      break;
      case 3:
        nphi = 3;

        getphase(HAZE_GARCIA_SIEWERT, 0.0, 0.9, nmom, 0, phfunc);

        dtau[0][0][0] = 1.;

        utau[0]  = 0.0;
        utau[1]  = 0.5;
        utau[2]  = 1.;

        umu[0]   = -1.;
        umu[1]   = -0.5;
        umu[2]   = -0.1;
        umu[3]   = 0.1;
        umu[4]   = 0.5;
        umu[5]   = 1.;

        phi[0]   =   0.0;
        phi[1]   =  90.0;
        phi[2]   = 180.0;

        umu0  =   0.5;
        if (VERBOSE) printf("\n\n%s 14-16\n",title);

        // Correct answers
        grfldir[0]=1.57080,    grfldir[1]=5.77864E-01,grfldir[2]=2.12584E-01;
        grfldn[0] =0.,         grfldn[1] =7.02764E-01,grfldn[2] =8.03294E-01;
        gflup[0]  =2.25487E-01,gflup[1]  =1.23848E-01,gflup[2]  =0.0;
        gdfdt[0]  =3.85003E-01,gdfdt[1]  =3.37317E-01,gdfdt[2]  =2.16403E-01;
        guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=8.70812E-01,guu[0][4][0][0]=2.24960E-01,guu[0][5][0][0]=2.27572E-02,
        guu[0][0][1][0]=4.77016E-02,guu[0][1][1][0]=3.02631E+00,guu[0][2][1][0]=1.41195E+00,guu[0][3][1][0]=6.97692E-01,guu[0][4][1][0]=1.09130E-01,guu[0][5][1][0]=9.32861E-03,
        guu[0][0][2][0]=8.38488E-02,guu[0][1][2][0]=2.70538E+00,guu[0][2][2][0]=8.76523E-01,guu[0][3][2][0]=0.,         guu[0][4][2][0]=0.,         guu[0][5][2][0]=0.,
        guu[0][0][0][1]=0.,         guu[0][1][0][1]=0.,         guu[0][2][0][1]=0.,         guu[0][3][0][1]=8.88117E-02,guu[0][4][0][1]=5.77411E-02,guu[0][5][0][1]=2.27572E-02,
        guu[0][0][1][1]=4.77016E-02,guu[0][1][1][1]=5.80971E-02,guu[0][2][1][1]=1.04502E-01,guu[0][3][1][1]=9.16071E-02,guu[0][4][1][1]=2.95842E-02,guu[0][5][1][1]=9.32861E-03,
        guu[0][0][2][1]=8.38488E-02,guu[0][1][2][1]=9.42187E-02,guu[0][2][2][1]=8.95457E-02,guu[0][3][2][1]=0.,         guu[0][4][2][1]=0.,         guu[0][5][2][1]=0.,
        guu[0][0][0][2]=0.,         guu[0][1][0][2]=0.,         guu[0][2][0][2]=0.,         guu[0][3][0][2]=6.98247E-02,guu[0][4][0][2]=5.02877E-02,guu[0][5][0][2]=2.27572E-02,
        guu[0][0][1][2]=4.77016E-02,guu[0][1][1][2]=2.58544E-02,guu[0][2][1][2]=6.25954E-02,guu[0][3][1][2]=5.91273E-02,guu[0][4][1][2]=2.47702E-02,guu[0][5][1][2]=9.32861E-03,
        guu[0][0][2][2]=8.38488E-02,guu[0][1][2][2]=3.99383E-02,guu[0][2][2][2]=4.67155E-02,guu[0][3][2][2]=0.,         guu[0][4][2][2]=0.,         guu[0][5][2][2]=0.0;
      break;
    }

    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
          &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
}


//========================== psgtest05() ==============================
//---------------------------------------------------------------------
//---  Test Problem 5:  Cloud C.1 Scattering, Beam Source          ----
//---  (Compare to Ref. GS, Tables 19-20)                          ----
//---------------------------------------------------------------------
void psgtest05(void) {

  int nstr   = 48;
  int nlyr   = 1;
  int ntau   = 3;
  int numu   = 6;
  int nphi   = 1;
  int nmom   = 299;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0=0.1;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, iod, iss;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double ****pmoms = array4D(nlam, nbin, nlyr, mmom+1, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu = NULL;
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field
  char title[100];

  dtau[0][0][0]  = 64.0;

  umu[0]    = -1.;
  umu[1]    = -0.5;
  umu[2]    = -0.1;
  umu[3]    =  0.1;
  umu[4]    =  0.5;
  umu[5]    =  1.;

  phi[0]    =  0.0;

  fstar[0]  =  1.0;
  umu0      =  1.0;
  fisot     =  0.0;
  tsrf[0]   =  0.0;
  esrf[0]   = 1.0 - tsrf[0];

  for (icas = 1; icas <= 2; icas++) {
    sprintf(title,"Test Case No. 5-%02d,  Cloud C.1 Scattering, Ref. GS, Table ",icas);
    switch(icas) {
      case 1:
        utau[0]  =  0.0;
        utau[1]  = 32.0;
        utau[2]  = 64.0;
        getphase(CLOUD_GARCIA_SIEWERT, 0.0, 1.0, nmom, 0, phfunc);
        if (VERBOSE) printf("\n\n%s 19\n",title);

        // Correct answers
        grfldir[0]=3.14159,grfldir[1]=3.97856E-14,grfldir[2]=5.03852E-28;
        grfldn[0] =0.,     grfldn[1] =2.24768,    grfldn[2] =4.79851E-01;
        gflup[0]  =2.66174,gflup[1]  =1.76783,    gflup[2]  =0.0;
        gdfdt[0]  =0.,     gdfdt[1]  =0.,         gdfdt[2]  =0.0;
        guu[0][0][0][0]=0.,         guu[0][1][0][0]=0.,         guu[0][2][0][0]=0.,         guu[0][3][0][0]=4.58927E-01,guu[0][4][0][0]=7.72983E-01,guu[0][5][0][0]=1.07196E+00,
        guu[0][0][1][0]=7.53662E-01,guu[0][1][1][0]=6.96362E-01,guu[0][2][1][0]=6.50541E-01,guu[0][3][1][0]=6.27631E-01,guu[0][4][1][0]=5.81809E-01,guu[0][5][1][0]=5.24532E-01,
        guu[0][0][2][0]=1.95230E-01,guu[0][1][2][0]=1.31990E-01,guu[0][2][2][0]=7.20655E-02,guu[0][3][2][0]=0.,         guu[0][4][2][0]=0.,         guu[0][5][2][0]=0.0;
      break;
      case 2:
        utau[0]  =  3.2;
        utau[1]  = 12.8;
        utau[2]  = 48.0;
        getphase(CLOUD_GARCIA_SIEWERT, 0.0, 0.9, nmom, 0, phfunc);
        if (VERBOSE) printf("\n\n%s 20\n",title);

        // Correct answers
        grfldir[0]=1.28058E-01,grfldir[1]=8.67322E-06,grfldir[2]=4.47729E-21;
        grfldn[0] =1.74767,    grfldn[1] =2.33975E-01,grfldn[2] =6.38345E-05;
        gflup[0]  =2.70485E-01,gflup[1]  =3.74252E-02,gflup[2]  =1.02904E-05;
        gdfdt[0]  =3.10129E-01,gdfdt[1]  =4.52671E-02,gdfdt[2]  =1.25021E-05;
        guu[0][0][0][0]=6.79623E+01,guu[0][1][0][0]=2.21027E-01,guu[0][2][0][0]=1.36619E-01,guu[0][3][0][0]=1.14084E-01,guu[0][4][0][0]=8.73870E-02,guu[0][5][0][0]=8.81626E-02,
        guu[0][0][1][0]=2.05706E-01,guu[0][1][1][0]=4.92736E-02,guu[0][2][1][0]=2.65449E-02,guu[0][3][1][0]=2.02154E-02,guu[0][4][1][0]=1.29661E-02,guu[0][5][1][0]=9.51334E-03,
        guu[0][0][2][0]=3.41286E-05,guu[0][1][2][0]=1.39916E-05,guu[0][2][2][0]=7.47039E-06,guu[0][3][2][0]=5.65602E-06,guu[0][4][2][0]=3.58245E-06,guu[0][5][2][0]=2.57858E-06;
        break;
      break;
    }

    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
          &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
  return;
}


//========================== psgtest06() ==============================
//---------------------------------------------------------------------*
//---  Test Problem 6:  No Scattering, Increasingly Complex Sources ----
//-------------------------------------------------------------------*--
void psgtest06(void) {

  int nstr   = 16;
  int nlyr   = 1;
  int ntau   = 3;
  int numu   = 4;
  int nphi   = 1;
  int nmom   = 0;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0=0.1;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, iod, iss;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);

  double *rfldir = array1D(ntau,1);// ......... Direct-beam flux (without delta-m scaling)
  double *rfldn = array1D(ntau,1);// .......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *flup = array1D(ntau,1);// ........... Diffuse up-flux
  double *fldn = array1D(ntau,1);// ........... Diffuse down-flux
  double *fldir = array1D(ntau,1);// .......... Direct flux
  double *dfdt = array1D(ntau,1);// ........... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****uu = array4D(nlam,numu,ntau,nphi,1);//.  Corrected intensity field

  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field
  char title[100];

  fstar[0]= 200.0/M_PI;
  umu0    =   0.5;
  fisot   =   0.0;
  temis   =   1.;

  for (icas = 1; icas <= 3; icas++) {
    sprintf(title,"Test Case No. 6-%02d: No Scattering; Source = Beam",icas);

    switch(icas) {
      case 1:
        // Transparent medium, beam source
        ntau = 2;

        umu[0]    = -1.;
        umu[1]    = -0.1;
        umu[2]    =  0.1;
        umu[3]    =  1.;

        phi[0]    = 90.0;

        utau[0]   = 0.0;
        utau[1]   = 0.0;

        dtau[0][0][0]  = 0.0;
        tsrf[0] = 0.0;
        esrf[0] = 1.0-tsrf[0];

        if (VERBOSE) printf("\n\n%s; Bottom Albedo = %.0f\n",title,tsrf[0]);

        // Correct answers
        grfldir[0]=100.,grfldir[1]=100.,
        grfldn[0] =0.,  grfldn[1] =0.,
        gflup[0]  =0.,  gflup[1]  =0.,
        gdfdt[0]  =200.,gdfdt[1]  =200.,
        guu[0][0][0][0]=0.,guu[0][1][0][0]=0.,guu[0][2][0][0]=0.,guu[0][3][0][0]=0.,
        guu[0][0][1][0]=0.,guu[0][1][1][0]=0.,guu[0][2][1][0]=0.,guu[0][3][1][0]=0.0;
      break;
      case 2:
        // Add some optical depth
        ntau = 3;

        umu[0]   = -1.;
        umu[1]   = -0.1;
        umu[2]   =  0.1;
        umu[3]   =  1.;

        phi[0]   = 90.0;

        utau[0]  = 0.0;
        utau[1]  = 0.5;
        utau[2]  = 1.;

        dtau[0][0][0] = 1.0;

        tsrf[0] = 0.0;
        esrf[0] = 1.0-tsrf[0];

        if (VERBOSE) printf("\n\n%s; Bottom Albedo = %.0f\n",title,tsrf[0]);

        // Correct answers
        grfldir[0]=100.,grfldir[1]=3.67879E+01,grfldir[2]=1.35335E+01;
        grfldn[0] =0.,  grfldn[1] =0.,         grfldn[2] =0.0;
        gflup[0]  =0.,  gflup[1]  =0.,         gflup[2]  =0.0;
        gdfdt[0]  =200.,gdfdt[1]  =7.35759E+01,gdfdt[2]  =2.70671E+01;
        guu[0][0][0][0]=0.,guu[0][1][0][0]=0.,guu[0][2][0][0]=0.,guu[0][3][0][0]=0.,
        guu[0][0][1][0]=0.,guu[0][1][1][0]=0.,guu[0][2][1][0]=0.,guu[0][3][1][0]=0.,
        guu[0][0][2][0]=0.,guu[0][1][2][0]=0.,guu[0][2][2][0]=0.,guu[0][3][2][0]=0.0;
      break;
      case 3:
        // Add some isotropic reflection
        ntau = 3;

        umu[0]    = -1.;
        umu[1]    = -0.1;
        umu[2]    =  0.1;
        umu[3]    =  1.;

        phi[0]    = 90.0;

        utau[0]   = 0.0;
        utau[1]   = 0.5;
        utau[2]   = 1.;

        dtau[0][0][0] = 1.0;
        tsrf[0] = 0.5;
        esrf[0] = 1.0-tsrf[0];

        if (VERBOSE) printf("\n\n%s; Bottom Albedo=%3.1f Lambert\n",title,tsrf[0]);

        // Correct answers
        grfldir[0]=100.,       grfldir[1]=3.67879E+01,grfldir[2]=1.35335E+01;
        grfldn[0] =0.,         grfldn[1] =0.,         grfldn[2] =0.0;
        gflup[0]  =1.48450E+00,gflup[1]  =2.99914E+00,gflup[2]  =6.76676E+00;
        gdfdt[0]  =2.02010E+02,gdfdt[1]  =7.79962E+01,gdfdt[2]  = 4.06006E+01;
        guu[0][0][0][0]=0.,guu[0][1][0][0]=0.,guu[0][2][0][0]=9.77882E-05,guu[0][3][0][0]=7.92386E-01,
        guu[0][0][1][0]=0.,guu[0][1][1][0]=0.,guu[0][2][1][0]=1.45131E-02,guu[0][3][1][0]=1.30642E+00,
        guu[0][0][2][0]=0.,guu[0][1][2][0]=0.,guu[0][2][2][0]=2.15393E+00,guu[0][3][2][0]=2.15393E+00;
      break;
    }

    // Force re-allocation of arrays
    dfdt=NULL, fldir=NULL, fldn=NULL, flup=NULL, rfldir=NULL, rfldn=NULL, uu = NULL;

    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
          &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
  return;
}


//========================== psgtest07() ==============================
//---------------------------------------------------------------------
//---  Test Problem 7:  Absorption + Scattering + All Possible     ----
//---  Sources, Various Surface Reflectivities ( One Layer )       ----
//--- (Compare 7a,f Fluxes and Intensities to Ref. KS, Tables I-II ----
//---------------------------------------------------------------------
void psgtest07(void) {

  // Maximum values for memory allocation
  int nlyr   = 1;
  int nphi   = 2;
  int ntau   = 3;
  int numu   = 4;
  int nmom   = 12;
  int nstr   = 16;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2;
  double *lamc  = array1D(1,1), dlam;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);

  double *rfldir = array1D(ntau,1);// ......... Direct-beam flux (without delta-m scaling)
  double *rfldn = array1D(ntau,1);// .......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *flup = array1D(ntau,1);// ........... Diffuse up-flux
  double *fldn = array1D(ntau,1);// ........... Diffuse down-flux
  double *fldir = array1D(ntau,1);// .......... Direct flux
  double *dfdt = array1D(ntau,1);// ........... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****uu = array4D(nlam,numu,ntau,nphi,1);//.  Corrected intensity field

  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field
  double gg; char title[300];

  for (icas = 1; icas <= 4; icas++) {
  //for (icas = 1; icas <= 1; icas++) {
    sprintf(title,"Test Case No. 7-%02d",icas);
    switch(icas) {
      case 1:
        nstr = 16;
        nmom = 16;
        nphi = 1;
        ntau = 2;
        numu = 2;

        dtau[0][0][0] = 1.;

        gg = 0.05;
        getphase(HENYEY_GREENSTEIN, gg, 0.1, nmom, 0, phfunc);

        temper[0]   = 200.0;
        temper[1]   = 300.0;

        lamc[0]     = 1e4/800.0;
        dlam        = 1e4/300.0 - 1e4/800.0;

        utau[0]     =   0.0;
        utau[1]     =   1.;

        umu[0]      =  -1.;
        umu[1]      =   1.;

        phi[0]      =   0.0;

        umu0        =   0.5;
        fstar[0]    =   0.0;
        fisot       =   0.0;
        tsrf[0]     =   0.0;
        esrf[0]     =   1.0 - tsrf[0];
        ttemp       =   0.0;
        btemp       =   0.0;
        temis       =   1.;

        if (VERBOSE) printf("\n\n%s: Absorption + Scattering, Internal Thermal Sources; Ref. KS, Table I, tau = %3.1f, a = %3.1f, g = %4.2f\n", title,dtau[0][0][0],phfunc[0][0][0][0],gg);

        // Correct answers
        grfldir[0]= 0.,         grfldir[1]= 0.,
        grfldn[0] = 0.,         grfldn[1] = 1.21204E+02,
        gflup[0]  = 8.62936E+01,gflup[1]  = 0.,
        gdfdt[0]  =-5.13731E+01,gdfdt[1]  =-5.41036E+02;
      break;
      case 2:

        ntau = 2;
        numu = 2;
        nphi = 1;

        dtau[0][0][0] = 100.0;

        gg = 0.75;
        getphase(HENYEY_GREENSTEIN, gg, 0.95, nmom, 0, phfunc);

        temper[0]   =  200.0;
        temper[1]   =  300.0;

        lamc[0]     = 1e4/2703.01;
        dlam        = 1e4/2702.99 - 1e4/2703.01;

        utau[0]     = 0.0;
        utau[1]     = 100.0;

        umu[0]      = -1.;
        umu[1]      =  1.;

        phi[0]      = 0.0;

        umu0        = 0.5;
        fstar[0]    = 0.0;
        fisot       = 0.0;
        tsrf[0]     = 0.0;
        esrf[0]     = 1.0 - tsrf[0];
        ttemp       = 0.0;
        btemp       = 0.0;
        temis       = 1.;

        if (VERBOSE) printf("\n\n%s: Absorption + Scattering, Internal Thermal Sources; Ref. KS, Table II, tau = %5.1f, a = %4.2f, g = %4.2f\n", title,dtau[0][0][0],phfunc[0][0][0][0],gg);

        // Correct answers
        grfldir[0]=0.,         grfldir[1]= 0.,
        grfldn[0] =0.,         grfldn[1] = 2.07786E-05,
        gflup[0]  =1.10949E-06,gflup[1]  = 0.,
        gdfdt[0]  =8.23219E-08,gdfdt[1]  =-5.06461E-06,
        guu[0][0][0][0]=0.,         guu[0][1][0][0]=4.65744E-07,
        guu[0][0][1][0]=7.52311E-06,guu[0][1][1][0]=0.0;
      break;
      case 3:

        nstr   =  12;
        nmom   =  12;
        ntau   =   3;
        numu   =   4;
        nphi   =   2;

        temper[0] = 300.0;
        temper[1] = 200.0;

        gg = 0.8;
        getphase(HENYEY_GREENSTEIN, gg, 0.5, nmom, 0, phfunc);

        dtau[0][0][0]    =     1.;

        lamc[0]     = 1e4/50000.0;
        dlam        = 1e4/1e-20 - 1e4/50000.0;

        utau[0]     =     0.0;
        utau[1]     =     0.5;
        utau[2]     =     1.;

        umu[0]      =    -1.;
        umu[1]      =    -0.1;
        umu[2]      =     0.1;
        umu[3]      =     1.;

        phi[0]      =     0.0;
        phi[1]      =    90.0;

        fstar[0]    =   200.0/M_PI;
        umu0        =     0.5;
        fisot       =   100.0;
        ttemp       =   100.0;
        btemp       =   320.0;
        temis       =     1.;
        tsrf[0]     =     0.0;
        esrf[0]     = 1.0 - tsrf[0];

        if (VERBOSE) printf("\n\n%s: Absorption + Henyey-Greenstein Scattering, All Sources, Bottom Albedo = %.0f\n",title,tsrf[0]);

        // Correct answers
        grfldir[0]= 100.,       grfldir[1]=3.67879E+01,grfldir[2]=1.35335E+01;
        grfldn[0] = 3.19830E+02,grfldn[1] =3.54099E+02,grfldn[2] =3.01334E+02;
        gflup[0]  = 4.29572E+02,gflup[1]  =4.47018E+02,gflup[2]  =5.94576E+02;
        gdfdt[0]  =-8.04270E+01,gdfdt[1]  =2.51589E+02,gdfdt[2]  =7.15964E+02;
        guu[0][0][0][0]=1.01805E+02,guu[0][1][0][0]=1.01805E+02,guu[0][2][0][0]=1.46775E+02,guu[0][3][0][0]=1.49033E+02,
        guu[0][0][1][0]=1.06583E+02,guu[0][1][1][0]=1.28565E+02,guu[0][2][1][0]=1.04464E+02,guu[0][3][1][0]=1.59054E+02,
        guu[0][0][2][0]=9.66519E+01,guu[0][1][2][0]=8.65854E+01,guu[0][2][2][0]=1.89259E+02,guu[0][3][2][0]=1.89259E+02,
        guu[0][0][0][1]=1.01805E+02,guu[0][1][0][1]=1.01805E+02,guu[0][2][0][1]=1.29641E+02,guu[0][3][0][1]=1.49033E+02,
        guu[0][0][1][1]=1.06583E+02,guu[0][1][1][1]=1.06408E+02,guu[0][2][1][1]=9.48418E+01,guu[0][3][1][1]=1.59054E+02,
        guu[0][0][2][1]=9.66519E+01,guu[0][1][2][1]=7.49310E+01,guu[0][2][2][1]=1.89259E+02,guu[0][3][2][1]=1.89259E+02;
      break;
      case 4:

        temper[0] = 300.0;
        temper[1] = 200.0;

        gg = 0.8;
        getphase(HENYEY_GREENSTEIN, gg, 0.5, nmom, 0, phfunc);

        dtau[0][0][0]    =     1.;

        lamc[0]     = 1e4/50000.0;
        dlam        = 1e4/1e-20 - 1e4/50000.0;

        utau[0]     =     0.0;
        utau[1]     =     0.5;
        utau[2]     =     1.;

        umu[0]      =    -1.;
        umu[1]      =    -0.1;
        umu[2]      =     0.1;
        umu[3]      =     1.;

        phi[0]      =     0.0;
        phi[1]      =    90.0;

        tsrf[0]     = 1.0;
        esrf[0]     = 1.0 - tsrf[0];

        if (VERBOSE) printf("\n\n%s: Absorption + Henyey-Greenstein Scattering, All Sources, Bottom Albedo = %.0f\n",title,tsrf[0]);

        // Correct answers
        grfldir[0]= 100.,       grfldir[1]=3.67879E+01,grfldir[2]=1.35335E+01;
        grfldn[0] = 3.19830E+02,grfldn[1] =3.50555E+02,grfldn[2] =2.92063E+02;
        gflup[0]  = 3.12563E+02,gflup[1]  =2.68126E+02,gflup[2]  =3.05596E+02;
        gdfdt[0]  =-1.68356E+02,gdfdt[1]  =1.01251E+02,gdfdt[2]  =4.09326E+02;
        guu[0][0][0][0]=1.01805E+02,guu[0][1][0][0]=1.01805E+02,guu[0][2][0][0]=1.40977E+02,guu[0][3][0][0]=9.62764E+01,
        guu[0][0][1][0]=1.06203E+02,guu[0][1][1][0]=1.23126E+02,guu[0][2][1][0]=9.19545E+01,guu[0][3][1][0]=8.89528E+01,
        guu[0][0][2][0]=9.56010E+01,guu[0][1][2][0]=7.25576E+01,guu[0][2][2][0]=9.72743E+01,guu[0][3][2][0]=9.72743E+01,
        guu[0][0][0][1]=1.01805E+02,guu[0][1][0][1]=1.01805E+02,guu[0][2][0][1]=1.23843E+02,guu[0][3][0][1]=9.62764E+01,
        guu[0][0][1][1]=1.06203E+02,guu[0][1][1][1]=1.00969E+02,guu[0][2][1][1]=8.23318E+01,guu[0][3][1][1]=8.89528E+01,
        guu[0][0][2][1]=9.56010E+01,guu[0][1][2][1]=6.09031E+01,guu[0][2][2][1]=9.72743E+01,guu[0][3][2][1]=9.72743E+01;
      break;
    }

    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
          &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
  return;
}


//========================== psgtest08() ==============================
//---------------------------------------------------------------------
//---  Test Problem 8:  Absorbing/Isotropic-Scattering Medium      ----
//---  With Two Computational Layers                               ----
//--- (Compare Fluxes To Ref. OS, Table 1)                         ----
//---------------------------------------------------------------------
void psgtest08(void) {

  int nstr =  8;
  int nlyr =  2;
  int ntau =  3;
  int numu =  4;
  int nphi =  1;
  int nmom = nstr;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, lc, nmax=nstr/2;
  double *lamc  = array1D(1,1), dlam;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu = NULL;

  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field
  double gg; char title[300];

  umu[0]    = -1.;
  umu[1]    = -0.2;
  umu[2]    =  0.2;
  umu[3]    =  1.;

  phi[0]    = 60.0;

  fstar[0]  =  0.0;
  fisot     =  1./M_PI;
  tsrf[0]   =  0.0;
  esrf[0]   =  1.0 - tsrf[0];
  umu0      = 0.5;

  for (icas = 1; icas <= 3; icas++) {
    switch(icas) {
      case 1:
        dtau[0][0][0] = .25;
        dtau[1][0][0] = .25;
        getphase(ISOTROPIC, 0.0, 0.5, nmom, 0, phfunc);
        getphase(ISOTROPIC, 0.0, 0.3, nmom, 1, phfunc);
        utau[0]  = .0;
        utau[1]  = .25;
        utau[2]  = .5;
        if (VERBOSE) printf("\n\nTest Case No. 8a:  Ref. OS, Table 1, Line 4 (Two Inhomogeneous Layers)\n");

        // Correct answers
        grfldir[0]=0.,         grfldir[1]=0.,         grfldir[2]=0.0;
        grfldn[0] =1.,         grfldn[1] =7.22235E-01,grfldn[2] =5.13132E-01;
        gflup[0]  =9.29633E-02,gflup[1]  =2.78952E-02,gflup[2]  =0.0;
        gdfdt[0]  =1.12474E+00,gdfdt[1]  =6.51821E-01,gdfdt[2]  =5.63361E-01;
        guu[0][0][0][0]=3.18310E-01,guu[0][1][0][0]=3.18310E-01,guu[0][2][0][0]=5.62566E-02,guu[0][3][0][0]=1.94423E-02,
        guu[0][0][1][0]=2.62711E-01,guu[0][1][1][0]=1.36952E-01,guu[0][2][1][0]=1.84909E-02,guu[0][3][1][0]=5.52188E-03,
        guu[0][0][2][0]=2.10014E-01,guu[0][1][2][0]=5.60376E-02,guu[0][2][2][0]=0.,         guu[0][3][2][0]=0.0;
      break;
      case 2:
        dtau[0][0][0] = .25;
        dtau[1][0][0] = .25;
        getphase(ISOTROPIC, 0.0, 0.80, nmom, 0, phfunc);
        getphase(ISOTROPIC, 0.0, 0.95, nmom, 1, phfunc);
        utau[0]  = .0;
        utau[1]  = .25;
        utau[2]  = .5;
        if (VERBOSE) printf("\n\nTest Case No. 8b:  Ref. OS, Table 1, Line 1 (Two Inhomogeneous Layers)\n");

        // Correct answers
        grfldir[0]=0.,         grfldir[1]=0.,         grfldir[2]=0.0;
        grfldn[0] =1.,         grfldn[1] =7.95332E-01,grfldn[2] =6.50417E-01;
        gflup[0]  =2.25136E-01,gflup[1]  =1.26349E-01,gflup[2]  =0.0;
        gdfdt[0]  =5.12692E-01,gdfdt[1]  =3.56655E-01,gdfdt[2]  =5.68095E-02;
        guu[0][0][0][0]=3.18310E-01,guu[0][1][0][0]=3.18310E-01,guu[0][2][0][0]=1.23687E-01,guu[0][3][0][0]=4.95581E-02,
        guu[0][0][1][0]=2.77499E-01,guu[0][1][1][0]=1.83950E-01,guu[0][2][1][0]=8.35695E-02,guu[0][3][1][0]=2.50575E-02,
        guu[0][0][2][0]=2.40731E-01,guu[0][1][2][0]=1.29291E-01,guu[0][2][2][0]=0.,         guu[0][3][2][0]=0.0;
      break;
      case 3:
        dtau[0][0][0] = 1.;
        dtau[1][0][0] = 2.;
        getphase(ISOTROPIC, 0.0, 0.80, nmom, 0, phfunc);
        getphase(ISOTROPIC, 0.0, 0.95, nmom, 1, phfunc);
        utau[0]  = 0.0;
        utau[1]  = 1.;
        utau[2]  = 3.;
        if (VERBOSE) printf("\n\nTest Case No. 8c:  Ref. OS, Table 1, Line 13 (Two Inhomogeneous Layers)\n");

        // Correct answers
        grfldir[0]=0.,         grfldir[1]=0.,         grfldir[2]=0.0;
        grfldn[0] =1.,         grfldn[1] =4.86157E-01,grfldn[2] =1.59984E-01;
        gflup[0]  =3.78578E-01,gflup[1]  =2.43397E-01,gflup[2]  =0.0;
        gdfdt[0]  =5.65095E-01,gdfdt[1]  =2.76697E-01,gdfdt[2]  =1.35679E-02;
        guu[0][0][0][0]=3.18310E-01,guu[0][1][0][0]=3.18310E-01,guu[0][2][0][0]=1.49335E-01,guu[0][3][0][0]=1.04766E-01,
        guu[0][0][1][0]=1.89020E-01,guu[0][1][1][0]=9.88158E-02,guu[0][2][1][0]=9.65192E-02,guu[0][3][1][0]=6.54445E-02,
        guu[0][0][2][0]=6.84762E-02,guu[0][1][2][0]=2.96698E-02,guu[0][2][2][0]=0.,         guu[0][3][2][0]=0.0;
      break;
    }

    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
          &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
  return;
}


//========================== psgtest09() ==============================
//---------------------------------------------------------------------
//---  Test Problem 9:  General Emitting/Absorbing/Scattering      ----
//---  Medium with Every Computational Layer Different.            ----
//--- (Compare 9a,b Fluxes to Ref. DGIS, Tables VI-VII, beta = 0)  ----
//---------------------------------------------------------------------
void psgtest09(void) {

  int nstr =  8;
  int nlyr =  6;
  int nmom =  8;
  int ntau =  5;
  int numu =  4;
  int nphi =  3;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, k, lc;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0, gg, albedo;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=array4D(nlam,numu,ntau,nphi,1);
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field

  fstar[0] = 0.0;
  fisot = 1./M_PI;
  umu0  = 0.5;

  for (lc = 1; lc <= nlyr; lc++) dtau[lc-1][0][0] = (double)lc;

  utau[0]  = 0.0;
  utau[1]  = 1.05;
  utau[2]  = 2.1;
  utau[3]  = 6.;
  utau[4]  = 21.;

  umu[0]   = -1.;
  umu[1]   = -0.2;
  umu[2]   =  0.2;
  umu[3]   =  1.;

  for (icas = 1; icas <= 3; icas++) {
    switch(icas) {
      case 1:

        nphi     = 1;
        phi[0]   = 60.0;

        for (lc = 1; lc <= nlyr; lc++) {
          albedo = 0.6+(double)lc*0.05;
          getphase(ISOTROPIC, 0.0, albedo, nmom, lc-1, phfunc);
        }

        tsrf[0]  = 0.0;
        esrf[0]  = 1.0 - tsrf[0];
        if (VERBOSE) printf("\n\nTest Case No. 9a:  Ref. DGIS, Tables VI-VII, beta=l=0 (multiple inhomogeneous layers)\n");

        // Correct answers
        grfldir[0]=0.,         grfldir[1]=0.,         grfldir[2]=0.,         grfldir[3]=0.,         grfldir[4]=0.0;
        grfldn[0] =1.,         grfldn[1] =3.55151E-01,grfldn[2] =1.44265E-01,grfldn[3] =6.71445E-03,grfldn[4] =6.16968E-07;
        gflup[0]  =2.27973E-01,gflup[1]  =8.75098E-02,gflup[2]  =3.61819E-02,gflup[3]  =2.19291E-03,gflup[4]  =0.0;
        gdfdt[0]  =8.82116E-01,gdfdt[1]  =2.32366E-01,gdfdt[2]  =9.33443E-02,gdfdt[3]  =3.92782E-03,gdfdt[4]  =1.02500E-07;
        guu[0][0][0][0]=3.18310E-01,guu[0][1][0][0]=3.18310E-01,guu[0][2][0][0]=9.98915E-02,guu[0][3][0][0]=5.91345E-02,
        guu[0][0][1][0]=1.53507E-01,guu[0][1][1][0]=5.09531E-02,guu[0][2][1][0]=3.67006E-02,guu[0][3][1][0]=2.31903E-02,
        guu[0][0][2][0]=7.06614E-02,guu[0][1][2][0]=2.09119E-02,guu[0][2][2][0]=1.48545E-02,guu[0][3][2][0]=9.72307E-03,
        guu[0][0][3][0]=3.72784E-03,guu[0][1][3][0]=1.08815E-03,guu[0][2][3][0]=8.83316E-04,guu[0][3][3][0]=5.94743E-04,
        guu[0][0][4][0]=2.87656E-07,guu[0][1][4][0]=1.05921E-07,guu[0][2][4][0]=0.,         guu[0][3][4][0]=0.0;
      break;
      case 2:

        utau[0]  = 0.0;
        utau[1]  = 1.05;
        utau[2]  = 2.1;
        utau[3]  = 6.;
        utau[4]  = 21.;

        umu[0]   = -1.;
        umu[1]   = -0.2;
        umu[2]   =  0.2;
        umu[3]   =  1.;

        nphi     = 1;
        phi[0]   = 60.0;

        double pmom[9][2];
        pmom[0][0] = 1.0;
        pmom[1][0] = 2.00916/3.;
        pmom[2][0] = 1.56339/5.;
        pmom[3][0] = 0.67407/7.;
        pmom[4][0] = 0.22215/9.;
        pmom[5][0] = 0.04725/11.;
        pmom[6][0] = 0.00671/13.;
        pmom[7][0] = 0.00068/15.;
        pmom[8][0] = 0.00005/17.;
        for (lc = 1; lc <= nlyr; lc++) {
          albedo = 0.6+(double)lc*0.05;
          for (k = 0; k <= nmom; k++) {
            phfunc[k][lc-1][0][0] = pmom[k][0]*albedo;
          }
        }

        if (VERBOSE) printf("\n\nTest Case No. 9b:  Ref. DGIS, Tables VI-VII, beta=0,l=%1d (multiple inhomogeneous layers)\n",nmom);

        // Correct answers
        grfldir[0]=0.,         grfldir[1]=0.,         grfldir[2]=0.,         grfldir[3]=0.,         grfldir[4]=0.0;
        grfldn[0] =1.,         grfldn[1] =4.52357E-01,grfldn[2] =2.36473E-01,grfldn[3] =2.76475E-02,grfldn[4] =7.41853E-05;
        gflup[0]  =1.00079E-01,gflup[1]  =4.52014E-02,gflup[2]  =2.41941E-02,gflup[3]  =4.16016E-03,gflup[4]  =0.0;
        gdfdt[0]  =8.04577E-01,gdfdt[1]  =2.55330E-01,gdfdt[2]  =1.30976E-01,gdfdt[3]  =1.36227E-02,gdfdt[4]  =1.22022E-05;
        guu[0][0][0][0]=3.18310E-01,guu[0][1][0][0]=3.18310E-01,guu[0][2][0][0]=7.39198E-02,guu[0][3][0][0]=1.32768E-02,
        guu[0][0][1][0]=1.96609E-01,guu[0][1][1][0]=5.92369E-02,guu[0][2][1][0]=3.00230E-02,guu[0][3][1][0]=7.05566E-03,
        guu[0][0][2][0]=1.15478E-01,guu[0][1][2][0]=3.01809E-02,guu[0][2][2][0]=1.52672E-02,guu[0][3][2][0]=4.06932E-03,
        guu[0][0][3][0]=1.46177E-02,guu[0][1][3][0]=3.85590E-03,guu[0][2][3][0]=2.38301E-03,guu[0][3][3][0]=7.77890E-04,
        guu[0][0][4][0]=3.37742E-05,guu[0][1][4][0]=1.20858E-05,guu[0][2][4][0]=0.,         guu[0][3][4][0]=0.0;
      break;
      case 3:

        temper[0] = 600.0;
        for (lc = 1; lc <= nlyr; lc++) {
          gg = (double)lc/7.;
          albedo = 0.6+(double)lc*0.05;
          getphase(HENYEY_GREENSTEIN, gg, albedo, nmom, lc-1, phfunc);
          temper[lc] = 600.+(double)lc*10.0;
        }

        nphi       =     3;
        phi[0]     =  60.0;
        phi[1]     = 120.0;
        phi[2]     = 180.0;

        lamc[0]    = 1e4/1000.0;
        dlam       = 1e4/999.0 - 1e4/1000.0;

        fstar[0] = 1.0;
        fisot    = 1.0;
        tsrf[0]  = 0.5;
        esrf[0]  = 1.0 - tsrf[0];
        btemp    = 700.0;
        ttemp    = 550.0;
        temis    = 1.0;
        if (VERBOSE) printf("\n\nTest Case No. 9c:  Generalization of 9A to include all possible complexity\n");

        // Correct answers
        grfldir[0]=1.57080E+00,grfldir[1]=1.92354E-01,grfldir[2]=2.35550E-02;grfldir[3]=9.65131E-06,grfldir[4]=9.03133E-19;
        grfldn[0] =6.09217E+00,grfldn[1] =4.97279E+00,grfldn[2] =4.46616E+00;grfldn[3] =4.22731E+00,grfldn[4] =4.73767E+00;
        gflup[0]  =4.68414E+00,gflup[1]  =4.24381E+00,gflup[2]  =4.16941E+00;gflup[3]  =4.30667E+00,gflup[4]  =5.11524E+00;
        gdfdt[0]  =3.49563E+00,gdfdt[1]  =8.81206E-01,gdfdt[2]  =3.50053E-01;gdfdt[3]  =1.93471E-02,gdfdt[4]  =7.15349E-02;
        guu[0][0][0][0]=1.93920,guu[0][1][0][0]=1.93920,guu[0][2][0][0]=1.61855,guu[0][3][0][0]=1.43872,
        guu[0][0][1][0]=1.66764,guu[0][1][1][0]=1.44453,guu[0][2][1][0]=1.38339,guu[0][3][1][0]=1.33890,
        guu[0][0][2][0]=1.48511,guu[0][1][2][0]=1.35009,guu[0][2][2][0]=1.33079,guu[0][3][2][0]=1.32794,
        guu[0][0][3][0]=1.34514,guu[0][1][3][0]=1.35131,guu[0][2][3][0]=1.35980,guu[0][3][3][0]=1.37918,
        guu[0][0][4][0]=1.48927,guu[0][1][4][0]=1.54270,guu[0][2][4][0]=1.62823,guu[0][3][4][0]=1.62823,
        guu[0][0][0][1]=1.93920,guu[0][1][0][1]=1.93920,guu[0][2][0][1]=1.57895,guu[0][3][0][1]=1.43872,
        guu[0][0][1][1]=1.66764,guu[0][1][1][1]=1.42925,guu[0][2][1][1]=1.37317,guu[0][3][1][1]=1.33890,
        guu[0][0][2][1]=1.48511,guu[0][1][2][1]=1.34587,guu[0][2][2][1]=1.32921,guu[0][3][2][1]=1.32794,
        guu[0][0][3][1]=1.34514,guu[0][1][3][1]=1.35129,guu[0][2][3][1]=1.35979,guu[0][3][3][1]=1.37918,
        guu[0][0][4][1]=1.48927,guu[0][1][4][1]=1.54270,guu[0][2][4][1]=1.62823,guu[0][3][4][1]=1.62823,
        guu[0][0][0][2]=1.93920,guu[0][1][0][2]=1.93920,guu[0][2][0][2]=1.56559,guu[0][3][0][2]=1.43872,
        guu[0][0][1][2]=1.66764,guu[0][1][1][2]=1.42444,guu[0][2][1][2]=1.37034,guu[0][3][1][2]=1.33890,
        guu[0][0][2][2]=1.48511,guu[0][1][2][2]=1.34469,guu[0][2][2][2]=1.32873,guu[0][3][2][2]=1.32794,
        guu[0][0][3][2]=1.34514,guu[0][1][3][2]=1.35128,guu[0][2][3][2]=1.35979,guu[0][3][3][2]=1.37918,
        guu[0][0][4][2]=1.48927,guu[0][1][4][2]=1.54270,guu[0][2][4][2]=1.62823,guu[0][3][4][2]=1.62823;
      break;
    }

    psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
          &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

    testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, grfldir, grfldn, gflup, gdfdt, guu);
  }
  return;
}

//========================== psgtest10() ==============================
//---------------------------------------------------------------------
//---  Test Problem 10: Compare ds.flag.usrang = TRUE vs. FALSE    ----
//---  Take Problem 9c (our most general case) but only 4 Streams  ----
//---------------------------------------------------------------------
void psgtest10(void) {

  int nstr =  4;
  int nlyr =  6;
  int nmom =  nstr;
  int ntau =  3;
  int numu =  4;
  int nphi =  2;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, k, lc;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0, gg, albedo;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=NULL;
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field

  fstar[0]  = 1.0;
  umu0      = 0.5;
  fisot     = 1.0;
  tsrf[0]   = 0.5;
  esrf[0]   = 1.0 - tsrf[0];
  btemp     = 700.0;
  ttemp     = 550.0;
  temis     = 1.0;

  lamc[0] = 1e4/1000.0;
  dlam    = 1e4/999.0 - 1e4/1000.0;

  temper[0] = 600.0;
  for (lc = 1; lc <= nlyr; lc++) {
    dtau[lc-1][0][0] = (double)lc;
    albedo = 0.6+(double)lc*.05;
    gg = (double)lc/(nlyr+1);
    getphase(HENYEY_GREENSTEIN, gg, albedo, nmom, lc-1, phfunc);
    temper[lc] = 600.+(double)lc*10.0;
  }

  utau[0] =  0.0;
  utau[1] =  2.1;
  utau[2] = 21.;

  phi[0] =  60.0;
  phi[1] = 120.0;

  umu[0] = -0.788675129;
  umu[1] = -0.211324871;
  umu[2] =  0.211324871;
  umu[3] =  0.788675129;

  if (VERBOSE) printf("\n\nTest Case No. 10a:  like 9c, ds.flag.usrang = TRUE\n");

  psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

  // Case 2
  double *dfdt2=NULL, *fldir2=NULL, *fldn2=NULL, *flup2=NULL, *rfldir2=NULL, *rfldn2=NULL, ****uu2=NULL;
  numu = 0;

  if (VERBOSE) printf("\n\nTest Case No. 10a:  like 9c, ds.flag.usrang = FALSE\n");

  psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt2, &fldir2, &fldn2, &flup2, &rfldir2, &rfldn2, &uu2);

  numu = nstr;
  testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, rfldir2, rfldn2, flup2, dfdt2, uu2);

  return;
}


//========================== psgtest11() ==============================
//---------------------------------------------------------------------
//---  Test Problem 11: Single-Layer vs. Multiple Layers           ----
//---  11a: Results at user levels for one computational layer     ----
//---  11b: Single layer of 11a subdivided into multiple           ----
//---       computational layers at the 11a user levels            ----
//---------------------------------------------------------------------
void psgtest11(void) {

  int nstr =  16;
  int nlyr =  3;
  int nmom =  nstr;
  int ntau =  4;
  int numu =  4;
  int nphi =  2;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, k, lc;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0, gg;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=NULL;
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field

  fstar[0]  = 1.0;
  umu0      = 0.5;
  fisot     = 0.5/M_PI;
  tsrf[0]   = 0.5;
  esrf[0]   = 1.0 - tsrf[0];

  nlyr = 1;
  dtau[0][0][0] = 1.0;
  getphase(ISOTROPIC, 0.0, 0.9, nmom, 0, phfunc);

  utau[0] =  0.00;
  utau[1] =  0.05;
  utau[2] =  0.50;
  utau[3] =  1.00;

  phi[0] =  0.0;
  phi[1] = 90.0;

  umu[0] = -1.0;
  umu[1] = -0.1;
  umu[2] =  0.1;
  umu[3] =  1.0;

  if (VERBOSE) printf("\n\nTest Case No. 11a: One Isotropic-Scattering Layer\n");

  psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

  // Case 2
  double *dfdt2=NULL, *fldir2=NULL, *fldn2=NULL, *flup2=NULL, *rfldir2=NULL, *rfldn2=NULL, ****uu2=NULL;
  ntau = 0;
  nlyr = 3;
  for (lc = 1; lc <= nlyr; lc++) {
    dtau[lc-1][0][0] = utau[lc]-utau[lc-1];
    getphase(ISOTROPIC, 0.0, 0.9, nmom, lc-1, phfunc);
  }

  if (VERBOSE) printf("\n\nTest Case No. 11b: Same as 11a but treated as multiple layers\n");

  psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt2, &fldir2, &fldn2, &flup2, &rfldir2, &rfldn2, &uu2);

  ntau = nlyr+1;
  testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, rfldir2, rfldn2, flup2, dfdt2, uu2);

  return;
}


//=========================== psgtest12() ==============================
//----------------------------------------------------------------------
//---  Test Problem 12: Test Absorption-Optical-Depth Shortcut      ----
//---  compares cases where the DISORT shortcut for absorption      ----
//---  optical depth > 10 is not used (12a), then is used (12b)     ----
//---  (this shortcut is only employed when ds.flag.planck = FALSE) ----
//----------------------------------------------------------------------
void psgtest12(void) {

  int nstr =  20;
  int nlyr =  4;
  int nmom =  nstr;
  int ntau =  4;
  int numu =  4;
  int nphi =  1;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, k, lc;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0, gg;
  double *alts = array1D(nlyr+1,1); for (i=0;i<=nlyr;i++) alts[i]=nlyr-i;
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double ****pmoms = array4D(nlam, nbin, nlyr, mmom+1, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=NULL;
  double *grfldir = array1D(ntau,1);// ........ Direct-beam flux (without delta-m scaling)
  double *grfldn = array1D(ntau,1);// ......... Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  double *gflup = array1D(ntau,1);// .......... Diffuse up-flux
  double *gdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  double ****guu = array4D(nlam,numu,ntau,nphi,1);//  Corrected intensity field

  nlyr = 1;
  ntau = 4;

  utau[0] =  0.0;
  utau[1] = 10.0;
  utau[2] = 19.9;
  utau[3] = 20.1;

  umu[0] = -1.0;
  umu[1] = -0.1;
  umu[2] =  0.1;
  umu[3] =  1.0;

  phi[0] =  0.0;

  dtau[0][0][0] = 20.1;
  getphase(HENYEY_GREENSTEIN, 0.9, 0.5, nmom, 0, phfunc);

  fstar[0]  = 1.0;
  umu0      = 1.0;
  fisot     = 0.0;
  tsrf[0]   = 1.0;
  esrf[0]   = 1.0 - tsrf[0];

  if (VERBOSE) printf("\n\nTest Case No. 12a:  Overhead Beam Striking Absorbing/Scattering Medium\n");

  psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

  // Case 2
  double *dfdt2=NULL, *fldir2=NULL, *fldn2=NULL, *flup2=NULL, *rfldir2=NULL, *rfldn2=NULL, ****uu2=NULL;
  nlyr = ntau-1;
  for (lc = 1; lc <= nlyr; lc++) {
    dtau[lc-1][0][0] = utau[lc]-utau[lc-1];
    getphase(HENYEY_GREENSTEIN, 0.9, 0.5, nmom, lc-1, phfunc);
  }

  if (VERBOSE) printf("\n\nTest Case No. 12b: Same as 12a but uses shortcut for absorption optical depth > 10\n");

  psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, numu, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt2, &fldir2, &fldn2, &flup2, &rfldir2, &rfldn2, &uu2);

  testprint(nphi, ntau, numu, phi, utau, umu, rfldir, rfldn, flup, dfdt, uu, rfldir2, rfldn2, flup2, dfdt2, uu2);
}


//========================== psgtest14() ============================
//----------------------------------------------------------------------
//---  Test Problem 14: Like Test 10b, but compare twostr() to     -----
//---  disort() with 4 streams.                                    -----
//----------------------------------------------------------------------
void psgtest14(void)
{

  int nstr =  2;
  int nlyr =  6;
  int nmom =  nstr;
  int ntau =  3;
  int numu =  nstr;
  int nphi =  2;
  int mmom   = (nmom>nstr ? nmom : nstr);

  double umu0;
  double radius=1e6;
  double *phi = array1D(nphi,1);
  double *umu  = array1D(numu,1);
  double *utau = array1D(ntau,1);
  long nlam=1, l1=0, l2=nlam; int nbin=1, i, icas, nmax=nstr/2, k, lc;
  double *lamc  = array1D(1,1), dlam=0.0;
  double *fstar = array1D(1,1);
  double *tsrf  = array1D(1,1);
  double *esrf  = array1D(1,1);
  double fisot=0.0, btemp=0.0, ttemp=0.0, temis=0.0, gg, alb;
  double *alts = array1D(nlyr+1,1);
  double *temper = array1D(nlyr+1,1);
  double ***dtau = array3D(nlyr, nbin, nlam, 1);
  double ****phfunc = array4D(mmom+1, nlyr, nbin, nlam, 1);
  double ****pmoms = array4D(nlam, nbin, nlyr, mmom+1, 1);
  double *dfdt=NULL, *fldir=NULL, *fldn=NULL, *flup=NULL, *rfldir=NULL, *rfldn=NULL, ****uu=NULL;

  fstar[0]  = 1.0;
  umu0      = 0.5;
  fisot     = 1.0;
  tsrf[0]   = 0.5;
  esrf[0]   = 1.0 - tsrf[0];
  btemp     = 700.0;
  ttemp     = 550.0;
  temis     = 1.0;
  lamc[0]   = 1e4/1000.0;
  dlam      = 1e4/999.0 - 1e4/1000.0;

  temper[0] = 600.0;
  for (lc = 0; lc < nlyr; lc++) {
    dtau[lc][0][0] = lc+1.0;
    gg = (double)(lc+1)/(nlyr+1.0);
    alb = 0.6 + (double)(lc+1.0)*0.05;
    getphase(HENYEY_GREENSTEIN, gg, alb, nmom, lc, phfunc);
    temper[lc+1] = 600.0 + (double)(lc+1)*10.0;
    alts[lc] = 10.*(double)(nlyr-lc+1);
  } alts[lc] = 10.0;

  utau[0] =  0.0;
  utau[1] =  2.1;
  utau[2] = 21.0;

  phi[0] =  60.0;
  phi[1] = 120.0;

  if (VERBOSE) printf("\n\nTest Case No. 14a: disort() as in 10b\n");

  psgdort(0, nlam, 1, nbin, nlyr, nmax, nmom, 0, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt, &fldir, &fldn, &flup, &rfldir, &rfldn, &uu);

  // Case 2 - Two streams approximation
  double *dfdt2=NULL, *fldir2=NULL, *fldn2=NULL, *flup2=NULL, *rfldir2=NULL, *rfldn2=NULL, ****uu2=NULL;
  //radius = 6378.0;

  psgdort(0, nlam, 1, nbin, nlyr, 0, nmom, 0, nphi, ntau, -1, umu0, phi, radius, alts, dtau, phfunc, fstar, tsrf, esrf, temis, NULL, fisot, lamc, dlam, temper, btemp, ttemp,
        &umu, &utau, &dfdt2, &fldir2, &fldn2, &flup2, &rfldir2, &rfldn2, &uu2);

  testprint(nphi, ntau, 2, phi, utau, umu, rfldir2, rfldn2, flup2, dfdt2, uu2, rfldir, rfldn, flup, dfdt, uu);

  return;
}


// ---------------------------------------------------------------------------------------
// Print PSGDORT results and, directly beneath them, their ratios to the correct answers
// Print number of non-unit ratios that occur but try
// to count just the cases where there is a real disagreement and not
// those where flux or intensity are down at their noise level (defined as
// 10^(-6) times their maximum value).  d(flux)/d(tau) is treated the
// same as fluxes in this noise estimation even though it is a different
// type of quantity (although with flux units).
// ---------------------------------------------------------------------------------------
void testprint(int nphi, int ntau, int numu, double *phi, double *utau, double *umu, double *rfldir, double *rfldn, double *flup, double *dfdt, double ****uu,
               double *grfldir, double *grfldn, double *gflup, double *gdfdt, double ****guu) {

  int lu, iu, j, numbad=0, onlyfl=(numu<=0);
  double flxmax=0.0, fnoise, umax=0.0, unoise, rat1, rat2, rat3, rat4, ratv[nphi];

  // Define the significance of the differences
  for (lu = 0; lu < ntau; lu++) {
    if (grfldir[lu]>flxmax) flxmax = grfldir[lu];
    if (grfldn[lu]>flxmax) flxmax = grfldn[lu];
    if (gflup[lu]>flxmax) flxmax = gflup[lu];
  } fnoise = 1.e-6*flxmax;

  if (VERBOSE) printf("\n\n                  <-------------- FLUXES -------------->\n"
         "      Optical       Downward       Downward         Upward    d(Net Flux)\n"
         "        Depth         Direct        Diffuse        Diffuse    / d(Op Dep)\n");

  for (lu = 0; lu < ntau; lu++) {
    if (VERBOSE) printf("NEW %9.4f%15.4e%15.4e%15.4e%15.4e\n",
                   utau[lu],rfldir[lu],rfldn[lu],flup[lu],dfdt[lu]);

    if (VERBOSE) printf("%13.4f%15.4e%15.4e%15.4e%15.4e\n",
                   utau[lu],grfldir[lu],grfldn[lu],gflup[lu],gdfdt[lu]);

    rat1 = (rfldir[lu]+fnoise)/(grfldir[lu]+fnoise);
    rat2 = (rfldn[lu]+fnoise)/(grfldn[lu]+fnoise);
    rat3 = (flup[lu]+fnoise)/(gflup[lu]+fnoise);
    rat4 = (dfdt[lu]+fnoise)/(gdfdt[lu]+fnoise);
    if (VERBOSE) printf("                 (%9.4f)    (%9.4f)    (%9.4f)    (%9.4f)\n",rat1,rat2,rat3,rat4);

    if (rat1<0.99 || rat1>1.01) numbad++;
    if (rat2<0.99 || rat2>1.01) numbad++;
    if (rat3<0.99 || rat3>1.01) numbad++;
    if (rat4<0.99 || rat4>1.01) numbad++;
  }

  if (!onlyfl) {

    for (j = 0; j < nphi; j++) {
      for (lu = 0; lu < ntau; lu++) {
        for (iu = 0; iu < numu; iu++) {
          if (guu[0][iu][lu][j]>umax) umax = guu[0][iu][lu][j];
        }
      }
    } unoise = 1e-6*umax;

    if (unoise>1e-30) {
      if (VERBOSE) printf("\n\n --------  I N T E N S I T I E S  --------*"
             "\n\n             Polar   Azimuthal Angles (Degrees)"
               "\n   Optical   Angle"
               "\n     Depth  Cosine");
      for (j = 0; j < nphi; j++) if (VERBOSE) printf("%10.1f    ",phi[j]);
      if (VERBOSE) printf("\n");

      for (lu = 0; lu < ntau; lu++) {
        for (iu = 0; iu < numu; iu++) {
          if (iu == 0 && VERBOSE) {
            printf("\n%10.3f%8.3f",utau[lu],umu[iu]);
            for (j = 0; j < nphi; j++) printf("%14.4e",uu[0][iu][lu][j]);
            printf("\n");
          }
          if (iu > 0 && VERBOSE) {
            printf("          %8.3f",umu[iu]);
            for(j = 0; j < nphi; j++) printf("%14.4e",uu[0][iu][lu][j]);
            printf("\n");
          }
          for (j = 0; j < nphi; j++) {
            ratv[j] = (uu[0][iu][lu][j]+unoise)/(guu[0][iu][lu][j]+unoise);
            if (ratv[j]<0.99 || ratv[j]>1.01) numbad++;
          }
          if (VERBOSE) printf("                  ");
          for (j = 0; j < nphi; j++) if (VERBOSE) printf("   (%9.4f)",ratv[j]);
          if (VERBOSE) printf("\n");
        }
      }
    }
  }
  if (numbad > 0) printf(" ====  %4d  SERIOUSLY NON-UNIT RATIOS    ====\n", numbad);
}


// ---------------------------------------------------------------------------------
// Calculate phase function Legendre expansion coefficients in various special cases
// ---------------------------------------------------------------------------------
void getphase(int iphas, double gg, double albedo, int nmom, int ilyr, double ****phfunc) {

  const double cldmom[299] = {
    2.544,3.883,4.568,5.235,5.887,6.457,7.177,7.859,8.494,9.286,9.856,10.615,11.229,11.851,12.503,
    13.058,13.626,14.209,14.660,15.231,15.641,16.126,16.539,16.934,17.325,17.673,17.999,18.329,18.588,
    18.885,19.103,19.345,19.537,19.721,19.884,20.024,20.145,20.251,20.330,20.401,20.444,20.477,20.489,
    20.483,20.467,20.427,20.382,20.310,20.236,20.136,20.036,19.909,19.785,19.632,19.486,19.311,19.145,
    18.949,18.764,18.551,18.348,18.119,17.901,17.659,17.428,17.174,16.931,16.668,16.415,16.144,15.883,
    15.606,15.338,15.058,14.784,14.501,14.225,13.941,13.662,13.378,13.098,12.816,12.536,12.257,11.978,
    11.703,11.427,11.156,10.884,10.618,10.350,10.090,9.827,9.574,9.318,9.072,8.822,8.584,8.340,8.110,
    7.874,7.652,7.424,7.211,6.990,6.785,6.573,6.377,6.173,5.986,5.790,5.612,5.424,5.255,5.075,4.915,
    4.744,4.592,4.429,4.285,4.130,3.994,3.847,3.719,3.580,3.459,3.327,3.214,3.090,2.983,2.866,2.766,
    2.656,2.562,2.459,2.372,2.274,2.193,2.102,2.025,1.940,1.869,1.790,1.723,1.649,1.588,1.518,1.461,
    1.397,1.344,1.284,1.235,1.179,1.134,1.082,1.040,0.992,0.954,0.909,0.873,0.832,0.799,0.762,0.731,
    0.696,0.668,0.636,0.610,0.581,0.557,0.530,0.508,0.483,0.463,0.440,0.422,0.401,0.384,0.364,0.349,
    0.331,0.317,0.301,0.288,0.273,0.262,0.248,0.238,0.225,0.215,0.204,0.195,0.185,0.177,0.167,0.160,
    0.151,0.145,0.137,0.131,0.124,0.118,0.112,0.107,0.101,0.097,0.091,0.087,0.082,0.079,0.074,0.071,
    0.067,0.064,0.060,0.057,0.054,0.052,0.049,0.047,0.044,0.042,0.039,0.038,0.035,0.034,0.032,0.030,
    0.029,0.027,0.026,0.024,0.023,0.022,0.021,0.020,0.018,0.018,0.017,0.016,0.015,0.014,0.013,0.013,
    0.012,0.011,0.011,0.010,0.009,0.009,0.008,0.008,0.008,0.007,0.007,0.006,0.006,0.006,0.005,0.005,
    0.005,0.005,0.004,0.004,0.004,0.004,0.003,0.003,0.003,0.003,0.003,0.003,0.002,0.002,0.002,0.002,
    0.002,0.002,0.002,0.002,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
    0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  const double hazelm[82] = {
    2.41260,3.23047,3.37296,3.23150,2.89350,2.49594,2.11361,1.74812,1.44692,1.17714,0.96643,0.78237,
    0.64114,0.51966,0.42563,0.34688,0.28351,0.23317,0.18963,0.15788,0.12739,0.10762,0.08597,0.07381,
    0.05828,0.05089,0.03971,0.03524,0.02720,0.02451,0.01874,0.01711,0.01298,0.01198,0.00904,0.00841,
    0.00634,0.00592,0.00446,0.00418,0.00316,0.00296,0.00225,0.00210,0.00160,0.00150,0.00115,0.00107,
    0.00082,0.00077,0.00059,0.00055,0.00043,0.00040,0.00031,0.00029,0.00023,0.00021,0.00017,0.00015,
    0.00012,0.00011,0.00009,0.00008,0.00006,0.00006,0.00005,0.00004,0.00004,0.00003,0.00003,0.00002,
    0.00002,0.00002,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001};
  register int k;

  phfunc[0][ilyr][0][0] = albedo;
  for (k = 1; k < nmom; k++) phfunc[k][ilyr][0][0] = 0.0;

  switch(iphas) {
    case RAYLEIGH:
      phfunc[2][ilyr][0][0] = 0.1*albedo;
    break;
    case HENYEY_GREENSTEIN:
      for(k = 1; k <= nmom; k++) phfunc[k][ilyr][0][0] = pow(gg,(double)k)*albedo;
    break;
    case HAZE_GARCIA_SIEWERT:
      // Haze-L phase function
      for (k = 1; k <= nmom; k++) phfunc[k][ilyr][0][0] = hazelm[k-1]/(double)(2*k+1)*albedo;
    break;
    case CLOUD_GARCIA_SIEWERT:
      // Cloud C.1 phase function
      for (k = 1; k <= nmom; k++) phfunc[k][ilyr][0][0] = cldmom[k-1]/(double)(2*k+1)*albedo;
    break;
  }
  return;
}
