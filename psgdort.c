// --------------------------------------------------------------------------------------------------------------------
// PSG Discrete Ordinate Radiative Transfer module, PSGDORT, NASA/GSFC, Villanueva, 2021
// Based on cdisort (rewritten in C from Fortran DISORT v2.1) by Dowling 2011
// Pseudospherical geometry implemented by Arve Kylling
//----------------------- References ----------------------------------------------------------------------------------
//    DGIS: Devaux C, Grandjean P, Ishiguro Y, Siewert CE, 1979, On Multi-Region Problems in Radiative Transfer, Astrophys. Space Sci. 62, 225-233
//      GS: Garcia RDM, Siewert CE, 1985, Benchmark Results in Radiative Transfer, Transport Theory and Statistical Physics 14, 437-483
//      DS: Dahlback, A. and K. Stamnes 1991: A new spherical model for computing the radiation field available for photolysis and heating at twilight, Planet. Space Sci. 39, 671-683.
//      KS: Kylling A, Stamnes K, 1992, Efficient yet accurate solution of the linear transport equation in the presence of internal sources: The exponential-linear-in-depth approximation, J. Comp. Phys., 102, 265-276
//       L: Lenoble J, ed, 1985:  Radiative Transfer in Absorbing and Scattering Atmospheres: Standard Computational Procedures, Deepak Publishing, Hampton, Virginia
//      NT: Nakajima T, Tanaka M, 1988,  Algorithms for Radiative Intensity Calculations in Moderately Thick Atmospheres Using a Truncation Approximation, J.Q.S.R.T. 40, 51-69
//      OS: Ozisik M, Shouman S, 1980,  Source Function Expansion Method for Radiative Transfer in a Two-Layer Slab, J.Q.S.R.T. 24, 441-449
//      SS: Stamnes K, Swanson R, 1981,  A New Look at the Discrete Ordinate Method for Radiative Transfer Calculations in Anisotropically Scattering Atmospheres, J. Atmos. Sci. 38, 387-399
//      SD: Stamnes K, Dale H, 1981, A New Look at the Discrete Ordinate Method for Radiative Transfer Calculations in Anisotropically Scattering Atmospheres. II: Intensity Computations, J. Atmos. Sci. 38, 2696-2706
//      S1: Stamnes K, 1982, On the Computation of Angular Distributions of Radiation in Planetary Atmospheres, J.Q.S.R.T. 28, 47-51
//      S2: Stamnes K, 1982, Reflection and Transmission by a Vertically Inhomogeneous Planetary Atmosphere, Planet. Space Sci. 30, 727-732
//      SC: Stamnes K, Conklin P, 1984, A New Multi-Layer Discrete Ordinate Approach to Radiative Transfer in Vertically Inhomogeneous Atmospheres, J.Q.S.R.T. 31, 273-282
//      SW: Sweigart A, 1970, Radiative Transfer in Atmospheres Scattering According to the Rayleigh Phase Function with Absorption, The Astrophysical Journal Supplement Series 22, 1-80
//    STWJ: Stamnes K, Tsay SC, Wiscombe W, Jayaweera K, 1988, A Numerically Stable Algorithm for Discrete-Ordinate-Method Radiative Transfer in Multiple Scattering and Emitting Layered Media, Appl. Opt. 27, 2502-2509
//    STWL: Stamnes K, Tsay SC, Wiscombe W, Laszlo I: A General-Purpose Numerically Stable Computer Code for Discrete-Ordinate-Method Radiative Transfer in Scattering and Emitting Layered Media, DISORT Report v1.1 (2000)
// VH1,VH2: Van de Hulst, H.C., 1980: Multiple Light Scattering, Tables, Formulas and Applications, Volumes 1 and 2, Academic Press, New York.
//       W: Wiscombe, W., 1977:  The Delta-M Method: Rapid Yet Accurate Radiative Flux Calculations, J. Atmos. Sci. 34, 1408-1422
//-----------------------------------------------------------------------------------------------------------------
// Function prototypes
void asymtx(double **aa, double **evec, double *eval, int m, int ia, int ievec, double *wk);
void sgeco(double **a, int lda, int n, int *ipvt, double *rcond, double *z);
void sgesl(double **a, int lda, int n, int *ipvt, double *b, int job);
void sgbco(double **abd, int lda, int n, int ml, int mu, int *ipvt, double *rcond, double *z);
void sgbfa(double **abd, int lda, int n, int ml, int mu, int *ipvt, int *info);
void sgbsl(double **abd, int lda, int n, int ml, int mu, int *ipvt, double *b, int job);
double fplanck(double lam, double dl, double T);

// Multiple-scattering radiative-transfer solver
void psgdort(
  long l1,// ................ Index of first wavelength bin to process (use 0 or one point)
  long l2,// ................ Index of last wavelength bin to process (not inclusive) (use 1 for one point)
  int nfine,// .............. Number of fine sub-divisions for each wavelength bin (use 1 for one point)
  int nbin,// ............... Number of correlated-k bins (use 1 for one point)
  int nlyr,// ............... Number of layers
  int nmax,// ............... Number of stream pairs
  int nmom,// ............... Maximum number of Legendre polynomials in the phase function
  int numu,// ............... Number of desired emission angles (use 1 for one observing angle)
  int nphi,// ............... Number of desired azimuths angles (use 1 for one observing angle)
  int ntau,// ............... Number of vertical/tau locations that fluxes are desired (use 1 for one observation location)
  int ulyr,// ............... Location of the observer in layer index (use -1 if tau values for the location of the observer are provided)
  double umu0,// ............ Cosine of the Solar incidence angle
  double *dphi,// ........... Solar azimuth angles (nphi values)
  double radius,// .......... Planetary radius (in the same units as the altitudes)
  double *alts,// ........... Altitude of the layer borders (nlyr+1 points)
  double ***dtau,// ......... Opacity at each layer, size is nlyr x nbin x [nfine x l2]
  double ****phfunc,// ...... Scattering phase function (L0=scatt_albedo, L1, ... Lnmom) size is [nmom+1] x nlyr x nbin x [nfine x l2]
  double *fstar,// .......... Flux of the star for each wavelength bin in the same units as srclyr or in [W/sr/um/m2] for dlam=0 or [W/sr/m2] otherwise (see fplanck)
  double *asurf,// .......... Albedo of the surface for each wavelength bin
  double *esurf,// .......... Emissivity of the surface for each wavelength bin
  double temis,// ........... Emissivity of the top of the atmosphere (0.0 is typical)
  double **srclyr,// ........ Source function for each layer border (allowing non-LTE sources), size is [nlyr+1] x [nfine x l2]
  double fisot,// ........... Isotropic flux, in same units as srclyr or fstar or fplanck
  double *lamc,// ........... Central wavelength [um] of each fine wavelength bin, only needed when thermal emission is required and no srclyr is provided
  double dlam,// ............ Width [um] of the of the fine wavelength bins, only needed when thermal emission is required and no srclyr is provided
  double *temper,// ......... Temperature [K] at the each layer, only needed when thermal emission is required and no srclyr is provided
  double btemp,// ........... Temperature [K] of the lower boundary (e.g., surface), only needed when thermal emission is required and no srclyr is provided
  double ttemp,// ........... Temperature [K] of the top boundary (e.g., TOA), only needed when thermal emission is required and no srclyr is provided
  double **pumu, double **putau, double **pdfdt, double **puavg, double **pfldir, double **pfldn, double **pflup, double **prfldir, double **prfldn, double *****puu) {

  // Allocate and initialize scalars
  int twostr=0, numu2=numu; if (nmax<=0) { twostr=1; nmax=1; numu=0; }
  int nstr=2*nmax, naz=nstr, mmom=(nstr>nmom ? nstr : nmom);
  long lf,li,lkf, nlam=(l2-l1)*nfine*nbin; int i, j, k, l, bi, iu, iq, jq, kq, it, ipnt, mazim, deltam, needm, ncut, nl, lc, lu, lyu, kconv, usrtau=1, wsrc, usrang=1;
  double delm0, sum, alpha, beta, gpmigm, gpplgm, rcond, psum0, psum1, fact, expa, diff, dirint, refflx;
  double tt, abstau, fm, fdntot, plsorc, azerr, cosphi, azterm, cratio, fbeam, tplanck=0.0, bplanck=0.0, big = sqrt(DBL_MAX)/1.e+10, q0a, q2a, deltat, q0, q2, albmax=1.0 - 1000.*DBL_EPSILON;
  double small = 1.e+30*DBL_MIN, little = 1.e+20*DBL_MIN, large = log(DBL_MAX)-20.0, val;
  if (numu<=0) { numu=nstr; usrang=0; }
  if (ulyr>=0) { ntau=1; }
  if (ntau<=0) { ntau=nlyr+1; usrtau=0; }
  if ((lamc!=NULL && temper!=NULL) || (srclyr!=NULL)) wsrc=1; else wsrc=0;

  // Allocate output arrays (if not previously allocated)
  if (*pdfdt==NULL) *pdfdt = array1D(ntau,1);// .......... Flux divergence d(net flux)/d(optical depth),where 'net flux' includes the direct beam (an exact result;  not from differencing fluxes)
  if (*puavg==NULL) *puavg = array1D(ntau,1);// .......... Mean intensity (including the direct beam)
  if (*pfldir==NULL) *pfldir = array1D(ntau,1);// ........ Direct-beam flux (delta-M scaled)
  if (*pfldn==NULL) *pfldn = array1D(ntau,1);// .......... Diffuse down-flux (delta-M scaled)
  if (*pflup==NULL) *pflup = array1D(ntau,1);// .......... Diffuse up-flux
  if (*prfldir==NULL) *prfldir = array1D(ntau,1);// ...... Direct-beam flux (without delta-m scaling)
  if (*prfldn==NULL) *prfldn = array1D(ntau,1);// ........ Diffuse down-flux (total minus direct-beam) (without delta-m scaling)
  if (*pumu==NULL) *pumu = array1D(numu,1);// ..,......... Zenith angles at desired intensities
  if (*puu==NULL) *puu = array4D(nlam,numu,ntau,nphi,1);// Corrected intensity field
  if (*putau==NULL) *putau = array1D(ntau,1);// .......... Optical depths of user output levels
  double *umu = *pumu, *utau=*putau, *dfdt=*pdfdt, *uavg=*puavg, *fldir=*pfldir, *fldn=*pfldn, *flup=*pflup, *rfldir=*prfldir, *rfldn=*prfldn, ****uu=*puu;

  // Allocate variables common to full discrete-ordinate solver and two-streams approximation
  double *b = array1D(nstr*nlyr,1);// .................... Right-hand side vector of eq. SC(5) going into SOLVE0,1, returns as solution vector vector  L, the constants of integration
  double **cband = array2D(9*nmax-2,nstr*nlyr,1);// ...... Matrix of left-hand side of the linear system eq. SC(5), scaled by eq. SC(12), in banded form required by LINPACK solution routines
  double *ch = array1D(nlyr,1);// ........................ The Chapman-factor to correct for pseudo-spherical geometry in the direct beam
  double *chtau = array1D(nlyr+1,1);// ................... The optical depth in spherical geometry
  double ***chscl = array3D(nlyr,nlyr,2,1);// ............ Chapman pseudo-spherical scalers
  double *cmu = array1D(nstr,1);// ....................... Computational polar angles (Gaussian)
  double *cwt = array1D(nstr,1);// ....................... Quadrature weights corresponding to CMU
  double *dtauc = array1D(nlyr,1);// ..................... Computational layer optical depth
  double *dtaucpr = array1D(nlyr,1);// ................... Computational-layer optical depths (delta-M-scaled if DELTAM = 1, otherwise equal to DTAUC)
  double *expbea = array1D(nlyr+1,1);// .................. Transmission of direct beam in delta-M optical depth coordinates
  double *flyr = array1D(nlyr,1);// ...................... Separated fraction in delta-M method
  double *ggprim = array1D(nlyr, 1);// ................... Phase asymmetry factor scaled by delta-M
  double **gl = array2D(nstr, nlyr, 1);// ................ Phase function Legendre polynomial expansion coefficients, calculated from PMOM by including single-scattering albedo, factor 2K+1, and (if DELTAM=1) the delta-M scaling
  double ***gu = array3D(numu,nstr,nlyr,1);// ............ Eigenvectors interpolated to user polar angles (g  in eqs. SC(3) and S1(8-9), i.e. g without the l factor)
  int *ipvt = malloc(nstr*nlyr*sizeof(int));// ........... Integer vector of pivot indices for LINPACK routines
  int *layru = malloc(ntau*sizeof(int));// ............... Computational layer in which user output level UTAU(LU) is located
  double **ll = array2D(nstr,nlyr,1);// .................. Constants of integration L in eq. SC(1), obtained by solving scaled version of eq. SC(5)
  double *lsrc = array1D(nlyr+1,1);// .................... Layer source for internal emission at the layer borders
  double *oprim = array1D(nlyr,1);// ..................... Single scattering albedo after delta-M scaling
  double **pmom = array2D(mmom+1, nlyr,1);// ............. Legendre polynomials of the phase function
  double *ssalb = array1D(nlyr,1);// ..................... Single scattering albedo of the layer
  double *tauc = array1D(nlyr+1,1);// .................... Cumulative optical depth (un-delta-M-scaled)
  double *taucpr = array1D(nlyr+1,1);// .................. Cumulative optical depth (delta-M-scaled if DELTAM = 1, otherwise equal to TAUC)
  double **u0c = array2D(nstr,ntau,1);// ................. Azimuthally-averaged intensity at quadrature angle
  double *utaupr = array1D(ntau,1);// .................... Optical depths of user output levels in delta-M coordinates; equal to UTAU(LU) if no delta-M
  double ***ylm0 = array3D(nstr,nstr,1,1);// ............. Normalized associated Legendre polynomial of subscript L at the beam angle
  double ***ylmc = array3D(nstr,nstr,nstr,1);// .......... Normalized associated Legendre polynomial of subscript L at the computational angles
  double ***ylmu = array3D(nstr,nstr,numu,1);// .......... Normalized associated Legendre polynomial of subscript L at the user angles

  // Allocate memory for the different methods
  double *diag=NULL, *kk2=NULL, *rr=NULL, *subd=NULL, *superd=NULL, *xp0=NULL, *xp1=NULL, *yb0d=NULL, *yb0u=NULL, *yb1d=NULL, *yb1u=NULL, *yp0d=NULL, *yp0u=NULL, *yp1d=NULL, *yp1u=NULL, *zba=NULL, *zpa=NULL;
  double **amb=NULL, **apb=NULL, **array=NULL, **bdr=NULL, *bem=NULL, **cc=NULL, *emu=NULL, *eval=NULL, **evecc=NULL, ***gc=NULL, **kk=NULL, *phasa=NULL, *phasm=NULL, *phast=NULL, *psi0=NULL, *psi1=NULL, **rmu=NULL, **uum=NULL, *wk=NULL;
  double *xba=NULL, **xb0=NULL, **xb1=NULL, *xr0=NULL, *xr1=NULL, *z=NULL, *z0=NULL, **z0u=NULL, *z1=NULL, **z1u=NULL, **zb0u=NULL, **zb1u=NULL, **zbau=NULL, *zbs0=NULL, *zbs1=NULL, *zbeama=NULL, **zbeam0=NULL, **zbeam1=NULL, *zj=NULL, *zju=NULL, **zplk0=NULL, **zplk1=NULL;
  if (twostr) {
    // Two-streams variables
    diag = array1D(2*nlyr,1);// ......................... Two-streams matrix element
    kk2 = array1D(nlyr,1);// ............................ Eigenvalues coeff. of two-streams solution
    rr = array1D(nlyr,1);// ............................. Eigenvectors of two-streams solution
    subd = array1D(2*nlyr,1);// ......................... Two-streams matrix element
    superd = array1D(2*nlyr,1);// ....................... Two-streams matrix element
    xp0 = array1D(nlyr,1);// ............................ Constant for two-streams solution
    xp1 = array1D(nlyr,1);// ............................ Constant for two-streams solution
    yb0d = array1D(nlyr,1);// ........................... Incident beam particular solution of two-streams solution
    yb0u = array1D(nlyr,1);// ........................... Incident beam particular solution of two-streams solution
    yb1d = array1D(nlyr,1);// ........................... Incident beam particular solution of two-streams solution
    yb1u = array1D(nlyr,1);// ........................... Incident beam particular solution of two-streams solution
    yp0d = array1D(nlyr,1);// ........................... Thermal particular solution of two-streams solution
    yp0u = array1D(nlyr,1);// ........................... Thermal particular solution of two-streams solution
    yp1d = array1D(nlyr,1);// ........................... Thermal particular solution of two-streams solution
    yp1u = array1D(nlyr,1);// ........................... Thermal particular solution of two-streams solution
    zba = array1D(nlyr,1);// ............................ Alpha coefficient of incident for two-streams solution
    zpa = array1D(nlyr,1);// ............................ Alpha coefficient of thermal for two-streams solution

  } else {
    // Full multiple scattering variables
    amb = array2D(nmax,nmax,1);// ....................... First matrix factor in reduced eigenvalue problem of eqs. SS(12), STWJ(8E), STWL(23f) (used only in solve_eigen);
    apb = array2D(nmax,nmax,1);// ....................... Second matrix factor in reduced eigenvalue problem of eqs. SS(12), STWJ(8E), STWL(23f) (used only in solve_eigen);
    array = array2D(nstr,nstr,1);// ..................... Scratch matrix for solve_eigen(), upbeam() and upisot()
    bdr = array2D(nmax, nmax,1);// ...................... Bottom-boundary bidirectional reflectivity for a given azimuthal component.  First index always, refers to a computational angle, second index: if zero, refers to incident beam angle UMU0, if non-zero, refers to a computational angle.
    bem = array1D(nmax, 0);// ........................... Bottom-boundary directional emissivity at computational angles.
    cc = array2D(nstr,nstr,1);// ........................ C-sub-IJ in eq. SS(5)
    emu = array1D(numu,1);// ............................ Emissivity of the surface at numu angles
    eval = array1D(nmax,1);// ........................... Temporary storage for eigenvalues of eq. SS(12)
    evecc = array2D(nstr,nstr,1);// ..................... Complete eigenvectors of SS(7) on return from solve_eigen; stored permanently in  GC
    gc = array3D(nstr,nstr,nlyr,1);// ................... Eigenvectors at polar quadrature angles, g in eq. SC(1)
    kk = array2D(nstr,nlyr,1);// ........................ Eigenvalues of coeff. matrix in eq. SS(7)
    phasa = array1D(nlyr,1);// .......................... Nakajima/Tanaka phase function
    phasm = array1D(nlyr,1);// .......................... Nakajima/Tanaka phase function
    phast = array1D(nlyr,1);// .......................... Nakajima/Tanaka phase function
    psi0 = array1D(nstr,1);// ........................... Sum just after square bracket in  eq. SD(9)
    psi1 = array1D(nstr,1);// ........................... Sum in  eq. STWL(31d)
    rmu = array2D(numu,nmax,1);// ....................... Bottom-boundary bidirectional reflectivity for a given azimuthal component.  First index always
    uum = array2D(numu,ntau,1);// ....................... Expansion coefficients when the intensity (u-super-M) is expanded in Fourier cosine series
    wk = array1D(nstr,1);// ............................. Scratch array
    xba = array1D(nlyr,1);// ............................ alfa in eq. KS(7)
    xb0 = array2D(nstr,nlyr,1);// ....................... x-sub-zero in KS(7)
    xb1 = array2D(nstr,nlyr,1);// ....................... x-sub-one in KS(7)
    xr0 = array1D(nlyr,1);// ............................ X-sub-zero in expansion of thermal source function preceding eq. SS(14)(has no mu-dependence); b-sub-zero in eq. STWL(24d)
    xr1 = array1D(nlyr,1);// ............................ X-sub-one in expansion of thermal source function; see eqs. SS(14-16); b-sub-one in STWL(24d)
    z = array1D(nstr*nlyr,1);// ......................... Scratch array used in solve0(), albtrans() to solve a linear system for the constants of integration
    z0 = array1D(nstr,1);// ............................. Solution vectors Z-sub-zero of eq. SS(16)
    z0u = array2D(numu,nlyr,1);// ....................... Z-sub-zero in eq. SS(16) interpolated to user angles from an equation derived from SS(16)
    z1 = array1D(nstr,1);// ............................. Solution vectors Z-sub-one  of eq. SS(16)
    z1u = array2D(numu,nlyr,1);// ....................... Z-sub-one  in eq. SS(16) interpolated to user angles from an equation derived from SS(16)
    zb0u = array2D(numu,nlyr,1);// ...................... x-sub-zero in KS(7) at user angles -umu-
    zb1u = array2D(numu,nlyr,1);// ...................... x-sub-one in KS(7) at user angles -umu-
    zbau = array2D(numu,nlyr,1);// ...................... x-sub-zero in KS(7) at user angles -umu-
    zbs0 = array1D(nstr,1);// ........................... solution vectors z-sub-zero of Eq. KS(10-11)
    zbs1 = array1D(nstr,1);// ........................... solution vectors z-sub-one  of Eq. KS(10-11)
    zbeama = array1D(nlyr,1);// ......................... Permanent storage for -zbs0,zbs1-, but rD-ordered
    zbeam0 = array2D(nstr,nlyr,1);// .................... Permanent storage for -zbs0,zbs1-, but rD-ordered
    zbeam1 = array2D(nstr,nlyr,1);// .................... Permanent storage for -zbs0,zbs1-, but rD-ordered
    zj = array1D(nstr,1);// ............................. Right-hand side vector  X-sub-zero in eq. SS(19), also the solution vector Z-sub-zero after solving that system
    zju = array1D(numu,1);// ............................ Right-hand side vector  X-sub-zero in eq. SS(19), also the solution vector Z-sub-zero after solving that system
    zplk0 = array2D(nstr,nlyr,1);// ..................... Permanent storage for the thermal source vectors plk[].zero obtained by solving eq. SS(16)
    zplk1 = array2D(nstr,nlyr,1);// ..................... Permanent storage for the thermal source vectors plk[].one  obtained by solving eq. SS(16)
  }


  // ---------------------------------------------------------------
  // Initialize parameters for the scattering calculation
  // ---------------------------------------------------------------
  // Compute gaussian quadrature angles and weights ----------------
  if (twostr && fstar[0]*umu0>little) { cmu[0]=0.5; cwt[0] = 1.0; }
  else if (nmax == 1) { cmu[0]=0.5; cwt[0] = 1.0; }
  else {
    double tol = 10.*1e-16, en = nmax, nnp1 = nmax*(nmax+1.0), cona = (nmax-1.0)/(8.0*nmax*nmax*nmax);
    double xi, p2pri, ppr, tmp, p, pm2, pm1, x, t, prod;
    int iter, k, nn, lim = nmax/2;

    for (k = 0; k < lim; k++) {
      // Initial guess for k-th root of Legendre polynomial, from Davis/Rabinowitz (2.7.3.3a)
      t = (4.0*(k+1.0) - 1.0)*M_PI/(4.0*nmax + 2.0);
      x = cos(t+cona/tan(t));

      // Upward recurrence for Legendre polynomials
      for (iter = 0; iter < 1000; iter++) {
        pm2 = 1.0; pm1 = x;
        for (nn = 2; nn <= nmax; nn++) {
          p   = ((2.0*nn - 1.0)*x*pm1 - (nn-1.0)*pm2)/nn;
          pm2 = pm1;
          pm1 = p;
        }
        // Newton Method
        tmp   = 1.0/(1.0 - x*x);
        ppr   = en*(pm2 - x*p)*tmp;
        p2pri = (2.0*x*ppr - nnp1*p)*tmp;
        xi    = x - p/ppr*(1.0 + p/ppr*p2pri/(2.0*ppr));

        // Check for convergence
        if (fabs(xi-x) <= tol) break; else x = xi;
      }

     // Iteration finished--calculate weights, abscissae for (-1,1)
      cmu[k] = -x;
      cwt[k] = 2.0/(tmp*(en*pm2)*(en*pm2));
      cmu[nmax-k-1] = -cmu[k];
      cwt[nmax-k-1] =  cwt[k];
    }

    // Set middle abscissa and weight for rules of odd order
    if (nmax % 2 != 0) {
      cmu[lim] = 0.0;
      prod     = 1.0;
      for (k = 3; k <= nmax; k+=2) prod *= (double)k/(k-1);
      cwt[lim] = 2.0/(prod*prod);
    }
    // Convert from (-1,1) to (0,1)
    for (k = 0; k < nmax; k++) {
      cmu[k] = 0.5*cmu[k]+0.5;
      cwt[k] = 0.5*cwt[k];
    }
  }

  // Downward (neg) angles and weights
  for (k = 0; k < nmax; k++) {
   cmu[k + nmax] = -cmu[k];
   cwt[k + nmax] =  cwt[k];
  }

  // Check if we need to define the observational angles to be the quadrature angles
  if (!usrang) {
    for (iu = 0;    iu < nmax; iu++) umu[iu] = -cmu[nmax-1-iu];
    for (iu = nmax; iu < nstr; iu++) umu[iu] =  cmu[iu-nmax];
  }

  // Dither observation angles if they match exactly the incidence angle or are close to tangential
  for (iu=0; iu<numu; iu++) {
    if (umu[iu]>1.0) umu[iu]=1.0;
    else if (umu[iu]<-1.0) umu[iu]=-1.0;
    else if (umu[iu]<1e-5 && umu[iu]>=0) umu[iu]=1e-5;
    else if (umu[iu]>-1e-5 && umu[iu]<0) umu[iu]=-1e-5;
    if (fabs(umu[iu]+umu0)<1e-6) { if ((umu[iu]+1e-6)>1.0) umu[iu] -= 1e-6; else umu[iu] += 1e-6; }
  }

  // Compute Legendre polynomials -----------------------------------------
  for (mazim = 0; mazim < nstr && !twostr; mazim++) {
    for (k = 0; k < nmax+2; k++) {
      double ***ylm, mus[nmax+numu], tmp1, tmp2; int nmu;
      if (k==0) { nmu = 1; mus[0] = -umu0; ylm = ylm0; }  // Solar beam Legendre polynomials
      else if (k==1) { nmu = numu; for (i=0; i<nmu; i++) mus[i] = umu[i]; ylm = ylmu; } // User angle Legendre polynomials
      else { nmu = nmax; for (i=0; i<nmu; i++) mus[i] = cmu[i]; ylm = ylmc; } // Gaussian angles Legendre polynomials

      if (mazim==0) {
        // Upward recurrence for ordinary Legendre polynomials
        for (i = 0; i < nmu; i++) {
          ylm[mazim][0][i] = 1.0;
          ylm[mazim][1][i] = mus[i];
        }
        for (l = 2; l < nstr; l++) {
          for (i = 0; i < nmu; i++) {
            ylm[mazim][l][i] = ((double)(2*l-1)*mus[i]*ylm[mazim][l-1][i] - (double)(l-1)*ylm[mazim][l-2][i])/l;
          }
        }
      }
      else {
        for (i = 0; i < nmu; i++) {
          for (l = 0; l < mazim; l++) ylm[mazim][l][i] = ylm[mazim-1][l][i];
          // Y-sub-m-super-m; derived from D/A eqs. (11,12), STWL(58c)
          ylm[mazim][mazim][i] = -sqrt((1.0-1.0/(2.0*mazim))*(1.0-pow(mus[i],2.0)))*ylm[mazim][mazim-1][i];
          // Y-sub-(m+1)-super-m; derived from D/A eqs.(13,14) using eqs.(11,12), STWL(58f)
          if (mazim<nstr-1) ylm[mazim][mazim+1][i] = sqrt(2.0*mazim + 1.0)*mus[i]*ylm[mazim][mazim][i];
        }
        // Upward recurrence; D/A eq.(10), STWL(58a)
        for (l = mazim+2; l < nstr; l++) {
          tmp1 = sqrt((l-mazim  )*(l+mazim  ));
          tmp2 = sqrt((l-mazim-1)*(l+mazim-1));
          for (i = 0; i < nmu; i++) {
            ylm[mazim][l][i] = ((double)(2*l-1)*mus[i]*ylm[mazim][l-1][i] - tmp2*ylm[mazim][l-2][i])/tmp1;
          }
        }
      }
    }

    // Get normalized associated Legendre polynomials with negative arguments from those with positive arguments; Dave/Armstrong eq. (15), STWL(59)
    double sgn = -1.0;
    if (mazim>0) for (l=0; l<mazim; l++) for (i=nmax; i<nstr; i++) ylmc[mazim][l][i] = ylmc[mazim-1][l][i];
    for (l=mazim; l<nstr; l++) {
      sgn *= -1.0;
      for (i=nmax; i<nstr; i++) ylmc[mazim][l][i] = sgn*ylmc[mazim][l][i-nmax];
    }
  }

  // Compute Chapman scalers for incidence angle
  int hinv=0; if (alts[1]>alts[0]) hinv=1; // Direction of layers, default(hinv:0) is TOA is layer=0, while hinv:1 bottom is layer=0
  double zd, zdm1, xp, xpsinz, rj, rjp1, dhj, dsj, fact2;
  double sunangle = acos(umu0)*180.0/M_PI;
  double sinz = sin(sunangle*M_PI/180.0);
  for (lc=0; lc<nlyr; lc++) {
    for (k=0; k<2; k++) {
      if (hinv) { zdm1 = alts[nlyr-lc]; zd = alts[nlyr-lc-1]; } else { zdm1 = alts[lc]; zd = alts[lc+1]; }
      xp = radius + zd + (zdm1 - zd)*(k==0 ? 0.0 : 0.5);
      xpsinz = xp*sinz;
      if (sunangle>90.0 && xpsinz<radius) {
        chscl[lc][lc][k] = 1e20;
        continue;
      }

      // Find index of layer in which the screening height lies
      i = lc;
      if (sunangle>90.0) {
        for (j=lc; j<nlyr; j++) {
          if (hinv) { zdm1 = alts[nlyr-j]; zd = alts[nlyr-j-1]; } else { zdm1 = alts[j+1]; zd = alts[j]; }
          if (xpsinz<(zdm1+radius) && xpsinz>=(zd+radius)) i=j;
        }
      }

      // Iterate through the top layers
      for (j=0; j<=i; j++) {
        // Include factor of 2 for zenang > 90., second sum in eq. B2 (DS)
        fact = 1.0; fact2 = 1.0;
        if (j>lc) fact=2.0;
        else if (j==lc && lc==i && sunangle>90.0) fact2=-1.0;
        if (hinv) { zdm1 = alts[nlyr-j]; zd = alts[nlyr-j-1]; } else { zdm1 = alts[j]; zd = alts[j+1]; }
        rj   = radius + zdm1;
        rjp1 = radius + zd;
        if (j==lc && i==lc) rjp1 = xp;
        dhj = zdm1 - zd;

        if (i>lc && j==i) {
          dsj = sqrt(rj*rj - xpsinz*xpsinz);
        } else {
          dsj = sqrt(rj*rj-xpsinz*xpsinz) - fact2*sqrt(rjp1*rjp1-xpsinz*xpsinz);
        }
        chscl[lc][j][k] = fact*dsj/dhj;
      }

      // Add third term in eq. B2 (DS)
      if (i>lc) {
        if (hinv) { zdm1 = alts[nlyr-lc]; zd = alts[nlyr-lc-1]; } else { zdm1 = alts[lc]; zd = alts[lc+1]; }
        rjp1 = radius + zd;
        dhj  = zdm1 - zd;
        dsj  = sqrt(xp*xp-xpsinz*xpsinz)-sqrt(rjp1*rjp1-xpsinz*xpsinz);
        chscl[lc][lc][k] += dsj/dhj;
      }
    }
  }

  // ---------------------------------------------------------------
  // Iterate across wavelength and correlated-k bins
  // ---------------------------------------------------------------
  for (li=l1; li<l2; li++) for (lf=li*nfine; lf<(li+1)*nfine; lf++) for (bi=0; bi<nbin;bi++) {
    // Assign layer values
    for (k=0, tt=0.0, deltam=1, needm=0, tauc[0]=0.0; k<nlyr; k++) {
      if (hinv) lc=nlyr-k-1; else lc=k;
      dtauc[k]  = dtau[lc][bi][lf];
      tauc[k+1] = tauc[k] + dtauc[k];
      tt += dtauc[k];
      ssalb[k] = phfunc[0][lc][bi][lf]; if (ssalb[k]>albmax) ssalb[k] = albmax;
      for (j=0; j<=nmom; j++) pmom[j][k] = phfunc[j][lc][bi][lf]/(phfunc[0][lc][bi][lf] + DBL_EPSILON);
      if (pmom[nstr][k]!=0.0) needm=1;
    }
    if (!usrtau) {
      for (k=0; k<ntau; k++) utau[k] = tauc[k];
    } else if (ulyr>=0) {
      if (hinv) lc=nlyr-ulyr; else lc=ulyr;
      utau[0] = tauc[lc];
    }
    if (!needm || umu0<0.1) deltam=0;
    lkf = (lf - l1*nfine)*nbin + bi;
    fbeam = fstar[li]*M_PI;

    // Apply delta-M scaling and move description of computational layers to local variables
    for (k = 0, ncut=0, abstau=0.0; k < nlyr; k++) {
      if (abstau < 10.0) ncut = k+1;
      abstau += (1.0-ssalb[k])*dtauc[k];
      if (!deltam) {
        // No delta-M transformation
        ggprim[k]   = pmom[1][k];
        oprim[k]    = ssalb[k];
        dtaucpr[k]  = dtauc[k];
        taucpr[k+1] = tauc[k+1];
        for (i = 0; i<nstr; i++) gl[i][k] = (2.0*i + 1.0)*oprim[k]*pmom[i][k];
        fm = 0.0;
      } else {
        // Do delta-M transformation
        fm          = pmom[nstr][k]; if (fm>=1.0-1e-16) fm = 1.0 - 1e-16;
        oprim[k]    = ssalb[k]*(1.0 - fm)/(1.0 - fm*ssalb[k]);
        dtaucpr[k]  = (1.0 - fm*ssalb[k])*dtauc[k];
        ggprim[k]   = (pmom[1][k] - fm)/(1.0 - fm);
        taucpr[k+1] = taucpr[k] + dtaucpr[k];
        for (i = 0; i < nstr; i++) gl[i][k] = (2.0*i + 1.0)*oprim[k]*(pmom[i][k] - fm)/(1.0 - fm);
      }
      flyr[k] = fm;
    }

    // If no thermal emission, cut off medium below absorption optical depth above 10
    if (!(abstau>=10.0 && !wsrc && nlyr>1)) ncut = nlyr;

    // Compute pseudo-spherical Chapman opacities
    double taup=0, chtau_tmp;
    for (lc = 0; lc < ncut; lc++) {
      if (sunangle>90) nl = nlyr; else nl = lc+1;
      taup = taucpr[lc] + dtaucpr[lc]/2.0;
      if (chscl[lc][lc][0]<1e20) for (i = 0, chtau[lc+1]=0.0; i < nl; i++) chtau[lc+1] += dtaucpr[i]*chscl[lc][i][0];
      else chtau[lc+1]=1e20;
      expbea[lc+1] = exp(-chtau[lc+1]);
      if (chscl[lc][lc][1]<1e20) for (i = 0, chtau_tmp=0.0;   i < nl; i++) chtau_tmp   += dtaucpr[i]*chscl[lc][i][1];
      else chtau_tmp=1e20;
      if (chtau_tmp>0.0) ch[lc] = taup/chtau_tmp; else ch[lc]=1.0/umu0;
    } chtau[0] = 0.0;
    if (sunangle>90) expbea[0] = expbea[1]; else expbea[0] = 1.0;

    // Set arrays defining location of user output levels within delta-M-scaled computational mesh
    for (lu = 0; lu < ntau; lu++) {
      for (lc = 0; lc < nlyr; lc++) if (utau[lu] >= tauc[lc] && utau[lu] <= tauc[lc+1]) break;
      if (deltam) {
        utaupr[lu] = taucpr[lc] + (1.0 - ssalb[lc]*flyr[lc])*(utau[lu] - tauc[lc]);
      } else {
        utaupr[lu] = utau[lu];
      }
      layru[lu] = lc;
    }

    // Calculate planck / emissions functions
    if (wsrc && bi==0) {
      if (srclyr!=NULL) {
        // User provided sources
        for (k = 0; k <= nlyr; k++) {
          if (hinv) lc=nlyr-k; else lc=k;
          if (k==0)    tplanck = srclyr[lc][lf]*temis;
          if (k==nlyr) bplanck = srclyr[lc][lf];
          lsrc[k] = srclyr[lc][lf];
        }
      } else {
        // Planck's thermal sources
        tplanck = fplanck(lamc[lf], dlam, ttemp)*temis;
        bplanck = fplanck(lamc[lf], dlam, btemp);
        for (k = 0; k <= nlyr; k++) {
          if (hinv) lc=nlyr-k; else lc=k;
          if (hinv && k==nlyr) lsrc[k] = bplanck;
          else lsrc[k] = fplanck(lamc[lf], dlam, temper[lc]);
        }
      }
    }

    // ---------------------------------------------------------------
    // Compute the two-streams approximation if asked
    // ---------------------------------------------------------------
    if (twostr) {
      double fact1, fact2, fact3, fact4, q_1, q_2, qq, q1a, q1, denomb, denomp, z0p, z0m, arg, sgn, wk0, wk1, wk2, rpp1_m, rpp1_p, rp_m, rp_p;
      double xb0d, xb0u, xb1d, xb1u;
      int nrow, irow, nloop, info;

      // ---------------------------------------------------------------
      // Analytically calculate the homogenous and particular solutions
      // ---------------------------------------------------------------
      // Iterate through computational layers
      for (lc = 0; lc < ncut; lc++) {
        // Calculate eigenvalues -kk- and eigenvector -rr-, eqs. KST(20-21)
        beta   = 0.5*(1.0 - 3.0*ggprim[lc]*cmu[0]*cmu[0]);
        fact1  = 1.0 - oprim[lc];
        fact2  = 1.0 - oprim[lc] + 2.0*oprim[lc]*beta;
        kk2[lc] = (1.0/cmu[0])*sqrt(fact1*fact2);
        rr[lc] = (sqrt(fact2)-sqrt(fact1))/(sqrt(fact2)+sqrt(fact1));

        if (fbeam > 0.0) {
          // Set coefficients in KST(22) for beam source
          q_1 = fbeam/(4.*M_PI)*oprim[lc]*(1.0 - 3.0*ggprim[lc]*cmu[0]*umu0);
          q_2 = fbeam/(4.*M_PI)*oprim[lc]*(1.0 + 3.0*ggprim[lc]*cmu[0]*umu0);

          if (umu0 >= 0.) qq = q_2; else qq = q_1;
          q0a = exp(-chtau[lc]);
          q0  = q0a*qq;
          if (q0 <= small) {
            q1a = 0.;
            q2a = 0.;
          } else {
            q1a = exp(-chtau[lc]);
            q2a = exp(-chtau[lc+1]);
          }
          q1 = q1a*qq;
          q2 = q2a*qq;

          // Calculate alpha coefficient
          deltat = taucpr[lc+1]-taucpr[lc];
          zba[lc] = 1.0/ch[lc];
          if (fabs(zba[lc]*taucpr[lc]) > large || fabs(zba[lc]*taucpr[lc+1]) > large) zba[lc] = 0.;

          // Dither alpha if it is close to an eigenvalue
          denomb = fact1*fact2-pow(zba[lc]*cmu[0],2.0);
          if (denomb < 1.e-03) zba[lc] = 1.02*zba[lc];
          q0 = q0a*q_1;
          q2 = q2a*q_1;

          // Set constants in eq. KST(22)
          if (deltat < 1.e-07) xb1d = 0.0; else xb1d = 1./deltat*(q2*exp(zba[lc]*taucpr[lc+1])-q0*exp(zba[lc]*taucpr[lc]));
          xb0d = q0*exp(zba[lc]*taucpr[lc])-xb1d*taucpr[lc];
          q0 = q0a*q_2;
          q2 = q2a*q_2;

          if (deltat < 1.e-07) xb1u = 0.; else xb1u = 1./deltat*(q2*exp(zba[lc]*taucpr[lc+1])-q0*exp(zba[lc]*taucpr[lc]));
          xb0u = q0*exp(zba[lc]*taucpr[lc])-xb1u*taucpr[lc];

          // Calculate particular solutions for incident beam source in pseudo-spherical geometry, eqs. KST(24-25)
          denomb    = fact1*fact2-pow(zba[lc]*cmu[0],2.0);
          yb1d[lc]  = (oprim[lc]*beta*xb1d+(1.-oprim[lc]*(1.-beta)+zba[lc]*cmu[0])*xb1u)/denomb;
          yb1u[lc]  = (oprim[lc]*beta*xb1u+(1.-oprim[lc]*(1.-beta)-zba[lc]*cmu[0])*xb1d)/denomb;
          z0p       = xb0u-cmu[0]*yb1d[lc];
          z0m       = xb0d+cmu[0]*yb1u[lc];
          yb0d[lc]  = (oprim[lc]*beta*z0m+(1.-oprim[lc]*(1.-beta)+zba[lc]*cmu[0])*z0p)/denomb;
          yb0u[lc]  = (oprim[lc]*beta*z0p+(1.-oprim[lc]*(1.-beta)-zba[lc]*cmu[0])*z0m)/denomb;
        } else {
          zba[lc] =0.0;
          yb1d[lc]=0.0; yb1u[lc]=0.0; yb0d[lc]=0.0; yb0u[lc]=0.0;
        }

        if (wsrc) {
          // Set coefficients in KST(22) for thermal / non-LTE sources
          // Calculate alpha coefficient
          q0     = (1.-oprim[lc])*lsrc[lc];
          q1     = (1.-oprim[lc])*(lsrc[lc+1] + lsrc[lc])/2.0;
          q2     = (1.-oprim[lc])*lsrc[lc+1];
          deltat = taucpr[lc+1]-taucpr[lc];

          if ((q2 < q0*1.e-02 || q2 <= little) && q1 > little && q0 > little) {
            // Case 1: source small at bottom layer; alpha eq. KS(50)
            val = 2.0/deltat*log(q0/q1);
            zpa[lc] = (val>big ? big : val);
            if (zpa[lc]*taucpr[lc] >= log(big)) xp0[lc] = big; else xp0[lc] = q0;
            xp1[lc] = 0.;

          } else if ((q2 <= q1*1.e-02 || q2 <= little) && (q1 <= q0*1.e-02 || q1 <= little) && q0 > little) {
            // Case 2: Source small at center and bottom of layer
            zpa[lc] = big/taucpr[ncut];
            xp0[lc] = q0;
            xp1[lc] = 0.;

          } else if (q2 <= little && q1 <= little && q0 <= little) {
            // Case 3: All sources zero
            zpa[lc] = 0.;
            xp0[lc] = 0.;
            xp1[lc] = 0.;

          } else if ( ( fabs((q2-q0)/q2) < 1.e-04 && fabs((q2-q1)/q2) < 1.e-04 ) || deltat < 1.e-04) {
            // Case 4: Sources same at center, bottom and top of layer or layer optically very thin
            zpa[lc] = 0.;
            xp0[lc] = q0;
            xp1[lc] = 0.;

          } else {
            // Case 5: Normal case
            val = pow(q1/q2,2.0)-q0/q2;
            arg = (val>0 ? val : 0);
            // alpha eq. (44). For source that has its maximum at the top of the layer, use negative solution
            sgn = 1.0; if (lsrc[lc] > lsrc[lc+1]) sgn = -1.;
            fact3 = log(q1/q2+sgn*sqrt(arg));

            // Be careful with log of numbers close to one
            if (fabs(fact3) <= 0.005) {
              // numbers close to one
              q1    = 0.99*q1;
              fact3 = log(q1/q2+sgn*sqrt(arg));
            }

            zpa[lc] = 2./deltat*fact3;
            if (fabs(zpa[lc]*taucpr[lc+1]) > log(DBL_MAX)-log(q0*100.)) zpa[lc] = 0.;

            // Dither alpha if it is close to an eigenvalue
            denomp = fact1*fact2-pow(zpa[lc]*cmu[0],2.0);
            if (denomp < 1.e-03) zpa[lc] *= 1.01;

            // Set constants in eqs. KST(22)
            if (deltat < 1.e-07) xp1[lc] = 0.0; else xp1[lc] = 1./deltat*(q2*exp(zpa[lc]*taucpr[lc+1])-q0*exp(zpa[lc]*taucpr[lc]));
            xp0[lc] = q0*exp(zpa[lc]*taucpr[lc])-xp1[lc]*taucpr[lc];
          }

          // Calculate particular solutions eqs. KST(24-25) for internal thermal so
          denomp   = fact1*fact2-pow(zpa[lc]*cmu[0],2.0);
          yp1d[lc] = (oprim[lc]*beta*xp1[lc]+(1.-oprim[lc]*(1.-beta)+zpa[lc]*cmu[0])*xp1[lc])/denomp;
          yp1u[lc] = (oprim[lc]*beta*xp1[lc]+(1.-oprim[lc]*(1.-beta)-zpa[lc]*cmu[0])*xp1[lc])/denomp;
          z0p      = xp0[lc]-cmu[0]*yp1d[lc];
          z0m      = xp0[lc]+cmu[0]*yp1u[lc];
          yp0d[lc] = (oprim[lc]*beta*z0m+(1.-oprim[lc]*(1.-beta)+zpa[lc]*cmu[0])*z0p)/denomp;
          yp0u[lc] = (oprim[lc]*beta*z0p+(1.-oprim[lc]*(1.-beta)-zpa[lc]*cmu[0])*z0m)/denomp;
        } else {
          zpa[lc] =0.0; xp0[lc]=0.0;  xp1[lc]=0.0;
          yp1d[lc]=0.0; yp1u[lc]=0.0; yp0d[lc]=0.0; yp0u[lc]=0.0;
        }
      }

      // ---------------------------------------------------------------
      // Solve for constants of integration in homogeneous solution (general boundary conditions)
      // ---------------------------------------------------------------
      // First top row, top boundary condition
      irow = 0; lc = 0;
      diag[irow]   = rr[lc]*exp(-kk2[lc]*taucpr[lc+1]);
      superd[irow] = 1.0;

      // next from layer no. 2 to nlyr-1
      nloop = ncut-1;
      for (lc = 0; lc < nloop; lc++) {
        irow++;
        wk0          = exp(-kk2[lc  ]*(taucpr[lc+1]-taucpr[lc]));
        wk1          = exp(-kk2[lc+1]*(taucpr[lc+2]-taucpr[lc+1]));
        subd[irow]   = 1.-rr[lc]*rr[lc+1];
        diag[irow]   = (rr[lc]-rr[lc+1])*wk0;
        superd[irow] = -(1.0 - pow(rr[lc+1],2.0))*wk1;
        irow++;
        subd[irow]   = (1.0 - pow(rr[lc],2.0))*wk0;
        diag[irow]   = (rr[lc]-rr[lc+1])*wk1;
        superd[irow] = -(1.-rr[lc+1]*rr[lc]);
      }

      // bottom layer
      irow++;
      lc = ncut-1;
      wk2 = exp(-kk2[lc]*(taucpr[lc+1]-taucpr[lc]));
      if (ncut<nlyr) {
        subd[irow] = 1.;
        diag[irow] = rr[lc]*wk2;
      } else {
        subd[irow] = 1.0 - 2.0*asurf[li]*cmu[0]*rr[lc];
        diag[irow] = (rr[lc] - 2.0*asurf[li]*cmu[0])*wk2;
      }

      // Construct -b-, for parallel beam + bottom reflection + thermal emission at top and/or bottom
      // Top boundary, right-hand-side of eq. KST(28)
      lc = 0; irow = 0;
      b[irow] = -yb0d[lc]-yp0d[lc]+fisot+tplanck;

      // Continuity condition for layer interfaces, right-hand-side of eq. KST(29)
      for (lc = 0; lc < nloop; lc++) {
        fact1     = exp(-zba[lc+1]*taucpr[lc+1]);
        fact2     = exp(-zpa[lc+1]*taucpr[lc+1]);
        fact3     = exp(-zba[lc  ]*taucpr[lc+1]);
        fact4     = exp(-zpa[lc  ]*taucpr[lc+1]);
        rpp1_m    = fact1*(yb0d[lc+1]+yb1d[lc+1]*taucpr[lc+1])+fact2*(yp0d[lc+1]+yp1d[lc+1]*taucpr[lc+1]);
        rp_m      = fact3*(yb0d[lc  ]+yb1d[lc  ]*taucpr[lc+1])+fact4*(yp0d[lc  ]+yp1d[lc  ]*taucpr[lc+1]);
        rpp1_p    = fact1*(yb0u[lc+1]+yb1u[lc+1]*taucpr[lc+1])+fact2*(yp0u[lc+1]+yp1u[lc+1]*taucpr[lc+1]);
        rp_p      = fact3*(yb0u[lc  ]+yb1u[lc  ]*taucpr[lc+1])+fact4*(yp0u[lc  ]+yp1u[lc  ]*taucpr[lc+1]);
        b[++irow] = rpp1_p-rp_p-rr[lc+1]*(rpp1_m-rp_m);
        b[++irow] = rpp1_m-rp_m-rr[lc  ]*(rpp1_p-rp_p);
      }

      // Bottom boundary
      lc = ncut-1;
      if (ncut<nlyr) { //lyrcut
        // Right-hand-side of eq. KST(30)
        b[++irow] = -exp(-zba[lc]*taucpr[lc])*(yb0u[lc]+yb1u[lc]*taucpr[lc+1])
                    -exp(-zpa[lc]*taucpr[lc])*(yp0u[lc]+yp1u[lc]*taucpr[lc+1]);
      } else {
        sum = cmu[0]*asurf[li]*(exp(-zba[lc]*taucpr[lc])*(yb0d[lc]+yb1d[lc]*taucpr[lc+1])
                               +exp(-zpa[lc]*taucpr[lc])*(yp0d[lc]+yp1d[lc]*taucpr[lc+1]));
        if (umu0 <= 0.) refflx = 0.; else refflx = 1.0;
        b[++irow] = 2.0*sum + asurf[li]*umu0*fbeam/M_PI*refflx*expbea[lc+1] + esurf[li]*bplanck
                   -exp(-zba[lc]*taucpr[lc])*(yb0u[lc]+yb1u[lc]*taucpr[lc+1])
                   -exp(-zpa[lc]*taucpr[lc])*(yp0u[lc]+yp1u[lc]*taucpr[lc+1]);
      }

      // Solve for constants of integration by inverting matrix KST(38-41)
      nrow = irow+1;
      for (irow = 0; irow < nrow; irow++) {
        cband[0][irow] = 0.;
        cband[2][irow] = diag[irow];
      }
      for (irow = 0; irow < nrow-1; irow++) cband[1][irow+1] = superd[irow];
      for (irow = 1; irow < nrow;   irow++) cband[3][irow-1] = subd[irow];

      // Solve the matrix
      sgbfa(cband,(9*nmax-2),nrow,1,1,ipvt,&info);
      sgbsl(cband,(9*nmax-2),nrow,1,1,ipvt,b,0);

      // Unpack the solution
      irow = 0;
      for (lc = 0; lc < ncut; lc++) {
        ll[0][lc] = b[irow++];  // downward direction
        ll[1][lc] = b[irow++];  // upward direction
      }

      // ---------------------------------------------------------------
      // Compute upward and downward fluxes, mean intensities and flux divergences
      // ---------------------------------------------------------------
      for (lu = 0; lu < ntau; lu++) { u0c[0][lu]=0.0; u0c[1][lu]=0.0; }
      if (wsrc) {
        for (lu = 0; lu < ntau; lu++) {
          lyu         = layru[lu];
          fact1       = exp(-zpa[lyu]*utaupr[lu]);
          u0c[0][lu] += fact1*(yp0d[lyu]+yp1d[lyu]*utaupr[lu]);
          u0c[1][lu] += fact1*(yp0u[lyu]+yp1u[lyu]*utaupr[lu]);
        }
      }

      // Loop over user levels
      for (lu = 0; lu < ntau; lu++) {
        fldn[lu]=0.0; flup[lu]=0.0; rfldn[lu]=0.0;
        uavg[lu]=0.0; dfdt[lu]=0.0; dirint=0.0;
        lyu = layru[lu];
        if (lyu<ncut) {
          if (fbeam > 0.0) {
            fact1      = exp(-zba[lyu]*utaupr[lu]);
            u0c[0][lu] += fact1*(yb0d[lyu]+yb1d[lyu]*utaupr[lu]);
            u0c[1][lu] += fact1*(yb0u[lyu]+yb1u[lyu]*utaupr[lu]);
            fact1      = fbeam*exp(-utaupr[lu]/ch[lyu]);
            dirint     = fact1;
            fldir[lu]  = fabs(umu0)*fact1;
            rfldir[lu] = fabs(umu0)*fbeam*exp(-utau[lu]/ch[lyu]);
          }
          fact1      = ll[0][lyu]*exp( kk2[lyu]*(utaupr[lu]-taucpr[lyu+1]));
          fact2      = ll[1][lyu]*exp(-kk2[lyu]*(utaupr[lu]-taucpr[lyu  ]));
          u0c[0][lu] += fact2+rr[lyu]*fact1;
          u0c[1][lu] += fact1+rr[lyu]*fact2;

          // Calculate fluxes and mean intensities; downward and upward fluxes from eq. KST(9)
          fact1     = 2.*M_PI*cmu[0];
          fldn[lu]  = fact1*u0c[0][lu];
          flup[lu]  = fact1*u0c[1][lu];
          fdntot    = fldn[lu]+fldir[lu];
          rfldn[lu] = fdntot-rfldir[lu];

          // Mean intensity from eq. KST(10)
          uavg[lu] = u0c[0][lu]+u0c[1][lu];
          uavg[lu] = (2.*M_PI*uavg[lu]+dirint)/(4.*M_PI);

          // Flux divergence from eqs. KST(11-12)
          plsorc   = 1./(1.-oprim[lyu])*exp(-zpa[lyu]*utaupr[lu])*(xp0[lyu]+xp1[lyu]*utaupr[lu]);
          dfdt[lu] = (1.-ssalb[lyu])*4.*M_PI*(uavg[lu]-plsorc);
        }

        // Define intensities common to all phi from the average fluxes
        for (j = 0; j < nphi; j++) {
          if (numu2>1 || numu2==0) {
            uu[lkf][0][lu][j] = rfldn[lu]/M_PI;
            uu[lkf][1][lu][j] = flup[lu]/M_PI;
          } else {
            uu[lkf][1][lu][j] = flup[lu]/M_PI;
          }
        }
      }
      continue;
    }

    // ---------------------------------------------------------------
    // Matrix discrete ordinate solutions
    // Iterate across azimuthal components (eq STWJ 5, STWL 6)
    // ---------------------------------------------------------------
    if (fbeam<=0.0 || fabs(1.0-umu0)<1e-5 || (fabs(1.0-umu[0])<1e-5 && numu<=1)) naz=1; else naz=nstr;

    for (mazim = 0, kconv=0; mazim < naz; mazim++) {
      if (mazim == 0) delm0 = 1.0; else delm0 = 0.0;

      // Assign Fourier expansion coefficient of surface bidirectional reflectance
      if (ncut==nlyr) {
        for (iq = 0; iq < nmax; iq++) {
          bem[iq] = esurf[li]*(mazim==0 ? 1.0 : 0.0);
          for (jq = 0; jq < nmax; jq++) bdr[iq][jq] = asurf[li]*(mazim==0 ? 1.0 : 0.0);
        }
        for (iu = 0; iu < numu; iu++) {
          for (iq = 0; iq < nmax; iq++) rmu[iu][iq] = asurf[li]*(mazim==0 ? 1.0 : 0.0);
          emu[iu] = esurf[li]*(mazim==0 ? 1.0 : 0.0);
        }
      }

      // --------------------------------------------------------------
      // Iterate across layers
      // --------------------------------------------------------------
      for (lc = 0; lc < ncut; lc++) {
        // Solve eigenfunction problem in eq. STWJ(8B), STWL(23f); return eigenvalues and eigenvectors
        // Calculate quantities in eqs. SS(5-6), STWL(8b,15,23f)
        for (iq = 0; iq < nmax; iq++) {
          for (jq = 0; jq < nstr; jq++) {
            for (l = mazim, sum=0.0; l < nstr; l++) sum += gl[l][lc]*ylmc[mazim][l][iq]*ylmc[mazim][l][jq];
            cc[iq][jq] = 0.5*sum*cwt[jq];
          }
          for (jq = 0; jq < nmax; jq++) {
            // Fill remainder of array using symmetry relations  C(-mui,muj) = C(mui,-muj) and C(-mui,-muj) = C(mui,muj)
            cc[iq+nmax][jq     ] = cc[iq][jq+nmax];
            cc[iq+nmax][jq+nmax] = cc[iq][jq     ];
            // Get factors of coeff. matrix of reduced eigenvalue problem
            alpha       = cc[iq][jq     ]/cmu[iq];
            beta        = cc[iq][jq+nmax]/cmu[iq];
            amb[iq][jq] = alpha-beta;
            apb[iq][jq] = alpha+beta;
          }
          amb[iq][iq] -= 1.0/cmu[iq];
          apb[iq][iq] -= 1.0/cmu[iq];
        }

        // Finish calculation of coefficient matrix of reduced eigenvalue problem:
        // get matrix product (alpha+beta)*(alpha-beta); SS(12),STWL(23f)
        for (iq = 0; iq < nmax; iq++) {
          for (jq = 0; jq < nmax; jq++) {
            for (kq = 0, sum=0.0; kq < nmax; kq++) sum += apb[iq][kq]*amb[kq][jq];
            array[iq][jq] = sum;
          }
        }

        // Find (real) eigenvalues and eigenvectors
        asymtx(array, evecc, eval, nmax, nmax, nstr, wk);

        for (iq = 0; iq < nmax; iq++) {
          eval[iq]     = sqrt(fabs(eval[iq]));
          kk[iq+nmax  ][lc] =  eval[iq];
          kk[nmax-1-iq][lc] = -eval[iq]; // Add negative eigenvalue
        }

        // Find eigenvectors (G+) + (G-) from SS(10) and store temporarily in APB array
        for (jq = 0; jq < nmax; jq++) {
          for (iq = 0; iq < nmax; iq++) {
            for (kq = 0, sum=0.0; kq < nmax; kq++) sum += amb[iq][kq]*evecc[kq][jq];
            apb[iq][jq] = sum/eval[jq];
          }
        }

        for (jq = 0; jq < nmax; jq++) {
          for (iq = 0; iq < nmax; iq++) {
            gpplgm = apb[iq][jq];
            gpmigm = evecc[iq][jq];
            // Recover eigenvectors G+,G- from their sum and difference; stack them to get eigenvectors of full system
            // SS(7) (JQ = eigenvector number)
            evecc[iq     ][jq] = 0.5*(gpplgm+gpmigm);
            evecc[iq+nmax][jq] = 0.5*(gpplgm-gpmigm);
            // Eigenvectors corresponding to negative eigenvalues (corresp. to reversing sign of 'k' in SS(10) )
            gpplgm *= -1;
            evecc[iq     ][jq+nmax] = 0.5*(gpplgm+gpmigm);
            evecc[iq+nmax][jq+nmax] = 0.5*(gpplgm-gpmigm);
            gc[nmax+iq  ][nmax+jq  ][lc] = evecc[iq     ][jq     ];
            gc[nmax-iq-1][nmax+jq  ][lc] = evecc[iq+nmax][jq     ];
            gc[nmax+iq  ][nmax-jq-1][lc] = evecc[iq     ][jq+nmax];
            gc[nmax-iq-1][nmax-jq-1][lc] = evecc[iq+nmax][jq+nmax];
          }
        }

        // Find the incident-beam particular solution of SS(18), STWL(24a).
        if (fbeam > 0.0) {
          // Calculate x-sub-zero in STWJ(6d)
          for (iq = 0; iq < nstr; iq++) {
            for (k = mazim, sum=0.0; k <= nstr-1; k++) sum += gl[k][lc]*ylmc[mazim][k][iq]*ylm0[mazim][k][0];
            zj[iq] = (2.0 - delm0)*fbeam*sum/(4.0*M_PI);
          }
          q0a = exp( -chtau[lc  ] );
          q2a = exp( -chtau[lc+1] );

          // Calculate alfa coefficient
          deltat = taucpr[lc+1] - taucpr[lc];
          if (deltat>1e-7) tt = 1.0/deltat; else tt=0.0;
          xba[lc] = 1.0/ch[lc];
          if (fabs(xba[lc]*taucpr[lc]) > large || fabs(xba[lc]*taucpr[lc+1]) > large) xba[lc] = 0.0;
          if (fabs(xba[lc]*taucpr[lc]) > 50 || fabs(xba[lc]*taucpr[lc+1]) > 50) tt = 0.0; // Disable 1st derivative for large opacities
          // Dither alfa if it is close to one of the quadrature angles
          if (fabs(xba[lc]) > 0.00001) {
            for (iq = 0; iq < nmax; iq++) {
              if (fabs((fabs(xba[lc])-1.0/cmu[iq])/xba[lc] ) < 0.05 ) xba[lc] = xba[lc] * 1.001;
            }
          }

          for (iq = 0; iq < nstr; iq++) {
            q0 = q0a * zj[iq];
            q2 = q2a * zj[iq];
            // x-sub-zero and x-sub-one in Eqs. KS(48-49)
            xb1[iq][lc] = tt * (q2*exp(xba[lc]*taucpr[lc+1]) - q0*exp(xba[lc]*taucpr[lc]));
            xb0[iq][lc] = q0 * exp(xba[lc]*taucpr[lc]) - xb1[iq][lc]*taucpr[lc];
          }

          // Get coefficients at umu for pseudo-spherical source
          // Calculate x-sub-zero in STWJ(6d)
          for (iu=0; iu<numu; iu++) {
            for (k = mazim, sum=0.0; k <= nstr-1; k++) sum += gl[k][lc]*ylmu[mazim][k][iu]*ylm0[mazim][k][0];
            zju[iu] = (2.0 - delm0)*fbeam*sum/(4.0*M_PI);
            q0 = q0a*zju[iu]*exp(xba[lc]*taucpr[lc  ]);
            q2 = q2a*zju[iu]*exp(xba[lc]*taucpr[lc+1]);
            // x-sub-zero and x-sub-one in Eqs. KS(48-49)
            zb1u[iu][lc] = tt * (q2 - q0);
            zb0u[iu][lc] = q0 - zb1u[iu][lc]*taucpr[lc];
          }

          for (iq = 0; iq < nstr; iq++) {
            for (jq = 0; jq < nstr; jq++) array[iq][jq] = -cc[iq][jq];
            array[iq][iq] += 1.0 + xba[lc]*cmu[iq];
            zbs1[iq] = xb1[iq][lc];
          }

          // Find L-U (lower/upper triangular) decomposition of ARRAY and see if it is nearly singular
          sgeco(array,nstr,nstr,ipvt,&rcond,wk);

          // Dither alpha if rcond to small
          if (rcond < 1.0e-4) {
            if (xba[lc]==0.0) xba[lc]=5e-9;
            xba[lc] = xba[lc] * (1.0+5e-8);
            for (iq = 0; iq < nstr; iq++) {
              for (jq = 0; jq < nstr; jq++) array[iq][jq] = -cc[iq][jq];
              array[iq][iq] += 1.0 + xba[lc]*cmu[iq];
              zbs1[iq] = xb1[iq][lc];
            }
            // Solve linear equations KS(10-11) with dithered alpha
            sgeco(array,nstr,nstr,ipvt,&rcond,wk);
          }

          for (iq = 0; iq < nstr; iq++) wk[iq] = zbs1[iq];
          sgesl(array,nstr,nstr,ipvt,wk,0);
          for (iq = 0; iq < nstr; iq++) {
            zbs1[iq] = wk[iq];
            zbs0[iq] = xb0[iq][lc] + cmu[iq]*zbs1[iq];
          }
          for (iq = 0; iq < nstr; iq++) wk[iq] = zbs0[iq];
          sgesl(array,nstr,nstr,ipvt,wk,0);
          for (iq = 0; iq < nstr; iq++) zbs0[iq] = wk[iq];

          // Save the zbeam array
          zbeama[lc] = xba[lc];
          for (iq = 0; iq < nmax; iq++) {
            zbeam0[iq+nmax  ][lc] = zbs0[iq     ];
            zbeam1[iq+nmax  ][lc] = zbs1[iq     ];
            zbeam0[nmax-iq-1][lc] = zbs0[iq+nmax];
            zbeam1[nmax-iq-1][lc] = zbs1[iq+nmax];
          }
        }

        // Calculate particular solutions of eq. SS(15), STWL(25) for thermal emission source
        if (wsrc && mazim==0) {
          xr1[lc] = 0.0; xr0[lc] = 0.0;
          if (dtaucpr[lc] > 1e-4) xr1[lc] = (lsrc[lc+1]-lsrc[lc])/dtaucpr[lc];
          xr0[lc] = lsrc[lc] - xr1[lc]*taucpr[lc];

          // Finds the particular solution of thermal radiation of STWL(25)
          for (iq = 0; iq < nstr; iq++) {
            for (jq = 0; jq < nstr; jq++) array[iq][jq] = -cc[iq][jq];
            array[iq][iq] += 1.0;
            z1[iq] = (1.0 - oprim[lc])*xr1[lc];
          }
          // Solve linear equations: same as in upbeam, except zj replaced by z1 and z0
          sgeco(array,nstr,nstr,ipvt,&rcond,wk);

          for (iq = 0; iq < nstr; iq++) wk[iq] = z1[iq];
          sgesl(array,nstr,nstr,ipvt,wk,0);
          for (iq = 0; iq < nstr; iq++) z1[iq] = wk[iq];
          for (iq = 0; iq < nstr; iq++) z0[iq] = (1.0 - oprim[lc])*xr0[lc] + cmu[iq]*z1[iq];
          for (iq = 0; iq < nstr; iq++) wk[iq] = z0[iq];
          sgesl(array,nstr,nstr,ipvt,wk,0);
          for (iq = 0; iq < nstr; iq++) z0[iq] = wk[iq];
          for (iq = 0; iq < nmax; iq++) {
            zplk0[nmax+iq  ][lc] = z0[iq     ];
            zplk1[nmax+iq  ][lc] = z1[iq     ];
            zplk0[nmax-iq-1][lc] = z0[iq+nmax];
            zplk1[nmax-iq-1][lc] = z1[iq+nmax];
          }
        }

        // Interpolate eigenvectors to user angle
        for (iq = 0; iq < nstr; iq++) {
          for (l = mazim; l < nstr; l++) {
            // Inner sum in SD(8) times all factors in outer sum but PLM(mu)
            for (jq = 0, sum = 0.0; jq < nstr; jq++) sum += cwt[jq]*ylmc[mazim][l][jq]*evecc[jq][iq];
            wk[l] = 0.5*gl[l][lc]*sum;
          }
          // Finish outer sum in SD(8) and store eigenvectors
          for (iu=0; iu < numu; iu++) {
            for (l = mazim, sum = 0.0; l < nstr; l++) sum += wk[l]*ylmu[mazim][l][iu];
            if (iq < nmax) gu[iu][nmax+iq][lc] = sum; else gu[iu][nstr-1-iq][lc] = sum;
          }
        }

        // Interpolate source terms to user angle
        // Solar light
        if (fbeam > 0.0) {
          // Beam source terms; eq. SD(9)
          for (iq = mazim; iq < nstr; iq++) {
            psum0 = 0.0; psum1 = 0.0;
            for (jq = 0; jq < nstr; jq++) {
              psum0 += cwt[jq]*ylmc[mazim][iq][jq]*zbs0[jq];
              psum1 += cwt[jq]*ylmc[mazim][iq][jq]*zbs1[jq];
            }
            psi0[iq] = 0.5*gl[iq][lc]*psum0;
            psi1[iq] = 0.5*gl[iq][lc]*psum1;
          }
          for (iu=0; iu < numu; iu++) {
            psum0 = 0.0; psum1 = 0.0;
            for (iq = mazim; iq < nstr; iq++) {
              psum0 += ylmu[mazim][iq][iu]*psi0[iq];
              psum1 += ylmu[mazim][iq][iu]*psi1[iq];
            }
            zb0u[iu][lc] = psum0 + zb0u[iu][lc];
            zb1u[iu][lc] = psum1 + zb1u[iu][lc];
            zbau[iu][lc] = zbeama[lc];
          }
        }

        // Thermal source terms, STWJ(27c), STWL(31c)
        if (wsrc && mazim == 0) {
          for (iq = mazim; iq < nstr; iq++) {
            psum0 = 0.0; psum1 = 0.0;
            for (jq = 0; jq < nstr; jq++) {
              psum0 += cwt[jq]*ylmc[mazim][iq][jq]*z0[jq];
              psum1 += cwt[jq]*ylmc[mazim][iq][jq]*z1[jq];
            }
            psi0[iq] = 0.5*gl[iq][lc]*psum0;
            psi1[iq] = 0.5*gl[iq][lc]*psum1;
          }
          for (iu=0; iu < numu; iu++) {
            psum0 = 0.0; psum1 = 0.0;
            for (iq = mazim; iq < nstr; iq++) {
              psum0 += ylmu[mazim][iq][iu]*psi0[iq];
              psum1 += ylmu[mazim][iq][iu]*psi1[iq];
            }
            z0u[iu][lc] = psum0 + (1.0 - oprim[lc])*xr0[lc];
            z1u[iu][lc] = psum1 + (1.0 - oprim[lc])*xr1[lc];
          }
        }
      }

      // Calculate coefficient matrix for the set of equations obtained from the
      // boundary conditions and the continuity-of-intensity-at-layer-interface equations.
      // Store in the special banded-matrix format required by LINPACK routines
      int ncd = 3*nmax-1, lda = 3*ncd+1, nshift = lda-2*nstr+1, ncol = 0, jcol, irow, nncol;
      for (iq = 0; iq < 9*nmax-2; iq++) memset(cband[iq],0,(nstr*nlyr)*sizeof(double));

      // Use continuity conditions of eq. STWJ(17) to form coefficient matrix in STWJ(20);
      // employ scaling transformation STWJ(22)
      for (lc = 0; lc < ncut; lc++) {
        jcol = 0;
        for (iq = 0; iq < nmax; iq++) {
          ncol += 1;
          irow  = nshift-jcol;
          wk[iq] = exp(kk[iq][lc]*dtaucpr[lc]);
          for (jq = 0; jq < nstr; jq++) {
            cband[irow+nstr-1][ncol-1] =  gc[jq][iq][lc];
            cband[irow     -1][ncol-1] = -gc[jq][iq][lc]*wk[iq];
            irow++;
          }
          jcol++;
        }

        for (iq = nmax; iq < nstr; iq++) {
          ncol += 1;
          irow = nshift-jcol;
          for (jq = 0; jq < nstr; jq++) {
            cband[irow+nstr-1][ncol-1] =  gc[jq][iq][lc]*wk[nstr-1-iq];
            cband[irow     -1][ncol-1] = -gc[jq][iq][lc];
            irow++;
          }
          jcol++;
        }
      }

      // Use top boundary condition of STWJ(20a) for first layer
      jcol = 0;
      for (iq = 0; iq < nmax; iq++) {
        expa = exp(kk[iq][0]*taucpr[1]);
        irow = nshift-jcol+nmax;
        for (jq = nmax-1; jq >= 0; jq--) {
          cband[irow-1][jcol] = gc[jq][iq][0]*expa;
          irow++;
        }
        jcol++;
      }

      for (iq = nmax; iq < nstr; iq++) {
        irow = nshift-jcol+nmax;
        for (jq = nmax-1; jq >= 0; jq--) {
          cband[irow-1][jcol] = gc[jq][iq][0];
          irow++;
        }
        jcol++;
      }

      // Use bottom boundary condition of STWJ(20c) for last layer
      nncol = ncol-nstr;
      jcol  = 0;
      for (iq = 0; iq < nmax; iq++) {
        nncol++;
        irow = nshift-jcol+nstr;
        for (jq = nmax; jq < nstr; jq++) {
          if (ncut<nlyr || delm0==0.0) {
            // No azimuthal-dependent intensity if Lambert surface;
            // no intensity component if truncated bottom layer
            cband[irow-1][nncol-1] = gc[jq][iq][ncut-1];
          }
          else {
            for (k = 0, sum=0.0; k < nmax; k++) sum += cwt[k]*cmu[k]*bdr[jq-nmax][k]*gc[nmax-1-k][iq][ncut-1];
            cband[irow-1][nncol-1] = gc[jq][iq][ncut-1] - (1.0+delm0)*sum;
          }
          irow++;
        }
        jcol++;
      }

      for (iq = nmax; iq < nstr; iq++) {
        nncol++;
        irow = nshift-jcol+nstr;
        expa = wk[nstr-1-iq];
        for (jq = nmax; jq < nstr; jq++) {
          if (ncut<nlyr || delm0==0.0) {
            cband[irow-1][nncol-1] = gc[jq][iq][ncut-1]*expa;
          }
          else {
            for (k = 0, sum=0.0; k < nmax; k++) sum += cwt[k]*cmu[k]*bdr[jq-nmax][k]*gc[nmax-1-k][iq][ncut-1];
            cband[irow-1][nncol-1] = (gc[jq][iq][ncut-1]-(1.0+delm0)*sum)*expa;
          }
          irow++;
        }
        jcol++;
      }

      // Solve for constants of integration in homogeneous solution (general boundary conditions)
      // Construct B, STWJ(20a,c) for parallel beam+bottom
      // reflection+thermal emission at top and/or bottom
      memset(b,0,nstr*nlyr*sizeof(double));
      if (mazim > 0 && fbeam > 0.0) {
        // Azimuth-dependent case
        for (iq = 0; iq < nmax; iq++) {
          b[iq] = -zbeam0[nmax-1-iq][0]; // Top boundary
          b[ncol-nmax+iq] = -exp(-zbeama[ncut-1]*taucpr[ncut])*(zbeam0[iq+nmax][ncut-1] + zbeam1[iq+nmax][ncut-1]*taucpr[ncut]); // Bottom boundary
        }

        // Continuity condition for layer interfaces of eq. STWJ(20b)
        it = nmax; diff = 0;
        for (lc = 0; lc < ncut-1; lc++) {
          for (iq = 0; iq < nstr; iq++, it++) {
            b[it] = exp(-zbeama[lc+1]*taucpr[lc+1])*(zbeam0[iq][lc+1] + zbeam1[iq][lc+1]*taucpr[lc+1])
                   -exp(-zbeama[lc  ]*taucpr[lc+1])*(zbeam0[iq][lc  ] + zbeam1[iq][lc  ]*taucpr[lc+1]) + diff;
          }
        }
      } else {
        // Azimuth-independent case
        if (fbeam == 0.0) {
          for (iq = 0; iq < nmax; iq++) b[iq] = -zplk0[nmax-1-iq][0] + fisot + tplanck; // Top boundary
          if ( ncut<nlyr) {
            // No intensity component for truncated bottom layer
            for (iq = 0; iq < nmax; iq++) b[ncol-nmax+iq] = -zplk0[iq+nmax][ncut-1]-zplk1[iq+nmax][ncut-1]*taucpr[ncut]; // Bottom boundary
          } else {
            for (iq = 0; iq < nmax; iq++) {
              for (jq = 0, sum = 0.0; jq < nmax; jq++) sum += cwt[jq]*cmu[jq]*bdr[iq][jq]*(zplk0[nmax-1-jq][ncut-1]+zplk1[nmax-1-jq][ncut-1]*taucpr[ncut]);
              b[ncol-nmax+iq] = 2.0*sum + bem[iq]*bplanck - zplk0[iq+nmax][ncut-1]-zplk1[iq+nmax][ncut-1]*taucpr[ncut];
            }
          }
          // Continuity condition for layer interfaces, STWJ(20b)
          it = nmax;
          for (lc = 0; lc < ncut-1; lc++) {
            for (iq = 0; iq < nstr; iq++, it++) {
              b[it] = zplk0[iq][lc+1] - zplk0[iq][lc] + (zplk1[iq][lc+1]-zplk1[iq][lc])*taucpr[lc+1];
            }
          }
        } else {
          for (iq = 0; iq < nmax; iq++) b[iq] = -zbeam0[nmax-1-iq][0]-zplk0[nmax-1-iq][0] + fisot + tplanck;
          if (ncut<nlyr) {
            for (iq = 0; iq < nmax; iq++) {
              b[ncol-nmax+iq] = -exp(-zbeama[ncut-1]*taucpr[ncut])*(zbeam0[iq+nmax][ncut-1] + zbeam1[iq+nmax][ncut-1]*taucpr[ncut]) - zplk0[iq+nmax][ncut-1] - zplk1[iq+nmax][ncut-1]*taucpr[ncut];
            }
          } else {
            for (iq = 0; iq < nmax; iq++) {
              for (jq = 0, sum=0.0; jq < nmax; jq++) {
                sum += cwt[jq]*cmu[jq]*bdr[iq][jq]*
                       ( exp(-zbeama[ncut-1]*taucpr[ncut])*
                       ( zbeam0[nmax-1-jq][ncut-1] + zbeam1[nmax-1-jq][ncut-1]*taucpr[ncut])
                       + zplk0[nmax-1-jq][ncut-1] + zplk1[nmax-1-jq][ncut-1]*taucpr[ncut]);
              }
              if (umu0 <= 0.) refflx = 0.; else refflx = 1.0;
              b[ncol-nmax+iq] = 2.0*sum +
                                ( bdr[iq][0]*umu0*refflx*fbeam/M_PI)*expbea[ncut]
                                -  exp(-zbeama[ncut-1]*taucpr[ncut])*
                                (zbeam0[iq+nmax][ncut-1]+zbeam1[iq+nmax][ncut-1]*taucpr[ncut])
                                + bem[iq]*bplanck
                                -zplk0[iq+nmax][ncut-1]-zplk1[iq+nmax][ncut-1]*taucpr[ncut];
            }
          }
          it = nmax;
          for (lc = 0; lc < ncut-1; lc++) {
            for (iq = 0; iq < nstr; iq++, it++) {
              b[it] = exp(-zbeama[lc+1]*taucpr[lc+1])*
                      (zbeam0[iq][lc+1]+zbeam1[iq][lc+1]*taucpr[lc+1])
                      -exp(-zbeama[lc]*taucpr[lc+1])*
                      (zbeam0[iq][lc]+zbeam1[iq][lc]*taucpr[lc+1])
                      +zplk0[iq][lc+1]-zplk0[iq][lc]+
                      (zplk1[iq][lc+1]-zplk1[iq][lc])*taucpr[lc+1];
            }
          }
        }
      }

      // Find L-U (lower/upper triangular) decomposition of band matrix
      // CBAND and test if it is nearly singular (note: CBAND is
      // destroyed) (CBAND is in LINPACK packed format)
      rcond = 0.;
      ncd   = 3*nmax-1;
      sgbco(cband,(9*nmax-2),ncol,ncd,ncd,ipvt,&rcond,z);

      // Solve linear system with coeff matrix CBAND and R.H. side(s) B
      // after CBAND has been L-U decomposed. Solution is returned in B.
      sgbsl(cband,(9*nmax-2),ncol,ncd,ncd,ipvt,b,0);

      for (lc = 0; lc < ncut; lc++) {
        ipnt = (lc+1)*nstr-nmax;
        for (iq = 0; iq < nmax; iq++) {
          ll[nmax-iq-1][lc] = b[ipnt-iq-1];
          ll[nmax+iq  ][lc] = b[ipnt+iq  ];
        }
      }

      // Compute upward and downward fluxes (only for the first wavelength/bin)
      if (mazim == 0 && lf==l1 && bi==0) {
        for (lu = 0; lu < ntau; lu++) {
          for (iq = 0; iq < nstr; iq++) u0c[iq][lu] = 0.0;
          dfdt[lu] = 0.0; uavg[lu] = 0.0;
          flup[lu] = 0.0; fldn[lu] = 0.0; rfldn[lu] = 0.0;
          rfldir[lu] = 0.0; dirint = 0.0; fldir[lu] = 0.0;
          lyu = layru[lu];
          if (ncut<nlyr && lyu >= ncut) continue; // No radiation reaches this level

          if (fbeam>0.0) {
            fact   = exp(-utaupr[lu]/ch[lyu]);
            dirint = fbeam*fact;
            if (umu0>0) {
              rfldir[lu] = umu0*fbeam*exp(-utau[lu]/ch[lyu]);
              fldir[lu]  = umu0*fbeam*fact;
            }
          }

          for (iq = 0; iq < nmax; iq++) {
            for (jq = 0, sum=0.0; jq < nmax; jq++) sum += gc[iq][jq][lyu]*ll[jq][lyu]*exp(-kk[jq][lyu]*(utaupr[lu]-taucpr[lyu+1]));
            for (jq = nmax; jq < nstr; jq++)       sum += gc[iq][jq][lyu]*ll[jq][lyu]*exp(-kk[jq][lyu]*(utaupr[lu]-taucpr[lyu  ]));
            u0c[iq][lu] = sum;
            if (fbeam > 0.0) u0c[iq][lu] += exp(-zbeama[lyu]*utaupr[lu])*( zbeam0[iq][lyu] + zbeam1[iq][lyu]*utaupr[lu] );
            u0c[iq][lu]+= zplk0[iq][lyu]+zplk1[iq][lyu]*utaupr[lu];
            uavg[lu]   += cwt[nmax-1-iq]*u0c[iq][lu];
            fldn[lu]   += cwt[nmax-1-iq]*u0c[iq][lu]*cmu[nmax-1-iq];
          }

          for (iq = nmax; iq < nstr; iq++) {
            for (jq = 0, sum=0.0; jq < nmax; jq++) sum += gc[iq][jq][lyu]*ll[jq][lyu]*exp(-kk[jq][lyu]*(utaupr[lu]-taucpr[lyu+1]));
            for (jq = nmax; jq < nstr; jq++)       sum += gc[iq][jq][lyu]*ll[jq][lyu]*exp(-kk[jq][lyu]*(utaupr[lu]-taucpr[lyu  ]));
            u0c[iq][lu] = sum;
            if (fbeam > 0.0) u0c[iq][lu] += exp(-zbeama[lyu]*utaupr[lu])*( zbeam0[iq][lyu] + zbeam1[iq][lyu]*utaupr[lu] );
            u0c[iq][lu]+= zplk0[iq][lyu] + zplk1[iq][lyu]*utaupr[lu];
            uavg[lu]   += cwt[iq-nmax]*u0c[iq][lu];
            flup[lu]   += cwt[iq-nmax]*u0c[iq][lu]*cmu[iq-nmax];
          }

          flup[lu]  *= 2.0*M_PI;
          fldn[lu]  *= 2.0*M_PI;
          fdntot     = fldn[lu]+fldir[lu];
          rfldn[lu]  = fdntot-rfldir[lu];
          uavg[lu]   = (2.0*M_PI*uavg[lu]+dirint)/(4.0*M_PI);
          plsorc     = xr0[lyu]+xr1[lyu]*utaupr[lu];
          dfdt[lu]   = (1.0-ssalb[lyu])*4.0*M_PI*(uavg[lu]-plsorc);
        }
      }

      // Compute azimuthal intensity components at output angles
      // Incorporate constants of integration into interpolated eigenvectors
      double alfa, sgn, palint, plkint, bndint, bnddfu, bnddir, dfuint, dtau, dtau1, dtau2, exp1=0.0, exp2=0.0, f0n, f1n, denom, expn; int negumu, lyrstr, lyrend;
      for (lc = 0; lc < ncut; lc++) {
        for (iq = 0; iq < nstr; iq++) {
          for (iu = 0; iu < numu; iu++) {
            gu[iu][iq][lc] *= ll[iq][lc];
          }
        }
      }

      // Loop over levels at which intensities are desired ('user output levels')
      for (iu=0; iu < numu; iu++) memset(uum[iu],0,ntau*sizeof(double));
      for (lu = 0; lu < ntau; lu++) {
        lyu = layru[lu];
        if (ncut<nlyr && lyu >= ncut) continue;
        for (iu=0; iu < numu; iu++) {
          negumu = (umu[iu] < 0.0);
          if (negumu) {
            lyrstr = 0;
            lyrend = lyu-1;
            sgn    = -1.0;
          } else {
            lyrstr = lyu+1;
            lyrend = ncut-1;
            sgn    = 1.0;
          }

          // For downward intensity, integrate from top to LYU-1 in eq. S1(8); for upward,
          // integrate from bottom to LYU+1 in S1(9)
          palint = 0.;
          plkint = 0.;
          for (lc = lyrstr; lc <= lyrend; lc++) {
            dtau = dtaucpr[lc];
            exp1 = exp((utaupr[lu]-taucpr[lc  ])/umu[iu]);
            exp2 = exp((utaupr[lu]-taucpr[lc+1])/umu[iu]);

            if (wsrc && mazim == 0) {
              // Eqs. STWL(36b,c, 37b,c)
              f0n     = sgn*(exp1-exp2);
              f1n     = sgn*((taucpr[lc  ]+umu[iu])*exp1
                            -(taucpr[lc+1]+umu[iu])*exp2);
              plkint += z0u[iu][lc]*f0n + z1u[iu][lc]*f1n;
            }

            if (fbeam > 0.0) {
              denom  =  sgn*1.0/(zbau[iu][lc]*umu[iu]+1.0);
              palint += (zb0u[iu][lc]*denom*( exp(-zbau[iu][lc]*taucpr[lc  ])*exp1
                                             -exp(-zbau[iu][lc]*taucpr[lc+1])*exp2 )
                       + zb1u[iu][lc]*denom*((taucpr[lc  ]+sgn*denom*umu[iu])
                                           *exp(-zbau[iu][lc]*taucpr[lc  ])*exp1
                                            -(taucpr[lc+1]+sgn*denom*umu[iu] )
                                           *exp(-zbau[iu][lc]*taucpr[lc+1])*exp2));
            }

            // KK is negative
            for (iq = 0; iq < nmax; iq++) {
              wk[iq] = exp(kk[iq][lc]*dtau);
              denom  = 1.0 + umu[iu]*kk[iq][lc];
              if (fabs(denom) < 0.0001) {
                // L'Hospital limit
                expn = (dtau/umu[iu])*exp2;
              } else {
                expn = sgn*(exp1*wk[iq]-exp2)/denom;
              }
              palint += gu[iu][iq][lc]*expn;
            }

            // KK is positive
            for (iq = nmax; iq < nstr; iq++) {
              denom = 1.0 + umu[iu]*kk[iq][lc];
              if (fabs(denom) < 0.0001) {
                // L'Hospital limit
                expn = -(dtau/umu[iu])*exp1;
              } else {
                expn = sgn*(exp1-exp2*wk[nstr-1-iq])/denom;
              }
              palint += gu[iu][iq][lc]*expn;
            }
          }

          // Calculate contribution from user output level to next computational level
          dtau1 = utaupr[lu]-taucpr[lyu  ];
          dtau2 = utaupr[lu]-taucpr[lyu+1];

          if ((fabs(dtau1) >= 1.e-6 || !negumu) && (fabs(dtau2) >= 1.e-6 ||  negumu)) {
            if (negumu) exp1 = exp(dtau1/umu[iu]); else exp2 = exp(dtau2/umu[iu]);
            if (fbeam > 0.0) {
              if ( negumu) {
                expn = exp1;
                alfa = zbau[iu][lyu];
                denom = (-1.0/(alfa*umu[iu]+1.0));
                palint += zb0u[iu][lyu]*denom*(     -exp(-alfa*utaupr[lu   ])
                                               +expn*exp(-alfa*taucpr[lyu  ]))
                         +zb1u[iu][lyu]*denom*( -(utaupr[lu]-umu[iu]*denom)*     exp(-alfa*utaupr[lu])
                                                +(taucpr[lyu  ]-umu[iu]*denom)*expn*exp(-alfa*taucpr[lyu ]));
              } else {
                expn = exp2;
                alfa = zbau[iu][lyu];
                denom = (1.0/(alfa*umu[iu]+1.0));
                palint += zb0u[iu][lyu]*denom*(     +exp(-alfa*utaupr[lu   ])
                                               -expn*exp(-alfa*taucpr[lyu+1]))
                         +zb1u[iu][lyu]*denom*(  (utaupr[lu]+umu[iu]*denom)*     exp(-alfa*utaupr[lu])
                                                 -(taucpr[lyu+1]+umu[iu]*denom)*expn*exp(-alfa*taucpr[lyu+1]));
              }
            }

            // KK is negative
            dtau = dtaucpr[lyu];
            for (iq = 0; iq < nmax; iq++) {
              denom = 1.0 + umu[iu]*kk[iq][lyu];
              if (fabs(denom) < 0.0001) expn = -dtau2/umu[iu]*exp2;
              else if (negumu) expn = (exp(-kk[iq][lyu]*dtau2)-exp(kk[iq][lyu]*dtau )*exp1)/denom;
              else expn = (exp(-kk[iq][lyu]*dtau2)-exp2)/denom;
              palint += gu[iu][iq][lyu]*expn;
            }

            // KK is positive
            for (iq = nmax; iq < nstr; iq++) {
              denom = 1.0+umu[iu]*kk[iq][lyu];
              if (fabs(denom) < 0.0001) expn = -(dtau1/umu[iu])*exp1;
              else if (negumu) expn = (exp(-kk[iq][lyu]*dtau1)-exp1)/denom;
              else expn = (exp(-kk[iq][lyu]*dtau1)-exp(-kk[iq][lyu]*dtau )*exp2)/denom;
              palint += gu[iu][iq][lyu]*expn;
            }

            if (wsrc && mazim == 0) {
              // Eqs. STWL (35-37) with tau-sub-n-1 replaced by tau for upward, and
              // tau-sub-n replaced by tau for downward directions
              if (negumu) {
                expn = exp1;
                fact = taucpr[lyu  ]+umu[iu];
              } else {
                expn = exp2;
                fact = taucpr[lyu+1]+umu[iu];
              }
              f0n     = 1.0-expn;
              f1n     = utaupr[lu]+umu[iu]-fact*expn;
              plkint += z0u[iu][lyu]*f0n + z1u[iu][lyu]*f1n;
            }
          }

          // Calculate intensity components attenuated at both boundaries.
          // Note: no azimuthal intensity component for isotropic surface
          bndint = 0.;
          if (negumu && mazim == 0) {
            bndint = (fisot+tplanck)*exp(utaupr[lu]/umu[iu]);
          } else if (!negumu) {
            if (ncut<nlyr || mazim>0) {
              uum[iu][lu] = palint + plkint;
              continue;
            }
            for (jq = nmax; jq < nstr; jq++) wk[jq] = exp(-kk[jq][nlyr-1]*dtaucpr[nlyr-1]);
            bnddfu = 0.;
            for (iq = nmax-1; iq >= 0; iq--) {
              dfuint = 0.;
              for (jq = 0;   jq < nmax; jq++) dfuint += gc[iq][jq][nlyr-1]*ll[jq][nlyr-1];
              for (jq= nmax; jq < nstr; jq++) dfuint += gc[iq][jq][nlyr-1]*ll[jq][nlyr-1]*wk[jq];
              if (fbeam > 0.) dfuint += exp(-zbeama[nlyr-1]*taucpr[nlyr]) * (zbeam0[iq][nlyr-1] + zbeam1[iq][nlyr-1]*taucpr[nlyr]);
              dfuint += delm0*(zplk0[iq][nlyr-1]+zplk1[iq][nlyr-1]*taucpr[nlyr]);
              bnddfu += (1.+delm0)*rmu[iu][nmax-1-iq]*cmu[nmax-1-iq]*cwt[nmax-1-iq]*dfuint;
            }
            bnddir = 0.;
            if (fbeam > 0. && umu0 >0.) bnddir = umu0*fbeam/M_PI*rmu[iu][0]*expbea[nlyr];
            bndint = (bnddfu+bnddir+delm0*emu[iu]*bplanck)*exp((utaupr[lu]-taucpr[nlyr])/umu[iu]);
          }
          uum[iu][lu] = palint+plkint+bndint;
        }
      }

      // Save azimuthally averaged intensities
      // And see if we have reached convergance
      if (mazim == 0) {
        for (lu = 0; lu < ntau; lu++) {
          for (iu = 0; iu < numu; iu++) {
            for (j = 0; j < nphi; j++) {
              uu[lkf][iu][lu][j] = uum[iu][lu];
            }
          }
        }
      } else {
        // Increment intensity by current azimuthal component (Fourier cosine series);  eq SD(2), STWL(6)
        azerr = 0.;
        for (j=0; j<nphi; j++) {
          cosphi = cos((double)mazim*dphi[j]*M_PI/180.0);
          for (lu = 0; lu < ntau; lu++) {
            for (iu = 0; iu < numu; iu++) {
              azterm = uum[iu][lu]*cosphi;
              uu[lkf][iu][lu][j] += azterm;
              cratio = fabs(azterm)/(fabs(uu[lkf][iu][lu][j]) + 1e-16);
              azerr  = (azerr > cratio ? azerr : cratio);
            }
          }
        }
        if (azerr <= 1e-5) kconv++;
        if (kconv >= 2) break;
      }
    }

    // Apply Nakajima/Tanaka intensity correction
    if (fbeam>0.0 && deltam) {
      double dtheta=10.0, theta0, thetap, ctheta, plm1, plm2, pl, wbar, fbar, stau;
      double ussndm=0.0, ussp=0.0, gbar, pspike, tmp, umu1, umu2, x1, ans, exp0, exp1=0.0, xi, duims, *phase, *omega; int lyr;
      theta0 = acos(-umu0)*180.0/M_PI;

      // Iterate across all zenith angles
      for (iu = 0; iu < numu; iu++) {
        // Iterate across all azimuthal angles
        for (j=0; j<nphi; j++) {
          // Calculate cosine of scattering angle, eq. STWL(4)
          ctheta = -umu0*umu[iu] + sqrt((1.0 - umu0*umu0)*(1.0 - umu[iu]*umu[iu]))*cos(dphi[j]*M_PI/180.0);

          // Initialize phase function
          for (lc = 0; lc < ncut; lc++) { phasa[lc] = 1.0; phasm[lc] = 1.0; }

          // Initialize Legendre poly. recurrence
          plm1 = 1.0; plm2 = 0.0;
          for (k = 1; k <= nmom; k++) {
            // Calculate Legendre polynomial of P-sub-l by upward recurrence
            pl   = ((double)(2.0*k - 1.0)*ctheta*plm1 - (double)(k-1.0)*plm2)/k;
            plm2 = plm1;
            plm1 = pl;

            // Calculate actual phase function
            for (lc = 0; lc < ncut; lc++) phasa[lc] += (double)(2.0*k + 1.0)*pl*pmom[k][lc];
            // Calculate delta-M transformed phase function
            if (k <= nstr-1) {
              for (lc = 0; lc < ncut; lc++) {
                phasm[lc] += (double)(2.0*k + 1.0)*pl*(pmom[k][lc] - flyr[lc])/(1.0 - flyr[lc]);
              }
            }
          }

          // Apply TMS method, eq. STWL(68)
          for (lc = 0; lc < ncut; lc++) {
            phast[lc] = phasa[lc]/(1.0 - flyr[lc]*ssalb[lc]);
          }
          for (lu = 0; lu < ntau; lu++) {
            if (ncut==nlyr || layru[lu] < ncut) {
              for (k=0;k<2;k++) {
                // Single scattering component
                if (k==0) { phase=phast; omega=ssalb; } else { phase=phasm; omega=oprim; }
                lyu  = layru[lu];
                tt = utaupr[lu];
                ans  = 0.0;
                exp0=exp(-tt/umu0);

                if (fabs(umu[iu]+umu0) <= 100.*DBL_EPSILON) {
                  // Calculate downward intensity when umu=umu0, eq. STWL (65e)
                  for (lyr = 0; lyr <= lyu-1; lyr++) ans += omega[lyr]*phase[lyr]*(taucpr[lyr+1]-taucpr[lyr]);
                  ans = fbeam/(4.0*M_PI*umu0)*exp0*(ans + omega[lyu]*phase[lyu]*(tt-taucpr[lyu]));
                } else {
                  if (umu[iu] > 0.0) {
                    // Upward intensity, eq. STWL (65b)
                    for (lyr = lyu; lyr < nlyr; lyr++) {
                      exp1  = exp(-((taucpr[lyr+1]-tt)/umu[iu]+taucpr[lyr+1]/umu0));
                      ans  += omega[lyr]*phase[lyr]*(exp0-exp1);
                      exp0  = exp1;
                    }
                  } else {
                    // Downward intensity, eq. STWL (65d)
                    for (lyr = lyu; lyr >= 0; lyr--) {
                      exp1  = exp(-((taucpr[lyr]-tt)/umu[iu]+taucpr[lyr]/umu0));
                      ans  += omega[lyr]*phase[lyr]*(exp0-exp1);
                      exp0  = exp1;
                    }
                  }
                  ans *= fbeam/(4.*M_PI*(1.0+umu[iu]/umu0));
                }
                if (k==0) ussndm=ans; else ussp=ans;
              }
              uu[lkf][iu][lu][j] += ussndm - ussp;
            }
          }

          thetap = acos(umu[iu])*180.0/M_PI;
          if (umu[iu] < 0.0 && fabs(theta0-thetap) <= dtheta) {
            // Emerging direction is in the aureole (theta0 +/- dtheta).
            // Apply IMS method for correction of secondary scattering below top level.
            lc = 0;
            if (utau[0] <= 100.*DBL_EPSILON) lc = 1;
            for (lu = lc; lu < ntau; lu++) {
              if (ncut==nlyr || layru[lu] < ncut) {
                // Calculate secondary scattering
                // Calculate vertically averaged value of single scattering albedo and separated fraction f, eq. STWL (A.15)
                lyu  = layru[lu];
                tt   = utau[lu]-tauc[lyu];
                wbar = ssalb[lyu]*tt;
                fbar = flyr[lyu]*wbar;
                stau = tt;
                for (lyr = 0; lyr <= lyu-1; lyr++) {
                  wbar += dtauc[lyr]*ssalb[lyr];
                  fbar += dtauc[lyr]*ssalb[lyr]*flyr[lyr];
                  stau += dtauc[lyr];
                }

                if (wbar <= 1e-4 || fbar <= 1e-4 || stau <= 1e-4 || fbeam <= 1e-4) continue;
                fbar /= wbar;
                wbar /= stau;

                // Calculate pspike = (2p"-p"*p")
                pspike = 1.0; gbar = 1.0; plm1 = 1.0; plm2 = 0.0;
                // pspike for l <= 2n-1
                for (k = 1; k < nstr; k++) {
                  pl      = ((double)(2.0*k-1.0)*ctheta*plm1 - (double)(k-1.0)*plm2)/k;
                  plm2    = plm1;
                  plm1    = pl;
                  pspike += gbar*(2.0-gbar)*(double)(2.0*k+1.0)*pl;
                }

                // pspike for l > 2n-1
                for (k = nstr; k <= nmom; k++) {
                  pl   = ((double)(2.0*k-1.0)*ctheta*plm1 - (double)(k-1.0)*plm2)/k;
                  plm2 = plm1;
                  plm1 = pl;
                  tt   = utau[lu]-tauc[lyu];
                  gbar = pmom[k][lyu]*ssalb[lyu]*tt;
                  for (lyr = 0; lyr <= lyu-1; lyr++) gbar += pmom[k][lyr]*ssalb[lyr]*dtauc[lyr];
                  tmp = fbar*wbar*stau;
                  if (tmp <= 1e-4) gbar = 0.0; else gbar /= tmp;
                  pspike += gbar*(2.0-gbar)*(double)(2.0*k+1.0)*pl;
                }
                umu1 = -umu[iu];
                umu2 = umu0/(1.0-fbar*wbar);
                tt   = utau[lu];
                x1   = (umu2-umu1)/(umu2*umu1);
                exp1 = exp(-tt/umu1);
                if (x1 != 0.0) {
                  xi = ((tt*x1-1.0)*exp(-tt/umu2)+exp1)/(x1*x1*umu1*umu2);
                } else {
                  xi = tt*tt*exp1/(2.0*umu1*umu2);
                }
                duims = fbeam/(4.*M_PI)*pow(fbar*wbar,2.0)/(1.0-fbar*wbar)*pspike*xi;

                // Calculate IMS correction term, eq. STWL (A.13)
                uu[lkf][iu][lu][j] -= duims;
              }
            }
          }
        }
      }
    }
  }

  // Free operational arrays
  free(b); free2D(cband,9*nmax-2); free(ch); free(chtau); free3D(chscl,nlyr,nlyr); free(cmu); free(cwt);
  free(dtauc); free(dtaucpr); free(expbea); free(flyr); free(ggprim); free2D(gl,nstr); free3D(gu,numu,nstr);
  free(ipvt); free(layru); free2D(ll,nstr); free(lsrc); free(oprim); free2D(pmom,nstr+1); free(ssalb);
  free(tauc); free(taucpr); free2D(u0c,nstr); free(utaupr); free3D(ylm0,nstr,nstr);
  free3D(ylmc,nstr,nstr); free3D(ylmu,nstr,nstr);

  if (twostr) {
    // Two-streams arrays
    free(diag); free(kk2); free(rr); free(subd); free(superd); free(xp0); free(xp1); free(yb0d); free(yb0u);
    free(yb1d); free(yb1u); free(yp0d); free(yp0u); free(yp1d); free(yp1u); free(zba); free(zpa);
  } else {
    // Full discrete-ordinates arrays
    free2D(amb,nmax); free2D(apb,nmax); free2D(array,nstr); free2D(bdr,nmax); free(bem); free2D(cc,nstr);
    free(emu); free(eval); free2D(evecc,nstr); free3D(gc,nstr,nstr); free2D(kk,nstr);
    free(phasa); free(phasm); free(phast); free(psi0); free(psi1); free2D(rmu,numu); free2D(uum,numu); free(wk);
    free(xba); free2D(xb0,nstr); free2D(xb1,nstr); free(xr0); free(xr1);  free(z); free(z0);
    free2D(z0u,numu); free(z1); free2D(z1u,numu); free2D(zb0u,numu); free2D(zb1u,numu); free2D(zbau,numu);
    free(zbs0); free(zbs1); free(zbeama); free2D(zbeam0,nstr); free2D(zbeam1,nstr); free(zj); free(zju);
    free2D(zplk0,nstr); free2D(zplk1,nstr);
  }

  return;
}


// -- Differential equations solver  ----------------------------------------
// Solves eigenfunction problem for real asymmetric matrix for which it
// is known a priori that the eigenvalues are real. This is an adaptation
// of a subroutine EIGRF in the IMSL library to use real instead of complex
// arithmetic, accounting for the known fact that the eigenvalues and
// eigenvectors in the discrete ordinate solution are real.

// EIGRF is based primarily on EISPACK routines.  The matrix is first
// balanced using the Parlett-Reinsch algorithm.  Then the Martin-Wilkinson
// algorithm is applied. There is a statement 'j = wk(i)' that converts a
// double precision variable to an integer variable; this seems dangerous
// to us in principle, but seems to work fine in practice.
// -------------------------------------------------------------------------
void asymtx(double **aa, double **evec, double *eval, int m, int ia, int ievec, double *wk) {

  const double c1 = 0.4375, c2 = 0.5, c3 = 0.75, c4 = 0.95, c5 = 16.0, c6 = 256.0;
  int noconv,notlas,i,ii,in,j,k,ka,kkk,l,lb=0,lll,n,n1,n2;
  double col,discri,f,g,h,p=0,q=0,r=0,repl,rnorm,row,s,scale,sgn,t,tol=DBL_EPSILON,uu,vv,w,x,y,z;

  // Handle 1x1 and 2x2 special cases
  if (m == 1) {
    eval[0]   = aa[0][0];
    evec[0][0] = 1.;
    return;

  } else if (m == 2) {
    discri = aa[0][0]-aa[1][1];
    discri = discri*discri + 4.0*aa[0][1]*aa[1][0];
    sgn = 1.0;
    if (aa[0][0] < aa[1][1]) sgn = -1.0;
    eval[0] = 0.5*(aa[0][0]+aa[1][1]+sgn*sqrt(discri));
    eval[1] = 0.5*(aa[0][0]+aa[1][1]-sgn*sqrt(discri));
    evec[0][0] = 1.;
    evec[1][1] = 1.;
    if (aa[0][0] == aa[1][1] && (aa[1][0] == 0. || aa[0][1] == 0.)) {
      rnorm     = fabs(aa[0][0])+fabs(aa[0][1])+fabs(aa[1][0])+fabs(aa[1][1]);
      w         = tol*rnorm;
      evec[1][0] =  aa[1][0]/w;
      evec[0][1] = -aa[0][1]/w;
    } else {
      evec[1][0] = aa[1][0]/(eval[0]-aa[1][1]);
      evec[0][1] = aa[0][1]/(eval[1]-aa[0][0]);
    }
    return;
  }

  // Initialize output variables
  memset(eval,0,m*sizeof(double));
  for (i = 1; i <= ievec; i++) memset(evec[i-1],0,ievec*sizeof(double));
  for (i = 1; i <= m; i++) evec[i-1][i-1] = 1.0;

  // Balance the input matrix and reduce its norm by diagonal similarity transformation stored in wk
  // then search for rows isolating an eigenvalue and push them down.
  rnorm = 0.; l = 1; k = m;

S50:
  kkk = k;
  for (j = kkk; j >= 1; j--) {
    row = 0.;
    for (i = 1; i <= k; i++) {
      if (i != j) row += fabs(aa[j-1][i-1]);
    }
    if (row == 0.) {
      wk[k-1] = (double)j;
      if (j != k) {
        for (i = 1; i <= k; i++) {
          repl    = aa[i-1][j-1];
          aa[i-1][j-1] = aa[i-1][k-1];
          aa[i-1][k-1] = repl;
        }
        for (i = l; i <= m; i++) {
          repl    = aa[j-1][i-1];
          aa[j-1][i-1] = aa[k-1][i-1];
          aa[k-1][i-1] = repl;
        }
      }
      k--;
      goto S50;
    }
  }

  // Search for columns isolating an eigenvalue and push them left
S100:
  lll = l;
  for (j = lll; j <= k; j++) {
    col = 0.;
    for (i = l; i <= k; i++) {
      if (i != j) col += fabs(aa[i-1][j-1]);
    }
    if (col == 0.) {
      wk[l-1] = (double)j;
      if (j != l) {
        for (i = 1; i <= k; i++) {
          repl    = aa[i-1][j-1];
          aa[i-1][j-1] = aa[i-1][l-1];
          aa[i-1][l-1] = repl;
        }
        for (i = l; i <= m; i++) {
          repl    = aa[j-1][i-1];
          aa[j-1][i-1] = aa[l-1][i-1];
          aa[l-1][i-1] = repl;
        }
      }
      l++;
      goto S100;
    }
  }

  // Balance the submatrix in rows L through K
  for (i = l; i <= k; i++) wk[i-1] = 1.;
  noconv = 1;
  while (noconv) {
    noconv = 0;
    for (i = l; i <= k; i++) {
      col = 0.; row = 0.;
      for (j = l; j <= k; j++) {
        if (j != i) {
          col += fabs(aa[j-1][i-1]);
          row += fabs(aa[i-1][j-1]);
        }
      }

      f = 1.;
      g = row/c5;
      h = col+row;
      while (col < g) {
        f   *= c5;
        col *= c6;
      }
      g = row*c5;
      while (col >= g) {
        f   /= c5;
        col /= c6;
      }

      // Now balance
      if ((col+row)/f < c4*h) {
        wk[i-1]  *= f;
        noconv  = 1;
        for (j = l; j <= m; j++) aa[i-1][j-1] /= f;
        for (j = 1; j <= k; j++) aa[j-1][i-1] *= f;
      }
    }
  }

  if (k-1 >= l+1) {
    // Transfer A to a Hessenberg form.
    for (n = l+1; n <= k-1; n++) {
      h       = 0.;
      wk[n+m-1] = 0.;
      scale   = 0.;
      // Scale column
      for (i = n; i <= k; i++) scale += fabs(aa[i-1][n-1-1]);
      if (scale != 0.) {
        for (i = k; i >= n; i--) {
          wk[i+m-1]  = aa[i-1][n-1-1]/scale;
          h       += pow(wk[i+m-1],2.0);
        }
        g        = -(wk[n+m-1]>=0.0 ? fabs(sqrt(h)) : -fabs(sqrt(h)));
        h       -= wk[n+m-1]*g;
        wk[n+m-1] -= g;
        // Form (I-(U*UT)/H)*A
        for (j = n; j <= m; j++) {
          f = 0.;
          for (i = k; i >= n; i--) f += wk[i+m-1]*aa[i-1][j-1];
          for (i = n; i <= k; i++) aa[i-1][j-1] -= wk[i+m-1]*f/h;
        }
        // Form (i-(u*ut)/h)*a*(i-(u*ut)/h)
        for (i = 1; i <= k; i++) {
          f = 0.;
          for (j = k; j >= n; j--) f += wk[j+m-1]*aa[i-1][j-1];
          for (j = n; j <= k; j++) aa[i-1][j-1] -= wk[j+m-1]*f/h;
        }
        wk[n+m-1]   *= scale;
        aa[n-1][n-1-1] = scale*g;
      }
    }

    for (n = k-2; n >= l; n--) {
      n1 = n+1;
      n2 = n+2;
      f = aa[n+1-1][n-1];
      if( f != 0.) {
        f *= wk[n+1+m-1];
        for (i = n+2; i <= k; i++) wk[i+m-1] = aa[i-1][n-1];
        if (n+1 <= k) {
          for (j = 1; j <= m; j++) {
            for (i = n+1, g=0.0; i <= k; i++) g += wk[i+m-1]*evec[i-1][j-1];
            for (i = n+1, g/=f;  i <= k; i++) evec[i-1][j-1] += g*wk[i+m-1];
          }
        }
      }
    }
  }

  n = 1;
  for (i = 1; i <= m; i++) {
    for (j = n; j <= m; j++) rnorm += fabs(aa[i-1][j-1]);
    n = i;
    if (i < l || i > k) eval[i-1] = aa[i-1][i-1];
  }

  n = k;
  t = 0.;

  // Search for next eigenvalues
S400:
  if (n < l) goto S550;
  in = 0;
  n1 = n-1;
  n2 = n-2;

  // Look for single small sub-diagonal element
S410:
  for (i = l; i <= n; i++) {
    lb = n+l-i;
    if (lb == l) break;
    s = fabs(aa[lb-1-1][lb-1-1])+fabs(aa[lb-1][lb-1]);
    if (s == 0.) s = rnorm;
    if (fabs(aa[lb-1][lb-1-1]) <= tol*s) break;
  }

  x = aa[n-1][n-1];
  if (lb == n) {
    // One eigenvalue found
    aa[n-1][n-1] = x+t;
    eval[n-1] = aa[n-1][n-1];
    n       = n1;
    goto S400;
  }

  y = aa[n1-1][n1-1];
  w = aa[n-1][n1-1]*aa[n1-1][n-1];

  if (lb == n1) {
    // Two eigenvalues found
    p         = (y-x)*c2;
    q         = p*p+w;
    z         = sqrt(fabs(q));
    aa[n-1][n-1]   = x+t;
    x         = aa[n-1][n-1];
    aa[n1-1][n1-1] = y+t;
    // Real pair
    z        = p + (p>=0.0 ? fabs(z) : -fabs(z));
    eval[n1-1] = x+z;
    eval[n-1]  = eval[n1-1];

    if (z != 0.0) eval[n-1] = x-w/z;
    x = aa[n-1][n1-1];
    // Employ scale factor in case X and Z are very small
    r = sqrt(x*x+z*z);
    p = x/r;
    q = z/r;
    // Row modification
    for (j = n1; j <= m; j++) {
      z        = aa[n1-1][j-1];
      aa[n1-1][j-1] =  q*z+p*aa[n-1][j-1];
      aa[n-1 ][j-1] = -p*z+q*aa[n-1][j-1];
    }
    // Column modification
    for (i = 1; i <= n; i++) {
      z        = aa[i-1][n1-1];
      aa[i-1][n1-1] =  q*z+p*aa[i-1][n-1];
      aa[i-1][ n-1] = -p*z+q*aa[i-1][n-1];
    }
    // Accumulate transformations
    for (i = l; i <= k; i++) {
      z          = evec[i-1][n1-1];
      evec[i-1][n1-1] =  q*z+p*evec[i-1][n-1];
      evec[i-1][ n-1] = -p*z+q*evec[i-1][n-1];
    }
    n = n2;
    goto S400;
  }

  if (in == 30) {
     // No convergence after 30 iterations.
    return;
  }

  // Form shift
  if (in == 10 || in == 20) {
    t += x;
    for (i = l; i <= n; i++) {
      aa[i-1][i-1] -= x;
    }
    s = fabs(aa[n-1][n1-1])+fabs(aa[n1-1][n2-1]);
    x = c3*s;
    y = x;
    w = -c1*s*s;
  }

  in++;

  // Look for two consecutive small sub-diagonal elements
  for (j = lb; j <= n2; j++) {
    i  = n2+lb-j;
    z  = aa[i-1][i-1];
    r  = x-z;
    s  = y-z;
    p  = (r*s-w)/aa[i+1-1][i-1]+aa[i-1][i+1-1];
    q  = aa[i+1-1][i+1-1]-z-r-s;
    r  = aa[i+2-1][i+1-1];
    s  = fabs(p)+fabs(q)+fabs(r);
    p /= s;
    q /= s;
    r /= s;

    if (i == lb) break;

    uu = fabs(aa[i-1][i-1-1])*(fabs(q)+fabs(r));
    vv = fabs(p)*(fabs(aa[i-1-1][i-1-1])+fabs(z)+fabs(aa[i+1-1][i+1-1]));

    if (uu <= tol*vv) break;
  }

  aa[i+2-1][i-1] = 0.;
  for (j = i+3; j <= n; j++) {
    aa[j-1][j-2-1] = 0.;
    aa[j-1][j-3-1] = 0.;
  }

  // Double QR step involving rows K to N and columns M to N
  for (ka = i; ka <= n1; ka++) {
    notlas = (ka != n1);
    if (ka == i) {
      s = (p>=0.0 ? fabs(sqrt(p*p+q*q+r*r)) : -fabs(sqrt(p*p+q*q+r*r)));
      if (lb != i) aa[ka-1][ka-1-1] *= -1;
    }
    else {
      p = aa[ka-1][ka-1-1];
      q = aa[ka+1-1][ka-1-1];
      r = 0.;
      if (notlas) r = aa[ka+2-1][ka-1-1];
      x = fabs(p)+fabs(q)+fabs(r);
      if (x == 0.) continue;
      p /= x;
      q /= x;
      r /= x;
      s  = (p>=0.0 ? fabs(sqrt(p*p+q*q+r*r)) : -fabs(sqrt(p*p+q*q+r*r)));
      aa[ka-1][ka-1-1] = -s*x;
    }

    p += s;
    x  = p/s;
    y  = q/s;
    z  = r/s;
    q /= p;
    r /= p;

    // Row modification
    for (j = ka; j <= m; j++) {
      p = aa[ka-1][j-1]+q*aa[ka+1-1][j-1];
      if (notlas) {
        p          += r*aa[ka+2-1][j-1];
        aa[ka+2-1][j-1] -= p*z;
      }
      aa[ka+1-1][j-1] -= p*y;
      aa[ka-1][j-1]   -= p*x;
    }

    // Column modification
    for (ii = 1; ii <= (n<ka+3?n:ka+3); ii++) {
      p = x*aa[ii-1][ka-1]+y*aa[ii-1][ka+1-1];
      if (notlas) {
        p           += z*aa[ii-1][ka+2-1];
        aa[ii-1][ka+2-1] -= p*r;
      }
      aa[ii-1][ka+1-1] -= p*q;
      aa[ii-1][ka-1]   -= p;
    }

    // Accumulate transformations
    for (ii = l; ii <= k; ii++) {
      p = x*evec[ii-1][ka-1]+y*evec[ii-1][ka+1-1];
      if (notlas) {
        p                  += z*evec[ii-1][ka+2-1];
        evec[ii-1][ka+2-1] -= p*r;
      }
      evec[ii-1][ka+1-1] -= p*q;
      evec[ii-1][ka-1]   -= p;
    }
  }

  goto S410;

  // All evals found, now backsubstitute real vector
S550:
  if (rnorm != 0.) {
    for (n = m; n >= 1; n--) {
      n2      = n;
      aa[n-1][n-1] = 1.;
      for (i = n-1; i >= 1; i--) {
        w = aa[i-1][i-1]-eval[n-1];
        if (w == 0.) w = tol*rnorm;
        r = aa[i-1][n-1];
        for (j = n2; j <= n-1; j++) r += aa[i-1][j-1]*aa[j-1][n-1];
        aa[i-1][n-1] = -r/w;
        n2      = i;
      }
    }
    // End backsubstitution vectors of isolated evals
    for (i = 1; i <= m; i++) {
      if (i < l || i > k) {
        for (j = i; j <= m; j++) {
          evec[i-1][j-1] = aa[i-1][j-1];
        }
      }
    }
    // Multiply by transformation matrix
    if (k != 0) {
      for (j = m; j >= l; j--) {
        for (i = l; i <= k; i++) {
          z = 0.;
          for (n = l; n <= (j<k?j:k); n++) z += evec[i-1][n-1]*aa[n-1][j-1];
          evec[i-1][j-1] = z;
        }
      }
    }
  }
  for (i = l; i <= k; i++) {
    for (j = 1; j <= m; j++) {
      evec[i-1][j-1] *= wk[i-1];
    }
  }

  // Interchange rows if permutations occurred
  for (i = l-1; i >= 1; i--) {
    j = wk[i-1];
    if (i != j) {
      for (n = 1; n <= m; n++) {
        repl      = evec[i-1][n-1];
        evec[i-1][n-1] = evec[j-1][n-1];
        evec[j-1][n-1] = repl;
      }
    }
  }
  for (i = k+1; i <= m; i++) {
    j = wk[i-1];
    if (i != j) {
      for (n = 1; n <= m; n++) {
        repl      = evec[i-1][n-1];
        evec[i-1][n-1] = evec[j-1][n-1];
        evec[j-1][n-1] = repl;
      }
    }
  }

  return;
}

// --------------------------------------------------------
// Matrix operations and inversions
// --------------------------------------------------------
int isamax(int n, double **a, int i0, int j0) {
  register int ans=0,i;
  double smax,xmag;

  if (n <= 0) {
    ans = 0;
  }
  else if (n == 1) {
    ans = 1;
  }
  else {
    smax = 0.;
    for (i = 1; i <= n; i++) {
      xmag = fabs(a[i0+i-1][j0]);
      if (smax < xmag) {
        smax = xmag;
        ans  = i;
      }
    }
  }
  return ans;
}

// Factors a real band matrix by elimination.
void sgbfa(double **abd, int lda, int n, int ml, int mu, int *ipvt, int *info) {
  register int i,i0,j,j0,j1,ju,jz,k,kp1,l,lm,m,mm,nm1,mn;
  double t;

  m = ml+mu+1;
  *info = 0;
  // Zero initial fill-in columns
  j0 = mu+2;
  j1 = (n<m ? n : m) - 1;
  for (jz = j0; jz <= j1; jz++) {
    i0 = m+1-jz;
    for (i=0;i<ml-i0+1;i++) abd[i0-1+i][jz-1]=0.0;
  }
  jz = j1;
  ju = 0;

  // Gaussian elimination with partial pivoting
  nm1 = n-1;
  for (k = 1; k <= nm1; k++) {
    kp1 = k+1;
   // Zero next fill-in column
    jz++;
    if (jz <= n) for (i=0;i<ml;i++) abd[i][jz-1]=0.0;

    // find L = pivot index
    lm      = (ml<n-k ? ml : n-k);
    l       = isamax(lm+1, abd, m-1, k-1)+m-1;
    ipvt[k-1] = l+k-m;
    if (abd[l-1][k-1] == 0.) {
      // zero pivot implies this column already triangularized
      *info = k;
    } else {
      // Interchange if necessary
      if (l != m) {
        t        = abd[l-1][k-1];
        abd[l-1][k-1] = abd[m-1][k-1];
        abd[m-1][k-1] = t;
      }
      // Compute multipliers
      t = -1./abd[m-1][k-1];
      for (i=0;i<lm;i++) abd[m+i][k-1]*=t;
      // Row elimination with column indexing
      mn = (ju>mu+ipvt[k-1] ? ju : mu+ipvt[k-1]);
      ju = (mn<n ? mn : n);
      mm = m;
      for (j = kp1; j <= ju; j++) {
        l--;
        mm--;
        t = abd[l-1][j-1];
        if (l != mm) {
          abd[l-1][j-1]  = abd[mm-1][j-1];
          abd[mm-1][j-1] = t;
        }
        for(i=0;i<lm;i++) abd[mm+i][j-1] += t*abd[m+i][k-1];
      }
    }
  }
  ipvt[n-1] = n;
  if (abd[m-1][n-1] == 0.) *info = n;
  return;
}

// Solves the real band system - A * X = B  or  transpose(A) * X = B
void sgbsl(double **abd, int lda, int n, int ml, int mu, int *ipvt, double *b, int job) {
  register int i,k,kb,l,la,lb,lm,m,nm1;
  double t;

  m   = mu+ml+1;
  nm1 = n-1;
  if (job == 0) {
   // Solve  A*X = B;  first solve L*Y = B
    if (ml != 0) {
      for (k = 1; k <= nm1; k++) {
        lm = (ml<n-k ? ml : n-k);
        l  = ipvt[k-1];
        t  = b[l-1];
        if (l != k) {
          b[l-1] = b[k-1];
          b[k-1] = t;
        }
        for(i=0;i<lm;i++) b[k+i] += t*abd[m+i][k-1];
      }
    }
    // Now solve  U*X = Y
    for (kb = 1; kb <= n; kb++) {
      k     = n+1-kb;
      b[k-1] /= abd[m-1][k-1];
      lm    = (k<m ? k : m) - 1;
      la    = m-lm;
      lb    = k-lm;
      t     = -b[k-1];
      for(i=0;i<lm;i++) b[lb-1+i] += t*abd[la-1+i][k-1];
    }
  }
  else {
    // Solve  trans(A)*X = B; first solve trans(U)*Y = B
    for (k = 1; k <= n; k++) {
      lm   = (k<m ? k : m) - 1;
      la   = m-lm;
      lb   = k-lm;
      for (i=0,t=0.0;i<lm;i++) t += abd[la-1+i][k-1]*b[lb-1+i];
      b[k-1] = (b[k-1]-t)/abd[m-1][k-1];
    }
    // Now solve trans(L)*X = Y
    if (ml != 0) {
      for (kb = 1; kb <= nm1; kb++) {
        k     = n-kb;
        lm    = (ml<n-k ? ml : n-k);
        for (i=0,t=0.0;i<lm;i++) t += abd[m+i][k-1]*b[k+i];
        l     = ipvt[k-1];
        if (l != k) {
          t    = b[l-1];
          b[l-1] = b[k-1];
          b[k-1] = t;
        }
      }
    }
  }
  return;
}

// Factors a real band matrix by Gaussian elimination and estimates the condition of the matrix.
void sgbco(double **abd, int lda, int n, int ml, int mu, int *ipvt, double *rcond, double *z) {
  int info;
  register int i,is,j,ju,k,kb,kp1,l,la,lm,lz,m,mm,mn;
  double anorm,ek,s,sm,t,wk,wkm,ynorm;

  // Compute 1-norm of A
  anorm = 0.;
  l  = ml+1;
  is = l+mu;
  for (j = 1; j <= n; j++) {
    for (i=0,sm=0.0;i<l;i++) sm+=fabs(abd[is-1+i][j-1]);
    anorm = (anorm>sm ? anorm:sm);
    if (is > ml+1) is--;
    if (j <= mu) l++;
    if (j >= n-ml) l--;
  }

  // Factor
  sgbfa(abd,lda,n,ml,mu,ipvt,&info);
  ek = 1.0;
  memset(z,0,n*sizeof(double));

  m  = ml+mu+1;
  ju = 0;
  for (k = 1; k <= n; k++) {
    if (z[k-1] != 0.) ek = (-z[k-1]>=0.0 ? fabs(ek) : -fabs(ek));
    if (fabs(ek-z[k-1]) > fabs(abd[m-1][k-1])) {
      s = fabs(abd[m-1][k-1])/fabs(ek-z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
      ek *= s;
    }
    wk  =  ek-z[k-1];
    wkm = -ek-z[k-1];
    s   = fabs(wk);
    sm  = fabs(wkm);
    if (abd[m-1][k-1] != 0.) {
      wk  /= abd[m-1][k-1];
      wkm /= abd[m-1][k-1];
    }
    else {
      wk  = 1.;
      wkm = 1.;
    }
    kp1 = k+1;
    mn  = (ju>mu+ipvt[k-1] ? ju : mu+ipvt[k-1]);
    ju  = (mn<n ? mn : n);
    mm  = m;
    if (kp1 <= ju) {
      for (j = kp1; j <= ju; j++) {
        mm--;
        sm     += fabs(z[j-1]+wkm*abd[mm-1][j-1]);
        z[j-1] += wk*abd[mm-1][j-1];
        s      += fabs(z[j-1]);
      }
      if (s < sm) {
        t  = wkm-wk;
        wk = wkm;
        mm = m;
        for (j = kp1; j <= ju; j++) {
          mm--;
          z[j-1] += t*abd[mm-1][j-1];
        }
      }
    }
    z[k-1] = wk;
  }

  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;

  // Solve trans(L)*Y = W
  for (kb = 1; kb <= n; kb++) {
    k  = n+1-kb;
    lm = (ml<n-k ? ml : n-k);
    if (k < n) for (i=0;i<lm;i++) z[k-1] += abd[m+i][k-1]*z[k+i];
    if (fabs(z[k-1]) > 1.) {
      s = 1./fabs(z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
    }
    l    = ipvt[k-1];
    t    = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;
  }

  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;

  ynorm = 1.0;
  // Solve L*V = Y
  for (k = 1; k <= n; k++) {
    l    = ipvt[k-1];
    t    = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;
    lm   = (ml<n-k ? ml : n-k);
    if (k < n) for(i=0;i<lm;i++) z[k+i] += t*abd[m+i][k-1];
    if (fabs(z[k-1]) > 1.) {
      s = 1.0/fabs(z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
      ynorm *= s;
    }
  }

  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;

  ynorm *= s;
  // Solve  U*Z = W
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (fabs(z[k-1]) > fabs(abd[m-1][k-1])) {
      s = fabs(abd[m-1][k-1])/fabs(z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
      ynorm *= s;
    }
    if (abd[m-1][k-1] != 0.) z[k-1] /= abd[m-1][k-1]; else z[k-1] = 1.;
    lm = (k<m ? k : m) - 1;
    la = m-lm;
    lz = k-lm;
    t  = -z[k-1];
    for(i=0;i<lm;i++) z[lz-1+i] += t*abd[la-1+i][k-1];
  }

  // Make znorm = 1.
  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;

  ynorm *= s;
  if (anorm != 0.) *rcond = ynorm/anorm; else *rcond = 0.0;
  return;
}

// Factors a real matrix by Gaussian elimination
void sgefa(double **a, int lda, int n, int *ipvt, int *info) {
  register int i,j,k,kp1,l,nm1;
  double t;

  // Gaussian elimination with partial pivoting
  *info = 0;
  nm1   = n-1;
  for (k = 1; k <= nm1; k++) {
    kp1 = k+1;
    // Find L = pivot index
    l = isamax(n-k+1, a, k-1, k-1)+k-1;
    ipvt[k-1] = l;
    if (a[l-1][k-1] == 0.) {
      // Zero pivot implies this column already triangularized
      *info = k;
    }
    else {
      // Interchange if necessary
      if (l != k) {
        t           = a[l-1][k-1];
        a[l-1][k-1] = a[k-1][k-1];
        a[k-1][k-1] = t;
      }
      // Compute multipliers
      t = -1./a[k-1][k-1];
      for (i=0;i<n-k;i++) a[k+i][k-1] *= t;

      // Row elimination with column indexing
      for (j = kp1; j <= n; j++) {
        t = a[l-1][j-1];
        if (l != k) {
          a[l-1][j-1] = a[k-1][j-1];
          a[k-1][j-1] = t;
        }
        for(i=0;i<n-k;i++) a[k+i][j-1] += t*a[k+i][k-1];
      }
    }
  }
  ipvt[n-1] = n;
  if (a[n-1][n-1] == 0.) *info = n;
  return;
}

//  Solves the real system, A*X = B  or  transpose(A)*X = B
void sgesl(double **a, int lda, int n, int *ipvt, double *b, int job) {
  register int i,k,kb,l,nm1;
  double t;

  nm1 = n-1;
  if (job == 0) {
    // Solve  A*X = B; first solve L*Y = B
    for (k = 1; k <= nm1; k++) {
      l = ipvt[k-1];
      t = b[l-1];
      if (l != k) {
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
      for(i=0;i<n-k;i++) b[k+i] += t*a[k+i][k-1];
    }
    // Now solve  U*X = Y
    for (kb = 1; kb <= n; kb++) {
      k     = n+1-kb;
      b[k-1] /= a[k-1][k-1];
      t     = -b[k-1];
      for(i=0;i<k-1;i++) b[i] += t*a[i][k-1];
    }
  }
  else {
    // Solve trans(A)*X = B; first solve trans(U)*Y = B
    for (k = 1; k <= n; k++) {
      for (i=0,t=0.0;i<k-1;i++) t += a[i][k-1]*b[i];
      b[k-1] = (b[k-1]-t)/a[k-1][k-1];
    }
    // Now solve  trans(l)*x = y
    for (kb = 1; kb <= nm1; kb++) {
      k = n-kb;
      for (i=0;i<n-k;i++) b[k-1] += a[k+i][k-1]*b[k+i];
      l     = ipvt[k-1];
      if (l != k) {
        t    = b[l-1];
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
    }
  }
  return;
}

// Factors a real matrix by Gaussian elimination and estimates the condition of the matrix.
void sgeco(double **a, int lda, int n, int *ipvt, double *rcond, double *z) {
  int info;
  register int i,j,k,kb,kp1,l;
  double anorm,ek,s,sm,t,wk,wkm,ynorm,mn;

  // Compute 1-norm of A
  anorm = 0.;
  for (j = 1; j <= n; j++) {
    for (i=0,mn=0.0;i<n;i++) mn+=fabs(a[i][j-1]);
    anorm = (anorm>mn ? anorm : mn);
  }

  // Factor
  sgefa(a,lda,n,ipvt,&info);

  ek = 1.0;
  memset(z,0,n*sizeof(double));

  for (k = 1; k <= n; k++) {
    if (z[k-1] != 0.) ek = (-z[k-1]>=0.0 ? fabs(ek) : -fabs(ek));
    if (fabs(ek-z[k-1]) > fabs(a[k-1][k-1])) {
      s = fabs(a[k-1][k-1])/fabs(ek-z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
      ek *= s;
    }
    wk  =  ek-z[k-1];
    wkm = -ek-z[k-1];
    s   = fabs(wk);
    sm  = fabs(wkm);
    if (a[k-1][k-1] != 0.) {
      wk  /= a[k-1][k-1];
      wkm /= a[k-1][k-1];
    }
    else {
      wk  = 1.;
      wkm = 1.;
    }
    kp1 = k+1;
    if (kp1 <= n) {
      for (j = kp1; j <= n; j++) {
        sm   += fabs(z[j-1]+wkm*a[k-1][j-1]);
        z[j-1] += wk*a[k-1][j-1];
        s    += fabs(z[j-1]);
      }
      if (s < sm) {
        t  = wkm-wk;
        wk = wkm;
        for (j = kp1; j <= n; j++) {
          z[j-1] += t*a[k-1][j-1];
        }
      }
    }
    z[k-1] = wk;
  }

  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;

  // Solve trans(L)*Y = W
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (k < n) for (i=0;i<n-k;i++) z[k-1] += a[k+i][k-1]*z[k+i];
    if (fabs(z[k-1]) > 1.) {
      s = 1./fabs(z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
    }
    l    = ipvt[k-1];
    t    = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;
  }
  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;
  // Solve L*V = Y
  ynorm = 1.0;
  for (k = 1; k <= n; k++) {
    l    = ipvt[k-1];
    t    = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;
    if (k < n) for(i=0;i<n-k;i++) z[k+i] += t*a[k+i][k-1];
    if (fabs(z[k-1]) > 1.) {
      s = 1./fabs(z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
      ynorm *= s;
    }
  }
  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;

  // Solve U*Z = V
  ynorm *= s;
  for (kb = 1; kb <= n; kb++) {
    k = n+1-kb;
    if (fabs(z[k-1]) > fabs(a[k-1][k-1])) {
      s = fabs(a[k-1][k-1])/fabs(z[k-1]);
      for (i=0;i<n;i++) z[i]*=s;
      ynorm *= s;
    }
    if (a[k-1][k-1] != 0.) z[k-1] /= a[k-1][k-1]; else z[k-1] = 1.;
    t = -z[k-1];
    for(i=0;i<k-1;i++) z[i] += t*a[i][k-1];
  }
  // Make znorm = 1.0
  for (i=0,s=0.0;i<n;i++) s+=fabs(z[i]);
  for (i=0;i<n;i++) z[i]/=s;
  ynorm *= s;
  if (anorm != 0.) *rcond = ynorm/anorm; else *rcond = 0.;

  return;
}

// ----------------------------------------------------------------
// Module to compute the Planck function integrated between two wavelengths
// Returns in units of (W/m2) if dl>0 or (W/m2/um) otherwise
// ----------------------------------------------------------------
double fplanck(double lam, double dl, double T) {

  const double A1=(1./3.), A2=(-1./8.), A3=(1./60.), A4=(-1./5040.);
  const double A5=(1./272160.), A6=(-1./13305600.), SIGMA=(5.67032E-8), VCUT=1.5;
  const double vcp[7] = {10.25,5.7,3.9,2.9,2.3,1.9,0.0};
  const double vmax=log(DBL_MAX),sigdpi=SIGMA/M_PI,conc=15./pow(M_PI,4.);
  double del,ex,exm,hh,mv,oldval,val,val0,vsq,d[2],p[2],v[2],ans,dvv;
  int i,k,m,mmax,n,smallv;

  if (T<1e-4 || lam<1e-4) return 0;
  if (dl==0.0) return (2e24*HP*CS*CS/pow(lam,5.0)) / (exp(1e6*HP*CS/(KB*T*lam)) - 1.0); // Return (W/m2/um)
  v[0] = C2*(1e4/(lam+dl))/T;
  v[1] = C2*(1e4/(lam))/T;

  // Wavenumbers are very close.  Get integral by iterating Simpson rule to convergence.
  if (v[0] > DBL_EPSILON && v[1] < vmax && dl/lam < 1.e-2) {
    hh     = v[1]-v[0];
    oldval = 0.;
    val0   = pow(v[0],3.0)/(exp(v[0])-1.0) + pow(v[1],3.0)/(exp(v[1])-1.0);
    for (n = 1; n <= 10; n++) {
      del = hh/(2*n);
      val = val0;
      for (k = 1; k <= 2*n-1; k++) {
        dvv = v[0]+(double)k*del;
        val += (double)(2*(1+k%2))*pow(dvv,3.0)/(exp(dvv)-1.0);
      }
      val *= del*A1;
      if (fabs((val-oldval)/val) <= 1.e-6) break;
      oldval = val;
    }
    return sigdpi*pow(T,4.0)*conc*val;
  }

  // General case
  smallv = 0;
  for (i = 1; i <= 2; i++) {
    if(v[i-1] < VCUT) { // Use power series
      smallv++;
      vsq    = v[i-1]*v[i-1];
      p[i-1] = conc*vsq*v[i-1]*(A1+v[i-1]*(A2+v[i-1]*(A3+vsq*(A4+vsq*(A5+vsq*A6)))));
    } else {            // Use exponential series
      mmax=1; while (v[i-1] < vcp[mmax-1]) mmax++;
      ex     = exp(-v[i-1]);
      exm    = 1.;
      d[i-1] = 0.;
      for (m = 1; m <= mmax; m++) {
        mv      = (double)m*v[i-1];
        exm     = ex*exm;
        d[i-1] += exm*(6.+mv*(6.+mv*(3.+mv)))/(1.0*m*m*m*m);
      }
      d[i-1] *= conc;
    }
  }

  // Handle ill-conditioning
  if (smallv == 2) {
    ans = p[1]-p[0];     // wnumlo and wnumhi both small
  } else if (smallv == 1) {
    ans = 1.0-p[0]-d[1]; // wnumlo small, wnumhi large
  } else {
    ans = d[0]-d[1];     // wnumlo and wnumhi both large
  }
  ans *= sigdpi*pow(T,4.0);
  return ans;
}
