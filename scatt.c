// ---------------------------------------------------------
// Surface scattering functions
// ---------------------------------------------------------
// Program to compute single-scattering phase functions
double scattphase(double phase, char *phasemodel, double *params) {
  double fphase=1.0, rd=PI/180.0, phi1,phi2,phi3,phiLS, f1,f2,fc, b,c;
  phase = fabs(phase);
  if (phase>180.0) phase=360-phase;

  // Henyey-Greenstein HG1
  if (strcmp(phasemodel,"HG1")==0) {
    b = params[0];
    fphase = (1.0-b*b) / pow(1.0 + 2.0*b*cos(phase*rd) + b*b,3.0/2.0);

  // Henyey-Greenstein HG2 / HGH
  } else if (strcmp(phasemodel,"HG2")==0 || strcmp(phasemodel,"HGH")==0) {
    b = params[0];
    if (strcmp(phasemodel,"HGH")==0) c=3.29*exp(-17.4*b*b) - 0.908; else c = params[1];
    fphase = ((1.0+c)/2.0) * (1.0-b*b) / pow(1.0 - 2.0*b*cos(phase*rd) + b*b,3.0/2.0);
    fphase+= ((1.0-c)/2.0) * (1.0-b*b) / pow(1.0 + 2.0*b*cos(phase*rd) + b*b,3.0/2.0);

  // Legendre polynomial LP2
  } else if (strcmp(phasemodel,"LP2")==0) {
    b = params[0]; c = params[1];
    fc = cos(phase*rd);
    fphase = 1.0 + b*fc + c*(1.5*fc*fc - 0.5);

  // Lumme-Bowel (HG) phase function
  } else if (strcmp(phasemodel,"HG")==0 || strcmp(phasemodel,"HM")==0) {
    if (phase==180.0) {
      fphase=0.0;
    } else {
      f1 = params[0];
      fc = tan(phase*rd/2.0);
      phi1 = exp(-3.33*pow(fc,0.63));
      phi2 = exp(-1.87*pow(fc,1.22));
      if (f1>1.0) f1=1.0;
      else if (f1<0.0) f1=0.0;
      fphase = (1.0-f1)*phi1 + f1*phi2;
      phiLS = 1.0 - (sin(phase*rd/2.0)*tan(phase*rd/2.0)*log(1.0/tan((phase*rd/4.0)+1e-6)));
      fphase /= phiLS;
    }

  // Muinonen H,G1,G2 phase function, with spline factors and H,G12 from Penttila+2016
  } else if (strcmp(phasemodel,"HG1G2")==0 || strcmp(phasemodel,"HG12")==0) {
    double x12[6]= {7.5,30.0,60.0,90.0,120.0,150.0};
    double y1[6] = {7.5e-1, 3.3486016e-1, 1.3410560e-1, 5.1104756e-2, 2.1465687e-2, 3.6396989e-3};
    double d1[6] = {-1.9098593, -5.5463432e-1, -2.4404599e-1, -9.4980438e-2, -2.1411424e-2, -9.1328612e-2};
    double y2[6] = {9.25e-1, 6.2884169e-1, 3.1755495e-1, 1.2716367e-1, 2.2373903e-2, 1.6505689e-4};
    double d2[6] = {-5.7295780e-1, -7.6705367e-1, -4.5665789e-1, -2.8071809e-1, -1.1173257e-1, -8.6573138e-8};
    double x3[9] = {0.0,0.3,1.0,2.0,4.0,8.0,12.0,20.0,30.0};
    double y3[9] = {1.0,8.3381185e-1, 5.7735424e-1, 4.2144772e-1, 2.3174230e-1, 1.0348178e-1, 6.1733473e-2, 1.6107006e-2, 0.0};
    double d3[9] = {-1.0630097e-1, -4.1180439e1, -1.0366915e1, -7.5784615, -3.6960950, -7.8605652e-1, -4.6527012e-1, -2.0459545e-1, 0.0};
    int i12=0,i3=0; double t,a,b;

    // Perform the splines
    if (phase>x12[5]) {
      phi1 = 0.0;
      phi2 = 0.0;
    } else if (phase<7.5) {
      phi1 = 1.0 - 1.90985931710274*phase*rd;
      phi2 = 1.0 - 0.572957795130823*phase*rd;
    } else {
      while (phase>x12[i12+1] && i12<5) i12++;
      t = (phase-x12[i12])/(x12[i12+1]-x12[i12]);
      a =  d1[i12  ]*rd*(x12[i12+1]-x12[i12]) - (y1[i12+1]-y1[i12]);
      b = -d1[i12+1]*rd*(x12[i12+1]-x12[i12]) + (y1[i12+1]-y1[i12]);
      phi1 = (1.0-t)*y1[i12] + t*y1[i12+1] + t*(1.0-t)*((1.0-t)*a + b*t);
      a =  d2[i12  ]*rd*(x12[i12+1]-x12[i12]) - (y2[i12+1]-y2[i12]);
      b = -d2[i12+1]*rd*(x12[i12+1]-x12[i12]) + (y2[i12+1]-y2[i12]);
      phi2 = (1.0-t)*y2[i12] + t*y2[i12+1] + t*(1.0-t)*((1.0-t)*a + b*t);
    }
    if (phase>x3[8]) {
      phi3=0.0;
    } else {
      while (phase>x3[i3+1] && i3<8) i3++;
      t = (phase-x3[i3])/(x3[i3+1]-x3[i3]);
      a =  d3[i3  ]*rd*(x3[i3+1]-x3[i3]) - (y3[i3+1]-y3[i3]);
      b = -d3[i3+1]*rd*(x3[i3+1]-x3[i3]) + (y3[i3+1]-y3[i3]);
      phi3 = ((1.0-t)*y3[i3]) + (t*y3[i3+1]) + (t*(1.0-t)*(((1.0-t)*a) + (b*t)));
    }

    if (strcmp(phasemodel,"HG1G2")==0) {
      f1=params[0]; f2=params[1];
    } else {
      fc=params[0];
      if (fc>1.0) fc=1.0;
      else if (fc<0.0) fc=0.0;
      f1 = 0.5351335*fc;
      f2 = 0.84293649 - 0.5351335*fc;
    }
    fphase = f1*phi1 + f2*phi2 + (1.0-f1-f2)*phi3;
    phiLS = 1.0 - (sin(phase*rd/2.0)*tan(phase*rd/2.0)*log(1.0/tan((phase*rd/4.0)+1e-6)));
    if (phiLS>0) fphase /= phiLS;

  } else if (strcmp(phasemodel,"EXP")==0) {
    // Exponential
    fc = phase; fphase = fc*params[0];
    fc*= phase; fphase+= fc*params[1];
    fc*= phase; fphase+= fc*params[2];
    fphase = exp(fphase);

  } else if (strcmp(phasemodel,"ROLO")==0) {
    // Lunar/ROLO phase
    fphase = params[0]*exp(-params[1]*phase);
    fc = phase; fphase+= fc*params[2];
    fc*= phase; fphase+= fc*params[3];
    fc*= phase; fphase+= fc*params[4];
    fc*= phase; fphase+= fc*params[5];
  }

  return fphase;
}

// Scattering function - It returns BRDF (bidirectional-reflectance distribution function) [sr-1]
double scattbrdf(double mu, double mu0, double phase, double fphase, double albedo, int model, double *params, double *emiss) {
  double rscat=0.0;
  if (albedo<0) albedo=0.0;
  if (albedo<1.0) *emiss = 1.0-albedo; else *emiss=0.0;

  // Calculate scattering functions
  if (model==0 || mu<1e-4 || mu0<1e-4) {
    // Lambert model
    rscat = albedo/PI*fphase;

  } else if (model==1) {
    // Hapke model
    double rd=PI/180.0, kphi, K=1.0, BSO, BS=0.0, BCO, BC=0.0, gamma, r0, x, He, Hu, Hu0, cosphi, sinphi2, phi, tanthe, denom, theta, fg, xtheta, fc, E1i, E1e, E2i, E2e, ni, ne;
    double mu0e=mu0, mue=mu, Si=1.0, tanhg=tan(phase*rd/2.0);
    if (albedo>=1.0) albedo=1.0-1e-6;
    gamma = sqrt(1.0-albedo);
    r0 = (1.0-gamma)/(1.0+gamma);

    // Porosity term
    kphi = 1.209*pow(params[5],2.0/3.0);
    if (kphi<1 && kphi>0) K = -log(1.0 - kphi)/kphi;

    // Shadow-hiding opposition effect (SHOE)
    BSO = params[2];
    if (params[3]>0) BS = 1.0/(1.0 + (1.0/params[3])*tanhg);

    // Coherent backscatter opposition effect (CHOE)
    BCO = params[6];
    if (params[7]>0) {
      x = (1.0/params[6])*tanhg;
      BC = 1.0/(1.0 + (1.3+K)*(x + x*x));
    }

    // Shadowing function
    if (params[4]>0) {
      double e=acos(mu), i=acos(mu0), sine=sin(e), sini=sin(i), cosg=cos(phase*rd);
      theta = (1.0-r0)*params[4]*rd;
      tanthe = tan(theta);
      cosphi = (cosg - mu0*mu)/(sini*sine + 1e-12);
      if (cosphi>1.0) cosphi=1.0;
      else if (cosphi<-1.0) cosphi=-1.0;
      phi = acos(cosphi);
      sinphi2 = sin(phi/2.0); sinphi2*=sinphi2;

      fc = 1.0/tan(theta)*1.0/tan(e);
      E1e = exp(-(2.0/PI)*fc);
      E2e = exp(-(1.0/PI)*fc*fc);
      fc = 1.0/tan(theta)*1.0/tan(i);
      E1i = exp(-(2.0/PI)*fc);
      E2i = exp(-(1.0/PI)*fc*fc);

      xtheta = 1.0/sqrt(1.0 + PI*tanthe*tanthe);
      fg = exp(-2.0*tan(phi/2.0));

      ne = xtheta*(mu  + sine*tanthe*(E2e/(2.0-E1e)));
      ni = xtheta*(mu0 + sini*tanthe*(E2i/(2.0-E1i)));

      if (i<e) {
        denom = 2.0 - E1e - (phi/PI)*E1i;
        mu0e = xtheta*(mu0 + sini*tanthe*(cosphi*E2e + sinphi2*E2i)/denom);
        mue  = xtheta*(mu  + sine*tanthe*(E2e - sinphi2*E2i)/denom);
        Si = (mue/ne)*(mu0/ni)*xtheta/(1.0 - fg + fg*xtheta*mu0/ni);
      } else {
        denom = 2.0 - E1i - (phi/PI)*E1e;
        mu0e = xtheta*(mu0 + sini*tanthe*(E2i - sinphi2*E2e)/denom);
        mue  = xtheta*(mu  + sine*tanthe*(cosphi*E2i + sinphi2*E2e)/denom);
        Si = (mue/ne)*(mu0/ni)*xtheta/(1.0 - fg + fg*xtheta*mu/ne);
      }
    }

    // Ambartsumian–Chandrasekhar scattering H-functions
    x = mu/K;   He  = 1.0/(1.0 - albedo*x*(r0 + (((1.0-2.0*r0*x)/2.0)*log((1.0+x)/x))));
    x = mue/K;  Hu  = 1.0/(1.0 - albedo*x*(r0 + (((1.0-2.0*r0*x)/2.0)*log((1.0+x)/x))));
    x = mu0e/K; Hu0 = 1.0/(1.0 - albedo*x*(r0 + (((1.0-2.0*r0*x)/2.0)*log((1.0+x)/x))));

    // Compute the Hapke BRDF
    rscat = K*albedo*Si/(4.0*PI*(mue+mu0e))*(mu0e/mu0);
    rscat*= fphase*(1.0 + BSO*BS) + Hu*Hu0 - 1.0;
    rscat*= (1.0 + BCO*BC);
    *emiss = gamma*He;

  } else if (model==2) {
    // Lommel-Seeliger
    rscat = albedo*fphase/(4.0*PI*(mu+mu0));

  } else if (model==3) {
    // Cox-Munk glitter/specular reflection model
    double sigma2, sigma, w, r, p, gp, tp, pdf, cotT,cotT2,lami=0.0,lame=0.0,sh;
    double tu0=tan(acos(mu0)), tu=tan(acos(mu)), rd=PI/180.0;

    // Fresnel equations (Jackson+Alpers 2010)
    sigma2 = 0.003 + 0.00512*params[0];
    w = phase*rd/2.0; if (w==0) w=0.01;
    r = asin(sin(w)/1.34);
    p = 0.5*(pow(sin(w-r)/sin(w+r),2.0) + pow(tan(w-r)/tan(w+r),2.0));

    // Adapted from Spurr+2004 and Jackson+Alpers 2010
    gp = (mu0+mu)/(2.0*cos(w)); if (gp>1.0) gp=1.0;
    tp = tan(acos(gp)); tp *= tp;
    pdf = exp(-tp/sigma2)/(PI*sigma2);

    // Shadowing effect (Ma+2015, http://dx.doi.org/10.1364/AO.54.009863)
    sigma = sqrt(sigma2);
    cotT = 1.0/(tu0+1e-6); cotT2 = cotT*cotT;
    lami = 0.5*((1.0/sqrt(PI))*(sigma/cotT)*exp(-cotT2/sigma2) - erfc(cotT/sigma));
    cotT = 1.0/(tu+1e-6); cotT2 = cotT*cotT;
    lame = 0.5*((1.0/sqrt(PI))*(sigma/cotT)*exp(-cotT2/sigma2) - erfc(cotT/sigma));
    sh = 1.0/(1.0 + lami + lame);

    rscat = (p/4.0)*sh*pdf*(1.0+tp)*(1.0+tp)/mu;  // Cox-Munk glint signal [sr-1]
    rscat+= albedo/PI*fphase;                     // Add non-glint surface signal [sr-1]

  } else {
    // Lambert model
    rscat = albedo/PI*fphase;
  }
  return rscat;
}
