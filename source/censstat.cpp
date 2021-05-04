/******************************************************************************/
/*  censstat.cpp                                                              */
/*        analysis of discrete levels with statistical model                  */
/******************************************************************************/

#include <iostream>
#include <cmath>

using namespace std;

#include "cens.h"
#include "polysq.h"

static int    LEVELCheckCompleteness(ENSDF *);
static double LEVELSpinCutoff(const int, ENSDF *);
static int    LEVELHighestSpin(const int, ENSDF *);
static int    LEVELMissing(const int, ENSDF *, double *, double *, double *);
static int    LEVELMissingCase1(const int, const int, double *, double *, double *);
static int    LEVELMissingCase2(const int, const int, double *, double *, double *);


/***********************************************************/
/*      Analysis of Nuclear Structure                       */
/***********************************************************/
void CENSStat(ENSDF *lib, StatProperty *stp)
{
  double *x = new double [lib->getNlevel()];
  double *y = new double [lib->getNlevel()];
  double a[2];


  /* check the highest level */
  /* up to which given information in ENSDF is complete */
  stp->ncomp = LEVELCheckCompleteness(lib);
  stp->ecomp = lib->getEnergy(stp->ncomp);

  if(lib->getNlevel() <= 3){
    stp->nmax = stp->ncomp;
    stp->emax = stp->ecomp;
    return;
  }

  /* up to which no-missing level can be assumed */
  stp->nmax = LEVELMissing(stp->ncomp,lib,x,y,a);
  stp->emax = lib->getEnergy(stp->nmax);

  if(stp->nmax == 0){
    stp->nmax = stp->ncomp;
    stp->emax = stp->ecomp;
    return;
  }

  /* spin cut-off parameter */
  stp->sigma2 = LEVELSpinCutoff(stp->ncomp,lib);

  /* highest spin */
  stp->jmax2 = LEVELHighestSpin(stp->ncomp,lib);

  /* constant temperature model */
  if(stp->nmax > 3){
    LSQPolynomial(stp->nmax-1,2,x,y,a);

    stp->temperature = 1.0 / a[1];
    stp->eshift = -stp->temperature * log(stp->nmax / (exp(stp->emax/stp->temperature) - 1.0));
  }
  else{
    stp->temperature = 0.0;
    stp->eshift = 0.0;
  }

  delete [] x;
  delete [] y;

}


/***********************************************************/
/*      Find Highest Level to Which Data Are Complete      */
/***********************************************************/
int LEVELCheckCompleteness(ENSDF *lib)
{
  int m = 0;

  for(int i=0 ; i<lib->getNlevel() ; i++){

    /* if more than one spins are given */
    if(lib->nspin[i] > 1){ m = i-1; break; }

    /* if spin is unknown */
    if((int)lib->spin[i][0].j < 0){ m = i-1; break; }

    /* if parity is unkown */
    if((int)lib->spin[i][0].p == 0){ m = i-1; break; }

    /* if no gamma-decay is given */
    if((i != 0) && (lib->gamma[i].getNgamma() == 0)){ m = i-1; break; }
  }

  return m;
}


/***********************************************************/
/* Calculate spin cutoff parameter from discrete levels    */
/***********************************************************/
double LEVELSpinCutoff(const int m, ENSDF *lib)
{
  double sig2 = 0.0;

  if(lib->getNlevel() <= 2 || m <= 2) return sig2;

  int mx = 0;
  for(int i=0 ; i<m ; i++){
    double s = 0.0;
    for(int j=0 ; j<lib->nspin[i] ; j++){
      if(lib->spin[i][j].j < (char)0) continue;

      s += (double)lib->spin[i][j].j / 2.0;
      mx ++;
    }

    s /= lib->nspin[i];
    sig2 += (s + 1.0) * s;
  }
  sig2 = sig2 / (2.0 * mx);

  return sig2;
}


/***********************************************************/
/* Highest spin of discrete levels                         */
/***********************************************************/
int LEVELHighestSpin(const int m, ENSDF *lib)
{
  int j2 = 0;

  if(lib->getA()%2 != 0) j2 = 1;

  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<lib->nspin[i] ; j++){
      if((int)lib->spin[i][j].j > j2) j2 = (int)lib->spin[i][j].j;
    }
  }

  return j2;
}


/***********************************************************/
/* Highest level upto which no missing level assmued       */
/***********************************************************/
int LEVELMissing(const int nc, ENSDF *lib, double *x, double *y, double *a)
{
  int nx = 20;
  if(nc > nx) nx = nc;

  /* exclude ground state */
  for(int i=0 ; i<lib->getNlevel()-1 ; i++){
    x[i] = lib->getEnergy(i+1);
    y[i] = log((double)i + 2.0);
  }

  int ncut = 0;
  /* when a large number of levels are given */
  if(lib->getNlevel() > nx){
    ncut = LEVELMissingCase1(nx,lib->getNlevel(),x,y,a);
    if(ncut == nx) ncut = LEVELMissingCase2(nc,nx,x,y,a);
  }
  /* small number case */
  else{
    ncut = LEVELMissingCase2(nc,lib->getNlevel(),x,y,a);
  }

  return ncut;
}


/***********************************************************/
/*  Highest level algorithm, Case 1                        */
/***********************************************************/
int LEVELMissingCase1(const int nx, const int nl, double *x, double *y, double *a)
{
  int ncut = 0;

  /* fit a linear function to the levels below I (>= Nx),
     and extrapolate into the higher region. */

  for(int i=nx ; i<nl ; i++){
    // fit the linear function a0 + a1 X[i] for the data [0:i]
    LSQPolynomial(i,2,x,y,a);

    /* If all the points in the higher region are below the fitted line,
       there might be missing levels. */
    int m = 0;
    for(int j=i+1 ; j<nl-1 ; j++){
      double z = a[0] + a[1]*x[j];
      if(y[j] > z) m ++;
    }

    ncut = i;
    if(m == 0) break;
  }

  return ncut;
}



/***********************************************************/
/*  Highest level algorithm, Case 2                        */
/***********************************************************/
int LEVELMissingCase2(const int nc, const int nx, double *x, double *y, double *a)
{
  int ncut = 0;

  LSQPolynomial(nx-1,2,x,y,a);

  double cmin = 1e+99;
  for(int i=nc ; i<nx-1 ; i++){
    int m = i - nc + 1;

    double c = 0.0;
    for(int j=0 ; j<=i ; j++){
      double z = a[0] + a[1]*x[j];
      c += (z - y[j]) * (z - y[j]);
    }
    c /= (double)m;

    if(c < cmin){
      cmin = c;
      ncut = i;
    }
  }

  return ncut;
}
