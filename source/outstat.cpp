/******************************************************************************/
/*  outstat.cpp                                                               */
/*        print statistical analysis results                                  */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "cens.h"
#include "polysq.h"

static void OUTLevelDensity(ENSDF *, StatProperty *);
static void OUTSpinDistribution(ENSDF *, StatProperty *);


/**********************************************************/
/*      Print Statistical Analysis Data                   */
/**********************************************************/
void OUTStatAnalysis(ENSDF *lib, StatProperty *stp)
{
//if(stp->sigma2 == 0.0) return;
  cout << setw(5) << lib->getZ();
  cout << setw(5) << lib->getA();
  cout << " Nc:";
  cout << setw(5) << stp->ncomp;
  cout << setprecision(4) << setw(12) << stp->ecomp;
  cout << " Nm:";
  cout << setw(5) << stp->nmax;
  cout << setprecision(4) << setw(12) << stp->emax;
  cout << " Temp:" << setprecision(5) << setw(12) << stp->temperature;
  cout << " E0: " << setprecision(5) << setw(12) << stp->eshift;
  cout << " Sig2:" << setprecision(5) << setw(12) << stp->sigma2;
  cout << " Jmax2:" << setw(5) << stp->jmax2;
  cout << endl;
}


/**********************************************************/
/*      Print Level Density Data                          */
/**********************************************************/
void OUTStatDensity(ENSDF *lib, StatProperty *stp)
{
  cout << "# ";
  cout << setw(5) << lib->getZ();
  cout << setw(5) << lib->getA() << endl;

  OUTLevelDensity(lib, stp);
  OUTSpinDistribution(lib, stp);
}


/**********************************************************/
/*      Print Level Density                               */
/**********************************************************/
void OUTLevelDensity(ENSDF *lib, StatProperty *stp)
{
  /* constant temperature model */
  cout << "# Cumulative Number of Levels" << endl;
  for(int i=0 ; i<lib->getNlevel()-1 ; i++){
    cout << setprecision(5) << setw(12) << lib->getEnergy(i);
    cout << setprecision(5) << setw(12) << i+1.0 << endl;
  }
  cout << endl;
  cout << endl;

  cout << "# Constant Temperature Model" << endl;
  double de = 0.1;
  for(int i=1 ; ; i++){
    double ex = i * de;
    if(ex > lib->getEnergy(lib->getNlevel()-1)) break;
    double nl = exp(-stp->eshift/stp->temperature) * (exp(ex/stp->temperature) - 1.0);
    cout << setprecision(5) << setw(12) << ex;
    cout << setprecision(5) << setw(12) << nl << endl;
  }
  cout << endl;
  cout << endl;

  cout << "# Highest Complete Level" << endl;
  cout << setprecision(5) << setw(12) << stp->ecomp;
  cout << setw(12) << stp->ncomp + 1 << endl;
  cout << endl;
  cout << endl;

  cout << "# Highest No Missing Level" << endl;
  cout << setprecision(5) << setw(12) << stp->emax;
  cout << setw(12) << stp->nmax + 1 << endl;
  cout << endl;
  cout << endl;
}


/**********************************************************/
/*      Print Spin Distribution                           */
/**********************************************************/
void OUTSpinDistribution(ENSDF *lib, StatProperty *stp)
{
  const int jmax = 21;
  double *sx = new double [jmax];

  for(int j=0 ; j<jmax ; j++) sx[j] = 0.0;

  /* histogram of spins given in ENSDF */
  for(int i=0 ; i<stp->nmax ; i++){

    /* when multiple candidates are given, average them */
    int n = 0;
    for(int j=0 ; j<lib->nspin[i] ; j++){
      if(lib->spin[i][j].j < (char)0) continue;
      n ++;
    }

    for(int j=0 ; j<lib->nspin[i] ; j++){
      if(lib->spin[i][j].j < (char)0) continue;

      int k = (int)lib->spin[i][j].j / 2.0;
      sx[k] += 1.0 / n;
    }
  }

  /* normalize distribution */
  double s = 0.0;
  for(int j=0 ; j<jmax ; j++) s += sx[j];
  if(s > 0){
    for(int j=0 ; j<jmax ; j++) sx[j] /= s;
  }

  double j0 = 0.0;
  if(lib->getA()%2 != 0) j0 = 0.5;

  cout << "# Spin Disribution from Levels" << endl;
  for(int j=0 ; j<jmax ; j++){
    cout << setprecision(5) << setw(12) << j + j0;
    cout << setprecision(5) << setw(12) << sx[j] << endl;
  }
  cout << endl;
  cout << endl;

  /* spin distribution by spin cut-off formula, not exactly normalized */
  cout << "# Spin Disribution" << endl;
  double dj = 0.1;
  for(int j=0 ; ; j++){
    double x = dj * j;
    double y = (x + 0.5)/stp->sigma2 * exp( -(x + 0.5) * (x + 0.5) /(2 * stp->sigma2) );
    cout << setprecision(5) << setw(12) << x;
    cout << setprecision(5) << setw(12) << y << endl;
    if(x >= (double) jmax) break;
  }

  delete [] sx;
}
