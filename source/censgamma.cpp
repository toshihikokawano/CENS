/******************************************************************************/
/*  censgamma.cpp                                                             */
/*        determine final states of gamma-decay and branching ratios          */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <cmath>

using namespace std;

#include "cens.h"
#include "terminate.h"

static void GAMFinalState(ENSDF *);
static void GAMNormalizeBranch(ENSDF *);

#undef DEBUG
#ifdef DEBUG
static void print(ENSDF *);
#endif


/***********************************************************/
/*      Clean Decay Matrix                                 */
/***********************************************************/
void CENSGamma(ENSDF *lib)
{
  GAMFinalState(lib);
  GAMNormalizeBranch(lib);

#ifdef DEBUG
  print(lib);
#endif
}


/***********************************************************/
/*      Find Final State and Fix Gamma-Ray Energy          */
/***********************************************************/
void GAMFinalState(ENSDF *lib)
{
  const double eps = 1e-3;

  /* initial state */
  for(int i0=1 ; i0<lib->getNlevel() ; i0++){
    double e0 = lib->getEnergy(i0);

    /* look at each gamma line */
    for(int j=0 ; j<lib->gamma[i0].getNgamma() ; j++){
      double eg = lib->gamma[i0].getEnergy(j);

      /* final state, determine by the least difference
         between level and gamma energies */
      int    k = 0;
      double z = abs((e0 - lib->getEnergy(k)) / eg - 1.0);
      for(int i1=1 ; i1<=i0-1 ; i1++){
        double e1 = lib->getEnergy(i1);
        double r  = abs((e0 - e1) / eg - 1.0);
        if(r < z){
          z = r;
          k = i1;
        }
      }
      if(z > eps){
        message << "energy mismatch detected: ";
        message << "initial " << e0 << ": final " << lib->getEnergy(k);
        message << ": dE = " << e0 - lib->getEnergy(k) << ": Egamma " << eg;
        message << ": delta " << abs((e0 - lib->getEnergy(k))/eg - 1.0);
        Notice("GAMFinalState");
      }

      /* set final state */
      lib->gamma[i0].fstate[j] = k;

      /* fix gamma-ray energy */
      lib->gamma[i0].energy[j] = e0 - lib->getEnergy(k);
    }
  }
}


/***********************************************************/
/*      Renormalize Branching Ratios (no fix)              */
/***********************************************************/
void GAMNormalizeBranch(ENSDF *lib)
{
  /* initial state */
  for(int i0=1 ; i0<lib->getNlevel() ; i0++){

    /* when only one gamma-ray */
    if(lib->gamma[i0].getNgamma() == 1) lib->gamma[i0].branch[0] = 1.0;

    else{
      /* for each gamma decay */
      double s = 0.0;
      for(int j=0 ; j<lib->gamma[i0].getNgamma() ; j++){
        s += lib->gamma[i0].getBranch(j);
      }
      /* if no branching ratios given, leave them */
      if(s > 0.0){
        s = 1.0 / s;
        for(int j=0 ; j<lib->gamma[i0].getNgamma() ; j++) lib->gamma[i0].branch[j] *= s;
      }
    }
  }
}


#ifdef DEBUG
/***********************************************************/
/*      Debugging Print                                    */
/***********************************************************/
void print(ENSDF *lib)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(5);

  for(int i = 0 ; i < lib->getNlevel() ; i++){

    cout << setw(5) << i;
    cout << setw(13) << lib->getEnergy(i) << endl;

    for(int j=0 ; j<lib->gamma[i].getNgamma() ; j++){
      cout << "        ";
      cout << setw(5)  << lib->gamma[i].getFstate(j);
      cout << setw(13) << lib->gamma[i].getEnergy(j);
      cout << setw(13) << lib->gamma[i].getBranch(j) << endl;
    }
    cout << endl;
  }
}
#endif
