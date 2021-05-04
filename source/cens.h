const int MaxDiscreteLevels = 10000;
const int MaxGammaLines     = 100;


#ifndef __ENSDF_H__
#define __ENSDF_H__
#include "ensdf.h"
#endif

//------------------------------------------------------------------------------
//     Class

/**********************************************************/
/*   Statistical Properties of Nuclear Structure          */
/**********************************************************/
class StatProperty{
 private:
 public:
  int    ncomp;       // number of levels to which complete information given
  int    nmax;        // no missing level assumed
  int    jmax2;       // highest spin (x2)
  double ecomp;       // energy at ncomp
  double emax;        // energy at nmiss;
  double sigma2;      // sigma2 from discrete levels
  double temperature; // temperature in the constant temperature model
  double eshift;      // energy shift, E0

  StatProperty(){
    ncomp = 0;
    nmax = 0;
    jmax2 = 0;
    ecomp = 0.0;
    emax = 0.0;
    sigma2 = 0.0;
    temperature = 0.0;
    eshift = 0.0;
  }
};


//------------------------------------------------------------------------------
//     Prototype Definitions

// cens.cpp
void NoticeMessage  (std::string);
void NoticeMessage  (std::string, const int);
void NoticeMessage  (std::string, const double);
void NoticeMessage  (std::string, std::string);
void WarningMessage (std::string);
void WarningMessage (std::string, const int);
void WarningMessage (std::string, const double);
void WarningMessage (std::string, std::string);

// censgamma.cpp
void CENSGamma (ENSDF *);

// censstat.cpp
void CENSStat (ENSDF *, StatProperty *);

// outxml.cpp
void OUTFxml (ENSDF *);

// outripl.cpp
void OUTFripl (const int, const int, ENSDF *);

// outstat.cpp
void OUTStatAnalysis (ENSDF *, StatProperty *);
void OUTStatDensity (ENSDF *, StatProperty *);
