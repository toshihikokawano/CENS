/******************************************************************************/
/*  outxml.cpp                                                                */
/*        print ENSDF data in XML                                             */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

#include "cens.h"
#include "xmltag.h"


/**********************************************************/
/*      Print Data in XML Format                          */
/**********************************************************/
void OUTFxml(ENSDF *lib)
{
  ostringstream attr;
  attr.str("");
  attr << "date=" << "\"" << setw(6) <<lib->date << "\"";
  attr << " Z="   << "\"" << lib->getZ() << "\"";
  attr << " A="   << "\"" << lib->getA() << "\"";
  attr << " energy_unit=" << "\"" << lib->getUnit() << "\"";
  XMLTagOpen("ENSDF",attr.str());

  for(int i = 0 ; i < lib->getNlevel() ; i++){

    attr.str(""); attr << "number=" << i;
    XMLTagOpen("LEVEL",attr.str());
    XMLTagVal("LevelEnergy",lib->getEnergy(i));

    if(lib->getThalf(i) < 0.0)
      XMLTagVal("LevelHalfLife", "stable");
    else
      XMLTagVal("LevelHalfLife", lib->getThalf(i));

    if(lib->nspin[i] == 1){
      double s = (double)lib->spin[i][0].j/2.0;
      if(s < 0.0)
        XMLTagVal("LevelSpin","unknown");
      else
        XMLTagVal("LevelSpin",s);

      int p = (int)lib->spin[i][0].p;
      if(p == 1)
        XMLTagVal("LevelParity", "+");
      else if(p == -1)
        XMLTagVal("LevelParity", "-");
      else
        XMLTagVal("LevelParity", "unknown");
    }
    else{
      for(int n=0 ; n<lib->nspin[i] ; n++){
        XMLTagOpen("SPINS");
        XMLTagVal("SpinCandidate",(double)lib->spin[i][n].j/2.0);
        int p = (int)lib->spin[i][n].p;
        if(p == 1)
          XMLTagVal("ParityCandidate", "+");
        else if(p == -1)
          XMLTagVal("ParityCandidate", "-");
        else
          XMLTagVal("ParityCandidate", "unknown");
        XMLTagClose("SPINS");
      }
    }

    for(int j=0 ; j<lib->gamma[i].getNgamma() ; j++){
      attr.str(""); attr << "number=" << j;
      XMLTagOpen("GAMMA",attr.str());
      XMLTagVal("GammaEnergy",lib->gamma[i].getEnergy(j));
      XMLTagVal("GammaBranch",lib->gamma[i].getBranch(j));
      XMLTagVal("GammaConversionCoefficient",lib->gamma[i].getCvcoef(j));
      XMLTagClose("GAMMA");
    }
    XMLTagClose("LEVEL");
  }

  XMLTagClose("ENSDF");
}


