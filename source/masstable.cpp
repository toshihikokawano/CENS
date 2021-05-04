/******************************************************************************/
/*  masstable.cpp                                                             */
/*        find nuclear mass from Mass Table                                   */
/******************************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>

using namespace std;

#include "masstable.h"

#include "masstable_audi2012_frdm2012.h"  // AW12 + FRDM2012

/*************************************************/
/*  Mass Excess, Non-Stop Mode                   */
/*************************************************/
double mass_excess(int z, int a, bool *found)
{
  double    mx  = 0.0;
  unsigned int za = z*1000+a;

  *found = false;
  for(int i=0 ; i<nMassTable ; i++){
    if(MassTable[i].za == za){
      *found = true;
      mx = MassTable[i].mass;
      break;
    }
  }

  /*** huge value returned, so that separation energy will be infinity */
  if(! *found) mx = 1e+10;

  return(mx);
}

