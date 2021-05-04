/******************************************************************************/
/*  riplread.cpp                                                              */
/*        read internal conversion from RIPL data file                        */
/******************************************************************************/

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

#include "cens.h"
#include "terminate.h"


/***********************************************************/
/*      Read RIPL                                          */
/***********************************************************/
int RIPLRead(string ripldir, ENSDF *lib)
{
  const double  eps = 1e-5;
  ostringstream os;
  ifstream      fp;
  string        file, str;

  os << setw(3) << setfill('0') << lib->getZ();
  file = os.str();
  file = "z" + file + ".dat";

  if(ripldir.length() > 0){
    /* remove if dir name includes a slash at the end */
    if(ripldir[ripldir.length() - 1] == '/') ripldir[ripldir.length() - 1] = '\0';
    file = ripldir + '/' + file;
  }

  message << "RIPL file name " << file;
  Notice("RIPLRead");

  fp.open(&file[0]);
  if(!fp){
    message << "RIPL file " << file << " cannot open";
    TerminateCode("RIPLRead");
  }

  double *ex, **ic;
  int *nl, **fs;
  ex = new double [lib->getNsize()];
  ic = new double * [lib->getNsize()];
  nl = new int [lib->getNsize()];
  fs = new int * [lib->getNsize()];
  for(int i=0 ; i<lib->getNsize() ; i++){
    ic[i] = new double [lib->gamma[0].getNsize()];
    fs[i] = new int [lib->gamma[0].getNsize()];
  }

  bool found = false;
  string d;
  int nlev = 0;
  while(getline(fp,str)){

    /*** search for Z and A entry in the file */
    d = str.substr( 5, 5);  int a   = atoi(&d[0]);
    d = str.substr(10, 5);  int z   = atoi(&d[0]);
    d = str.substr(15, 5);  nlev    = atoi(&d[0]);

    if((z == lib->getZ()) && a == lib->getA()) found = true;

    /* for all discrete levels */
    for(int i1=0 ; i1<nlev ; i1++){
      getline(fp,str);
      d = str.substr( 4,10);  double e = atof(&d[0]);
      d = str.substr(34, 3);  int    n = atoi(&d[0]);

      if(found){
        ex[i1] = e;
        nl[i1] = n;
      }

      /* for gamma-rays */
      for(int j1=0 ; j1<n ; j1++){
        getline(fp,str);
        d = str.substr(39, 4);  int    m = atoi(&d[0]);
        d = str.substr(77,10);  double c = atof(&d[0]);

        if(found){
          fs[i1][j1] = m - 1;
          ic[i1][j1] = c;
        }
      }
    }
    if(found) break;
  }
  fp.close();

  /* initial and final energies for gamma transition in ENSDF */
  for(int i0 = 1 ; i0<lib->getNlevel() ; i0++){
    double e00 = lib->getEnergy(i0);
    for(int j0=0 ; j0<lib->gamma[i0].getNgamma() ; j0++){
      double e01 = lib->getEnergy( lib->gamma[i0].getFstate(j0) );

      /* see if IC is zero in ENSDF */
      if(lib->gamma[i0].getCvcoef(j0) == 0.0){

        bool found = false;
        /* look for the same gamma transition in RIPL */
        for(int i1 = 1 ; i1 < nlev ; i1++){
          double e10 = ex[i1];
          for(int j1 = 0 ; j1 < nl[i1] ; j1++){
            double e11 = ex[ fs[i1][j1] ];

            /* if two energies are close enough, this is it */
            double d0 = abs(e00 - e10);
            double d1 = abs(e01 - e11);
            if( (d0 <= eps) && (d1 <= eps) ){
              if(ic[i1][j1] > 0.0) lib->gamma[i0].cvcoef[j0] = ic[i1][j1];
              found = true;
              break;
            }
          }
          if(found) break;
        }
      }
    }
  }

  for(int i=0 ; i<lib->getNsize() ; i++){
    delete [] ic[i];
    delete [] fs[i];
  }

  delete [] ex;
  delete [] ic;
  delete [] fs;
  delete [] nl;

  return 0;
}
