/******************************************************************************/
/*  outripl.cpp                                                               */
/*        print discrete level data in RIPL format                            */
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "cens.h"
#include "elements.h"
#include "masstable.h"
#include "physicalconstant.h"


/**********************************************************/
/*      Print Data in RIPL Format                         */
/**********************************************************/
void OUTFripl(const int nc, const int nm, ENSDF *lib)
{
  /* count total gamma-rays */
  int nog  = 0;
  for(int i = 0 ; i < lib->getNlevel() ; i++) nog += lib->gamma[i].getNgamma();

  /* calculate separation energies */
  bool f1, f2;
  double sn = mass_excess(lib->getZ(),lib->getA()-1,&f2) + ENEUTRON
            - mass_excess(lib->getZ(),lib->getA()  ,&f1);

  double sp = mass_excess(lib->getZ()-1,lib->getA(),&f2) + EPROTON
            - mass_excess(lib->getZ()  ,lib->getA(),&f1);

  /* print header line */
  cout << setw(3) << lib->getA() << left << setw(2) << element_name[lib->getZ()];
  cout << right;
  cout << setw(5) << lib->getA() << setw(5) << lib->getZ();
  cout << setw(5) << lib->getNlevel() << setw(5) << nog;
  cout << setw(5) << nm << setw(5) << nc;

  cout.setf(ios::fixed, ios::floatfield);
  cout << setprecision(6);
  cout << setw(12) << sn;
  cout << setw(12) << sp;
  cout << endl;

  cout.setf(ios::scientific);
  cout << setprecision(4);

  /* print all levels */
  for(int i = 0 ; i < lib->getNlevel() ; i++){

    cout << setw(3) << i+1;

    /* level energy */
    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(6) << setw(11) << lib->getEnergy(i);

    /* spin and pairy */
    double s = (double)lib->spin[i][0].j/2.0;
    int    p = (int)lib->spin[i][0].p;
    if(s < 0.0) s = -1.0;

    cout << setprecision(1) << setw(6) << s;
    cout << setw(3) << p;

    /* half-life */
    cout.setf(ios::scientific, ios::floatfield);
    if(lib->getThalf(i) == 0.0) cout << "           ";
    else cout << setprecision(2) << setw(11) << lib->getThalf(i);

    /* number of gamma-rays */
    cout << setw(3) << lib->gamma[i].getNgamma();
    cout << endl;

    /* print gamma-rays */
    for(int j=0 ; j<lib->gamma[i].getNgamma() ; j++){
      cout << "                                      ";
      cout << setw(5) << lib->gamma[i].getFstate(j)+1;

      cout.setf(ios::fixed, ios::floatfield);
      cout << setprecision(3) << setw(11) << lib->gamma[i].getEnergy(j);

      cout.setf(ios::scientific, ios::floatfield);
      double r = 1.0 / (1.0 + lib->gamma[i].getCvcoef(j));
      cout << setprecision(3) << setw(11) << lib->gamma[i].getBranch(j) * r;
      cout << setprecision(3) << setw(11) << lib->gamma[i].getBranch(j);
      cout << setprecision(3) << setw(11) << lib->gamma[i].getCvcoef(j);
      cout << endl;
    }
  }
}
