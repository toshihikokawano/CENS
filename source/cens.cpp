/******************************************************************************/
/**                                                                          **/
/**   CENS : Convert ENSDF into Hauser-Feshbach code library                 **/
/**                                               Version 0.3  Jul. 2022     **/
/**                                                            T. Kawano     **/
/**                                       Los Alamos National Laboratory     **/
/******************************************************************************/

#include <iostream>
#include <unistd.h>

using namespace std;

#include "cens.h"
#include "terminate.h"
#include "elements.h"
#include "cfgread.h"

static string version = "0.3 (Jul. 2022)";

static void CENSHelp(void);
static void CENSAllocMemory(void);
static void CENSFreeMemory(void);

static bool verbflag = false;
static ENSDF lib;
static char cfgdat[WORD_LENGTH];

/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/
#define CENS_TOPLEVEL
ostringstream message;


/**********************************************************/
/*      CENS Main                                         */
/**********************************************************/
int main (int argc, char *argv[])
{
  cout.setf(ios::scientific, ios::floatfield);
  cerr.setf(ios::scientific, ios::floatfield);

  string   libname_in = "",  libname_out = "", elem = "";
  int      anum = 0, znum = 0, popt = 0;

  /*** command line options */
  int p;
  while((p = getopt(argc,argv,"o:z:a:e:p:vh")) != -1){
    switch(p){
    case 'o': libname_out = optarg;   break;
    case 'z': elem = optarg;
              znum = element_getZ(&elem[0]);
              if(znum == 0){
                message << "unknown nuclide name " << elem;
                TerminateCode("main");
              }                        break;
    case 'a':  anum = atoi(optarg);    break;
    case 'p':  popt = atoi(optarg);    break;
    case 'v':  verbflag = true;        break;
    case 'h':  CENSHelp();             break;
    default:                           break;
    }
  }
  if(optind < argc) libname_in = argv[optind];

  if( (znum > 0 && anum == 0) || (znum == 0 && anum > 0) ){
    message << "invalid Z and A numbers, Z = " << znum << " A = " << anum;
    TerminateCode("main");
  }
  ZAnumber za(znum,anum);

  /* allocate ENSDF memory */
  CENSAllocMemory();

  /* set energy unit */
  lib.setUnit("MeV");
  if(CFGRead("EnergyUnit",cfgdat)){
    if((string)cfgdat == "MeV") lib.setUnit("MeV");
    else if((string)cfgdat == "keV") lib.setUnit("keV");
    else{
      message << "unknown energy unit " << cfgdat;
      TerminateCode("main");
    }

    message << "Energy Unit changed into " << cfgdat;
    Notice("main");
  }

  /* read ENSDF file and store information in ENSDF object */
  string ensdfdir = "";
  if(CFGRead("ENSDFDirectory",cfgdat)){
    message << "ENSDF directory changed into " << cfgdat;
    Notice("main");
    ensdfdir = (string)cfgdat;
  }

  /* determine RIPL directory */
  string ripldir = "";
  if(CFGRead("RIPLDirectory",cfgdat)){
    message << "RIPL directory changed into " << cfgdat;
    Notice("main");
    ripldir = (string)cfgdat;
  }

  /* read ENSDF data file */
  ENSDFRead(za,ensdfdir,libname_in,&lib);

  /* print raw data */
  if(popt == 1) OUTFxml(&lib);

  else{
    /* adjust gamma-ray energies and minimum fix of branching ratios */
    CENSGamma(&lib);

    /* read RIPL file for internal conversion coefficents if not given in ENSDF */
    if(ripldir.length() > 0) RIPLRead(ripldir,&lib);

    if(popt == 2) OUTFxml(&lib);

    else{
      /* statistical model analysis */
      StatProperty stp;
      CENSStat(&lib, &stp);

      if(popt == 3) OUTStatAnalysis(&lib, &stp);

      else if(popt == 4) OUTStatDensity(&lib, &stp);

      else if(popt == 0){
        /* print results in RIPL format */
        OUTFripl(stp.ncomp,stp.nmax,&lib);
      }
    }
  }

  /* free allocated */
  CENSFreeMemory();

  return 0;
}


/**********************************************************/
/*      Allocate / Free Memory                            */
/**********************************************************/
void CENSAllocMemory()
{
  int mlevel = MaxDiscreteLevels;
  int mgamma = MaxGammaLines;

  /* total number of levels */
  if(CFGRead("MaxDiscreteLevels",cfgdat)){
    mlevel = atoi(cfgdat);
    message << "maximum number of levels changed from " << MaxDiscreteLevels;
    message << " to " << mlevel << " by configuration";
    Notice("CENSAllocMemory");
  }

  /* total number of gammas */
  if(CFGRead("MaxGammaLines",cfgdat)){
    mgamma = atoi(cfgdat);
    message << "maximum number of gamma lines changed from " << MaxGammaLines;
    message << " to " << mgamma << " by configuration";
    Notice("CENSAllocMemory");
  }

  /* allocate heap memory */
  lib.memalloc(mlevel, mgamma);
}

void CENSFreeMemory()
{
  lib.memfree();
}


/**********************************************************/
/*      Help                                              */
/**********************************************************/
void CENSHelp()
{
  cout << "cens: ver." << version << endl;
  cout <<
    " command line options\n"
    " % cens -z Znum -a Anum -p N\n"
    "      atomic and mass numbers are given in the command line options\n"
    "      cens looks for default location for the ENSDF file\n"
    " % cens -p N ENSDF_file\n"
    "      read given ENSDF file\n"
    "     -p output option\n"
    "        N = 0 (or no -p option): print RIPL format\n"
    "          = 1: raw ENSDF data in XML\n"
    "          = 2: fixed ENSDF data in XML\n"
    "          = 3: print level density information\n";
  cout << endl;
  exit(0);
}


/**********************************************************/
/*     Warning Notice Message                             */
/**********************************************************/
void WarningMessage()
{
  cerr << "WARNING   : " << message.str()      << endl;
  message.str("");
}

void Notice(string module){
  if(module == "NOTE"){
    cerr << " (._.) " << message.str() << endl;
  }
  else{
    if(verbflag) cerr << " (@_@) [" << module << "] " << message.str() << endl;
  }
  message.str("");
}


/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int TerminateCode(string module)
{
  CENSFreeMemory();
  cerr << "ERROR     :[" << module << "] " << message.str() << endl;
  exit(-1);
}


