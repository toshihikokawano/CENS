/******************************************************************************/
/*  ensdfread.cpp                                                             */
/*        read ENSDF data file, extract levels and gamma branching ratios     */
/******************************************************************************/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "cens.h"
#include "terminate.h"
#include "elements.h"
#include "physicalconstant.h"

static int      ENSDFReadIdentification(const string);
static ZAnumber ENSDFReadZA(const string);
static int      ENSDFSeekNextRecord(const char, const int, const int);
static void     ENSDFParseLevelLine(const string, ENSDF *, const double);
static void     ENSDFParseGammaLine(const string, Gamma *, const double);
static int      ENSDFParseSpinParity(const string, int *, int *);
static double   ENSDFParseHalfLife(const string);

static inline bool isNumeric(const char c)
{
 if( (0x30 <= c && c <= 0x39) ||  // numeric
     (c == 0x20) || (c == 0x2b) || (c == 0x2d) || (c == 0x2e) || (c == 0x45) || (c == 0x65) // space + - . E e
   ) return true;
 return false;
}

#undef DEBUG
#ifdef DEBUG
static void print(ENSDF *);
#endif

static string *dbase;
static int     nline = 0;

/***********************************************************/
/*      Read ENSDF                                         */
/***********************************************************/
int ENSDFRead(ZAnumber za, string ensdfdir, string libname, ENSDF *lib)
{
  ostringstream os;
  ifstream      fp;
  string        file, str;

  /* if ZA number is given, look for the ENSDF file in the current directory */
  if((za.getZ() > 0) && (za.getA() > 0)){
    os << setw(3) << setfill('0') << za.getZ() << setw(3) << setfill('0') << za.getA();
    file = os.str();
    file = "ENSDF" + file + ".dat";
  }
  /* when file name is given */
  else file = libname;

  /* when directory is given by config.dat and file name does not contain directory */
  if((ensdfdir.length() > 0) && !strchr(libname.c_str(),'/')){
    /* remove if dir name includes a slash at the end */
    if(ensdfdir[ensdfdir.length() - 1] == '/') ensdfdir[ensdfdir.length() - 1] = '\0';
    file = ensdfdir + '/' + file;
  }

  message << "ENSDF file name " << file;
  Notice("ENSDFRead");

  fp.open(&file[0]);
  if(!fp){
    message << "ENSDF file " << file << " cannot open";
    TerminateCode("ENSDFRead");
  }

  /* scan the file length */
  nline = 0;
  while(!fp.eof()){
    getline(fp,str);
    if(str.length() > 0) nline ++;
  }
  message << "ENSDF file length " << nline << " lines";
  Notice("ENSDFRead");

  /* grab the entire ENSDF file */
  dbase = new string [nline];

  fp.clear();
  fp.seekg(0,ios::beg);
  for(int i=0 ; i<nline; i++){
    getline(fp,str);
    if(str.length() > 0) dbase[i] = str;
  }
  fp.close();


  /* read first line in ENSDF datafile */
  int c0 = 0; // main counter
  if(lib->getZ() == 0){
    ZAnumber za = ENSDFReadZA(dbase[c0]);
    lib->setZA(za.getZ(),za.getA());
  }
  lib->date = ENSDFReadIdentification(dbase[c0++]);

  /* read L records */
  int *cl = new int [lib->getNsize() + 1]; // index of L record
  while(c0 < nline){
    c0 = ENSDFSeekNextRecord('l',c0,nline); if(c0 < 0) break;
    /* remember the current L card location */
    cl[lib->getNlevel()] = c0;

    ENSDFParseLevelLine(dbase[c0++],lib,lib->getUnit());
  }
  /* insert the last line */
  cl[lib->getNlevel()] = nline;

  message << "total number of given levels " << lib->getNlevel();
  Notice("ENSDFRead");

  /* read G records */
  for(int i = 1 ; i < lib->getNlevel()-1 ; i++){
    /* look for G records between two L records */
    int p0 = cl[i];
    int p1 = cl[i+1];

    /* scan all gamma-rays between p0 and p1 */
    while(p0 > 0){
      p0 = ENSDFSeekNextRecord('g',p0,p1); if(p0 < 0) break;
      ENSDFParseGammaLine(dbase[p0++],&lib->gamma[i],lib->getUnit());
    }
  }

#ifdef DEBUG
  print(lib);
#endif

  delete [] dbase;
  delete [] cl;

  return(0);
}


/***********************************************************/
/*      First Line (Header) in ENSDF                       */
/***********************************************************/
int ENSDFReadIdentification(const string s)
{
  int date = atoi(s.substr(74, 6).c_str());

  message << "NUCID nuclide identification  " << s.substr( 0, 5).c_str();
  Notice("ENSDFReadIdentification");

  message << "DSID data set identification  " << s.substr( 9,30).c_str();
  Notice("ENSDFReadIdentification");

  message << "PUB publication information   " << s.substr(65, 9).c_str();
  Notice("ENSDFReadIdentification");

  message << "DATE entered in ENSDF         " << date;
  Notice("ENSDFReadIdentification");

  return date;
}


/***********************************************************/
/*      Determine Z and A from ENSDF if not Provided       */
/***********************************************************/
ZAnumber ENSDFReadZA(const string s)
{
  int    a = atoi(s.substr(0,3).c_str());
  string e = s.substr(3,2);
  if(e[1] == ' ') e[1] = '\0';
  int    z = element_getZ(&e[0]);

  ZAnumber za(z,a);
  return za;
}


/***********************************************************/
/*      Move to L or G Card in Database                    */
/***********************************************************/
int ENSDFSeekNextRecord(const char c, const int p0, const int p1)
{
  int p = 0;
  for(p=p0 ; p<p1 ; p++){
    string str = dbase[p];

    char c5 = tolower(str[5]);
    char c6 = tolower(str[6]);
    char c7 = tolower(str[7]);

    if(c5 != ' ') continue;
    if(c6 == 'c' || c7 == 'd') continue;
    if(c6 == ' ' && c7 == c) break;
  }
  if(p >= p1) return -1;
  else return p;
}


/***********************************************************/
/*      Parse L Record in ENSDF                            */
/***********************************************************/
void ENSDFParseLevelLine(const string line, ENSDF *lib, const double u)
{
  int j[Candidate_Spin], p[Candidate_Spin];

  /* check if non-numeric letters are given in the energy field */
  for(int i=9 ; i<19 ; i++) if(!isNumeric(line[i])) return;

  /* discrete level energy, given in keV */
  double e = atof(line.substr( 9,10).c_str()) * 1e+3 / u;

  /* spin and parity, their candidates */
  int nc = ENSDFParseSpinParity(line,j,p);
  
  /* half-life */
  double t = ENSDFParseHalfLife(line);
  
  /* copy data to object */
  lib->setLevel(e,t,nc,j,p);
}


/***********************************************************/
/*      Parse G Record in ENSDF                            */
/***********************************************************/
void ENSDFParseGammaLine(const string line, Gamma *gam, const double u)
{
  /* gamma-ray energy */
  double g = atof(line.substr( 9,10).c_str()) * 1e+3 / u;

  /* gamma-ray intensity */
  double r = atof(line.substr(21, 8).c_str());

  /* conversion coefficient */
  double c = atof(line.substr(55, 7).c_str());

  /* copy data to object */
  gam->setGamma(g,r,c);
}


/***********************************************************/
/*      Extract Candidate Spins And Parities               */
/***********************************************************/
int ENSDFParseSpinParity(const string line, int *j, int *p)
{
  /* parse spin and parity */
  bool blnk = true;
  string sdata = "";
  const int slength = 18;
  for(int i=0 ; i<slength ; i++) sdata += ' '; // pad 18 blanks

  for(int i=0 ; i<slength ; i++){
    sdata[i] = line[i+21];                // copy column 21 to 38
    if(sdata[i] != ' ') blnk = false;     // check if spin/parity field is empty
  }

  /* when mulitple spins are given, they are not certain */
  bool uniq = (strchr(sdata.c_str(),',')) ? false : true;
  if(blnk) uniq = false;

  /* check if spin is a half-integer */
  bool half = (strchr(sdata.c_str(),'/')) ? true : false;
  
  int nc = 0;

  /* no information given */
  if(blnk){
    j[nc] = -1;
    p[nc] =  0;
    nc ++;
  }

  else{
    /* only one value given even if it is in parenthises */
    if(uniq){
      /* parity */
      p[nc] = 0;
      for(int i=0 ; i<slength ; i++){
        if(     sdata[i] == '+'){ p[nc] =  1; break; }
        else if(sdata[i] == '-'){ p[nc] = -1; break; }
      }

      /* spin */
      int i0 = 0; // numeric char start
      for(int i=0 ; i<slength ; i++){
        if(isdigit(sdata[i])){ i0 = i; break; }
      }

      int i1 = 0; // numeric char end
      for(int i=i0 ; i<slength ; i++){
        if(!isdigit(sdata[i])){ i1 = i; break; }
      }

      j[nc++] = atoi(sdata.substr(i0,i1-i0).c_str()) * ((half) ? 1 : 2);
    }

    /* when multiple candidates are given */
    else{
      nc = 1;
      for(int i=0 ; i<slength ; i++){
        if(sdata[i] == ',') nc++;  // count number of candidates
      }
      if(nc > Candidate_Spin){
        message << "too many spin candidates " << nc;
        TerminateCode("ENSDFParseSpinParity");
      }

      int i0 = 0, i1 = 0, i2 = 0, n = 0;
      while(n < nc){
        /* pointers i1: numeric text starts, i2: numeric text ends */
        for(int i=i0 ; i<slength ; i++){
          if(isdigit(sdata[i])){ i1 = i; break; }
        }
        for(int i=i1 ; i<slength ; i++){
          if(!isdigit(sdata[i])){ i2 = i; break; }
        }
        /* move i0 to the next testing field */
        for(int i=i2 ; i<slength ; i++){
          if(sdata[i] == ','){ i0 = i; break; }
        }

        /* parity could be at the end of numeric, if given */
        p[n] = 0;
        if(     sdata[i2] == '+') p[n] =  1;
        else if(sdata[i2] == '-') p[n] = -1;

        /* spin is given between i1 and i2 */
        j[n++] = atoi(sdata.substr(i1,i2-i1).c_str()) * ((half) ? 1 : 2);
      }
    }
  }

  return nc;
}


/***********************************************************/
/*      Extract Half Life                                  */
/***********************************************************/
double ENSDFParseHalfLife(const string line)
{
  double t = 0.0;

  /* parse half-life */
  bool blnk = true;
  string sdata = "";
  const int slength = 10;
  for(int i=0 ; i<slength ; i++) sdata += ' ';

  for(int i=0 ; i<slength ; i++){
    sdata[i] = tolower(line[i+39]);       // copy column 39 to 48, make it lowercase
    if(sdata[i] != ' ') blnk = false;     // check if half-life field is empty
  }
  /* when half-life is not given, return */
  if(blnk) return t;

  /* when STABLE, return a negative value */
  if(strstr(sdata.c_str(),"stable" )){ return -1.0; }

  /* convert unit, I'm not sure if one year is defined as exact 365 days or not */
  double f = 1.0;
  if     (strstr(sdata.c_str(),"ms" )){f = 1.0e-03;}
  else if(strstr(sdata.c_str(),"us" )){f = 1.0e-06;}
  else if(strstr(sdata.c_str(),"ns" )){f = 1.0e-09;}
  else if(strstr(sdata.c_str(),"ps" )){f = 1.0e-12;}
  else if(strstr(sdata.c_str(),"fs" )){f = 1.0e-15;}
  else if(strstr(sdata.c_str(),"as" )){f = 1.0e-18;}
  else if(strstr(sdata.c_str(),"ev" )){f = -1.0e-6;} // negative for unit in energy in MeV
  else if(strstr(sdata.c_str(),"kev")){f = -1.0e-3;}
  else if(strstr(sdata.c_str(),"mev")){f = -1.0e+0;}
  /* these one letter cases should be evaluated after */
  else if(strstr(sdata.c_str(),"y"  )){f = 365.2422 * 24.0 * 60.0 * 60.0;}
  else if(strstr(sdata.c_str(),"d"  )){f =            24.0 * 60.0 * 60.0;}
  else if(strstr(sdata.c_str(),"h"  )){f =                   60.0 * 60.0;}
  else if(strstr(sdata.c_str(),"m"  )){f =                          60.0;}


  int i0 = 0; // numeric char start, including period
  for(int i=0 ; i<slength ; i++){
    if(isdigit(sdata[i]) || (sdata[i] == '.')){
      i0 = i;
      break;
    }
  }

  int i1 = 0; // numeric char end
  for(int i=i0 ; i<slength ; i++){
    /* skip "EV" case, maybe not happen */
    if((sdata[i] == 'e') && (i <= slength-2)){
      if(sdata[i+1] == 'v'){
        i1 = i-1;
        break;
      }
    }

    if((!isdigit(sdata[i])) && (sdata[i] != '.') && (sdata[i] != 'e')){
      i1 = i;
      break;
    }
  }

  t = atof(sdata.substr(i0,i1-i0).c_str()) * f;

  /* when time is given in eV */
  if(t < 0.0){
    double g = -t;
    t = log(2.0) * HBAR / g; // life-time = h-bar / Gamma [MeV s / MeV]
  }

  return t;
}


#ifdef DEBUG
/***********************************************************/
/*      Debugging Print                                    */
/***********************************************************/
void print(ENSDF *lib)
{
  cout << setprecision(5);
  for(int i = 0 ; i < lib->getNlevel() ; i++){

    cout << setw(5) << i;
    cout << setw(13) << lib->getEnergy(i);
    cout << setw(13) << lib->getThalf(i) << endl;

    for(int n=0 ; n<lib->nspin[i] ; n++){
      cout << "        ";
      cout << setw(3) << n;
      cout << setw(3) << (int)lib->spin[i][n].j;
      cout << setw(3) << (int)lib->spin[i][n].p << endl;
    }

    for(int j=0 ; j<lib->gamma[i].getNgamma() ; j++){
      cout << "                ";
      cout << setw(13) << lib->gamma[i].getEnergy(j);
      cout << setw(13) << lib->gamma[i].getBranch(j);
      cout << setw(13) << lib->gamma[i].getCvcoef(j) << endl;
    }
    cout << endl;
  }
}
#endif
