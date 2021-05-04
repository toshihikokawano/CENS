const int Record_Length = 80;
const int Candidate_Spin = 5;

//------------------------------------------------------------------------------
//     Class

/**********************************************************/
/*   Z and A numbers of Nucleus                           */
/**********************************************************/
class ZAnumber{ 
 private:
  unsigned int Z;
  unsigned int A;
 public:
  ZAnumber(){
    Z = 0;
    A = 0;
  }
  ZAnumber(int z, int a){
    Z = z;
    A = a;
  }
  void setZA(int z, int a){
    Z = z;
    A = a;
  }
  unsigned int getZ(){ return (Z); }
  unsigned int getA(){ return (A); }
  unsigned int getN(){ return (A-Z); }
  ZAnumber operator+(ZAnumber x){
    ZAnumber y;
    y.Z = Z + x.Z;
    y.A = A + x.A;
    return y;
  }
  ZAnumber operator-(ZAnumber x){
    ZAnumber y;
    y.Z = Z - x.Z;
    y.A = A - x.A;
    return y;
  }
  bool operator==(ZAnumber x){
    bool z = false;
    if( (Z == x.Z) && (A == x.A) ) z = true;
    return z;
  }
};


/**********************************************************/
/*   Spin and Parity                                      */
/**********************************************************/
class Spin{
 public:
  char j;  // spin of level, doubled number
  char p;  // parity -1 for odd, +1 even, 0 for unknown

  Spin(){
    init();
  }

  void set(int a, int b){
    j = (char)a;
    p = (char)b;
  }

  void init(){
    j = -1; // negative for unknown
    p = 0;
  }
};


/**********************************************************/
/*   Gamma-Ray Branch                                     */
/**********************************************************/
class Gamma{
 private:
  int      nsize;     // allocated memory size
  bool     allocated; // flag to know if heap memory is allocated
 public:
  int      ngamma;    // number of gamma-rays
  int      *fstate;   // level index of final state
  double   *energy;   // gamma-ray energy
  double   *branch;   // relative intensity or branching ratio
  double   *cvcoef;   // conversion coefficient

  Gamma(){
    nsize = 0;
    ngamma = 0;
    allocated = false;
  }

  ~Gamma(){
    memfree();
  }

  void memalloc(int n){
    if(!allocated){
      nsize = n;
      fstate = new int [n];
      energy = new double [n];
      branch = new double [n];
      cvcoef = new double [n];
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      nsize = 0;
      delete [] fstate;
      delete [] energy;
      delete [] branch;
      delete [] cvcoef;
      allocated = false;
    }
  }

  void init(){
    if(allocated){
      for(int i=0 ; i<nsize ; i++){
        fstate[i] = 0;
        energy[i] = branch[i] = cvcoef[i] = 0.0;
      }
    }
  }

  bool setGamma(double a, double b, double c){
    if(ngamma >= nsize-1) return false;
    else{
      energy[ngamma] = a;
      branch[ngamma] = b;
      cvcoef[ngamma] = c;
      ngamma ++;
      return true;
    }
  }

  int getNgamma(void){ return ngamma; }
  int getNsize(void){ return nsize; }

  int getFstate(int i){
    int k = -1;
    if(0 <= i && i < ngamma) k = fstate[i];
    return k;
  }
  
  double getEnergy(int i){
    double e = -1.0;
    if(0 <= i && i < ngamma) e = energy[i];
    return e;
  }

  double getBranch(int i){
    double b = -1.0;
    if(0 <= i && i < ngamma) b = branch[i];
    return b;
  }

  double getCvcoef(int i){
    double c = -1.0;
    if(0 <= i && i < ngamma) c = cvcoef[i];
    return c;
  }
};


/**********************************************************/
/*   ENSDF Database                                       */
/**********************************************************/
class ENSDF{ 
 private:
  ZAnumber za;        // Z and A numbers
  int      nsize;     // allocated memory size for levels
  bool     allocated; // flag to know if heap memory is allocated
  int      nlevel;    // number of actual discrete levels
  double   ebase;     // energy unit, 1 for eV, 1000 for keV, etc.
 public:
  int      date;      // file created date
  double   *energy;   // excitation energy in MeV
  double   *thalf;    // half-life in second
  int      *nspin;    // number of candidate spins
  Spin     **spin;    // spin and parity
  Gamma    *gamma;    // gamma-rays

  ENSDF(){
    za.setZA(0,0);
    nsize = 0;
    nlevel = 0;
    ebase = 1.0;      // default energy unit = eV
    allocated = false;
  }

  ~ENSDF(){
    memfree();
  }

  void memalloc(int n, int m){
    if(!allocated){
      nsize = n;
      energy = new double [nsize];
      thalf = new double [nsize];
      nspin = new int [nsize];
      spin = new Spin * [nsize];
      for(int i=0 ; i<nsize ; i++) spin[i] = new Spin [Candidate_Spin];

      gamma = new Gamma [nsize];
      for(int i=0 ; i<nsize ; i++) gamma[i].memalloc(m);

      allocated = true;
    }

    init();
  }

  void memfree(){
    if(allocated){
      nsize = 0;
      delete [] energy;
      delete [] thalf;
      delete [] nspin;
      for(int i=0 ; i<nsize ; i++) delete [] spin[i];
      delete [] spin;
      for(int i=0 ; i<nsize ; i++) gamma[i].memfree();
      delete [] gamma;

      allocated = false;
    }
  }

  void init(){
    if(allocated){
      for(int i=0 ; i<nsize ; i++){
        energy[i] = 0.0;
        thalf[i]  = 0.0;
        nspin[i]  = 0;
        for(int j=0 ; j<Candidate_Spin ; j++) spin[i][j].init();
        gamma[i].init();
      }
    }
  }

  void setZA(unsigned int z, unsigned int a){
    za.setZA(z,a);
  }

  bool setLevel(double e, double t, int n, int *j, int *p){
    if(nlevel >= nsize-1) return false;
    else{
      energy[nlevel] = e;
      thalf[nlevel]  = t;
      nspin[nlevel]  = n;
      for(int i=0 ; i<n ; i++) spin[nlevel][i].set(j[i],p[i]);
      nlevel ++;
      return true;
    }
  }

  void setUnit(std::string u){
    if(u == "eV") ebase = 1.0;
    else if(u == "keV") ebase = 1.0e+3;
    else if(u == "MeV") ebase = 1.0e+6;
  }

  int getNlevel(void){ return nlevel; }
  int getNsize(void){ return nsize; }

  double getEnergy(int i){
    double e = -1.0;
    if(0 <= i && i < nlevel) e = energy[i];
    return e;
  }

  double getThalf(int i){
    double t = -1.0;
    if(0 <= i && i < nlevel) t = thalf[i];
    return t;
  }

  double getUnit(){ return ebase; }

  int getZ(){ return(za.getZ()); }
  int getA(){ return(za.getA()); }

  ZAnumber getZA() { return za; }
};



//------------------------------------------------------------------------------
//     Prototype Definition

// ensdfread.cpp
int  ENSDFRead(ZAnumber, std::string, std::string, ENSDF *);

// riplread.cpp
int  RIPLRead(std::string, ENSDF *);
