/*
   masstable.h : 
        define mass excess data file location
        prototype of function to read mass table
 */

class MassExcess{
 public:
  unsigned int za;    // Z*1000 + A
  float        mass;  // mass excess
};


/**************************************/
/*      masstable.cpp                 */
/**************************************/
double  mass_excess           (int, int, bool *);

