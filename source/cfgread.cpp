//------------------------------------------------------------------------------
// Configuration Data Read
//------------------------------------------------------------------------------

#include <cstring>
#include <iostream>
#include <fstream>

#include "cfgread.h"
static char cfgread_key[WORD_LENGTH], cfgread_dat[WORD_LENGTH];

bool CFGRead(const std::string key, char *dat)
{
  std::ifstream  fp;
  std::string line;
  std::string str = CONFIG_FILE;

  dat[0] = '\0';

  fp.open(&str[0]);
  if(!fp) return false;

  bool found = false;
  while(getline(fp,line)){
    if(line[0] == '#') continue;

    unsigned int n = line.length();
    if(n <= 1) continue;

    // find = sign
    unsigned int im = 0;
    for(im=0 ; im<n ; im++){
      if(line[im] == ' ' || line[im] == '\t') continue;
      if(line[im] == '=') break;
    }

    // item before =
    unsigned int k = 0;
    for(unsigned int i=0 ; i<im ; i++){
      if(line[i] == ' ' || line[i] == '\t') continue;
      cfgread_key[k++] = line[i];
      if(k == sizeof(cfgread_key)-1) break;
    }
    cfgread_key[k] = '\0';

    // item after =
    k = 0;
    for(unsigned int i=im+1 ; i<n ; i++){
      if(line[i] == ' ' || line[i] == '\t') continue;
      cfgread_dat[k++] = line[i];
      if(k == sizeof(cfgread_dat)-1) break;
    }
    cfgread_dat[k] = '\0';

    // check key if this is the same as the given keyword
    if(key == (std::string)cfgread_key){
      strncpy(dat,cfgread_dat,WORD_LENGTH);
      found = true;
      break;
    }
  }
  fp.close();

  return found;
}


