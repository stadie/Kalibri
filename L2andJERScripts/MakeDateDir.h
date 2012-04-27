#ifndef MakeDateDir_h
#define MakeDateDir_h


#include <iostream>
#include "TString.h"
#include "TDatime.h"

TString GetDateDir(){

  TDatime CurrentTime;
  //  CurrentTime.Print();
  
  TString helper_number="";
  TString today="";
  today+=CurrentTime.GetYear();
  today+="_";
  helper_number += CurrentTime.GetMonth();
  while(helper_number.Length()<2)helper_number="0"+helper_number;
  today+=helper_number;
  today+="_";
  helper_number="";
  helper_number+=CurrentTime.GetDay();
  while(helper_number.Length()<2)helper_number="0"+helper_number;
  today+=helper_number;
  today+="_plots";

  //  std::cout << today << std::endl;
  return today;
}

TString MakeDateDir(){
  TString today=GetDateDir();
  if(chdir(today) != 0){ 
    mkdir(today, S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir(today); 
  } 
  std::cout << "creating folder " << today << std::endl;
 return today;
}


#endif
