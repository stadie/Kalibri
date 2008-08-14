//
// Original Author:  Hartmut Stadie
//         Created:  Mon Jun 30 11:00:00 CEST 2008
// $Id: toy.cc,v 1.6 2008/07/25 15:34:01 stadie Exp $
//
#include "ToyMC.h"

#include <iostream>

#include "ConfigFile.h"

int main(int argc, char* argv[]) {
  if (argc!= 2) {
    std::cout << "Usage:\n    " << argv[0] << " <config.file>\n";
    return 1;
  }
  ToyMC* mc = new ToyMC();
  mc->init(argv[1]);
  mc->print();
  ConfigFile config(argv[1]);
  std::string file = config.read<std::string>("ToyMC output file","toy.root");
  int nevents = config.read<int>("ToyMC events",10);
  int mode = config.read<int>("ToyMC mode",1);
  if(mode == 1) 
    mc->makePhotonJet(file.c_str(),nevents);
  else if(mode == 2) 
    mc->makeDiJet(file.c_str(),nevents);
  else {
    std::cerr << "unknown mode:" << mode << '\n';
    delete mc;
    return 1;
  }
  delete mc;
  return 0;
}

