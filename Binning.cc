//!
//!  \brief provides bin ids
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: DiJetReader.h,v 1.18 2010/05/19 13:34:48 stadie Exp $
// ----------------------------------------------------------------   
#include "Binning.h"
#include "JetBin.h"
#include "ConfigFile.h"

#include <iostream>

Binning::Binning(const ConfigFile* config) : borders_(4)
{   
  const std::string name = "jet binning";
  std::vector<std::string> vars = bag_of_string(config->read<std::string>(name+" variables",""));
  std::cout << "Binning:  variable  #bins\n";
  for(std::vector<std::string>::const_iterator i = vars.begin() ; 
      i != vars.end() ; ++i) {
    std::vector<double> vd =  bag_of<double>(config->read<std::string>(name + " " + *i +" bins","0"));
    std::cout << "              " << *i << "    " << vd.size() <<  '\n';
    borders_[i - vars.begin()] = std::set<double>(vd.begin(),vd.end());
  } 
}


int Binning::findBin(double x1, double x2, double x3, double x4) const
{
  BorderSets::const_iterator i = borders_.begin();
  int id1 =  std::distance(i->begin(),i->lower_bound(x1));
  ++i;
  int id2 =  std::distance(i->begin(),i->lower_bound(x2));
  ++i;
  int id3 =  std::distance(i->begin(),i->lower_bound(x3));
  ++i;
  int id4 =  std::distance(i->begin(),i->lower_bound(x4));
  //std::cout << "x:" << x1 << "," << x2 << "," << x3 << "," << x4 << std::endl;
  //std::cout << "id:" <<  id1 + id2 * 100 + id3 * 10000 + id4 * 1000000
  //	    << std::endl;
  return id1 + id2 * 100 + id3 * 10000 + id4 * 1000000;
}
