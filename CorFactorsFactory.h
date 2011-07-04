//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactorsFactory.h,v 1.3 2011/02/15 12:53:14 stadie Exp $
//   
#ifndef CORFACTORSFACTORY_H
#define CORFACTORSFACTORY_H

#include <map>
#include <string>

class Jet;
class CorFactors;

class CorFactorsFactory
{
 public:
  CorFactorsFactory(const std::string& name);
  virtual ~CorFactorsFactory();
  virtual CorFactors* create(const Jet* j,int nPV = 1, double rho=-9999., double jetA=0.) = 0;
  virtual CorFactorsFactory* clone() const = 0;

  static CorFactorsFactory* get(std::string name) {
    CorFactorsFactory* ccf = map[name];
    return ccf ? ccf->clone() : 0;
  }
 private:
  std::string name_;
  static std::map<std::string,CorFactorsFactory*> map;
  class Cleaner
  {
  public:
    Cleaner() {}
    ~Cleaner()
      {
	for( std::map<std::string,CorFactorsFactory*>::iterator i = map.begin();
	     i != map.end() ; ++i) delete i->second;
	map.clear();
      }
  };
  friend class Cleaner;
};

#endif
