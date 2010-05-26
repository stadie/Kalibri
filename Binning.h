#ifndef BINNING_H
#define BINNING_H
//!
//!  \brief provides bin ids
//!
//!  \author Hartmut Stadie
//!  \date 2008/12/12
//!  $Id: DiJetReader.h,v 1.18 2010/05/19 13:34:48 stadie Exp $
// ----------------------------------------------------------------   

#include <string>
#include <map>
#include <set>
#include <vector>

class JetBin;
class ConfigFile;

class Binning 
{
 public:
  Binning(const ConfigFile* configfile);
  JetBin* operator()(double x1, double x2 = 0, double x3 = 0, double x4 = 0) const
  {
    BinMap::const_iterator jbi = jetbins_.find(findBin(x1,x2,x3,x4));
    if(jbi == jetbins_.end()) return 0;
    return jbi->second;
  }
  JetBin* setBin(JetBin* jb, double x1, double x2 = 0, double x3 = 0, double x4 = 0) {
    jetbins_.insert(std::make_pair(findBin(x1,x2,x3,x4),jb));
    return (*this)(x1,x2,x3,x4);
  }
  typedef  std::map<int,JetBin*> BinMap;
  const BinMap& bins() const { return jetbins_;}
  void clear() {
    jetbins_.clear();
  }
 private:
  int findBin(double x1, double x2, double x3, double x4) const;
  typedef  std::vector<std::set<double> > BorderSets;
  BorderSets borders_;
  BinMap jetbins_;
};

#endif


