//
//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactors.h,v 1.1 2009/11/25 13:07:45 stadie Exp $
// 
#include "CorFactorsFactory.h"  
#include "CorFactors.h"

#include <iostream>

std::map<std::string,CorFactorsFactory*> CorFactorsFactory::map;


CorFactorsFactory::CorFactorsFactory(const std::string& name) : name_(name)
{
  static Cleaner cleaner;
  std::cout << "creating CorFactorsFactory: " << name_ << '\n';

  map[name] = this;
}

CorFactorsFactory::~CorFactorsFactory()
{
  std::cout << "deleting CorFactorsFactory: " << name_ << '\n';
}
