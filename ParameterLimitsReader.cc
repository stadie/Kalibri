//
//    Reader for Parameter Limits
//
//    This class add user defined parameter limits
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: PhotonJetReader.h,v 1.1 2008/12/12 13:43:15 stadie Exp $
//   
#include "ParameterLimitsReader.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"

#include <iostream>

ParameterLimitsReader::ParameterLimitsReader(const std::string& configfile, TParameters* p) :
  EventReader(configfile,p)
{
  vector<double> limits = bag_of<double>(config->read<string>( "Jet Parameter Limits",""));
  
  if(limits.size() % 4 == 0) {
    for(unsigned int i = 0 ; i < limits.size() ; i += 4) {
      int index = (int)limits[i];
      for(int j = p->GetNumberOfTowerParameters() + index; 
	  j <  p->GetNumberOfParameters() ; 
	  j += p->GetNumberOfJetParametersPerBin()) {
	par_limits.push_back(ParameterLimit(j,limits[i+1],limits[i+2],
					    limits[i+3]));
      }
    }
  } else if(limits.size() > 1) {
    std::cout << "wrong number of arguments for Jet Parameter Limits:" 
	      << limits.size() << '\n';
  }
  limits.clear();
  limits = bag_of<double>(config->read<string>( "Tower Parameter Limits",""));
  
  if(limits.size() % 4 == 0) {
    for(unsigned int i = 0 ; i < limits.size() ; i += 4) {
      int index = (int)limits[i];
      for(int j = index; 
	  j <  p->GetNumberOfTowerParameters() ; 
	  j += p->GetNumberOfTowerParametersPerBin()) {
	par_limits.push_back(ParameterLimit(j,limits[i+1],limits[i+2],
					    limits[i+3]));
      }
    }
  } else if(limits.size() > 1) {
    std::cout << "wrong number of arguments for Tower Parameter Limits:" 
	      << limits.size() << '\n';
  }
  
  delete config;
  config = 0;
}

ParameterLimitsReader::~ParameterLimitsReader()
{
}

int ParameterLimitsReader::readEvents(std::vector<TData*>& data)
{
  for(std::vector<ParameterLimit>::const_iterator pl = par_limits.begin() ;
      pl != par_limits.end() ; ++pl) {
    std::cout << "adding limit for parameter " << pl->index+1 << " min:" 
	      << pl->min << " max:" << pl->max << " k:" << pl->k << '\n';
    TMeasurement* limitp  = new TMeasurement;
    limitp->pt  = pl->min;//@@Add new TMeasurement derivative
    limitp->eta = pl->max;
    TData_ParLimit * parlim = new TData_ParLimit(pl->index,limitp,pl->k,p->GetPars()+pl->index,p->parameter_limit);
    data.push_back(parlim);
  }
  return par_limits.size();
}
