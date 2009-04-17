#ifndef EVENTREADER_H
#define EVENTREADER_H

class TData;
class TParameters;
class ConfigFile;
class TMeasurement;

#include <vector>
#include <string>


//!
//!    \brief Abstract base class for all event readers
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/12
//!    $Id: EventReader.h,v 1.2 2009/01/16 08:46:40 stadie Exp $
// ----------------------------------------------------------------   
class EventReader
{
 public:
  EventReader(const std::string& configfile, TParameters* p);
  virtual ~EventReader();
  virtual int readEvents(std::vector<TData*>& data) = 0;
 protected:
  ConfigFile* config;
  TParameters* p;
  bool useTracks;
  double (*tower_error_param)(const double *x, const TMeasurement *xorig, double err);
  double (*jet_error_param)  (const double *x, const TMeasurement *xorig, double err);
  double (*track_error_param)(const double *x, const TMeasurement *xorig, double err);
};


#endif
