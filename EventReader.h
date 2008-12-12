//
//    Base class for all event readers
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: caliber.h,v 1.33 2008/11/20 16:38:03 stadie Exp $
//   
#ifndef EVENTREADER_H
#define EVENTREADER_H

class TData;
class TParameters;
class ConfigFile;
class TMeasurement;

#include <vector>
#include <string>

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
  double const (*tower_error_param)(double *const x, TMeasurement *const xorig, double const err);
  double const (*jet_error_param)  (double *const x, TMeasurement *const xorig, double const err);
  double const (*track_error_param)(double *const x, TMeasurement *const xorig, double const err);
};


#endif
