//
// $Id: EventReader.h,v 1.6 2009/11/25 13:07:45 stadie Exp $
//
#ifndef EVENTREADER_H
#define EVENTREADER_H

class Event;
class TParameters;
class ConfigFile;
class Measurement;
class CorFactors;
class CorFactorsFactory;

#include <vector>
#include <string>

class EventReader
{
 public:
  static unsigned int numberOfEventReaders_;   //!< Number of initialized event readers

  EventReader(const std::string& configfile, TParameters* p);
  virtual ~EventReader();
  virtual int readEvents(std::vector<Event*>& data) = 0;

 protected:
  // read CorFactors from Ntuple
  virtual CorFactors* createCorFactors(int jetid) const { return 0;}
 
  ConfigFile* config_;   //!< The configfile
  TParameters* par_;     //!< The parametrization
  bool useTracks_;       //!< True, if tracks are used in calibration
  CorFactorsFactory* corFactorsFactory_; //! Factory class for external source of CorFactors;

  double (*tower_error_param)(const double *x, const Measurement *xorig, double err);
  double (*jet_error_param)  (const double *x, const Measurement *xorig, double err);
  double (*track_error_param)(const double *x, const Measurement *xorig, double err);
};


#endif
