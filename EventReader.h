//
// $Id: EventReader.h,v 1.5 2009/11/24 16:52:59 stadie Exp $
//
#ifndef EVENTREADER_H
#define EVENTREADER_H

class Event;
class TParameters;
class ConfigFile;
class Measurement;
class CorFactors;

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
  virtual CorFactors* createCorFactors(int jetid) const { return 0;}

  ConfigFile* config_;   //!< The configfile
  TParameters* par_;     //!< The parametrization
  bool useTracks_;       //!< True, if tracks are used in calibration

  double (*tower_error_param)(const double *x, const Measurement *xorig, double err);
  double (*jet_error_param)  (const double *x, const Measurement *xorig, double err);
  double (*track_error_param)(const double *x, const Measurement *xorig, double err);
};


#endif
