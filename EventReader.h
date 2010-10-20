//
// $Id: EventReader.h,v 1.14 2010/10/12 08:38:59 stadie Exp $
//
#ifndef EVENTREADER_H
#define EVENTREADER_H

class Event;
class Parameters;
class ConfigFile;
class Measurement;
class CorFactors;
class CorFactorsFactory;
class TTree;
class JetConstraintEvent;
class Parametrization;
class Binning;

#include <vector>
#include <string>

class EventReader
{
 public:
  static unsigned int numberOfEventReaders_;   //!< Number of initialized event readers

  EventReader(const std::string& configfile, Parameters* p);
  virtual ~EventReader();
  virtual int readEvents(std::vector<Event*>& data) = 0;
  virtual int readControlEvents(std::vector<Event*>& control, int id) { 
    return 0;
  }

  static int addConstraints(std::vector<Event*>& data);

 protected:
  //! Read CorFactors from Ntuple
  virtual CorFactors* createCorFactors(int jetid) const { return 0;}
  //! Create TTree with data files
  TTree* createTree(const std::string& name) const;

  ConfigFile* config_;   //!< The configfile
  Parameters* par_;     //!< The parametrization
  bool useTracks_;       //!< True, if tracks are used in calibration
  CorFactorsFactory* corFactorsFactory_; //! Factory class for external source of CorFactors;
  //! Correct jets to L3 i.e. with L1*L2*L3
  bool correctToL3_;
  //!< Correct jets with L2*L3 corrections
  bool correctL2L3_;
  double weightRelToNtuple_;


  float (*tower_error_param)(const float *x, const Measurement *xorig, float err);
  float (*jet_error_param)  (const float *x, const Measurement *xorig, float err);
  float (*track_error_param)(const float *x, const Measurement *xorig, float err);  

  static std::vector<JetConstraintEvent*> constraints_;
  Parametrization *cp_;
  static Binning* binning_;
};


#endif
