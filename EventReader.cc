//
//    Base class for all event readers
//
//    first version: Hartmut Stadie 2008/12/12
//    $Id: caliber.h,v 1.33 2008/11/20 16:38:03 stadie Exp $
//   
#include "EventReader.h"

#include "ConfigFile.h"
#include "Parameters.h" 


EventReader::EventReader(const std::string& configfile, TParameters* param) 
  : config(0),p(param)
{
  config = new ConfigFile(configfile.c_str());
  
  useTracks = config->read<bool>("use Tracks",true);
  if(p->GetNumberOfTrackParameters() < 1) useTracks = false;
  if(useTracks) cout<<"Tracks are used to calibrate jets"<<endl;
  else cout<<"Only Calorimeter information is used"<<endl;

  
  //Error Parametrization...
  //...for tracks:
  track_error_param = p->track_error_parametrization;
  //...for tower:
  string te = config->read<string>("tower error parametrization","standard"); 
  if (te=="standard")
    tower_error_param = p->tower_error_parametrization;
  else if (te=="fast")
    tower_error_param = p->fast_error_parametrization;
  else if (te=="Jans E parametrization")
    tower_error_param = p->jans_E_tower_error_parametrization;
  else if(te=="const")
    tower_error_param = p->const_error_parametrization;
  else if(te=="toy")
    tower_error_param = p->toy_tower_error_parametrization;
  else if(te=="jet")
    tower_error_param = p->jet_only_tower_error_parametrization;
  else  
    tower_error_param = p->tower_error_parametrization;
  //...for jets:
  string je = config->read<string>("jet error parametrization","standard");
  if (je=="standard")
    jet_error_param   = p->jet_error_parametrization;
  else if (je=="fast")
    jet_error_param   = p->fast_error_parametrization;
  else if (je=="dummy")
    jet_error_param   = p->dummy_error_parametrization;
  else if(je=="const")
    jet_error_param = p->const_error_parametrization;
  else if(je=="toy")
    jet_error_param   = p->toy_jet_error_parametrization;
  else if(je=="jet et")
    jet_error_param   = p->jet_only_jet_error_parametrization_et;
  else if(je=="jet energy")
    jet_error_param   = p->jet_only_jet_error_parametrization_energy;
  else  
    jet_error_param   = p->jet_error_parametrization;
}

EventReader::~EventReader()
{
  delete config;
}
