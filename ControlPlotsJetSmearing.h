// $Id: ControlPlotsJetSmearing.h,v 1.1 2009/06/11 17:34:05 mschrode Exp $

#ifndef JS_CONTROLPLOTS_JETSMEARING_H
#define JS_CONTROLPLOTS_JETSMEARING_H

#include <string>

#include "TCanvas.h"
#include "TH1F.h"
#include "TObject.h"
#include "TPostScript.h"

#include "CalibData.h"
#include "ConfigFile.h"
#include "Parameters.h"



//!  \brief Generates validation plots for jet-smearing method
//!  \author Matthias Schroeder
//!  \date Thu May  7 11:30:28 CEST 2009 
//!  $Id: ControlPlotsJetSmearing.h,v 1.1 2009/06/11 17:34:05 mschrode Exp $
// --------------------------------------------------
class ControlPlotsJetSmearing {
 public:
  ControlPlotsJetSmearing(const std::string& configfile,const std::vector<TData*> * data, TParameters * param);
  ~ControlPlotsJetSmearing() {};
  
  void plotDijets() const;
  void plotResponse() const;
  void plotParameterScan(const std::vector<unsigned int>& pars) const;
  void setBinningResp(int nbins, double min, double max) { respNBins_ = nbins; respMin_ = min; respMax_ = max;}


 private:
  typedef std::vector<TData*>::const_iterator DataIt;

  const std::vector<TData*> * data_;
  const ConfigFile          * config_;
  mutable TParameters       * param_;
  
  int          respNBins_;
  double       respMin_;
  double       respMax_;
  std::string  dir_;

  void drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, std::string option, bool logy) const;
  void drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, std::string option, bool logy) const;
  void normHist(TH1F * h, std::string option = "") const
    { if( h->Integral(option.c_str()) ) h->Scale(1./h->Integral(option.c_str())); }
  void setGStyle() const;

};
#endif
