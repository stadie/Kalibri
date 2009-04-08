//!  \mainpage
//!
//!  \section label_sec_intro Introduction
//!  Package for data driven calibration using a global fit (see also the related
//!  <A HREF="https://twiki.cern.ch/twiki/bin/view/CMS/Calibration">
//!  Twiki Page</A>).
//!  \image html kalibri_workflow.jpg
//!
//!  \section label_sec_src Source Code
//!  The source code can be found
//!  <A HREF="http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Bromo/Calibration/CalibCore/">here</A>.
//!
//!  \section label_sec_relinfo Related information
//!  - C. Autermann:
//!    <A HREF="http://indico.cern.ch/getFile.py/access?contribId=7&resId=0&materialId=slides&confId=22705">
//!    A global fit approach to HCAL/jet calibration</A>,
//!    JetMET Meeting, 18th October 2007
//!  - R. Wolf:
//!    <A HREF="http://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=29582">
//!    Data-Driven Calorimeter Calibration Exploiting a Global-Fit Ansatz</A>,
//!    JetMET Meeting, 26th February 2008
//!  - T. Schum:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=1&amp;resId=0&amp;materialId=slides&amp;confId=627">
//!    HCAL Calibration using a Global Fit Ansatz</A>,
//!    Hamburg CMS Meeting, 12th March 2008
//!  - S. Naumann-Emme: Top as a Calibration Tool,
//!    2nd "Physics at the Terascale" Workshop, 27th November 2008
//!  - H. Stadie:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=1&amp;resId=0&amp;materialId=slides&amp;confId=1683">
//!    CMS Calorimeter and Jet Calibration</A>,
//!    CMS Hamburg Meeting, 28th January 2009
//!  - M. Schr&ouml;der:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=2&amp;resId=0&amp;materialId=slides&amp;confId=1683">
//!    Conceptual Studies for a Jet Energy Correction</A>,
//!    CMS Hamburg Meeting, 28th January 2009
//!  - S. Naumann-Emme:
//!    <A HREF="https://indico.desy.de/getFile.py/access?contribId=3&amp;resId=0&amp;materialId=slides&amp;confId=1683">
//!    Jet Energy Corrections from Top Quark Decays</A>,
//!    CMS Hamburg Meeting, 28th January 2009
//!
//!  \section label_sec_authors Authors
//!  - Christian Autermann
//!  - Ulla Gebbert
//!  - Robert Klanner
//!  - Bj&ouml;rn Kolodzey
//!  - Sebastian Naumann-Emme
//!  - Christian Sander
//!  - Matthias Schr&ouml;der
//!  - Torben Schum
//!  - Hartmut Stadie
//!  - Jan Thomsen
//!  - Roger Wolf
 

#ifndef caliber_h
#define caliber_h

#include <vector>
#include <string>

class TParameters;
class TControlPlots;
class TData;
class TMeasurement;


//!  \author Christian Autermann
//!  \date Wed Jul 18 13:54:50 CEST 2007
//!  $Id: caliber.h,v 1.40 2009/02/01 16:38:19 stadie Exp $
class TCaliber {
public :
  TCaliber(const std::string& f)
  : configfile(f),p(0),plots(0),deriv_step(1e-03),eps(1e-02),
  wlf1(1e-04),wlf2(0.9),print_parnderiv(false)
 {};
  ~TCaliber(){};

  void Init();
  void Run();
  void Done();
  const char * GetOutputFile(){ return output_file.c_str(); };

protected:  
  //internal functions
  void Run_Lvmini();

private:
  //internal variables
  int fit_method, n_gammajet_events, n_dijet_events;
  int n_trijet_events,n_trackcluster_events, n_zjet_events, n_top_events;
  std::string configfile, output_file;              //input/output
  //int use_GammaJetTowerMethod,use_DisplayMethod;    //plots
  //bool useMassConstraintW;
  //bool useMassConstraintTop;
 

  std::vector<int> _residualScalingScheme;          // Iteration scheme of scaling of residuals
  double OutlierChi2Cut;                            // Cut on outlier when no scaling is chosen
  int nthreads;
  bool flatten_spectra;
  std::vector<TData*> data;
  
  TParameters * p;    //fit parameters, depend on number of bins & geometry

  TControlPlots * plots;  //the control plots
  // control parameters of fit
  double deriv_step;
  float eps,wlf1,wlf2;
  bool print_parnderiv;
  std::vector<int> globaljetpars;
  std::vector<int> fixedpars;
};

#endif
