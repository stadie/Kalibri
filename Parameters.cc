// $Id: Parameters.cc,v 1.63 2011/06/23 07:53:36 stadie Exp $

#include <fstream>
#include <cassert>
#include <pwd.h>
#include <unistd.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <algorithm>

#include "Parameters.h"

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"

using namespace std;

Parameters* Parameters::instance_ = 0;
std::vector<Parameters*> Parameters::clones_;
Parametrization* Parameters::p_ = 0;
ResolutionParametrization* Parameters::resParam_ = 0;
long long Parameters::ncalls_ = 0;
long long Parameters::ntries_ = 0;
long long Parameters::nfails_ = 0;
long long Parameters::nwarns_ = 0;

// -----------------------------------------------------------------
Parametrization* Parameters::createParametrization(const std::string& name, const ConfigFile& config) {
  if(name == "StepParametrization") {
    return new StepParametrization();
  } else if(name == "StepParametrizationEnergy") {
    return new StepParametrizationEnergy();
  } else if(name == "StepEfracParametrization") {
    return new StepEfracParametrization();
  } else if(name == "StepJetParametrization") {
    return new StepJetParametrization();
  } else if(name == "MyParametrization") {
    return new MyParametrization();
  }  else if(name == "JetMETParametrization") {
    return new JetMETParametrization();
  }  else if(name == "GlobalScaleFactorParametrization") {
    return new GlobalScaleFactorParametrization();
  }  else if(name == "SimpleParametrization") {
    return new SimpleParametrization();
  }  else if(name == "ToyParametrization") {
    return new ToyParametrization();
  }  else if(name == "ToyJetParametrization") {
    return new ToyJetParametrization(); 
  }  else if(name == "ToyStepParametrization") {
    return new ToyStepParametrization();
  }  else if(name == "ToyStepJetParametrization") {
    return new ToyStepJetParametrization();
  } else if(name == "TrackParametrization") {
    return new TrackParametrization();
  } else if(name == "L2L3JetParametrization") {
    return new L2L3JetParametrization();
  } else if(name == "L2L3JetParametrization2") {
    return new L2L3JetParametrization2();
  } else if(name == "L2L3JetTrackParametrization") {
    return new L2L3JetTrackParametrization();  
  } else if(name == "ResidualJetParametrization") {
    return new ResidualJetParametrization();
  } else if(name == "ToySimpleInverseParametrization") {
    return new ToySimpleInverseParametrization();
  } else if(name == "GroomParametrization") {
    return new GroomParametrization();
  } else if(name == "EtaEtaParametrization") {
    return new EtaEtaParametrization();
  } else if(name == "PhiPhiParametrization") {
    return new PhiPhiParametrization();
  } else if(name == "BinnedEMFParametrization") {
    return new BinnedEMFParametrization();
  } else if(name == "BinnedPhiPhiParametrization") {
    return new BinnedPhiPhiParametrization();
  } else if(name == "BinnedPhiPhiParametrization2") {
    return new BinnedPhiPhiParametrization2();
  } else if(name == "BinnedScaledPhiPhiParametrization") {
    return new BinnedScaledPhiPhiParametrization();
  } else if(name == "BinnedScaledEtaEtaParametrization") {
    return new BinnedScaledEtaEtaParametrization();
  } else if(name == "SimplePhiPhiParametrization") {
    return new SimplePhiPhiParametrization();
  } else if(name == "MeanWidthParametrization") {
    return new MeanWidthParametrization();
  } else if(name == "MeanWidthParametrization_old") {
    return new MeanWidthParametrization_old();
  }
  return 0;
}


// -----------------------------------------------------------------
ResolutionParametrization* Parameters::createResolutionParametrization(const std::string& name, const ConfigFile& config) {
  std::vector<double> ptBinEdges = bag_of<double>(config.read<std::string>("PtAve bin edges","-2 -1"));
  if( ptBinEdges.front() == -2. ) {
    ptBinEdges.at(0) = config.read<double>("Et min cut on dijet",-1.);
    ptBinEdges.at(1) = config.read<double>("Et max cut on dijet",-1.);
  }
  assert( ptBinEdges.size() > 1 && ptBinEdges.front() >= 0. );
  unsigned int nPtBins = ptBinEdges.size()-1;
  ResolutionParametrization *param = 0;
  if(name == "ResolutionGaussAvePt") {
    param = new ResolutionGaussAvePt(nPtBins);
  } else if(name == "ResolutionGauss") {
    // Read spectrum from file and split into pdfs for different pt bins
    TH1 *hPtGen = 0;
    std::string spectrum = config.read<string>("jet spectrum","default.root");
    std::cout << "Getting spectrum from file '" << spectrum << "'\n";
    TFile file(spectrum.c_str(),"READ");
    file.GetObject("hPtGen",hPtGen);
    if( !hPtGen ) {
      std::cerr << "ERROR: No histogram 'hPtGen' found in file '" << file.GetName() << "'\n";
      exit(1);
    }
    hPtGen->SetDirectory(0);
    file.Close();

    // Split into spectra per bin
    std::vector<TH1*> spectra(nPtBins);
    std::cout << "Creating underlying pdfs per bin\n";
    for(unsigned int i = 0; i < nPtBins; ++i) {
      double ptTrueMin = 0.4*ptBinEdges[i];
      double ptTrueMax = 1.8*ptBinEdges[i+1];
      if( nPtBins == 1 ) {
	ptTrueMin = config.read<double>("Et genJet min",ptTrueMin);
	ptTrueMax = config.read<double>("Et genJet max",ptTrueMax);
      }
      int binMin = hPtGen->FindBin(ptTrueMin);
      int binMax = hPtGen->FindBin(ptTrueMax);
      TString name = "spectrum";
      name += i;
      spectra[i] = new TH1D(name,"",1+binMax-binMin,hPtGen->GetXaxis()->GetBinLowEdge(binMin),hPtGen->GetXaxis()->GetBinUpEdge(binMax));
      for(int xBin = 1; xBin <= spectra[i]->GetNbinsX(); ++xBin) {
	spectra[i]->SetBinContent(xBin,hPtGen->GetBinContent(binMin+xBin-1));
      }
      if( spectra[i]->Integral("width") ) spectra[i]->Scale(1./spectra[i]->Integral("width"));
    }

    param = new ResolutionGauss(nPtBins,ptBinEdges,spectra);

  } else {
    param = new ResolutionEmpty(1);
  }

  return param;
}


// -----------------------------------------------------------------
Parameters* Parameters::createParameters(const ConfigFile& config) 
{
  static Cleaner cleanup;
  if(  instance_ != 0  )
  {
    delete instance_; 
    instance_ = 0; 
  }  
  
  string parclass = config.read<string>("Parametrization Class","");
  //create Parameters
  if(parclass == "TStepParameters") {
    parclass = "StepParametrization";
  } else if(parclass == "TMyParameters") {
    parclass = "MyParametrization";
  } else if(parclass == "TStepParametersEnergy") {
    parclass = "StepParametrizationEnergy";
  } else if(parclass == "TStepEfracParameters") {
    parclass = "StepEfracParametrization";
  } else if(parclass == "TJetMEParameters") {
    parclass = "JetMETParametrization";
  }  else if(parclass == "TGlobalScaleFactorParameters") {
    parclass = "GlobalScaleFactorParametrization";
  }  else if(parclass == "TSimpleParameters") {
    parclass = "SimpleParametrization";
  }  else if(parclass == "TToyParameters") {
    parclass = "ToyParametrization";
  }  else if(parclass == "TToyJetParameters") {
    parclass = "ToyJetParametrization";
  }  else if(parclass == "TToyStepParametersEnergy") {
    parclass = "ToyStepParametrizationEnergy";
  } else if(parclass == "StepJetParametrization") {
    parclass = "StepJetParametrization";
  } else if(parclass == "TTrackParameters") {
    parclass = "TrackParametrization";
  }
  
  int mode = config.read<int>("Mode",0);
  if( mode == 0 ) {
    Parametrization *param = createParametrization(parclass,config);
    if(! param) {
      cerr << "Parameters::createParameters: could not instantiate class " << parclass << '\n';
      exit(1);
    }
    instance_ = new Parameters(param);
  } else if( mode == 1 ) {
    ResolutionParametrization *resParam = createResolutionParametrization(parclass,config);
    if( !resParam ) {
      cerr << "Parameters::createParameters: could not instantiate class " << parclass << '\n';
      exit(1);
    }
    instance_ = new Parameters(resParam);
  } else {
    std::cerr << "Parameters::createParameters: unknown mode '" << mode << "'\n";
    exit(1);
  }
  
  instance_->init(config);

  return instance_;
}


// -----------------------------------------------------------------
Parameters::Parameters() 
  : k_(0),parErrors_(0),parGCorr_(0),parCov_(0),trackEff_(0),
    s_(gsl_root_fsolver_alloc(gsl_root_fsolver_brent)) {
  gsl_set_error_handler_off();
}


// -----------------------------------------------------------------
Parameters::Parameters(Parametrization* p) 
  : k_(0),parErrors_(0),parGCorr_(0),parCov_(0),trackEff_(0),
    s_(gsl_root_fsolver_alloc(gsl_root_fsolver_brent)) {
  if(p) p_ = p;
  gsl_set_error_handler_off();
  if( !resParam_ ) resParam_ = new ResolutionEmpty(0);
}

// -----------------------------------------------------------------
Parameters::Parameters(ResolutionParametrization *resParam) 
  : k_(0),parErrors_(0),parGCorr_(0),parCov_(0),trackEff_(0),
    s_(gsl_root_fsolver_alloc(gsl_root_fsolver_brent)) {
  if( resParam ) resParam_ = resParam;
  if( !p_ ) p_ = new EmptyParametrization(resParam_->nPtBins()*resParam_->nParPerPtBin());
  gsl_set_error_handler_off();
}

// -----------------------------------------------------------------
Parameters& Parameters::operator=(const Parameters& p) {
    return *instance_;
}



// -----------------------------------------------------------------
void Parameters::init(const ConfigFile& config)
{
  eta_ntwr_used_   = config.read<unsigned>("maximum eta twr used",82); 
  eta_granularity_ = config.read<unsigned>("granularity in eta",1); 
  phi_granularity_ = config.read<unsigned>("granularity in phi",1); 
  eta_symmetry_    = config.read<bool>("symmetry in eta",false);
  eta_granularity_jet_ = config.read<unsigned>("jet granularity in eta",1); 
  phi_granularity_jet_ = config.read<unsigned>("jet granularity in phi",1); 
  eta_granularity_track_ = config.read<unsigned>("track granularity in eta",1); 
  phi_granularity_track_ = config.read<unsigned>("track granularity in phi",1); 

  if (eta_ntwr_used_%2 !=0){
    cerr << "WARNING: Use even number of eta towers! Forced exit."<< endl;    
    exit(1);
  }

  if (phi_ntwr_%phi_granularity_!=0) {
    cerr << "WARNING: Check phi granularity! Forced exit."<< endl;
    exit(1);
  }
      
  if (eta_symmetry_ && (eta_granularity_!=1 && eta_granularity_!=3 &&eta_granularity_!=4 &&eta_granularity_!=5&&eta_granularity_!=11&&
      eta_granularity_!=21 && eta_granularity_!=41 )){
    cerr << "WARNING: Check eta granularity! Should be 1, 3, 4, 5, 11, 21, or 41: Forced exit."<< endl;
    exit(1);
  }

  if (!eta_symmetry_ && (eta_granularity_!=2 && eta_granularity_!=6 &&eta_granularity_!=10&&eta_granularity_!=22&&
      eta_granularity_!=42 && eta_granularity_!=82 )){
    cerr << "WARNING: Check eta granularity! Should be 2, 6, 10, 22, 42 or 82: Forced exit."<< endl;
    exit(1);
  }

  start_values_ = bag_of<double>(config.read<string>("start values","")); 
  if ( start_values_.size()< p_->nTowerPars()){
    cerr<< "ERROR: Number of start values and free parameters does not match!"<<endl
        << "       There must be at least " << p_->nTowerPars() << " parameters!" << endl;
    exit(2);    
  }

  jet_start_values_ = bag_of<double>(config.read<string>("jet start values",""));
  if ( jet_start_values_.size()< p_->nJetPars()){
    if( jet_start_values_.size() == resParam_->nParPerPtBin() ) { // for resolution parametrisation; add better condition
      for(unsigned int i = 1; i < resParam_->nPtBins(); ++i) {
	for(unsigned int j = 0; j < resParam_->nParPerPtBin(); ++j) {
	  jet_start_values_.push_back(jet_start_values_[0+j]);
	}
      }
    } else {
      cerr<< "ERROR: Number of jet start values and free jet parameters does not match!"<<endl
	  << "       There must be at least " << p_->nJetPars() << " parameters!" << endl;
      exit(3);
    }
  }
  track_start_values_ = bag_of<double>(config.read<string>("track start values","")); 
  if ( track_start_values_.size()< p_->nTrackPars()){
    cerr<< "ERROR: Number of track start values and free track parameters does not match!"<<endl
        << "       There must be at least " << p_->nTrackPars() << " parameters!" << endl;
    exit(3);
  }
  global_jet_start_values_ = bag_of<double>(config.read<string>("global jet start values","")); 
  if( global_jet_start_values_.size() < p_->nGlobalJetPars() ) {
    cerr<< "ERROR: Number of global jet start values and free global jet parameters does not match!"<<endl
        << "       There must be at least " << p_->nGlobalJetPars() << " parameters!" << endl;
    exit(3);
  }

  // Initialize storage for parameter values and errors
  k_ = new double[numberOfParameters()];
  parErrors_ = new double[numberOfParameters()];
  parGCorr_ = new double[numberOfParameters()];
  parCov_ = new double[numberOfCovCoeffs()];
  for(int i = 0; i < numberOfParameters(); i++) {
    k_[i] = 0.;
    parErrors_[i] = 0.;
    parGCorr_[i] = 0.;
  }
  for(int i = 0; i < numberOfCovCoeffs(); ++i) {
    parCov_[i] = 0.;
  }
  trackEff_ = new double[169];
  isFixedPar_ = std::vector<bool>(numberOfParameters(),false); 

  for (unsigned int bin=0; bin<eta_granularity_*phi_granularity_; ++bin){
    for (unsigned int tp=0; tp < p_->nTowerPars(); ++tp){
      k_[ bin*p_->nTowerPars() + tp ] = start_values_[ tp ];
    }
  }
  for (unsigned int bin=0; bin<eta_granularity_jet_*phi_granularity_jet_; ++bin){
    for (unsigned int jp=0; jp < p_->nJetPars(); ++jp){
      int i = numberOfTowerParameters() + bin*p_->nJetPars() + jp;   
      k_[i] = jet_start_values_[jp];
    }
  }
  for (unsigned int bin=0; bin<eta_granularity_track_*phi_granularity_track_; ++bin){
    for (unsigned int trp=0; trp < p_->nTrackPars(); ++trp){
      int i = numberOfTowerParameters() + numberOfJetParameters() + bin*p_->nTrackPars() + trp;   
      k_[i] = track_start_values_[trp];
    }
  }

  for(int etabin=0; etabin<13; ++etabin)
    {
      for(int ptbin=0; ptbin<13; ++ptbin)
	{
	  trackEff_[13*etabin+ptbin] = 1;
	}
    }

  for (unsigned int gjp = 0 ; gjp < p_->nGlobalJetPars() ; ++gjp){
    int i = numberOfTowerParameters() + numberOfJetParameters() + numberOfTrackParameters() + gjp;   
    k_[i] = global_jet_start_values_[gjp];
  }

  // read predefined calibration contants from cfi
  // or txt file depending on the ending of the name
  std::vector<std::string> inputCalibration = bag_of_string(config.read<string>("input calibration",";"));

  // Check whether start values for calibration are to be read
  // from file
  if( inputCalibration.empty() ) {
    cout << "Using calibration start values from config file.\n";
  } else {
    cout << "Using calibration start values from file ";

    // Check whether calibration constants are given in one of the
    // Kalibri formats or in official JetMET format and prepare
    // vector for further processing
    std::string inputFormat = "UNKNOWN";
    if( inputCalibration.front() == "Kalibri" ) {
      inputFormat = "Kalibri";
      inputCalibration.erase(inputCalibration.begin());
    }
    else if( inputCalibration.front() == "JetMET" ) {
      inputFormat = "JetMET";
      inputCalibration.erase(inputCalibration.begin());
    }
    else if( inputCalibration.size() == 1 ) { // For backward compatibility
      inputFormat = "Kalibri";
    }

    // Read calibration constants
    if( inputCalibration.size() == 0 ) {
      std::cerr << "\nWARNING: No file name specified to read start values from.\n"; 
      std::cerr << "         Using start values from config file.\n";
    }
    else {
      if( inputFormat == "Kalibri" ) {
	// Read predefined calibration contants from cfi
	// or txt file depending on the ending of the name
	std::string inputFileName = inputCalibration.front();
	if( !inputFileName.substr(inputFileName.rfind(".")+1).compare("cfi") ) {
	  cout << inputFileName << std::endl;
	  readCalibrationCfi(inputFileName); 
	} 
	else if( !inputFileName.substr(inputFileName.rfind(".")+1).compare("txt") ) {
	  cout << inputFileName << std::endl;
	  readCalibrationTxt(inputFileName); 
	}
	else {
	  cerr << "\nERROR: Unknown file format: '" ;
	  cerr << inputFileName.substr(inputFileName.rfind(".")) << "'\n";
	  cerr << "       Using start values from config file.\n";
	}
      }
      else if( inputFormat == "JetMET" ) {
	cout << ":\n";
	readCalibrationJetMET(inputCalibration); 
      }
      else {
	std::cerr << "\nWARNING: Unknown input format.\n"; 
	std::cerr << "         Using start values from config file.\n";
      }
    }
  }

  std::string trackEffFileName = config.read<string>("track efficiency","");
  if(!trackEffFileName.empty()){
  cout << "Reading Track Efficiency from file '" << trackEffFileName << endl;
  readTrackEffTxt(trackEffFileName);
  }


  // Specific to resolution fit
  
  // PtBinning
  ptBinEdges_ = bag_of<double>(config.read<std::string>("PtAve bin edges","-2 -1"));
  if( ptBinEdges_.front() == -2. ) {
    ptBinEdges_.at(0) = config.read<double>("Et min cut on dijet",-1.);
    ptBinEdges_.at(1) = config.read<double>("Et max cut on dijet",-1.);
  }
  ptBinCenters_ = std::vector<double>(nPtBins());
  ptTrueMin_ = std::vector<double>(nPtBins());
  ptTrueMax_ = std::vector<double>(nPtBins());
  for(unsigned int i = 0; i < nPtBins(); i++) { 
    ptBinCenters_[i] = 0.5*( ptBinEdges_[i] + ptBinEdges_[i+1] ); // Just to have a scale, use mean
    ptTrueMin_[i] = 0.4*ptBinEdges_[i];
    ptTrueMax_[i] = 1.8*ptBinEdges_[i+1];
  }
  if( nPtBins() == 1 ) {
    ptTrueMin_[0] = config.read<double>("Et genJet min",ptTrueMin_[0]);
    ptTrueMax_[0] = config.read<double>("Et genJet max",ptTrueMax_[0]);
  }

  //TODO: Should be made available via config
  parNames_ = std::vector<std::string>(numberOfParameters(),"");
}

Parameters::~Parameters() {
  delete [] k_;
  delete [] parErrors_;
  delete [] trackEff_;
  delete [] parGCorr_;
  delete [] parCov_;
  for(FunctionMap::const_iterator i = funcmap_.begin() ; 
      i != funcmap_.end() ; ++i) {
    delete i->second;
  } 
  gsl_root_fsolver_free(s_); 
  if((this == instance_) && ncalls_) {
    std::cout << "Inversion statistics for Jet::expectedEt:\n";
    std::cout << "calls: " << ncalls_ << " average number of iterations:"
	      << (double)ntries_/ncalls_ << " failures:" << (double)nfails_/ncalls_*100
	      << "% warnings:" << (double)nwarns_/ntries_*100 << "%" <<std::endl;
  }
}
 
Parameters* Parameters::clone() const {
  assert(this == instance_);
  Parameters* c = new Parameters();
  c->eta_ntwr_used_ = eta_ntwr_used_;
  c->eta_symmetry_ = eta_symmetry_;
  c->eta_granularity_ = eta_granularity_; 
  c->phi_granularity_ = phi_granularity_;
  c->eta_granularity_jet_ = eta_granularity_jet_;
  c->phi_granularity_jet_ = phi_granularity_jet_; 
  c->eta_granularity_track_ = eta_granularity_track_; 
  c->phi_granularity_track_ = phi_granularity_track_;
  // do not copy entries that are not needed...
  //std::vector<double> start_values_, jet_start_values_, track_start_values_, global_jet_start_values_;
  //std::vector<std::string> parNames_;
  int n = numberOfParameters();
  c->k_ = new double[n];
  c->parErrors_ = new double[n];
  c->parGCorr_ = new double[n];
  for(int i = 0 ; i < n ; ++i) {
    c->k_[i] = k_[i]; 
    c->parErrors_[i] = parErrors_[i];
    c->parGCorr_[i] = parGCorr_[i];
  }
  c->parCov_ = new double[numberOfCovCoeffs()];
  for(int i = 0; i < numberOfCovCoeffs(); ++i) {
    c->parCov_[i] = parCov_[i];
  }
  c->trackEff_ = new double[169];

  c->isFixedPar_ = isFixedPar_; 
  c->ptBinEdges_ = ptBinEdges_;
  c->ptBinCenters_ = ptBinCenters_;
  c->ptTrueMin_ = ptTrueMin_;
  c->ptTrueMax_ = ptTrueMax_;

  //clone function map
  for(FunctionMap::const_iterator i = funcmap_.begin() ; 
      i != funcmap_.end() ; ++i) {
    Function* f = i->second->clone();
    //    std::cout << "cloning..." << i->second << ", " << f << '\n';
    f->changeParBase(k_,c->k_);
    c->funcmap_[i->first] = f;
  }
  clones_.push_back(c);

  return c;
}

void Parameters::removeClone(Parameters* p) {
  std::vector<Parameters*>::iterator i = find(clones_.begin(), 
					      clones_.end(),p);
  if(i != clones_.end()) {
    delete *i;
    clones_.erase(i);
  }
}
// -----------------------------------------------------------------
std::string Parameters::trim(std::string const& source, char const* delims) 
{
  std::string result(source);
  std::string::size_type index = result.find_last_not_of(delims);
  if(index != std::string::npos)
    result.erase(++index);

  index = result.find_first_not_of(delims);
  if(index != std::string::npos)
    result.erase(0, index);
  else
    result.erase();
  
  //replace all "," by " "  :
  std::string::size_type  pos = result.find(",");
  while(pos != string::npos) {
    result.replace(pos,1," ");
    pos = result.find(",",pos);
  }
    
  return result;
}



//!  \brief Read predefined calibration constants from txt file 
//!
//! fills start parameters for fit when read from txt file; expects 
//! 72 lines for 72 bins in phi for each eta bin ranging from -41
//! to 41 (skipping the 0) and the following parameter format:
//! maxEta minEta nPar towerParameters jetParameters separated by
//! blanks
// ---------------------------------------------------------------
void Parameters::readCalibrationTxt(std::string const& configFile)
{
  std::ifstream file(configFile.c_str());
  std::string line; // buffer line

  int      ietaBin=-42;
  //  unsigned iLines=  0;
  while( std::getline(file,line) ){
    unsigned int iphiBin = 1; // Only 1 phibin is written to the file
    // determine phi bin on the fly
    //    phiBin=(iLines%72)+1;      // phi counts from 1...72 for each eta bin
//     ++iLines;
    // determine eta bin on the fly
    if(iphiBin==1) ++ietaBin; // increas etaValue by for the first phi bin
    if(ietaBin==0) ++ietaBin; // and don't forget to skip the 0

    //    cout << "etaBin: " << etaBin << " :: " << "phiBin: " << phiBin << endl;

    // buffers for input parameters
    unsigned nPar=0; //this is not needed but read out for control reasons 
    double etaMax=0; //this is not needed but read out for control reasons  
    double etaMin=0; //this is not needed but read out for control reasons 
    double etMin=0;
    double etMax=0;
    std::vector<double> twrPars, jetPars, trkPars,globaljetPars;
    unsigned entry=0; // controls which parameter is to filled
    while( line.length()>line.substr(0, line.find(" ")).size() ){
      if( 0<line.find(" ")){
   	// extract value
	switch(++entry){
	case 1 : etaMin = std::atof( line.substr(0, line.find(" ")).c_str() ); 
	  break;
	case 2 : etaMax = std::atof( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	case 3 : nPar   = std::atoi( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	case 4 : etMin  = std::atoi( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	case 5 : etMax  = std::atoi( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	default:
	  if((entry-5)<=p_->nTowerPars()){
	    twrPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  else if((entry-5)<=p_->nTowerPars()+p_->nJetPars()){
	    jetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  else if((entry-5)<=p_->nTowerPars()+p_->nJetPars()+p_->nTrackPars()) {
	    trkPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  } else {
	    globaljetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  break;
	}
	// cut string
	line = line.substr(line.find(" "));
      }
      else{
	//cut string
	if(line.find(" ")<std::string::npos){
	  line = line.substr(line.find(" ")+1);
	}
      }
    }
    // catch last character
    //trkPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
    if((entry-5)<=p_->nTowerPars()){
      twrPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
    }
    else if((entry-5)<=p_->nTowerPars()+p_->nJetPars()){
      jetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
    }
    else if((entry-5)<=p_->nTowerPars()+p_->nJetPars()+p_->nTrackPars()) {
      trkPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
    } else {
      globaljetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
    }

    // fill parameters
    for(iphiBin = 1; iphiBin <= 72; iphiBin++) {
      int towerIdx = bin(etaBin(ietaBin),phiBin(iphiBin));
      if( towerIdx<0 ) continue;
      for (unsigned n=0; n< p_->nTowerPars(); ++n) {
	k_[towerIdx*p_->nTowerPars()+n] = twrPars[n];
	//e[towerIdx*p_->nTowerPars()+n] = NOT_READ_OUT;
      }
      int jetIdx = jetBin(jetEtaBin(ietaBin),jetPhiBin(iphiBin));
      if( jetIdx<0 ) continue;
      for (unsigned n=0; n<p_->nJetPars(); ++n) {
	k_[numberOfTowerParameters()+jetIdx*p_->nJetPars()+n] = jetPars[n];
	//e[GetNumberOfTowerParameters()+jetIdx*p_->nJetPars()+n] = NOT_READ_OUT;
      }
      int trackIdx = trackBin(trackEtaBin(ietaBin),trackPhiBin(iphiBin));
      if( trackIdx<0 ) continue;
      for (unsigned n=0; n<p_->nTrackPars(); ++n) {
	k_[numberOfTowerParameters()+numberOfJetParameters()+trackIdx*p_->nTrackPars()+n] = trkPars[n];
	//e[numberOfTowerParameters()+jetIdx*p_->nJetPars()+n] = NOT_READ_OUT;
      }
      for (unsigned n=0; n<p_->nGlobalJetPars(); ++n) {
	k_[numberOfTowerParameters()+numberOfJetParameters()+numberOfTrackParameters() +n] = globaljetPars[n];
      }
    }
  }
}



//!  \brief Read predefined calibration constants from cfi file 
// -----------------------------------------------------------------
void Parameters::readCalibrationCfi(std::string const& configFile)
{
  std::ifstream file(configFile.c_str());

  std::string line, name;
  char * dummy = new char[28];
  std::vector<int> eta, phi;
  std::vector<double> param[p_->nTowerPars()], error[p_->nTowerPars()];
  std::vector<int> eta_jet, phi_jet;
  std::vector<double> param_jet[p_->nJetPars()], error_jet[p_->nJetPars()];
  std::vector<int> eta_track, phi_track;
  std::vector<double> param_track[p_->nTrackPars()], error_track[p_->nTrackPars()];
  std::vector<double> param_globaljet, error_globaljet;
  int posEqual;
  while (std::getline(file,line)) {
    if (! line.length()) continue;
    if( line.find("#") != string::npos) continue;
    
    //Read Tower Calibration: ---------------------------------------------------
    //if ( line.find("module ccctm = CalibratedCaloTowerMaker") != string::npos ) {
    if ( line.find("block TowerCalibConstants = {") != string::npos ) {
      while( std::getline(file,line) ) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	if( name.find("TowMapEta") != string::npos) 
  	  eta = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("TowMapPhi") != string::npos) 
  	  phi = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i < p_->nTowerPars() ; ++i) {
	  sprintf(dummy,"TowerParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"TowerError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
    //Read Jet Calibration --------------------------------------------------------
    if ( line.find("block JetCalibConstants = {") != string::npos ) {
      while (std::getline(file,line)) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	//std::cout << name << ".\n";
	if( name.find("JetMapEta") != string::npos) 
	  eta_jet = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("JetMapPhi") != string::npos) 
  	  phi_jet = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i<p_->nJetPars(); ++i) {
	  sprintf(dummy,"JetParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"JetError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
	if( name.find("GlobalJetParams") != string::npos) 
	  param_globaljet = bag_of<double>(trim(line.substr(posEqual+1)));
	if( name.find("GlobalJetErrors") != string::npos) 
	  error_globaljet = bag_of<double>(trim(line.substr(posEqual+1)));
      }
    }
    //Read Track Calibration --------------------------------------------------------
    if ( line.find("block TrackCalibConstants = {") != string::npos ) {
      while (std::getline(file,line)) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	//std::cout << name << ".\n";
	if( name.find("TrackMapEta") != string::npos) 
	  eta_track = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("TrackMapPhi") != string::npos) 
  	  phi_track = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i<p_->nTrackPars(); ++i) {
	  sprintf(dummy,"TrackParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param_track[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"TrackError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error_track[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
  }
  //check if the read calibration is ok:
  bool ok=eta.size()==phi.size();
  for (unsigned i=0; i < p_->nTowerPars(); ++i){
    ok &= eta.size()==param[i].size();
    ok &= eta.size()==error[i].size();
  }
  //fill tower parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta.size(); ++i){
      int index = bin(etaBin(eta[i]),phiBin(phi[i]));
      if (index<0) continue;
      for (unsigned n=0; n< p_->nTowerPars(); ++n) {
        k_[index*p_->nTowerPars()+n] = param[n][i];
        parErrors_[index*p_->nTowerPars()+n] = error[n][i];
      }
    }
  }

  //check if the read calibration is ok:
  ok=eta_jet.size()==phi_jet.size();
  for (unsigned i=0; i<p_->nJetPars(); ++i){
    ok &= eta_jet.size()==param_jet[i].size();
    ok &= eta_jet.size()==error_jet[i].size();
  } 
  //fill Jet parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta_jet.size(); ++i){
      int index = jetBin(jetEtaBin(eta_jet[i]),jetPhiBin(phi_jet[i]));
      if (index<0) continue;
      for (unsigned n=0; n<p_->nJetPars(); ++n) {
        k_[numberOfTowerParameters() + index*p_->nJetPars()+n] = param_jet[n][i];
        parErrors_[numberOfTowerParameters() + index*p_->nJetPars()+n] = error_jet[n][i];
      }
    }
  }

  //check if the read calibration is ok:
  ok=eta_track.size()==phi_track.size();
  for (unsigned i=0; i<p_->nTrackPars(); ++i){
    ok &= eta_track.size()==param_track[i].size();
    ok &= eta_track.size()==error_track[i].size();
  } 
  //fill Track parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta_track.size(); ++i){
      int index = trackBin(trackEtaBin(eta_track[i]),trackPhiBin(phi_track[i]));
      if (index<0) continue;
      for (unsigned n=0; n<p_->nTrackPars(); ++n) {
        k_[numberOfTowerParameters() + numberOfJetParameters () + index*p_->nTrackPars()+n] = param_track[n][i];
        parErrors_[numberOfTowerParameters() + numberOfJetParameters () + index*p_->nTrackPars()+n] = error_track[n][i];
      }
    }
  }
  ok =  (param_globaljet.size() == p_->nGlobalJetPars());
  //fill global Jet parameters and errors:  
  if (ok) {
    for (unsigned n=0; n<p_->nGlobalJetPars(); ++n) {
      k_[numberOfTowerParameters() + numberOfJetParameters () + numberOfTrackParameters()+n] = param_globaljet[n];
      parErrors_[numberOfTowerParameters() + numberOfJetParameters () + numberOfTrackParameters()+n] = error_globaljet[n];
    }
  }
  delete[] dummy;
}



//!  \brief Read correction factors in CondDB format
// -----------------------------------------------------------------
void Parameters::readCalibrationJetMET(const std::vector<std::string>& inputFileNames) {
  std::string corrL2FileName;
  std::string corrL3FileName;
  std::string corrLResFileName;
  for(size_t i = 0; i < inputFileNames.size(); i++) {
    if( inputFileNames.at(i).find("L2Relative") != std::string::npos ) 
      corrL2FileName = inputFileNames.at(i);
    else if( inputFileNames.at(i).find("L3Absolute") != std::string::npos ) 
      corrL3FileName = inputFileNames.at(i);
    else if( inputFileNames.at(i).find("L2L3Residual") != std::string::npos ) 
      corrLResFileName = inputFileNames.at(i);
  }

  if( !corrL2FileName.empty() ) {
    std::cout << "  L2: " << corrL2FileName << std::endl;
    readCalibrationJetMETL2(corrL2FileName);
  }
  if( !corrL3FileName.empty() ) {
    std::cout << "  L3: " << corrL3FileName << std::endl;
    readCalibrationJetMETL3(corrL3FileName);
  }  
  if( !corrLResFileName.empty() ) {
    std::cout << "  LRes: " << corrLResFileName << std::endl;
    readCalibrationJetMETLRes(corrLResFileName);
  }
}



//!  \brief Read L2 correction factors in CondDB format
//!
//!  Read parameters of L2 correction from
//!  txt file in CondDB format i.e.
//!  <tt>etaMin etaMax nPar EtMin EtMax Par1 Par2 Par3 Par4 Par5 Par6</tt>.
//!  If there are more than 6 L2 parameters in the file,
//!  they are ignored.
//!
//!  The pt ranges of validity are not considered.
//!
//!  In case some eta bins are missing, default parameter
//!  values 1 0 0 are assumed.
//!
//!  \note There are scaling factors for the L2 parameters
//!  in the \p L2L3JetParametrization . The parameter
//!  values read by this method are scaled accordingly.
// -----------------------------------------------------------------
void Parameters::readCalibrationJetMETL2(const std::string& inputFileName) {
  // There are scaling factors in "L2L3JetParametrization"
  std::vector<double> scale(6,1.);
  scale.at(1) = 10.;
  scale.at(2) = 100.;
  scale.at(3) = 100.;
  scale.at(4) = 100.;

  std::vector< std::vector<double> > parL2;

  std::ifstream file;
  file.open(inputFileName.c_str());

  int    etaBin = -41;
  float  etaMin = 0.;
  float  etaMax = 0.;
  int      nPar = 0;
  float     val = -1.;  // Needs float precision as etaEdge() has
                        // float precision; do we need it there?

  if( file.is_open() ) {
    while( !file.eof() && etaBin < 42 ) {
      if( etaBin == 0 ) etaBin++;          // No bin index 0
      file >> etaMin;                      // Eta min
      file >> etaMax;                      // Eta max
      val = 0.;
      file >> val;                         // Number of values following
      if( val != 0 ) {                     // Avoid reading of empty last line
	nPar = static_cast<int>(val - 2);  // Number of L2 parameters in file
	std::vector<double> par(numberOfJetParametersPerBin(),0.);  // Storage of the L2 parameters in this bin
	file >> val;                       // Et min
	file >> val;                       // Et max
	// Store L2 parameters
	for(int i = 0; i < numberOfJetParametersPerBin(); i++) {
	  file >> val;
	  par.at(i) = scale.at(i) * val;
	}
	// In case of different numbers of parameters in
	// JetMET and Kalibri L2 parametrization
	for(int i = 0; i < nPar - numberOfJetParametersPerBin(); i++) {
	  file >> val;
	}      
	// In case some eta bin is missing,
	// add default parameters
	while( etaBin < 42 && 
	       etaMin != etaLowerEdge(etaBin) && 
	       etaMax != etaUpperEdge(etaBin)    ) {
	  std::cout << "    WARNING: No parameters for eta bin " << etaBin;
	  std::cout << "; using default parameters instead.\n";

	  std::vector<double> defaultPar(numberOfJetParametersPerBin(),0.);
	  for(int i = 0; i < numberOfJetParametersPerBin(); i++) {
	    if( i == 0 ) defaultPar.at(i) = 1.;
	    else         defaultPar.at(i) = 0.;
	  }
	  parL2.push_back(defaultPar);
	  etaBin++;
	}
      
	if( etaBin < 42 ) {
	  parL2.push_back(par);
	}
	etaBin++;
      }
    }
  }
  file.close();

  // In case last eta bins are missing,
  // add default parameters  
  while( etaBin < 42 ) {
    if( etaBin == 0 ) etaBin++;          // No bin index 0

    std::cout << "    WARNING: No parameters for eta bin " << etaBin;
    std::cout << "; using default parameters instead.\n";

    std::vector<double> defaultPar(numberOfJetParametersPerBin(),0.);
    for(int i = 0; i < numberOfJetParametersPerBin(); i++) {
      if( i == 0 ) defaultPar.at(i) = 1.;
      else         defaultPar.at(i) = 0.;
    }
    parL2.push_back(defaultPar);
    etaBin++;
  }

  // Write read constants to array
  int etaIdx = 0;
  for(etaBin = -41; etaBin <= 41; etaBin++, etaIdx++) {
    if( etaBin == 0 ) etaBin++;
    for(int phiBin = 1; phiBin <= 72; phiBin++) {
      int jetIdx = jetBin(jetEtaBin(etaBin),jetPhiBin(phiBin));
      if( jetIdx<0 ) continue;
      for(int i = 0; i < numberOfJetParametersPerBin(); i++) {
	k_[numberOfTowerParameters() +
	  jetIdx*numberOfJetParametersPerBin() + i] = parL2.at(etaIdx).at(i);
      }
    }
  }
}



//!  \brief Read L3 correction factors in CondDB format
//!
//!  Read parameters of L3 correction from
//!  txt file in CondDB format i.e.
//!  <tt>etaMin etaMax nPar EtMin EtMax Par1 Par2 Par3 Par4</tt>
//!
//!  The pt ranges of validity are not considered.
// -----------------------------------------------------------------
void Parameters::readCalibrationJetMETL3(const std::string& inputFileName) {
  std::ifstream file;
  file.open(inputFileName.c_str());

  std::vector<double> parL3;   // Storage for the L3 parameters
  double val  = -1.;
  int n = 0;

  if( file.is_open() ) {
    file >> val;                       // Eta min
    file >> val;                       // Eta max
    file >> val;                       // Number of values following
    n = static_cast<int>(val - 2);     // Number of L3 parameters
    file >> val;                       // Et min
    file >> val;                       // Et max
    for(int i = 0; i < n; i++) {       // Store L3 parameters
      file >> val;
      parL3.push_back(val);
    }
  }
  file.close();

  if( n == numberOfGlobalJetParameters() ) {
    for(int i = 0; i < numberOfGlobalJetParameters(); i++) {
      k_[numberOfTowerParameters() + 
	numberOfJetParameters() + 
	numberOfTrackParameters() + i] = parL3[i];
    }
    
  } else {
    std::cerr << "ERROR: Number of read global jet parameters too small.\n";
    std::cerr << "       Using start values from config file.\n";
  }
}

//!  \brief Read Lres correction factors in CondDB format
//!
//!  Read parameters of Lres correction from
//!  txt file in CondDB format i.e.
//!  <tt>etaMin etaMax nPar EtMin EtMax Par1</tt>.
//!
//!  The pt ranges of validity are not considered.
//!
// -----------------------------------------------------------------
void Parameters::readCalibrationJetMETLRes(const std::string& inputFileName) {
  std::vector<double> parLres;

  std::ifstream file;
  file.open(inputFileName.c_str());

  for(int etaBin = -41; etaBin <= 41; ++etaBin) {
    if(etaBin == 0 ) etaBin++;
    for(int phiBin = 1; phiBin <= 72; ++phiBin) {
      int jetIdx = jetBin(jetEtaBin(etaBin),jetPhiBin(phiBin));
      if( jetIdx<0 ) continue;
      k_[numberOfTowerParameters() + jetIdx*numberOfJetParametersPerBin()] = 1.0;
    }
  }
  if(! file.is_open() ) return;
  std::string header;
  std::getline(file,header);
  while(!file.eof()) {
    float  etaMin = 0.;
    float  etaMax = 0.;
    file >> etaMin;                      // Eta min
    file >> etaMax;                      // Eta max
    float val = 0.;
    file >> val;                         // Number of values following
    if( val != 0 ) {                     // Avoid reading of empty last line
      file >> val;                       // Et min
      file >> val;                       // Et max
      file >> val;                       // cor factor
      //std::cout << etaMin << ", " << etaMax << ", " << val << '\n';
      for(int etaBin = -41; etaBin <= 41; ++etaBin) {
	if(etaBin == 0 ) etaBin++;
	//std::cout << etaLowerEdge(etaBin) << ", " << etaUpperEdge(etaBin) << ":" <<
	// ((etaLowerEdge(etaBin) >= etaMin) && (etaUpperEdge(etaBin) <= etaMax)) << '\n';
	if((etaLowerEdge(etaBin) >= etaMin) && (etaUpperEdge(etaBin) <= etaMax)) {    
	  for(int phiBin = 1; phiBin <= 72; ++phiBin) {
	    int jetIdx = jetBin(jetEtaBin(etaBin),jetPhiBin(phiBin));
	    if( jetIdx<0 ) continue;
	    k_[numberOfTowerParameters() + jetIdx*numberOfJetParametersPerBin()] = val;
	  }
	}
      }
    }
  }
  file.close();
  /*
  for(int etaBin = -41; etaBin <= 41; ++etaBin) {
    if(etaBin == 0 ) etaBin++;
    int jetIdx = jetBin(jetEtaBin(etaBin),jetPhiBin(1));
    std::cout << k_[numberOfTowerParameters() + jetIdx*numberOfJetParametersPerBin()] << "\n";
  }
  */
}


//!  \brief Read track efficiency from txt file
// -----------------------------------------------------------------
void Parameters::readTrackEffTxt(std::string const& configFile)
{
  // ---------------------------------------------------------------
  //read Track Efficiency as used in JPT Algorithm
  // ---------------------------------------------------------------
  std::ifstream file(configFile.c_str());
  std::string line; // buffer line
  int      etaBin=  -1;
  int      ptBin=  0;
  unsigned iLines=  0; 
  if(! file) cout<<configFile.c_str()<<" does not exist"<<endl;
  while( std::getline(file,line) ){
    // determine pt bin on the fly
    ptBin=(iLines%13);      // pt counts from 0..12 for each eta bin
    // determine eta bin on the fly
    if(iLines%13==0) ++etaBin; // increas etaValue by 1 each 13 lines...
    ++iLines;

    //cout << "etaBin: " << etaBin << " :: " << "ptBin: " << ptBin << endl;

    unsigned entry=0; // controls which parameter is to filled
    while( line.length()>line.substr(0, line.find(" ")).size() ){
      if( 0<line.find(" ")){
   	// extract value
	switch(++entry){
	case 1 : if( std::atof( line.substr(0, line.find(" ")).c_str() ) != etaBin) cout<<"error in Track efficiency file"<<endl;
	  break;
	case 2 : if( std::atof( line.substr(0, line.find(" ")).c_str() ) != ptBin) cout<<"error in Track efficiency file"<<endl;
	  break; 
	case 3 :  break; 
	default:
	  trackEff_[iLines-1] = std::atof( line.substr(line.find(" ")).c_str() );
	  //cout<<iLines-1<<"  "<<etaBin<<"  "<<ptBin<<"  :  "<< trackEff[iLines-1]<<endl;
	  break;
	}
	// cut string
	line = line.substr(line.find(" "));
      }
      else{
	//cut string
	if(line.find(" ")<std::string::npos){
	  line = line.substr(line.find(" ")+1);
	}
      }
    }
  }
}



// -----------------------------------------------------------------
int Parameters::etaBin(int eta_id, int etagranu, int phigranu, bool etasym) const
{  
  assert(eta_id != 0);
  assert(eta_id <= 41);
  assert(eta_id >= -41);
  //This function knows the number of wanted eta-bins and returns 
  //in which eta-bin the tower with eta-ID "eta_id" is located.
  //Case 1 bin:
  //cout << "eta="<<eta_id<<", etagranu:"<< etagranu<< ", eta_ntwr_used:"<< eta_ntwr_used<< "eta sym:"
  //     << etasym << endl;
  if(etagranu<=1) return 0;
  if(etagranu==2) return (eta_id < 0) ? 0 : 1;

  //check if tower is within wanted etarange:
  if( etasym && std::abs(eta_id)*2>(int)eta_ntwr_used_)   return -2; 
  //calculate an index:
  unsigned index=(unsigned)(41+eta_id);
  if (eta_id>0) --index;
  if (index<0 || index>81) return -3;
  //std::cout << "eta id:" << eta_id << " index:" << index << '\n';
  //Lookup tables for some binning scenarios:                      eta 1   ->|<-  eta 2                                                |<- eta 3                         
  static const unsigned ta_42[82]={ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20,21,21,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41};
  static const unsigned ta_22[82]={ 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,15,15,15,15,16,16,16,16,17,17,17,17,18,18,18,18,19,19,19,19,20,20,20,21,21};
  static const unsigned ta_10[82]={ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9};
  static const unsigned ta_6[82]= { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5};
  
  static const unsigned ts_21[41]={ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20};
  static const unsigned ts_11[41]={ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9,10,10};
  static const unsigned ts_5[41]= { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4};
  static const unsigned ts_4[41]= { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  //                                                                         15(eta 1.305)                       27    29(2.964)                           41(5.191)
  static const unsigned ts_3[41]= { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  if (!etasym) {
    if (etagranu==82) return index;
    else if (etagranu==42) return ta_42[index];
    else if (etagranu==22) return ta_22[index];
    else if (etagranu==10) return ta_10[index];
    else if (etagranu==6) return ta_6[index];
  } else {
    if (etagranu==41) return std::abs(eta_id)-1;
    else if (etagranu==21) return ts_21[std::abs(eta_id)-1];
    else if (etagranu==11) return ts_11[abs(eta_id)-1];
    else if (etagranu== 5) return ts_5[ abs(eta_id)-1];
    else if (etagranu== 4) return ts_4[ abs(eta_id)-1];
    else if (etagranu== 3) return ts_3[ abs(eta_id)-1];
  }
  
  //Default value, should never be returned!
  return -4;
}


//!  \brief Return upper or lower eta edge
//!  \param etaBin Index of eta bin,
//!         \f$ -41 \leq \texttt{etaBin} \leq 41 \f$
//!  \param lowerEdge If true, return value of lower
//!                   edge else of upper edge
// -----------------------------------------------------------------
float Parameters::etaEdge(int const etaBin, bool lowerEdge)
{
  // return eta bin - eta edge mapping
  switch(etaBin){
  case -41: return (lowerEdge ? -5.191 : -4.889); break;
  case -40: return (lowerEdge ? -4.889 : -4.716); break;
  case -39: return (lowerEdge ? -4.716 : -4.538); break;
  case -38: return (lowerEdge ? -4.538 : -4.363); break;
  case -37: return (lowerEdge ? -4.363 : -4.191); break;
  case -36: return (lowerEdge ? -4.191 : -4.013); break;
  case -35: return (lowerEdge ? -4.013 : -3.839); break;
  case -34: return (lowerEdge ? -3.839 : -3.664); break;
  case -33: return (lowerEdge ? -3.664 : -3.489); break;
  case -32: return (lowerEdge ? -3.489 : -3.314); break;
  case -31: return (lowerEdge ? -3.314 : -3.139); break;
  case -30: return (lowerEdge ? -3.139 : -2.964); break;
  case -29: return (lowerEdge ? -2.964 : -2.853); break; 
  case -28: return (lowerEdge ? -2.853 :  -2.65); break;
  case -27: return (lowerEdge ?  -2.65 :   -2.5); break;
  case -26: return (lowerEdge ?   -2.5 : -2.322); break;
  case -25: return (lowerEdge ? -2.322 : -2.172); break;
  case -24: return (lowerEdge ? -2.172 : -2.043); break;
  case -23: return (lowerEdge ? -2.043 :  -1.93); break;
  case -22: return (lowerEdge ?  -1.93 :  -1.83); break;
  case -21: return (lowerEdge ?  -1.83 :  -1.74); break;
  case -20: return (lowerEdge ?  -1.74 : -1.653); break;
  case -19: return (lowerEdge ? -1.653 : -1.566); break;
  case -18: return (lowerEdge ? -1.566 : -1.479); break;
  case -17: return (lowerEdge ? -1.479 : -1.392); break;
  case -16: return (lowerEdge ? -1.392 : -1.305); break;
  case -15: return (lowerEdge ? -1.305 : -1.218); break;
  case -14: return (lowerEdge ? -1.218 : -1.131); break;
  case -13: return (lowerEdge ? -1.131 : -1.044); break;
  case -12: return (lowerEdge ? -1.044 : -0.957); break;
  case -11: return (lowerEdge ? -0.957 : -0.879); break;
  case -10: return (lowerEdge ? -0.879 : -0.783); break;
  case  -9: return (lowerEdge ? -0.783 : -0.696); break;
  case  -8: return (lowerEdge ? -0.696 : -0.609); break;
  case  -7: return (lowerEdge ? -0.609 : -0.522); break;
  case  -6: return (lowerEdge ? -0.522 : -0.435); break;
  case  -5: return (lowerEdge ? -0.435 : -0.348); break;
  case  -4: return (lowerEdge ? -0.348 : -0.261); break;
  case  -3: return (lowerEdge ? -0.261 : -0.174); break;
  case  -2: return (lowerEdge ? -0.174 : -0.087); break;
  case  -1: return (lowerEdge ? -0.087 :      0); break;
  case  +1: return (lowerEdge ?      0 :  0.087); break;
  case  +2: return (lowerEdge ?  0.087 :  0.174); break;
  case  +3: return (lowerEdge ?  0.174 :  0.261); break;
  case  +4: return (lowerEdge ?  0.261 :  0.348); break;
  case  +5: return (lowerEdge ?  0.348 :  0.435); break;
  case  +6: return (lowerEdge ?  0.435 :  0.522); break;
  case  +7: return (lowerEdge ?  0.522 :  0.609); break;
  case  +8: return (lowerEdge ?  0.609 :  0.696); break;
  case  +9: return (lowerEdge ?  0.696 :  0.783); break;
  case +10: return (lowerEdge ?  0.783 :  0.879); break;
  case +11: return (lowerEdge ?  0.879 :  0.957); break;
  case +12: return (lowerEdge ?  0.957 :  1.044); break;
  case +13: return (lowerEdge ?  1.044 :  1.131); break;
  case +14: return (lowerEdge ?  1.131 :  1.218); break;
  case +15: return (lowerEdge ?  1.218 :  1.305); break;
  case +16: return (lowerEdge ?  1.305 :  1.392); break;
  case +17: return (lowerEdge ?  1.392 :  1.479); break;
  case +18: return (lowerEdge ?  1.479 :  1.566); break;
  case +19: return (lowerEdge ?  1.566 :  1.653); break;
  case +20: return (lowerEdge ?  1.653 :   1.74); break;
  case +21: return (lowerEdge ?   1.74 :   1.83); break;
  case +22: return (lowerEdge ?   1.83 :   1.93); break;
  case +23: return (lowerEdge ?   1.93 :  2.043); break;
  case +24: return (lowerEdge ?  2.043 :  2.172); break;
  case +25: return (lowerEdge ?  2.172 :  2.322); break;
  case +26: return (lowerEdge ?  2.322 :    2.5); break;
  case +27: return (lowerEdge ?    2.5 :   2.65); break;
  case +28: return (lowerEdge ?   2.65 :  2.853); break;
  case +29: return (lowerEdge ?  2.853 :  2.964); break;
  case +30: return (lowerEdge ?  2.964 :  3.139); break;
  case +31: return (lowerEdge ?  3.139 :  3.314); break;
  case +32: return (lowerEdge ?  3.314 :  3.489); break;
  case +33: return (lowerEdge ?  3.489 :  3.664); break;
  case +34: return (lowerEdge ?  3.664 :  3.839); break;
  case +35: return (lowerEdge ?  3.839 :  4.013); break;
  case +36: return (lowerEdge ?  4.013 :  4.191); break;
  case +37: return (lowerEdge ?  4.191 :  4.363); break;
  case +38: return (lowerEdge ?  4.363 :  4.538); break;
  case +39: return (lowerEdge ?  4.538 :  4.716); break;
  case +40: return (lowerEdge ?  4.716 :  4.889); break;
  case +41: return (lowerEdge ?  4.889 :  5.191); break;
    //something went wrong;
  default : return -1; break;
  }
}



// -----------------------------------------------------------------
int Parameters::phiBin(int phi_id, int phigranu) const
//This function knows the number of wanted phi-bins and returns 
//in which phi-bin the tower with eta-ID "phi_id" is located.
{
  assert(phi_id >0);
  assert(phi_id <= 72);
  return (phi_id-1)*phigranu/phi_ntwr_;
}



// -----------------------------------------------------------------
void Parameters::print() const
{
  std::cout  << p_->name() << " resulting in:\n "
	     << eta_granularity_ << " x " << phi_granularity_ << " tower bins with " 
	     << numberOfTowerParametersPerBin() << " free parameters each, or " 
	     << numberOfTowerParameters() << " in total, and\n "
	     << eta_granularity_jet_ << " x " << phi_granularity_jet_ << " JES bins with " 
	     << numberOfJetParametersPerBin() << " free parameters each, or " 
	     << numberOfJetParameters() << " in total \n "
	     << eta_granularity_track_ << " x " << phi_granularity_track_ << " track bins with " 
	     << numberOfTrackParametersPerBin() << " free parameters each, or " 
	     << numberOfTrackParameters() << " in total \n "
	     << "and " << numberOfGlobalJetParameters() << " global jet parameters\n"; 
}



// -----------------------------------------------------------------
void Parameters::writeCalibrationTxt(const char* name)
{
  cout << "Writing calibration to file '" << name << "'" << endl;

  ofstream file(name, ofstream::binary);
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned iphi=1; iphi<=1; ++iphi){ // No phi binning
      int towerIdx = bin(etaBin(ieta),phiBin(iphi));
      int jetIdx = jetBin(jetEtaBin(ieta),jetPhiBin(iphi));
      int trackIdx = trackBin(trackEtaBin(ieta),trackPhiBin(iphi));
      if(towerIdx<0 || jetIdx<0 || trackIdx<0) continue;
      // write: lower eta | upper eta | nparameters, for
      // each eta id of the tower and n times for n phi bins
      file << std::setw(10) << etaLowerEdge(ieta) 
	   << std::setw(10) << etaUpperEdge(ieta)  
	   << std::setw(10) << 2 + p_->nTowerPars()+p_->nJetPars()+p_->nTrackPars()+p_->nGlobalJetPars();
      // Dummy: pt range of validity
      file << std::setw(8) << std::setprecision(4) << 4;
      file << std::setw(8) << std::setprecision(4) << 2000;
      // write: each tower parameter
      for(unsigned itower=0; itower<p_->nTowerPars(); ++itower){
	file << std::setw(12) << std::setprecision(4) << k_[towerIdx*p_->nTowerPars()+itower];
      }
      // write: each jet parameter
      for(unsigned ijet=0; ijet<p_->nJetPars(); ++ijet){
	file << std::setw(12) << std::setprecision(4) << k_[numberOfTowerParameters()+jetIdx*p_->nJetPars()+ijet];
      }
      // write: each track parameter
      for(unsigned itrack=0; itrack<p_->nTrackPars(); ++itrack){
	file << std::setw(12) << std::setprecision(4) << k_[numberOfTowerParameters()+numberOfJetParameters()+trackIdx*p_->nTrackPars()+itrack];
      }
      for(unsigned igjet=0; igjet<p_->nGlobalJetPars(); ++igjet){
	file << std::setw(12) << std::setprecision(4) << k_[numberOfTowerParameters()+numberOfJetParameters()+numberOfTrackParameters()+igjet];
      }
      // complete line
      file << std::endl;
    }
  }
  file.close();
}

//-----------------------------------------------------
void Parameters::writeCalibrationTex(const char* name, const ConfigFile& config)
{
  cout << "Writing calibration to file '" << name << "'" << endl;

  // Getting fitted jet parameter values from Parameters
  std::vector<double> pJetFit(numberOfJetParameters());
  std::vector<double> pJetError(numberOfJetParameters());
  std::vector<double> pJetGCorr(numberOfJetParameters());
  for(int i = 0; i < numberOfJetParameters(); i++) {
    int jetbin = 0;
    pJetFit[i] = jetParRef(jetbin)[i];
    pJetError[i] = jetParErrorRef(jetbin)[i];
    pJetGCorr[i] = jetParGlobalCorrCoeffRef(jetbin)[i];
  }
  
  // Getting scales from config file
  std::vector<double> pJetScale = bag_of<double>(config.read<string>("jet parameter scales","")); 

  // Getting fitted global parameter values from Parameters
  std::vector<double> pGlobalJetFit(numberOfGlobalJetParameters());
  std::vector<double> pGlobalJetParError(numberOfGlobalJetParameters());
  std::vector<double> pGlobalJetParGCorr(numberOfGlobalJetParameters());
  for(int i = 0; i < numberOfGlobalJetParameters(); i++) {
    pGlobalJetFit[i] = globalJetParRef()[i];
    pGlobalJetParError[i] = globalJetParErrorRef()[i];
    pGlobalJetParGCorr[i] = globalJetParGlobalCorrCoeffRef()[i];
  }

  // Write to tex file
  ofstream outfile(name, ofstream::binary);
  if(outfile.is_open()) {
    outfile << "\\documentclass{article}\n";
    outfile << "\\usepackage{color}\n";
    outfile << "\\definecolor{gray}{rgb}{0.5,0.5,0.5}\n";

    outfile << "\\begin{document}\n";

    // Parametriaztion
    outfile << "Parametrization class: \\texttt{";
    outfile << config.read<string>("Parametrization Class","");
    outfile << "}\n";

    // BFGS parameters
    outfile << "\\begin{flushleft}\n\\begin{tabular}{lcl}\n";
    outfile << texTabularLine<double>(config,"BFGS derivative step");
    outfile << texTabularLine<int>(config,"BFGS mvec");
    outfile << texTabularLine<int>(config,"BFGS niter");
    outfile << texTabularLine<double>(config,"BFGS eps");
    outfile << texTabularLine<double>(config,"BFGS 1st wolfe parameter");
    outfile << texTabularLine<double>(config,"BFGS 2nd wolfe parameter");
    outfile << "\\end{tabular}\\end{flushleft}\n";

    // Event selection cuts
    outfile << "\\begin{flushleft}\n\\begin{tabular}{lcl}\n";
    outfile << texTabularLine<double>(config,"Et genJet min");
    outfile << texTabularLine<double>(config,"Et genJet max");
    outfile << texTabularLine<double>(config,"DeltaR cut on jet matching");
    outfile << texTabularLine<int>(config,"Et cut on jet");
    outfile << texTabularLine<int>(config,"Eta cut on jet");
    outfile << texTabularLine<double>(config,"Min had fraction");
    outfile << texTabularLine<double>(config,"Max had fraction");
    outfile << texTabularLine<double>(config,"Relative n+1 Jet Et Cut");
    outfile << "\\end{tabular}\\end{flushleft}\n";

    // Dijet integration parameters
    outfile << "\\begin{flushleft}\n\\begin{tabular}{lcl}\n";
    outfile << texTabularLine<int>(config,"DiJet integration number of iterations");
    outfile << texTabularLine<double>(config,"DiJet integration epsilon");
    outfile << texTabularLine<double>(config,"DiJet integration min");
    outfile << texTabularLine<double>(config,"DiJet integration max");
    outfile << texTabularLine<double>(config,"DiJet integration pt bin edges");
    outfile << "\\end{tabular}\\end{flushleft}\n";

    // Start and fitted parameters
    outfile << "\\begin{center}\n";
    outfile << "\\begin{tabular}{ccccc}\n";
    outfile << "\\hline\n\\hline\n";
    outfile << "Index & Scale & Start value & Fitted value & Global correlation \\\\ \n\\hline \n";
    for(unsigned int i = 0; i < jet_start_values_.size() && i < pJetFit.size(); i++) {
      if( isFixedPar(i) ) {
	outfile << "\\textcolor{gray}{$" << i << "$} & \\textcolor{gray}{$ ";
	outfile << (i < pJetScale.size() ? pJetScale[i] : 1.0) << " $} & \\textcolor{gray}{$ ";
 	outfile << jet_start_values_.at(i) << "$ } & \\textcolor{gray}{$ ";
 	outfile << pJetFit[i] << " \\pm ";
 	outfile << pJetError[i] << " $} & \\textcolor{gray}{$ ";
	outfile << pJetGCorr[i] << " $} \\\\ \n";
      } else {
	outfile << "$" << i << "$ & $";
	outfile << (i < pJetScale.size() ? pJetScale[i] : 1.0) << "$ & $";
	outfile << jet_start_values_.at(i) << "$ & $";
	outfile << pJetFit[i] << " \\pm ";
	outfile << pJetError[i] << "$ & $";
	outfile << pJetGCorr[i] << "$ \\\\ \n";
      }
    }
    for(unsigned int i = 0; i < pGlobalJetFit.size(); i++) {
      if( i == 0 ) {
	outfile << "\\hline\n";
      }
      outfile << "$" << i << "$ & $";
      outfile << 1. << "$ & $";
      outfile << global_jet_start_values_.at(i) << "$ & $";
      outfile << pGlobalJetFit[i] << " \\pm ";
      outfile << pGlobalJetParError[i] << "$ & $";
      outfile << pGlobalJetParGCorr[i] << "$ \\\\ \n";
    }
    outfile << "\\hline\n\\hline\n";
    outfile << "\\end{tabular}\n";
    outfile << "\\end{center}\n";
    outfile << "\\end{document}\n";
    outfile.close();
  }
}



// -----------------------------------------------------------------
int Parameters::trackEffBin(double pt, double eta)
{
  int bin, etabin, ptbin;
  etabin = (int)(fabs(eta)/0.2);
  if (etabin > 12) etabin = 12;
  if (pt < 5) ptbin = (int)(pt); //bin 0-4
  else{
    if(pt < 30) ptbin = (int)(5+(pt-5)/5);  //bin 5-9
    else{
      if(pt < 40) ptbin = 10;
      else{
	if(pt<50) ptbin = 11;
	else ptbin = 12;
      }
    }
  }
  bin = 13 * etabin + ptbin;
  return bin;
}

const Function& Parameters::tower_function(int etaid, int phiid) {
  int id = bin(etaBin(etaid),phiBin(phiid));
  if (id <0) { 
    std::cerr<<"WARNING: Parameters::tower_function::index = " << id << endl; 
    exit(-2);  
  }
  int parIndex = id*numberOfTowerParametersPerBin(); 
  FunctionID fid(Tower, parIndex);
  Function* f = funcmap_[fid];
  if(! f) {
    f = new Function(&Parametrization::correctedTowerEt,0,
		     towerParRef(id),parIndex,numberOfTowerParametersPerBin(),p_);
    funcmap_[fid] = f;
  }
  return *f;
}

const Function& Parameters::jet_function(int etaid, int phiid) {
  int id = jetBin(jetEtaBin(etaid),jetPhiBin(phiid));
  if (id <0) { 
    std::cerr<<"WARNING: Parameters::jet_function::index = " << id << endl; 
    exit(-2);  
  }
  int parIndex = id * numberOfJetParametersPerBin() + numberOfTowerParameters();
  FunctionID fid(Jet, parIndex);
  Function* f = funcmap_[fid];
  if(! f) {
    f = new Function(&Parametrization::correctedJetEt,
		     p_->hasInvertedCorrection() ? &Parametrization::inverseJetCorrection : 0,
		     jetParRef(id),parIndex,numberOfJetParametersPerBin(),p_);
    
    funcmap_[fid] = f;
  }
  return *f;
}

const Function& Parameters::track_function(int etaid, int phiid) {
  int id = (etaid == 0) && (phiid == 0) ? 0: trackBin(trackEtaBin(etaid),trackPhiBin(phiid));
  if (id <0) { 
    std::cerr<<"WARNING: Parameters::track_function::index = " << id << endl; 
    exit(-2);  
  }
  int parIndex = id * numberOfTrackParametersPerBin() +
    numberOfTowerParameters() + numberOfJetParameters();
  FunctionID fid(Track, parIndex);
  Function* f = funcmap_[fid];
  if(! f) {
    f = new Function(&Parametrization::expectedResponse,0,
		     trackParRef(id),parIndex,numberOfTrackParametersPerBin(),p_);
    
    funcmap_[fid] = f;
  }
  return *f;
}

const Function& Parameters::global_jet_function() {
  int parIndex = numberOfTowerParameters()+numberOfJetParameters()+numberOfTrackParameters();
  FunctionID id(Global, parIndex);
  Function* f = funcmap_[id];
  if(! f) {
    f = new Function(&Parametrization::correctedGlobalJetEt,0,
		     globalJetParRef(),parIndex,numberOfGlobalJetParameters(),
		     p_);
    
    funcmap_[id] = f;
  }
  return *f;
}

const Function& Parameters::function(const Function& f) {
  FunctionID id(f.parFunc(), f.parIndex());
  return *(funcmap_[id]);
}


const ResolutionFunction& Parameters::function(const ResolutionFunction& f) {
  FunctionID id(Resolution, f.parIndex());
  return static_cast<const ResolutionFunction&>(*(funcmap_[id]));
}


const ResolutionFunction& Parameters::resolutionFitPDF(unsigned int ptBin, int etaid, int phiid) {
  int etaPhiBin = jetBin(jetEtaBin(etaid),jetPhiBin(phiid));
  if( etaPhiBin < 0 ) { 
    std::cerr<<"WARNING: Parameters::resolutionFitPDF::index = " << etaPhiBin << endl; 
    exit(-2);  
  }
  unsigned short int firstParIdx = numberOfTowerParameters() + etaPhiBin*numberOfJetParametersPerBin() + ptBin*resParam_->nParPerPtBin();
  double *firstParRef = k_ + firstParIdx;

//    std::cout << "\nptBin " << ptBin << std::endl;
//    std::cout << "jetIdx " << firstParIdx << std::endl;
//    std::cout << "jetParRef " << firstParRef << std::endl;
  FunctionID id(Resolution,firstParIdx);
  Function* f = funcmap_[id];
  if(! f) {
    f = new ResolutionFunction(resParam_->nParPerPtBin(),ptBin,firstParIdx,firstParRef,
			       findParStatus(firstParIdx,resParam_->nParPerPtBin()),
			       &ResolutionParametrization::pdfPtMeas,
			       &ResolutionParametrization::pdfPtTrue,
			       &ResolutionParametrization::pdfResp,
			       &ResolutionParametrization::pdfDijetAsym,
			       &ResolutionParametrization::dMeasMax,
			       resParam_);
    funcmap_[id] = f;
  }
  return static_cast<const ResolutionFunction&>(*(funcmap_[id]));
}

float Parameters::toy_tower_error_parametrization(const float *x, const Measurement *xorig, float errorig)
{
  float hadet = x[0] - xorig->EMF - xorig->OutF;
  if(hadet < 0.001) hadet = 0.001;
  float hade = hadet * xorig->E / xorig->pt; 
  //std::cout << "had Et:" << hadet << " , " << "had E:" << hade << '\n';
  
  float a = 4.44;
  float b = 1.11;
  float c = 0.03;
  
  float var = a*a/hade/hade + b*b/hade + c*c;
  //truncate variance accordingly
  float truncvar = - sqrt(var) * exp(-0.5/var) * sqrt(2/M_PI) + var * TMath::Erf(1/(sqrt(2 * var)));
  return sqrt(truncvar) * hadet;
}

float Parameters::toy_jet_error_parametrization(const float *x, const Measurement *xorig, float errorig)
{
  float a = 4.44;
  float b = 1.11;
  float c = 0.03;
  
  //return sqrt(a*a/x[0]/x[0] + b*b/x[0] + c*c)*x[0];
  
  float e   = x[0] * xorig->E / xorig->pt;
  float var = a*a/e/e + b*b/e + c*c;
  //truncate variance accordingly
  float truncvar = - sqrt(var) * exp(-0.5/var) * sqrt(2/M_PI) + var * TMath::Erf(1/(sqrt(2 * var)));
  return sqrt(truncvar) * x[0];
}


//!  The line is the following: "\texttt{<fieldname>} & = & <value> \\"
//!
//!  \param config Configfile
//!  \param fieldname Name of parameter in config file
//!  \return The line for the LaTeX tabular
// --------------------------------------------------
template<class T> std::string Parameters::texTabularLine(const ConfigFile& config, const std::string& fieldname) const {
  std::stringstream line;
  line << "\\texttt{" << fieldname << "} & = & $";
  line << config.read<T>(fieldname.c_str(),-1) << "$ \\\\ \n";
  
  return line.str();
}


//! The returned vector stores the indices of the
//! submatrix elements w.r.t. full covariance matrix
//! \p parCov_ .
// --------------------------------------------------
std::vector<int> Parameters::findCovIndices(int firstPar, int nPar) const {
  // Dimension of the submatrix
  int nCov = (nPar*nPar+nPar)/2;
  std::vector<int> indices(nCov);

  int idx = ((firstPar+1)*(firstPar+1)+firstPar+1)/2 - 1;
  int rowIdx = firstPar;
  int i = 0;
  while( i < nCov ) {
    indices[i] = idx;

    int max = ((rowIdx+1)*(rowIdx+1)+rowIdx+1)/2 - 1;
    if( idx == max ) {
      rowIdx++;
      idx += firstPar;
    }
    i++;
    idx++;
  }

  int maxNCov = numberOfCovCoeffs();
  assert( indices.back() < maxNCov );

  return indices;
}


// --------------------------------------------------
std::vector<bool> Parameters::findParStatus(int firstPar, int nPar) const {
  std::vector<bool> isFixed(nPar);
  for(int i = 0; i < nPar; i++) {
    isFixed[i] = isFixedPar(firstPar+i);
  }

  return isFixed;
}


// --------------------------------------------------
bool Parameters::findBin(double x, const std::vector<double> &binEdges, unsigned int &bin) const {
  bin = 0;
  bool inRange = false;
  if( x >= binEdges.front() && x <= binEdges.back() ) {
    inRange = true;
    for(unsigned int i = 0; i < (binEdges.size()-1); ++i) {
      if( x > binEdges[i] ) bin = i;
      else break;
    }
  }

  return inRange;
}


// --------------------------------------------------
bool Parameters::findRoot(double (* f) (double x, void * params), 
			  void* params, double& x1, double& x2, double eps) 
{
  ++ncalls_;
  F_.function = f;
  F_.params = params;
  if(gsl_root_fsolver_set(s_,&F_,x1,x2)) {
    std::cout << "Warning: root not bracketed\n";
    ++nfails_;
    return false;
  }
  int status, iter = 0;
  double r;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s_);
    r = gsl_root_fsolver_root(s_);
    x1 = gsl_root_fsolver_x_lower(s_);
    x2 = gsl_root_fsolver_x_upper(s_);
    status = gsl_root_test_interval(x1,x2,0, eps);
  }
  while(status == GSL_CONTINUE && iter < 100);
  ntries_ += iter;
  if(status != GSL_SUCCESS) {
    std::cout << "inversion failed:" << x1 << ", " << x2 << std::endl;
    ++nfails_;
    return false;
  }
  x2 = r;
  return true;
} 
