//
// $Id: EventReader.cc,v 1.15 2010/10/12 08:37:40 stadie Exp $
//
#include "EventReader.h"

#include "ConfigFile.h"
#include "Parameters.h" 
#include "Parametrization.h"
#include "CorFactorsFactory.h"
#include "JetConstraintEvent.h"
#include "Binning.h" 
#include "TChain.h"
#include "ToyMC.h"
#include "TTree.h"

unsigned int EventReader::numberOfEventReaders_ = 0;
std::vector<JetConstraintEvent*> EventReader::constraints_;
Binning* EventReader::binning_ = 0;

EventReader::EventReader(const std::string& configfile, TParameters* param) 
  : config_(0),par_(param),corFactorsFactory_(0),cp_(new ConstParametrization())
{
  numberOfEventReaders_++;

  config_ = new ConfigFile(configfile.c_str());
  
  useTracks_ = config_->read<bool>("use Tracks",true);
  if(par_->GetNumberOfTrackParameters() < 1) useTracks_ = false;
  weightRelToNtuple_ = 1.;

  //Error Parametrization...
  //...for tracks:
  track_error_param = par_->track_error_parametrization;
  //...for tower:
  string te = config_->read<string>("tower error parametrization","standard"); 
  if (te=="standard")
    tower_error_param = par_->tower_error_parametrization;
  else if (te=="fast")
    tower_error_param = par_->fast_error_parametrization;
  else if (te=="Jans E parametrization")
    tower_error_param = par_->jans_E_tower_error_parametrization;
  else if(te=="const")
    tower_error_param = par_->const_error_parametrization;
  else if(te=="toy")
    tower_error_param = par_->toy_tower_error_parametrization;
  else if(te=="jet")
    tower_error_param = par_->jet_only_tower_error_parametrization;
  else  
    tower_error_param = par_->tower_error_parametrization;
  //...for jets:
   std::string je = config_->read<string>("jet error parametrization","standard");
  if (je=="standard")
    jet_error_param   = par_->jet_error_parametrization;
  else if (je=="fast")
    jet_error_param   = par_->fast_error_parametrization;
  else if (je=="dummy")
    jet_error_param   = par_->dummy_error_parametrization;
  else if(je=="const")
    jet_error_param = par_->const_error_parametrization;
  else if(je=="toy")
    jet_error_param   = par_->toy_jet_error_parametrization;
  else if(je=="jet et")
    jet_error_param   = par_->jet_only_jet_error_parametrization_et;
  else if(je=="jet energy")
    jet_error_param   = par_->jet_only_jet_error_parametrization_energy;
  else  
    jet_error_param   = par_->jet_error_parametrization;

  
  std::string jcn = config_->read<string>("jet correction name","");
  
  corFactorsFactory_ = CorFactorsFactory::get(jcn);
  if(jcn !="" && (! corFactorsFactory_)) {
    std::cerr << "Failed to apply correction " << jcn << std::endl;
    exit(-1);
  } 
  if(corFactorsFactory_) {
    std::cout << "Jet corrections will be overwritten with " << jcn << " from " << std::endl; 
  }
  correctToL3_ = config_->read<bool>("correct jets to L3",false);
  correctL2L3_ = config_->read<bool>("correct jets L2L3",false);
  if( correctToL3_ && correctL2L3_ ) {
    std::cerr << "WARNING: Jets are corrected twice (to L3 and L2L3).\n" << std::endl;
    exit(-9);
  }

  // Print info only once for all readers
  if( numberOfEventReaders_ == 1 ) {
    // Track usage
    if(useTracks_)  std::cout<<"Tracks are used to calibrate jets"<< std::endl;
    else std::cout<<"Only Calorimeter information is used"<< std::endl;
    // Correction of jets
    if(correctToL3_) {
      std::cout << "Jets will be corrected to Level 3 (i.e. with L1 * L2 * L3)" << std::endl;
    } else if(correctL2L3_) {
      std::cout << "Jets will be corrected with L2 * L3" << std::endl;
    } 
  }

  if(! constraints_.size() ) {
    std::vector<double> jet_constraint = bag_of<double>(config_->read<std::string>( "jet constraints",""));
    if(jet_constraint.size() % 5 == 0) {
      for(unsigned int i = 0 ; i < jet_constraint.size() ; i += 5) {
	constraints_.push_back(new JetConstraintEvent(jet_constraint[i],jet_constraint[i+1],jet_constraint[i+2],jet_constraint[i+3],jet_constraint[i+4]));
      } 
    } else if(jet_constraint.size() > 1) {
      std::cout << "wrong number of arguments for jet constraint:" << jet_constraint.size() << '\n';
    }
    for(unsigned int i = 0 ; i < constraints_.size() ; ++i) {
      const JetConstraintEvent* jce = constraints_[i];
      std::cout << "adding constraint for jets with " << jce->minEta() << " < |eta| <  " 
		<< jce->maxEta() << " and " << jce->minPt() << " < pt < " << jce->maxPt() 
		<< " with weight " << jce->weight() << "\n";
    }
  }
  if(! binning_) {
    binning_ = new Binning(config_);
  }
}

EventReader::~EventReader()
{
  delete config_;
  for(unsigned int i = 0 ; i < constraints_.size() ; ++i) {
    delete constraints_[i];
  }
  constraints_.clear();
  delete binning_;
  binning_ = 0;
  delete corFactorsFactory_;
  delete cp_;
}

TTree * EventReader::createTree(const std::string& name) const {
  std::string treeName = config_->read<string>(name+" tree","CalibTree");
  std::vector<std::string> inputFileNames = bag_of_string(config_->read<std::string>(name+" input file","input/dijet.root"));  
  int nEvts = config_->read<int>("use "+name+" events",-1);

 
  std::string fileEnding = "";
  if( inputFileNames[0].size() > 5 ) {
    fileEnding = inputFileNames[0].substr(inputFileNames[0].size()-5,inputFileNames[0].size());
    if( fileEnding != ".root" && inputFileNames.size() == 1 ) {
      std::ifstream filelist;
      filelist.open(inputFileNames[0].c_str());
      inputFileNames.clear();
      if( filelist.is_open() ) {
	std::string name = "";
	while( !filelist.eof() ) {
	  filelist >> name;
	  if( filelist.eof() ) break;
	  inputFileNames.push_back(name);
	}
	filelist.close();
      } else {
	std::cerr << "ERROR opening file '" << inputFileNames[0] << "'\n";
	exit(1);
      }
    }
  }
  
  if( inputFileNames[0] == "toy" ) { 
    // Generate Toy MC sample
    std::cout << "\n" << name << " reader: generating ToyMC events\n";
    ToyMC* mc = new ToyMC();
    mc->init(config_);
    mc->print();
    TTree* tree = new TTree(treeName.c_str(),name.c_str());
    if( name == "Di-Jet" ) {
      mc->generateDiJetTree(tree,nEvts);
    } else if( name == "Gamma-Jet") {
      mc->generatePhotonJetTree(tree,nEvts);
    } else if( name == "Top") {
      mc->generateTopTree(tree,nEvts);
    }
    delete mc;
    return tree;
  } 
  
  TChain* chain = new TChain(treeName.c_str()); 
  std::cout << "\n" << name  << " reader: opening up to " << inputFileNames.size() << " files\n";
  for(unsigned int i = 0 ;  i < inputFileNames.size() ; ++i) {
    //std::cout << i << " " << inputFileNames[i] << std::endl;
    chain->Add(inputFileNames[i].c_str());
  }  
  return chain;
}
 
int EventReader::addConstraints(std::vector<Event*>& data) {
  unsigned int n = constraints_.size();
  for(unsigned int i = 0 ; i < n ; ++i) { 
    JetConstraintEvent* jce = constraints_[i];
    std::cout << "added constraint for jets with " << jce->minEta() << " < |eta| <  " 
	      << jce->maxEta() << " and " << jce->minPt() << " < pt < " << jce->maxPt() 
	      << " with weight " << jce->weight() << " and " << jce->nJets() 
	      << " jets " << "\n";
    data.push_back(jce);
  }
  constraints_.clear();
  return n;
}
