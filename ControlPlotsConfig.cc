// $Id: ControlPlotsConfig.cc,v 1.36 2012/11/21 18:50:14 kirschen Exp $

#include "ControlPlotsConfig.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>

#include "TFile.h"
#include "TSystem.h"

#include "ConfigFile.h"


//! \brief Constructor
//!
//! Reads the parameters for a profile control plot
//! of name \p name from the configuration file
//! \p configFile.
// --------------------------------------------------
ControlPlotsConfig::ControlPlotsConfig(const ConfigFile *configFile, const std::string &name)
  : config_(configFile), name_(name) {
  init();
}



//! This is the name used e.g. for histogram names and
//! consists of "binVariable()+binIdx".
// --------------------------------------------------
std::string ControlPlotsConfig::binName(int binIdx) const { 
  std::string name = binVariable();
  name += toString(binIdx);

  return name;
}



//! This is used for the profile histogram title
//! and consists of "min < varTitle(binning) < max (unit)".
// --------------------------------------------------
std::string ControlPlotsConfig::binTitle(double min, double max, bool showCut) const {
  std::string title = toString(min);
  title += " #leq ";
  title += varTitle(binVariable());
  title += " < ";
  title += toString(max);

  std::string unit = unitTitle(binVariable());
  if( unit != "" ) {
    title += " " + unit;
  }

  if(showCut){
    if(cutVariable() != "") {
      title += "  and  ";
      title += toString(cutMin());
      title += " #leq ";
      title += varTitle(cutVariable());
      title += " < ";
      title += toString(cutMax());
      std::string unit = unitTitle(cutVariable());
      if( unit != "" ) {
	title += " " + unit;
      }
    }
  }
  
  return title;
}



//! This is the name used e.g. for histogram names and
//! consists of "xVariable()+xBinIdx".
// --------------------------------------------------
std::string ControlPlotsConfig::xBinName(int xBinIdx) const { 
  std::string name = xVariable();
  name += toString(xBinIdx);

  return name;
}



//! This is a label specifying the y distribution
//! in x bin \p xBinIdx and with a value of the
//! binning variable between \p binMin and \p binMax.
//! It is drawn on the histogram and consists of
//! "min < varTitle() < max (unit), xBinMin < varTitle(x) < xBinMax (unit)"
// --------------------------------------------------
std::string ControlPlotsConfig::xBinTitle(int xBinIdx, double binMin, double binMax, bool showCut) const {
  std::string title = binTitle(binMin,binMax,showCut);
  if( xBinIdx >= 0 && xBinIdx < nXBins() ) {
    title += ",  ";
    title += toString(round(xBinEdges_.at(xBinIdx)));
    title += " #leq " + varTitle(xVariable()) + " #leq ";
    title += toString(round(xBinEdges_.at(xBinIdx+1)));
    std::string unit = unitTitle(xVariable());
    if( unit != "" ) {
      title += " " + unit;
    }
  }

  return title;
}



// --------------------------------------------------
std::string ControlPlotsConfig::yProfileTitle(ProfileType type) const {
  std::string title = "";

  if( type == Mean )
    title = "<" + yTitle() + ">";
  else if( type == StandardDeviation ) {
    if((yVariable()).find("Response") != std::string::npos)
      title = "#sigma( " + yTitle() + ") / <" + yTitle() + ">";
    else {
      title = "#sigma( " + yTitle() + ")";
    }
  }
  else if( type == GaussFitMean )
    title = "GaussFit < " + yTitle() + ">";
  else if( type == GaussFitWidth ) {
    if((yVariable()).find("Response") != std::string::npos)
      title = "GaussFit #sigma( " + yTitle() + ") / <" + yTitle() + ">";
    else 
      title = "GaussFit #sigma( " + yTitle() + ")";
  }
  else if( type == Median )
    title = "Median " + yTitle();
  else if( type == IQMean )
    title = "IQ Mean " + yTitle();
  else if( type == Chi2 )
    title = "#chi^{2} / ndof " + yTitle();
  else if( type == Probability )
    title = "Probability " + yTitle();
  else if( type == Quantiles )
    title = "Quantiles " + yTitle();
  else if( type == RatioOfMeans ) {
    //title = "(1 + <" + yTitle() + ">)/(1 - <" + yTitle() + ">)";
    title = "Relative Response";
  } 
  else if( type == RatioOfGaussFitMeans ) {
    //title = "GaussFit (1 + <" + yTitle() + ">)/(1 - <" + yTitle() + ">)";
    title = "GaussFit Relative Response";
  }  
  else if( type == RatioOfIQMeans ) {
    //title = "(1 + <" + yTitle() + ">)/(1 - <" + yTitle() + ">)";
    title = "IQMean Relative Response";
  } 
  else
    std::cerr << "WARNING: Undefined ProfileType '" << type << "'\n";    

  return title;
}



// --------------------------------------------------
int ControlPlotsConfig::color(const InputTag& tag) const {
  int color = 1;
  std::map<InputTag,int>::const_iterator it = colors_.find(tag);
  if( it != colors_.end() ) color = it->second;

  return color;
}



// --------------------------------------------------
int ControlPlotsConfig::markerStyle(const InputTag& tag) const {
  int markerStyle = 7;
  std::map<InputTag,int>::const_iterator it = markerStyles_.find(tag);
  if( it != markerStyles_.end() ) markerStyle = it->second;

  return markerStyle; 
}



// --------------------------------------------------
std::string ControlPlotsConfig::legendLabel(const InputTag& tag) const {
  std::string label = "DEFAULT";
  std::map<InputTag,std::string>::const_iterator it = legendLabels_.find(tag);
  if( it != legendLabels_.end() ) label = it->second;

  return label;
}



//! Possible names \p typeName are
//! - "Uncorrected": \p CorrectionType::Uncorrected
//! - "Kalibri": \p CorrectionType::Kalibri
//! - "L2L3"   : \p CorrectionType::L2L3
//! - "L2L3Res"   : \p CorrectionType::L2L3Res
//! - "L2L3L4" : \p CorrectionType::L2L3L4
//! - "L2L3ResL4" : \p CorrectionType::L2L3ResL4
//! - "L1L2L3" : \p CorrectionType::L1L2L3
//! - "L5" : \p CorrectionType::L5
// --------------------------------------------------
ControlPlotsConfig::CorrectionType ControlPlotsConfig::correctionType(const std::string &typeName) const {
  CorrectionType type = Uncorrected;

  if( typeName == "Uncorrected" )
    type = Uncorrected;
  else if( typeName == "Kalibri" )
    type = Kalibri;
  else if( typeName == "L2L3" )
    type = L2L3;
  else if( typeName == "L2L3Res" )
    type = L2L3Res;
  else if( typeName == "L2L3L4" )
    type = L2L3L4;
  else if( typeName == "L2L3ResL4" )
    type = L2L3ResL4;
  else if( typeName == "L1L2L3" )
    type = L1L2L3;
  else if( typeName == "L5" )
    type = L5;
  else
    std::cerr << "WARNING: Undefined CorrectionType '" << typeName << "'\n";

  return type;
}



//! Possible correction types \p corrType are
//! - \p CorrectionType::Uncorrected: "Uncorrected"
//! - \p CorrectionType::Kalibri: "Kalibri"
//! - \p CorrectionType::L2L3   : "L2L3"
//! - \p CorrectionType::L2L3Res   : "L2L3res"
//! - \p CorrectionType::L2L3L4 : "L2L3L4"
//! - \p CorrectionType::L2L3ResL4 : "L2L3ResL4"
//! - \p CorrectionType::L1L2L3   : "L1L2L3"
//! - \p CorrectionType::L5   : "L5"
// --------------------------------------------------
std::string ControlPlotsConfig::correctionTypeName(CorrectionType corrType) const {
  std::string name = "corrTypeName";

  if( corrType == Uncorrected )
    name = "Uncorrected";
  else if( corrType == Kalibri )
    name = "Kalibri";
  else if( corrType == L2L3 )
    name = "L2L3"; 
  else if( corrType == L2L3Res )
    name = "L2L3res";
  else if( corrType == L2L3L4 )
    name = "L2L3L4";
  else if( corrType == L2L3ResL4 )
    name = "L2L3ResL4";
  else if( corrType == L1L2L3 )
    name = "L1L2L3";
  else if( corrType == L5 )
    name = "L5";
  else
    std::cerr << "WARNING: Undefined CorrectionType '" << corrType << "'\n";

  return name;
}

//! Possible correction types \p corrType are
//!
// --------------------------------------------------
std::string  ControlPlotsConfig::sampleName(int sample) const
{
  std::map<int,std::string>::const_iterator i = sampleNames_.find(sample);
  return i->second;
}



//! Possible names \p typeName are
//! - "Mean": \p ProfileType::Mean
//! - "StandardDeviation": \p ProfileType::StandardDeviation
//! - "GaussFitMean": \p ProfileType::GaussFitMean
//! - "GaussFitWidth": \p ProfileType::GaussFitWidth
//! - "Median": \p ProfileType::Median
//! - "IQMean": \p ProfileType::IQMean - Interquartile mean
//! - "Chi2": \p ProfileType::Chi2
//! - "Probability": \p ProfileType::Probability
//! - "Quantiles": \p ProfileType::Quantiles
// --------------------------------------------------
ControlPlotsConfig::ProfileType ControlPlotsConfig::profileType(const std::string &typeName) const {
  ProfileType type = GaussFitMean;

  if( typeName == "Mean" )
    type = Mean;
  else if( typeName == "StandardDeviation" )
    type = StandardDeviation;
  else if( typeName == "RMS" )
    type = StandardDeviation;
  else if( typeName == "GaussFitMean" )
    type = GaussFitMean;
  else if( typeName == "GaussFitWidth" )
    type = GaussFitWidth;
  else if( typeName == "Median" )
    type = Median;
  else if( typeName == "IQMean" )
    type = IQMean;
  else if( typeName == "Chi2" )
    type = Chi2;
  else if( typeName == "Probability" )
    type = Probability;
  else if( typeName == "Quantiles" )
    type = Quantiles;  
  else if( typeName == "RatioOfMeans" )
    type = RatioOfMeans;
  else if( typeName == "RatioOfGaussFitMeans" )
    type = RatioOfGaussFitMeans;
  else if( typeName == "RatioOfIQMeans" )
    type = RatioOfIQMeans;
  else
    std::cerr << "WARNING: Undefined ProfileType '" << typeName << "'\n";

  return type;
}



//! Possible profile types \p profType are
//! - \p ProfileType::Mean: "Mean" 
//! - \p ProfileType::StandardDeviation: "StandardDeviation" 
//! - \p ProfileType::GaussFitMean: "GaussFitMean" 
//! - \p ProfileType::GaussFitWidth: "GaussFitWidth" 
//! - \p ProfileType::Median: "Median" 
//! - \p ProfileType::IQMean: "IQMean" Interquartile mean 
//! - \p ProfileType::Chi2: "Chi2" 
//! - \p ProfileType::Probability: "Probability" 
//! - \p ProfileType::Quantiles: "Quantiles" 
//! - \p ProfileType::RatioOfMeans: "RatioOfMeans" 
//! - \p ProfileType::RatioOfIQMeans: "RatioOfIQMeans" 
//! - \p ProfileType::RatioOfGaussFitMeans: "RatioOfGaussFitMeans" 
// --------------------------------------------------
std::string ControlPlotsConfig::profileTypeName(ProfileType profType) const {
  std::string name = "ProfileTypeName";

  if( profType == Mean )
    name = "Mean";
  else if( profType == StandardDeviation )
    name = "StandardDeviation";
  else if( profType == GaussFitMean )
    name = "GaussFitMean";
  else if( profType == GaussFitWidth )
    name = "GaussFitWidth";
  else if( profType == Median )
    name = "Median";
  else if( profType == IQMean )
    name = "IQMean";
  else if( profType == Chi2 )
    name = "Chi2";
  else if( profType == Probability )
    name = "Probability";
  else if( profType == Quantiles )
    name = "Quantiles";  
  else if( profType == RatioOfMeans )
    name = "RatioOfMeans";
  else if( profType == RatioOfGaussFitMeans )
    name = "RatioOfGaussFitMeans";
  else if( profType == RatioOfIQMeans )
    name = "RatioOfIQMeans";
  else
    std::cerr << "WARNING: Undefined ProfileType '" << profType << "'\n";    

  return name;
}


//!  Open the file \p outFileName_ for writing
//!  
//---------------------------------------------------------------
void ControlPlotsConfig::openRootFile() {
  // Create / open ROOT file for output
  TDirectory* old_directory = gDirectory;//store old gdirectory in order not to associate new histograms to the following TFile at creation time (default ROOT-behaviour)
  outFile_ = new TFile((outDirName()+"/"+outRootFileName_).c_str()
			     ,"UPDATE","Kalibri control plots");
  std::cout << "opened ROOT-file for writing: " << (outDirName()+"/"+outRootFileName_).c_str()<< std::endl;
  old_directory->cd();//->cd(old_directory);
}

//!  Close the file \p outFileName_ for writing
//!  
//---------------------------------------------------------------
void ControlPlotsConfig::closeRootFile() {
  // Close ROOT file for output
  outFile_->Close();
  delete outFile_;
  std::cout << "closing file" << std::endl;
}


//!  Write \p obj to the file \p outFileName_ into the
//!  directory \p outDirName(). Inside the ROOT file,
//!  \p obj is written into the directory \p outFile_:/name().
//!  If it does not exist, it is created first.
//---------------------------------------------------------------
void ControlPlotsConfig::toRootFile(TObject *obj) const {
  TDirectory* old_directory = gDirectory;
  std::string directory = outFile_->GetName();
  directory += ":";
  gDirectory->cd(directory.c_str());
  directory += "/" + name();
  bool dirExists = gDirectory->GetDirectory(directory.c_str());
  if( !dirExists ) {
    gDirectory->mkdir(name().c_str());
  }
  gDirectory->cd(directory.c_str());
  if( !(gDirectory->WriteTObject(obj)) ) {
    std::cerr << "ERROR writing object '" << obj->GetName() << "' to ROOT file." << std::endl;
  }
  old_directory->cd();

}

//!  Write \p obj to the file \p outFile_ into the
//!  directory \p outDirName(). Inside the ROOT file,
//!  \p obj is written into the directory \p outFile_:/name().
//!  If it does not exist, it is created first.
//!  Can be used at any place but has low performance due to 
//!  multiple file-operations (open/close for each object)
//---------------------------------------------------------------
void ControlPlotsConfig::safelyToRootFile(TObject *obj) const {
  // Create / open ROOT file for output
  TFile *outFile = new TFile((outDirName()+"/"+outRootFileName_).c_str()
			     ,"UPDATE","Kalibri control plots");

  std::string directory = outFile->GetName();
  directory += ":";
  gDirectory->cd(directory.c_str());
  directory += "/" + name();
  bool dirExists = gDirectory->GetDirectory(directory.c_str());
  if( !dirExists ) {
    gDirectory->mkdir(name().c_str());
  }
  gDirectory->cd(directory.c_str());
  if( !(gDirectory->WriteTObject(obj)) ) {
    std::cerr << "ERROR writing object '" << obj->GetName() << "' to ROOT file." << std::endl;
  }

  outFile->Close();
  delete outFile;
}



// --------------------------------------------------
void ControlPlotsConfig::init() {
  // Read binning
  std::vector<std::string> strVar = bag_of_string(config_->read<std::string>(name_+" bin variable","Eta"));
  binVar_ = strVar.at(0);
  binEdges_ = bag_of<double>(config_->read<std::string>(name_+" bin edges","0. 1."));
  assert( binEdges_.size() > 1 );
  for(size_t i = 1; i < binEdges_.size(); i++) {
    assert( binEdges_.at(i) > binEdges_.at(i-1) );
  }
  nBins_ = static_cast<int>(binEdges_.size())-1;

  // Read x axis
  strVar = bag_of_string(config_->read<std::string>(name_+" x variable","GenJetPt"));
  xVar_ = strVar.at(0);
  logX_ = false;
  if( strVar.size() == 2 ) {
    if( strVar.at(1) == "log" ) {
      logX_ = true;
    }
  }
  bool useDirectBinEdges=false;
  double min = 20.;
  double max = 100.;
  std::vector<double> var = bag_of<double>(config_->read<std::string>(name_+" x edges","15 10 1000"));
  if( var.size() == 3 ) {
    nXBins_ = static_cast<int>(var.at(0));
    min = var.at(1);
    max = var.at(2);
  } else {
    std::cerr << "WARNING: Wrong number of arguments in config line '" << name_ << " x edges'\n";
    std::cerr << "         Trying to take arguments as bin edges, instead.\n";
    xBinEdges_ = bag_of<double>(config_->read<std::string>(name_+" x edges","0. 1."));
    assert( xBinEdges_.size() > 1 );
    for(size_t i = 1; i < xBinEdges_.size(); i++) {
      assert( xBinEdges_.at(i) > xBinEdges_.at(i-1) );
    }
    nXBins_ = static_cast<int>(xBinEdges_.size())-1;
    if(nXBins_>2)useDirectBinEdges=true;
    else nXBins_ = 5;
  }
  if(!useDirectBinEdges){
    xBinEdges_ = std::vector<double>(nXBins_+1);
    if( logX_ ) {
      if( !equidistLogBins(xBinEdges_,nXBins_,min,max) )
	std::cerr << "ERROR creating equidistant logarithmic binning.\n";
    } else {
      double width = (max - min) / nXBins_;
      for(int i = 0; i < nXBins_+1; i++) {
	xBinEdges_.at(i) = min + width*i;
      }
    }
    for(int i = 0; i < nXBins_; i++) {
      assert( xBinEdges_.at(i) < xBinEdges_.at(i+1) );
    }
  }

  cutVar_ = config_->read<std::string>(name_+" cut variable","");
  var = bag_of<double>(config_->read<std::string>(name_+" cut edges","0 0"));
  if( var.size() == 2 ) {
    cutEdges_ = std::make_pair(var[0],var[1]);
  } else {
     std::cerr << "WARNING: Wrong number of arguments in config line '" << name_ << " cut edges'\n";
     cutEdges_ = std::make_pair(0.0,0.0);
  }
  
  // Store which profile types are to be drawn
  std::vector<std::string> profTypesStr = bag_of_string(config_->read<std::string>(name_+" profile types","Uncorrected"));
  for(std::vector<std::string>::const_iterator profTypesIt = profTypesStr.begin();
      profTypesIt != profTypesStr.end(); profTypesIt++) {
    profTypes_.push_back(profileType(*profTypesIt));
  }

  // Read y axis
  strVar = bag_of_string(config_->read<std::string>(name_+" y variable","GenJetResponse"));
  yVar_ = strVar[0];
  logY_ = false;
  if( strVar.size() == 2 ) {
    if( strVar.at(1) == "log" ) {
      logY_ = true;
    }
  }

  min = 0.;
  max = 2.;
  var = bag_of<double>(config_->read<std::string>(name_+" y edges","51 0 2 0.5 1.5"));
  if( var.size() == 3 + 2 *  profTypes_.size()) {
    nYBins_ = static_cast<int>(var[0]);
    min = var[1];
    max = var[2];
    for(ProfileTypeIt pi = profTypes_.begin() ; pi != profTypes_.end() ; ++pi) { 
      unsigned int index = 2 *(pi-profTypes_.begin() +1) + 1;
      yMinZoom_[*pi] = var[index];
      yMaxZoom_[*pi] = var[index+1];
    }
  } else {
    std::cerr << "WARNING: Wrong number of arguments in config line '" << name_ << " y edges'\n";
    nYBins_ = 51;
    min = 0.;
    max = 2.;
    for(ProfileTypeIt pi = profTypes_.begin() ; pi != profTypes_.end() ; ++pi) { 
      yMinZoom_[*pi] = 0.5;
      yMaxZoom_[*pi] = 1.5;
    }
  }
  yBinEdges_ = std::vector<double>(nYBins_+1);
  if( logY_ ) {
    if( !equidistLogBins(yBinEdges_,nYBins_,min,max) )
      std::cerr << "ERROR creating equidistant logarithmic binning.\n";
  } else {
    double width = (max - min) / nYBins_;
    for(int i = 0; i < nYBins_+1; i++) {
      yBinEdges_.at(i) = min + width*i;
    }
  }
  for(int i = 0; i < nYBins_; i++) {
    assert( yBinEdges_.at(i) < yBinEdges_.at(i+1) );
    assert( yMinZoom_ < yMaxZoom_ );
  }

  // Store which inputs are to be drawn
  // in the profile plots
  std::vector<std::string> samplesStr = bag_of_string(config_->read<std::string>(name_+" input samples","0:"));
  std::vector<int> sampleIds;
  //std::cout << samplesStr.size() << '\n';
  for(std::vector<std::string>::const_iterator s = samplesStr.begin() ;
      s != samplesStr.end() ; ++s) {
    size_t pos = s->find(':');
    //std::cout << "string: " << *s << '\n';
    std::string nums,sname = "";
    if(pos == std::string::npos) {
      nums = *s;
    } else {
      nums = s->substr(0,pos);
      if(pos+2 < s->size()) {
	sname = s->substr(pos+1);
      } 
    }
    std::istringstream iss(nums);
    int id;
    iss >> id;
    //std::cout << "sample id:" << id << " name:" << sname << '\n';
    sampleIds.push_back(id);
    sampleNames_[id] = sname;
  }
  
  std::vector<std::string> corrTypesStr[3];
  std::string corstrs = config_->read<std::string>(name_+" correction types","Uncorrected");
  corrTypesStr[0] = bag_of_string(corstrs);
  corrTypesStr[1] = bag_of_string(config_->read<std::string>(name_+" 1 correction types",corstrs));
  corrTypesStr[2] = bag_of_string(config_->read<std::string>(name_+" 2 correction types",corstrs));
  
  for(std::vector<int>::const_iterator samplesIt = sampleIds.begin() ;
      samplesIt != sampleIds.end() ; samplesIt++) {
    for(std::vector<std::string>::const_iterator corrTypesIt = corrTypesStr[*samplesIt].begin(); corrTypesIt != corrTypesStr[*samplesIt].end(); corrTypesIt++) {
      inputTags_.push_back(std::make_pair(*samplesIt,correctionType(*corrTypesIt)));
    }
  }
  
  // Store which input tags are to be drawn
  // in the distributions
  std::vector<std::string> corrTypesStr2[3];
  std::string corstrs2 = config_->read<std::string>(name_+" distributions",";");
  corrTypesStr2[0] = bag_of_string(corstrs2);
  corrTypesStr2[1] = bag_of_string(config_->read<std::string>(name_+" 1 distributions",corstrs2));
  corrTypesStr2[2] = bag_of_string(config_->read<std::string>(name_+" 2 distributions",corstrs2));
  //  std::vector<std::string> corrTypesStr2 = bag_of_string(config_->read<std::string>(name_+" distributions",";")); 
  for(std::vector<int>::const_iterator samplesIt = sampleIds.begin() ;
      samplesIt != sampleIds.end() ; samplesIt++) {
    for(std::vector<std::string>::const_iterator corrTypesIt = corrTypesStr2[*samplesIt].begin(); corrTypesIt != corrTypesStr2[*samplesIt].end(); corrTypesIt++) {
      inputTagsDistributions_.push_back(std::make_pair(*samplesIt,correctionType(*corrTypesIt)));
    }
  }

  // Store directory name for output
  outDirName_ = config_->read<std::string>("plots output directory","controlPlots");
  gSystem->MakeDirectory(outDirName_.c_str()); 

  // Store whether plots are only exported to root-file
  outOnlyRoot_= config_->read<bool>("plots only to root-file",0);

  //Initialize root-file suffix for export
  outRootFileName_="";

  // Store whether XY projections of all 2D histos are exported as .eps (or other plots format)
  outXYProjections_= config_->read<bool>("export all XY projections",0);
 
  // Store whether all Y projections of all 2D-histos are created, plotted with the Gauss Fits and written out.
  outAllGaussFitsforProfiles_= config_->read<bool>("export all fitProfileHistos",0);
 
  outFileType_ = config_->read<std::string>("plots format","eps");
  // Define style for different correction types
  // This should become configurable via config file
  for(std::vector<int>::const_iterator samplesIt = sampleIds.begin() ;
      samplesIt != sampleIds.end() ; samplesIt++) {
    colors_[std::make_pair(*samplesIt,Uncorrected)] = 1;
    colors_[std::make_pair(*samplesIt,Kalibri)] = 2;
    colors_[std::make_pair(*samplesIt,L2L3)] = kBlue;
    colors_[std::make_pair(*samplesIt,L2L3Res)] = kBlue+2;
    colors_[std::make_pair(*samplesIt,L2L3L4)] = 8;
    colors_[std::make_pair(*samplesIt,L2L3ResL4)] = 1;
    colors_[std::make_pair(*samplesIt,L1L2L3)] = 2;
    colors_[std::make_pair(*samplesIt,L5)] = 8;
    if(samplesIt - sampleIds.begin() == 0) {
      markerStyles_[std::make_pair(*samplesIt,Uncorrected)] = 20;
      markerStyles_[std::make_pair(*samplesIt,Kalibri)] = 21;
      markerStyles_[std::make_pair(*samplesIt,L2L3)] = 24;
      markerStyles_[std::make_pair(*samplesIt,L2L3Res)] = 27;
      markerStyles_[std::make_pair(*samplesIt,L2L3L4)] = 28;
      markerStyles_[std::make_pair(*samplesIt,L2L3ResL4)] = 20;
      markerStyles_[std::make_pair(*samplesIt,L1L2L3)] = 21;
      markerStyles_[std::make_pair(*samplesIt,L5)] = 28;
    } else {
      int style = -(samplesIt - sampleIds.begin());
      markerStyles_[std::make_pair(*samplesIt,Uncorrected)] = style;
      markerStyles_[std::make_pair(*samplesIt,Kalibri)] = style;
      markerStyles_[std::make_pair(*samplesIt,L2L3)] = style;
      markerStyles_[std::make_pair(*samplesIt,L2L3Res)] = style;
      markerStyles_[std::make_pair(*samplesIt,L2L3L4)] = style;
      markerStyles_[std::make_pair(*samplesIt,L2L3ResL4)] = style;
      markerStyles_[std::make_pair(*samplesIt,L1L2L3)] = style;
      markerStyles_[std::make_pair(*samplesIt,L5)] = style;
    }
    // Define default legend labels for the different corrections
    std::string name = sampleName(*samplesIt);
    legendLabels_[std::make_pair(*samplesIt,Uncorrected)] = name + " Uncorrected";
    legendLabels_[std::make_pair(*samplesIt,Kalibri)] = name + " Kalibri";
    legendLabels_[std::make_pair(*samplesIt,L2L3)] = name + " L2L3";
    legendLabels_[std::make_pair(*samplesIt,L2L3Res)] = name + " L2L3res";
    legendLabels_[std::make_pair(*samplesIt,L2L3L4)] = name +" L2L3L4";
    legendLabels_[std::make_pair(*samplesIt,L2L3ResL4)] = name +" L2L3ResL4";
    legendLabels_[std::make_pair(*samplesIt,L1L2L3)] = name + " L1L2L3";
    legendLabels_[std::make_pair(*samplesIt,L5)] = name + " L5";

    // Read optional legend labels
    std::vector<std::string> legLabelStr = bag_of_string(config_->read<std::string>(name_+" legend label",";"));
    for(std::vector<std::string>::const_iterator legLabelIt = legLabelStr.begin(); legLabelIt != legLabelStr.end(); legLabelIt++) {
      size_t pos = legLabelIt->find(":");
      if( pos != std::string::npos ) {
	legendLabels_[std::make_pair(*samplesIt,correctionType(legLabelIt->substr(0,pos)))] = name + " " + legLabelIt->substr(pos+1);
      }
    }
  }
}



//!  Filling \p bins with borders of \p nBins bins between \p first
//!  and \p last that are equidistant when viewed in log scale,
//!  so \p bins must have length \p nBins+1. If \p first, \p last
//!  or \p nBins are not positive, failure is reported.
// -------------------------------------------------------------
bool ControlPlotsConfig::equidistLogBins(std::vector<double>& bins, int nBins, double first, double last) const {
  if( nBins < 1 || first <= 0. || last <= 0. || first >= last ) return false;

  bins[0]     = first;
  bins[nBins] = last;
  const double firstLog = log10(bins[0]);
  const double lastLog  = log10(bins[nBins]);
  for (int i = 1; i < nBins; ++i) {
    bins[i] = pow(10., firstLog + i*(lastLog-firstLog)/(nBins));
  }

  return true;
}



//! The title consists of "varTitle(varName) (unit)"
// --------------------------------------------------
std::string ControlPlotsConfig::axisTitle(const std::string &varName) const {
  std::string title = varTitle(varName);
  std::string unit = unitTitle(varName);
  if( unit != "" ) {
    title += " (" + unit + ")";
  }

  return title;
}



// --------------------------------------------------
std::string ControlPlotsConfig::unitTitle(const std::string &varName) const {
  std::string title = "";

  //if( varName == "GenJetPt" )
  //  title = "GeV";

  return title;
}



// --------------------------------------------------
std::string ControlPlotsConfig::varTitle(const std::string &varName) const {
  std::string title = "";
  
  if( varName == "Phi" )
    title = "#varphi";
  else if( varName == "Eta" )
    title = "#eta";
  else if( varName == "AbsEta" )
    title = "|#eta|";
  else if( varName == "JetPt" )
    title = "p_{T} [GeV]";
  else if( varName == "Pt" )
    title = "probe jet p_{T} [GeV]";
  else if( varName == "Jet2Pt" )
    title = "barrel jet p_{T} [GeV]";
  else if( varName == "JetLead2Pt" )
    title = "leading barrel jet p_{T} [GeV]";
  else if( varName == "JetLeadPt" )
    title = "leading jet p_{T} [GeV]";
  else if( varName == "ThirdJetPt" )
    title = "third jet p_{T} [GeV]";
  else if( varName == "MeanPt" )
    //   title = "barrel jet p_{T}^ [GeV]";
     title = "#bar p_{T} [GeV]";
  else if( varName == "momentPhiPhi" ) 
    title = "#sigma_{#phi#phi}";
  else if( varName == "scaledPhiPhi" ) 
    title = "#sigma_{#phi#phi} ln(p_{T}/GeV)";
  else if( varName == "scaledEtaEta" ) 
    title = "#sigma_{#eta#eta} ln(p_{T}/GeV)";
  else if( varName == "momentEtaEta" ) 
    title = "#sigma_{#eta#eta}";
  else if( varName == "meanMoment" )  
    title = "jet width";
  //    title = "(#sigma_{#phi#phi} + #sigma_{#eta#eta})/2";
  else if( varName == "Flavor" )
    title = "Flavor gluon = 0, uds = 1";
  else if ( varName == "EMF" ) 
    title = "emf";
  else if( varName == "GenJetPt" )
    title = "p^{gen}_{T} [GeV]";
  else if( varName == "GenJetResponse" )
    title = "p_{T} / p^{gen}_{T}";
  else if( varName == "Asymmetry") 
    //    title = "(p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})";
    title = "Asymmetry";
  else if( varName == "MPFResponse") 
    title = "MPF response";
  else if( varName == "ThirdJetFraction") 
    title = "p_{3}^{proj.}/#bar p_{T}";
  else if( varName == "ThirdJetFractionPlain") 
    title = "p_{3}/#bar p_{T}";
  else if( varName == "NPU")
    title = "n_{PU}^{MC}";
  else if( varName == "MCNPUVtx")
    title = "n_{PU}^{MC}";
  else if( varName == "MCNPUTruth")
    title = "n_{PU, true}^{MC}";
  else if( varName == "VtxN")
    title = "Reconstructed vertices";
  else if( varName == "PF_CH_Fraction")
    title = "Charged hadrons";
  //    title = "PF charged hadron fraction";
  else if( varName == "PF_NH_Fraction")
    title = "Neutral hadrons";
  //    title = "PF neutral hadron fraction";
  else if( varName == "PF_PH_Fraction")
    title = "Photons";
  ///    title = "PF photon fraction";
  else if( varName == "PF_EL_Fraction")
    title = "Electrons";
  //    title = "PF electron fraction";
  else if( varName == "PF_HFHad_Fraction")
    title = "HF Hadronic";
  else if( varName == "PF_HFEm_Fraction")
    title = "HF Em";
  else if( varName == "DeltaPhi")
    title = "#Delta #varphi_{1,2}";
  else if( varName == "RunNumber")
    title = "Run number";
  else if( varName == "MCTruthResponse")
    title = "MC Response";

  return title;
}


// --------------------------------------------------
template <class T> std::string ControlPlotsConfig::toString(const T& t) const {
  std::stringstream ss;
  ss << t;
  return ss.str();
}


// --------------------------------------------------
void ControlPlotsConfig::setOutRootFileName(std::string outRootFileName) {
  outRootFileName_=outRootFileName;
}

// --------------------------------------------------
void ControlPlotsConfig::setOutDirName(std::string outDirName) {
  outDirName_=outDirName;
}

// --------------------------------------------------
void ControlPlotsConfig::setOutPlotSuffix(std::string outPlotSuffix) {
  outPlotSuffix_=outPlotSuffix;
}


// --------------------------------------------------
void ControlPlotsConfig::determineOutPlotSuffix(std::string name) {
  std::string outPlotSuffix;
  if(name=="EventWeightProcessor")outPlotSuffix="_EW";
  else if(name=="DiJetEventWeighting")outPlotSuffix="_DJ";
  else if(name=="PUTruthReweighting")outPlotSuffix="_TruthPU";
  else if(name=="PU weighting")outPlotSuffix="_PU";
  else if(name=="Event binning")outPlotSuffix="_EB";
  else if(name=="DiJetEventCuts")outPlotSuffix="_CU";
  this->setOutPlotSuffix(outPlotSuffix);
}
