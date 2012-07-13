// $Id: ControlPlotsProfile.cc,v 1.25 2012/05/24 21:46:28 kirschen Exp $

#include "ControlPlotsProfile.h"

#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TROOT.h"
#include "TImage.h"
#include "TStyle.h"
#include "TGraph.h"

#include "CalibData.h"
#include "ControlPlotsFunction.h"

#include <cstdlib>

//! \param config Configuration file
//! \param function The functions which return the x, y, and binning
//!                 variables from an event (see \p ControlPlotsFunction)
// ----------------------------------------------------------------   
ControlPlotsProfile::ControlPlotsProfile(const ControlPlotsConfig *config, const ControlPlotsFunction *function)
  : config_(config), function_(function) {
  if( !function_->isInit() ) {
    std::cerr << "ERROR: Initialization of ControlPlotsProfile with non-initialized Function.\n";
  }

  // Set up Bin objects
  bins_ = std::vector<Bin*>(config_->nBins());
  for(int i = 0; i < config_->nBins(); i++) {
    bins_.at(i) = new Bin(i,config_->binEdges()->at(i),
			  config_->binEdges()->at(i+1),
			  config_);
  }
  for(size_t i = 0; i < bins_.size(); i++) {
    assert( bins_.size() > 0 && bins_.at(i)->min() < bins_.at(i)->max() );
  }
  hXSpectrum_[0] = 0;
  hXSpectrum_[1] = 0;
  hXSpectrum_[2] = 0;


  // Create histogram for x spectrum 
  for(ControlPlotsConfig::InputTagsIterator i = config_->inputTagsBegin();
      i != config_->inputTagsEnd(); ++i) {
    int id = i->first;
    if(hXSpectrum_[id]) continue;
    std::string name = config_->name() + "_" + config_->sampleName(i->first) + "_" + config_->xVariable() + "Spectrum";
    hXSpectrum_[id] = new TH1D(name.c_str(),"",
			       config_->nXBins(),&(config_->xBinEdges()->front()));
    hXSpectrum_[id]->SetMarkerStyle(20);
    hXSpectrum_[id]->SetMarkerColor(config_->color(*i));
    hXSpectrum_[id]->SetLineColor(config_->color(*i));
    hXSpectrum_[id]->SetXTitle((config_->xTitle()).c_str());
    hXSpectrum_[id]->SetYTitle("Number of events");
    hXSpectrum_[id]->Sumw2();
    //hXSpectrum_[id]->GetXaxis()->SetNdivisions(505);
  }
}



// ----------------------------------------------------------------   
ControlPlotsProfile::~ControlPlotsProfile() {
  for(std::vector<Bin*>::iterator it = bins_.begin(); 
      it != bins_.end(); it++) {
    delete *it;
  }
  bins_.clear();
  delete hXSpectrum_[0];
  delete hXSpectrum_[1];
  delete hXSpectrum_[2];
}



//! Draws the 2D histograms, the (zoomed) profile histograms, and if
//! specified the y distributions and writes them to .eps file. The
//! histograms are also stored in a .root file.
// ----------------------------------------------------------------   
void ControlPlotsProfile::draw() {
  TCanvas *c1 = new TCanvas("c1","",500,500);
  TPad *p1 = new TPad("i1", "i1", 0.84, 0.86,0.99,0.99);
  p1->SetFillStyle(4000);  
  TImage* img = TImage::Open("kalibriLogoSmall.gif");
  //img->Scale(img->GetWidth(),img->GetHeight()); 
  p1->cd();
  img->Draw("XZ");  
  TPad *p2 = new TPad("i1", "i1", 0.79, 0.71,0.94,0.84);
  p2->SetFillStyle(4000);  
  //img->Scale(img->GetWidth(),img->GetHeight()); 
  p2->cd();
  img->Draw("XZ");
  c1->cd();
  std::string fileName;

  Bool_t only_to_root = config_->outOnlyRoot();
  Bool_t export_XY_projections = config_->outXYProjections();
  Bool_t outAllGaussFitsforProfiles = config_->outAllGaussFitsforProfiles();
  //  std::cout << "only to root: " << only_to_root << std::endl; 

  TLegend *leg = bins_.front()->createLegend();

  // Draw 2D histograms
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) {
    std::vector<TH1D*> X_projections_;
    std::vector<TH1D*> Y_projections_;
    for(ControlPlotsConfig::InputTagsIterator it = config_->inputTagsBegin() ; 
	it != config_->inputTagsEnd(); ++it) {
      c1->Clear();
      TH2D *h = (*binIt)->hYvsX(*it);
      if(export_XY_projections){
	X_projections_.push_back(h->ProjectionX());
	Y_projections_.push_back(h->ProjectionY());
	X_projections_.back()->SetMarkerStyle(config_->markerStyle(*it));
	//	std::cout << "markerstyle set to: " << config_->markerStyle(*it) << " and name is: "<< (*binIt)->hist2DFileName(*it) << std::endl;
	X_projections_.back()->SetMarkerColor(config_->color(*it));
	X_projections_.back()->SetLineColor(config_->color(*it));
	Y_projections_.back()->SetMarkerStyle(config_->markerStyle(*it));
	Y_projections_.back()->SetMarkerColor(config_->color(*it));
	Y_projections_.back()->SetLineColor(config_->color(*it));

	X_projections_.back()->SetName((*binIt)->hist2DFileName(*it)+(TString)"_X");
	Y_projections_.back()->SetName((*binIt)->hist2DFileName(*it)+(TString)"_Y");

      }	
      h->Draw("COLZ");
      if( config_->logX() ) c1->SetLogx(1);
      c1->RedrawAxis();      
      p1->DrawClone();
      config_->toRootFile(h);
      fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
      fileName += (*binIt)->hist2DFileName(*it) + "." + config_->outFileType();
      if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    }
    //config_->outPlotSuffix()
    if(export_XY_projections){
      c1->Clear();
      c1->cd();
      bool firstHist = true;
      for (unsigned int X_i=0 ; X_i < X_projections_.size(); X_i++ ){
	//	cout << " " << *it;
	X_projections_.at(X_i)->SetYTitle("Number of Events");
	X_projections_.at(X_i)->GetYaxis()->SetRangeUser(0,X_projections_.at(X_i)->GetMaximum()*2);
	TH1D* h=X_projections_.at(X_i);
	if( config_->logX() ) c1->SetLogx(1);

	if( firstHist ) {
	  //	if(X_i==0){
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HIST");
	  } else {
	    h->Draw("PE1X0");
	    //	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	  firstHist = false;
	} else { 
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HISTSAME");
	  } else {
	    h->Draw("PE1X0same"); 
	    //	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	}
	//	delete h;
      }
      p1->DrawClone();
      leg->Draw("same");
      fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
      fileName += X_projections_.at(0)->GetName() + (TString)"." + config_->outFileType();
      if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());


      c1->Clear();
      c1->cd();
      c1->SetLogx(0);
      firstHist = true;
      for (unsigned int Y_i=0 ; Y_i < Y_projections_.size(); Y_i++ ){
	//	cout << " " << *it;
	Y_projections_.at(Y_i)->SetYTitle("Number of Events");
	Y_projections_.at(Y_i)->GetYaxis()->SetRangeUser(0,Y_projections_.at(Y_i)->GetMaximum()*2);
	TH1D* h=Y_projections_.at(Y_i);

	if( firstHist ) {
	  //	if(Y_i==0){
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HIST");
	  } else {
	    h->Draw("PE1X0");
	    //	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	  firstHist = false;
	} else { 
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HISTSAME");
	  } else {
	    h->Draw("PE1X0same"); 
	    //	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	}
	//	delete h;
      }
      p1->DrawClone();
      leg->Draw("same");
      fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
      fileName += Y_projections_.at(0)->GetName() + (TString)"." + config_->outFileType();
      if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    }


    //ASK HARTMUT... DO I HAVE TO DELETE THESE?
    for (unsigned int Y_i=0 ; Y_i < Y_projections_.size(); Y_i++ ){
      delete Y_projections_.at(Y_i);
      delete X_projections_.at(Y_i);
    }
    
    
  }//end bin loop
  
  

  // Draw profile histograms
  TLine *hLine = bins_.front()->createHorizontalLine(); 
  TLine *hLine2 = new TLine(config_->xMin(),1.0,config_->xMax(),1.0);
  hLine2->SetLineStyle(2);
  hLine2->SetLineColor(1);
  
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.13);
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) {
    ControlPlotsConfig::ProfileTypeIt profTypeIt = config_->profileTypesBegin();
    for(; profTypeIt != config_->profileTypesEnd(); profTypeIt++) {
      c1->Clear();
      bool firstHist = true;
      for(ControlPlotsConfig::InputTagsIterator tagsIt = config_->inputTagsBegin() ; tagsIt != config_->inputTagsEnd() ; ++tagsIt) {
	TH1D *h = (*binIt)->hXProfile(*tagsIt,*profTypeIt);
	if( firstHist ) {
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HIST");
	  } else {
	    h->Draw("PE1X0");
	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	  h->GetXaxis()->SetMoreLogLabels();
	  h->GetXaxis()->SetNoExponent();
	  firstHist = false;
	} else { 
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HISTSAME");
	  } else {
	    h->Draw("PE1X0same"); 
	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	}
	config_->toRootFile(h);
      } 
      if( *profTypeIt == ControlPlotsConfig::RatioOfMeans 
	  || *profTypeIt == ControlPlotsConfig::RatioOfGaussFitMeans || *profTypeIt == ControlPlotsConfig::RatioOfIQMeans) hLine2->Draw("same");
      else if( *profTypeIt == ControlPlotsConfig::Mean
	       || *profTypeIt == ControlPlotsConfig::GaussFitMean  ) hLine->Draw("same");
      leg->Draw("same");
      
      if( config_->logX() ) c1->SetLogx(1);
      c1->RedrawAxis();
      p2->DrawClone();
      fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
      fileName += (*binIt)->profileFileName(*profTypeIt) + "." + config_->outFileType();
      if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    }
  }  

  // Draw profile histograms (zoom on y axis)
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) {
    ControlPlotsConfig::ProfileTypeIt profTypeIt = config_->profileTypesBegin();
    for(; profTypeIt != config_->profileTypesEnd(); profTypeIt++) {
      c1->Clear();
      bool firstHist = true;
      for(ControlPlotsConfig::InputTagsIterator tagsIt = config_->inputTagsBegin() ; tagsIt != config_->inputTagsEnd(); tagsIt++) {
	TH1D *h = (*binIt)->hXProfile(*tagsIt,*profTypeIt);
	h->GetYaxis()->SetRangeUser(config_->yMinZoom(*profTypeIt),config_->yMaxZoom(*profTypeIt));
	if( firstHist ) {  
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HIST");
	  } else {
	    h->Draw("PE1X0"); 
	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	  firstHist = false;
	  h->GetXaxis()->SetMoreLogLabels();
	  h->GetXaxis()->SetNoExponent();
	} else {
	  if(h->GetMarkerStyle() < 0) {
	    h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	    h->Draw("HISTSAME");
	  } else {
	    h->Draw("PE1X0same"); 
	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	  }
	}
      }  
      if( *profTypeIt == ControlPlotsConfig::RatioOfMeans 
	  || *profTypeIt == ControlPlotsConfig::RatioOfGaussFitMeans || *profTypeIt == ControlPlotsConfig::RatioOfIQMeans) hLine2->Draw("same"); 
      else if( *profTypeIt == ControlPlotsConfig::Mean
	       || *profTypeIt == ControlPlotsConfig::GaussFitMean  ) hLine->Draw("same");
      leg->Draw("same");

      if( config_->logX() ) c1->SetLogx(1);
      c1->RedrawAxis();
      p2->DrawClone();
      fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
      fileName += (*binIt)->profileFileName(*profTypeIt) + "_zoom." + config_->outFileType();
      if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    }
  }
  delete leg;
  // Draw Mikko's profile histograms (zoom on y axis)
  for(std::vector<Bin*>::iterator binIt = bins_.begin();
      binIt != bins_.end(); binIt++) { 
    c1->Clear();
    bool firstHist = true;
    leg = new TLegend(0.3,0.85-2*0.06,0.8,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    for(ControlPlotsConfig::InputTagsIterator tagsIt = config_->inputTagsBegin() ; tagsIt != config_->inputTagsEnd(); tagsIt++) {
      if( tagsIt->second == ControlPlotsConfig::Uncorrected) continue;
      
      TH1D *hg = (*binIt)->hXProfile(*tagsIt,ControlPlotsConfig::GaussFitMean);
      TH1D *hm = (*binIt)->hXProfile(*tagsIt,ControlPlotsConfig::Mean);
      if( (!hg) || (!hm)) continue;
      hg->GetYaxis()->SetRangeUser(config_->yMinZoom(ControlPlotsConfig::GaussFitMean),config_->yMaxZoom(ControlPlotsConfig::GaussFitMean));
      if( firstHist ) {  
	if(hg->GetMarkerStyle() < 0) {
	  hg->SetYTitle("Response");
	  hg->SetLineStyle(std::abs(hg->GetMarkerStyle()));
	  hg->Draw("HIST");
	} else { 
	  hg->SetYTitle("Response");
	  hg->Draw("PE1X0"); 
	  if(TGraph* g = makeGraph(hg)) g->Draw("LSAME");
	}
	hg->GetXaxis()->SetMoreLogLabels();
	hg->GetXaxis()->SetNoExponent();
	firstHist = false;
      } else {
	if(hg->GetMarkerStyle() < 0) {
	  hg->SetLineStyle(std::abs(hg->GetMarkerStyle()));
	  hg->Draw("HISTSAME");
	} else {
	  hg->Draw("PE1X0same"); 
	  if(TGraph* g = makeGraph(hg)) g->Draw("LSAME");
	}
      }
      hm->SetLineStyle(7);
      if(hm->GetMarkerStyle() < 0) {
	hm->SetLineStyle(std::abs(hm->GetMarkerStyle()-1));
	hm->Draw("HISTSAME");
      } else {
	hm->SetMarkerStyle(hm->GetMarkerStyle()+5);
	hm->Draw("PE1X0same"); 
	if(TGraph* g = makeGraph(hm)) g->Draw("CSAME");
      } 
      leg->AddEntry(hg,(config_->legendLabel(*tagsIt)+ " Gauss").c_str(),
		    "PL");
      leg->AddEntry(hm,(config_->legendLabel(*tagsIt)+ " Mean").c_str(),
		    "PL");
    }
    leg->Draw("same");
    hLine->Draw("same");
    if( config_->logX() ) c1->SetLogx(1);
    c1->RedrawAxis();
    p2->DrawClone();
    fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
    fileName += config_->name() +"_Special_"+config_->binName((*binIt)->id())+ "_zoom." + config_->outFileType();
    if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
    delete leg;
  }  










  leg = bins_.front()->createLegend();

  std::cout << "outAllGaussFitsforProfiles: " << outAllGaussFitsforProfiles << std::endl; 
  // Draw distributions
  if(outAllGaussFitsforProfiles){
    for(std::vector<Bin*>::iterator binIt = bins_.begin();
	binIt != bins_.end(); binIt++) {
      for(int n = 0; n < (*binIt)->nDistributions(); n++) {
	c1->Clear();
	c1->cd();
	bool firstHist = true;
	
	for( ControlPlotsConfig::InputTagsIterator tagsIt = config_->distributionInputTagsBegin() ;  tagsIt != config_->distributionInputTagsEnd(); ++tagsIt) {
	  TH1D *h = (*binIt)->hYDistribution(n,*tagsIt);
	  h->SetMarkerStyle(config_->markerStyle(*tagsIt));
	  h->SetMarkerColor(config_->color(*tagsIt));
	  h->SetLineColor(config_->color(*tagsIt));

	  h->GetXaxis()->SetRange(0,-1);
	  h->GetYaxis()->SetRangeUser(0,h->GetMaximum()*2);

	  //	  h->Dump();
	  if( firstHist ) {
	    //	if(X_i==0){
	    if(h->GetMarkerStyle() < 0) {
	      h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	      h->Draw("HIST");
	      std::cout << "should draw something" << std::endl;
	    } else {
	      h->Draw("PE1X0");
	      //	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	    }
	    firstHist = false;
	  } else { 
	    if(h->GetMarkerStyle() < 0) {
	      h->SetLineStyle(std::abs(h->GetMarkerStyle()));
	      h->Draw("HISTSAME");
	    } else {
	      h->Draw("PE1X0same"); 
	      //	    if(TGraph* g = makeGraph(h)) g->Draw("LSAME");
	    }
	  }
	  config_->toRootFile(h);
	  fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
	  fileName += (*binIt)->distributionFileName(n,*tagsIt) + "_Distribution." + config_->outFileType();
	}
 	c1->SetLogx(0);
	c1->RedrawAxis();
	leg->Draw();
	p2->DrawClone();
	if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
      }
    }
  }

  // Draw x spectrum
  c1->Clear();
  c1->cd();
  hXSpectrum_[0]->Draw("PE1");
  hXSpectrum_[0]->GetXaxis()->SetMoreLogLabels();
  hXSpectrum_[0]->GetXaxis()->SetNoExponent();
  if(hXSpectrum_[1]) hXSpectrum_[1]->Draw("HISTSAME");
  if(hXSpectrum_[2]) hXSpectrum_[2]->Draw("HISTSAME");
  if( config_->logX() )c1->SetLogx(1);
  c1->SetLogy(1);
  p2->DrawClone();
  fileName = config_->outDirName() + "/" + config_->outPlotSuffix() + "_";
  fileName += hXSpectrum_[0]->GetName();
  fileName += "." + config_->outFileType();
  if(!only_to_root)c1->SaveAs(fileName.c_str(),(config_->outFileType()).c_str());
  config_->toRootFile(hXSpectrum_[0]);
  if(hXSpectrum_[1]) config_->toRootFile(hXSpectrum_[1]);
  if(hXSpectrum_[2]) config_->toRootFile(hXSpectrum_[2]);
  // Clean up
  delete p1;
  delete c1;
  delete hLine;
  delete hLine2;
}



//! This method finds the bin into which the event \p evt falls
//! from \p ControlPlotsFunction::binValue(evt). The 2D histograms
//! in this bin are filled for the different correction types as
//! y vs x with the weight w. The values of x and y are given by
//! \p ControlPlotsFunction::xValue(evt) and
//! \p ControlPlotsFunction::yValue(evt).
//! \sa Bin::fill(), findBin()
// ----------------------------------------------------------------   
void ControlPlotsProfile::fill(const Event * evt, int id) {
  double cutv = function_->cutValue(evt);
  
  if((config_->cutMax() != config_->cutMin()) && ((cutv < config_->cutMin()) || (cutv >= config_->cutMax()))) return;
  
  double x = function_->xValue(evt);
  hXSpectrum_[id]->Fill(x,evt->weight());

  int bin = findBin(evt);
  if( bin >= 0 ) {
    for(ControlPlotsConfig::InputTagsIterator it = config_->inputTagsBegin(); 
	it != config_->inputTagsEnd(); ++it) {
      if(id != it->first) continue;
      double y = function_->yValue(evt,it->second);
      if( bins_[bin]->fill(x,y,evt->weight(),*it) ) {
	std::cerr << "ERROR when filling YvsX histograms of CorrectionType '" << it->first << ", " << it->second << "'\n";
      }
    }
  }
}



//! see also \p ControlPlotsProfile::Bin::fitProfiles()
// ----------------------------------------------------------------   
void ControlPlotsProfile::fitProfiles() {
  for(std::vector<Bin*>::iterator it = bins_.begin(); it != bins_.end(); it++) {
    (*it)->fitProfiles();
  }
}



//! The value of the binning variable for this event \p evt
//! is given by \p ControlPlotsFunction::binValue(evt).
// ----------------------------------------------------------------   
int ControlPlotsProfile::findBin(const Event * evt) const {
  int bin = -1;
  double binValue = function_->binValue(evt);
  if( binValue >= config_->min() && binValue < config_->max() ) {
    for(int i = 0, nbins = bins_.size(); i < nbins ; ++i) {
      bin = i;
      if( binValue <= bins_.at(i)->max() ) break;
    }
  }
  return bin;
}

//! create Graph for Hist
// ----------------------------------------------------------------   
TGraph* ControlPlotsProfile::makeGraph(TH1* h) {
  TGraph* gr = new TGraph(h);	
  for(int i = 0 ; i < gr->GetN() ; ++i) {
    //std::cout << i << ":" << gr->GetY()[i]  << '\n';
    if(gr->GetY()[i] == 0) {
      gr->RemovePoint(i);
      --i;
      //std::cout << "after:" << gr->GetY()[i]  << '\n';
    }
  }
  if(gr->GetN() == 0) return 0;
  return gr;
}

//! \param binIdx Index of this bin
//! \param min Minimum value of the binning variable in this bin
//! \param max Maximum value of the binning variable in this bin
//! \param config Configuration parameters
// ----------------------------------------------------------------   
ControlPlotsProfile::Bin::Bin(int binIdx, double min, double max, const ControlPlotsConfig *config)
  : idx_(binIdx), min_(min), max_(max), config_(config), nCallsFitProfiles_(0) {

  // Create one 2D histogram per correction type
  for(ControlPlotsConfig::InputTagsIterator i = config_->inputTagsBegin();
      i != config_->inputTagsEnd(); ++i) {
    std::string name = config_->name();
    name += "_"+config_->yVariable()+"Vs"+config_->xVariable();
    name += "_"+config_->sampleName(i->first);
    name += "_"+config_->correctionTypeName(i->second);
    name += "_"+config_->binName(idx_);

    TH2D * h = new TH2D(name.c_str(),"",
			config_->nXBins(),&(config_->xBinEdges()->front()),
			config_->nYBins(),config_->yMin(),config_->yMax());
    h->SetTitle((config_->binTitle(min,max)).c_str());
    h->SetXTitle((config_->xTitle()).c_str());
    h->SetYTitle((config_->yTitle()).c_str());
    h->SetMarkerColor(config_->color(*i));
    hYvxX_[*i] = h;
    h->Sumw2();
  }
}



// ----------------------------------------------------------------   
ControlPlotsProfile::Bin::~Bin() {
  for(std::map<ControlPlotsConfig::InputTag,TH2D*>::iterator it = hYvxX_.begin();
      it != hYvxX_.end(); it++) {
    delete it->second;
  }
  for(std::map< ControlPlotsConfig::InputTag, std::map< ControlPlotsConfig::ProfileType, TH1D*> >::iterator corrIt = hXProfile_.begin() ; corrIt != hXProfile_.end(); corrIt++) {
    for(std::map< ControlPlotsConfig::ProfileType, TH1D*>::iterator profIt = corrIt->second.begin() ; profIt != corrIt->second.end(); profIt++) {
      delete profIt->second;
    }
    for(std::vector<TH1D*>::iterator distIt = hYDistributions_[corrIt->first].begin() ; distIt != hYDistributions_[corrIt->first].end(); distIt++) {
      delete *distIt;
    }
  }
}



// ----------------------------------------------------------------   
TH2D *ControlPlotsProfile::Bin::hYvsX(const ControlPlotsConfig::InputTag& tag) {
  TH2D * h = 0;
  std::map<ControlPlotsConfig::InputTag,TH2D*>::iterator it = hYvxX_.find(tag);
  if( it != hYvxX_.end() ) {
    h = it->second;
  }

  return h;
}



// ----------------------------------------------------------------   
TH1D *ControlPlotsProfile::Bin::hXProfile(const ControlPlotsConfig::InputTag& tag, ControlPlotsConfig::ProfileType profType) {
  TH1D * h = 0;
  std::map< ControlPlotsConfig::InputTag, std::map< ControlPlotsConfig::ProfileType, TH1D* > >::iterator corrIt = hXProfile_.find(tag);
  if( corrIt != hXProfile_.end() ) {
    std::map< ControlPlotsConfig::ProfileType, TH1D* >::iterator profIt = corrIt->second.find(profType);
    if( profIt != corrIt->second.end() ) {
      h = profIt->second;
    }
  }

  return h;
}



// ----------------------------------------------------------------   
TH1D *ControlPlotsProfile::Bin::hYDistribution(int n,const  ControlPlotsConfig::InputTag& tag) 
{
  TH1D * h = 0;
  std::map< ControlPlotsConfig::InputTag, std::vector< TH1D* > >::iterator corrIt = hYDistributions_.find(tag);
  if( corrIt != hYDistributions_.end() ) {
    if( n >=0 && n < nDistributions() ) {
      h = corrIt->second.at(n);
    }
  }

  return h;
}



// ----------------------------------------------------------------   
int ControlPlotsProfile::Bin::fill(double x, double y, double w,const ControlPlotsConfig::InputTag& tag) 
{
  int status = -1;
  std::map<ControlPlotsConfig::InputTag,TH2D*>::iterator it = hYvxX_.find(tag);
  if( it != hYvxX_.end() ) {
    it->second->Fill(x,y,w);
    status = 0;
  }
    
  return status;
}



//!  Creates the y distributions for different x bins from the
//!  2D histograms. Different quantities are calculated from 
//!  these distributions and stored in the profiles vs x:
//!  - Mean
//!  - StandardDeviation
//!  - GaussFitMean
//!  - GaussFitWidth
//!  - Median
//!  - IQMean  - Interquartile Mean
//!  - Chi2
//!  - Probability
//!  - Quantiles
//!  - RatioOfMeans
//!  - RatioOfGaussFitMeans
//!  - RatioOfIQMeans
//---------------------------------------------------------------
int ControlPlotsProfile::Bin::fitProfiles() {
  std::cout << "Fitting profiles for " << config_->name() << std::endl;
  nCallsFitProfiles_++;
  if( nCallsFitProfiles_ > 1 ) {
    std::cerr << "WARNING: FitProfiles has been called already\n";
    // Deleting profiles
    for(std::map< ControlPlotsConfig::InputTag, std::map< ControlPlotsConfig::ProfileType, TH1D*> >::iterator corrIt = hXProfile_.begin() ; corrIt != hXProfile_.end(); corrIt++) {
      for(std::map< ControlPlotsConfig::ProfileType, TH1D*>::iterator profIt = corrIt->second.begin(); profIt != corrIt->second.end(); profIt++) {
	delete profIt->second;
      }
      for( std::vector<TH1D*>::iterator distIt = hYDistributions_[corrIt->first].begin() ; distIt != hYDistributions_[corrIt->first].end(); distIt++) {
	delete *distIt;
      }
      corrIt->second.clear();
    }
    hXProfile_.clear();
    hYDistributions_.clear();
  }
  int status = 0; // No error handling implement yet

  bool isResponse = ((config_->yVariable()).find("Response") != std::string::npos); 
  // Create profile histograms
  // Loop over CorrectionTypes
  for(std::map<ControlPlotsConfig::InputTag,TH2D*>::iterator corrIt = hYvxX_.begin();
      corrIt != hYvxX_.end(); corrIt++) {

    hYDistributions_[corrIt->first] = std::vector<TH1D*>(corrIt->second->GetNbinsX());
    std::map< ControlPlotsConfig::ProfileType, TH1D*> profMap;
    // Loop over *all possible* profile quantities
    for(int i = 0; i < ControlPlotsConfig::nProfileTypes; i++) {
      ControlPlotsConfig::ProfileType profType = static_cast<ControlPlotsConfig::ProfileType>(i);
      std::string name = corrIt->second->GetName();
      name += "_"+config_->profileTypeName(profType);
      TH1D *h = new TH1D(name.c_str(),"",
			 corrIt->second->GetNbinsX(),
			 corrIt->second->GetXaxis()->GetXbins()->GetArray());
      h->GetYaxis()->SetRangeUser(config_->yMin(),config_->yMax());
      h->SetTitle((config_->binTitle(min(),max())).c_str());
      h->SetXTitle((config_->xTitle()).c_str());
      h->SetYTitle((config_->yProfileTitle(profType)).c_str());
      h->SetMarkerStyle(config_->markerStyle(corrIt->first));
      h->SetMarkerColor(config_->color(corrIt->first));
      h->SetLineColor(config_->color(corrIt->first));	    
      h->Sumw2();
      // Add histogram to map
      profMap[profType] = h;
    } // End of loop over profile quantities

    hXProfile_[corrIt->first] = profMap;
  } // End of loop over CorrectionTypes

    // Fill profile histograms
    // Loop over CorrectionTypes
  for(std::map<ControlPlotsConfig::InputTag,TH2D*>::iterator corrIt = hYvxX_.begin();
      corrIt != hYvxX_.end(); corrIt++) {
    // Temp histogram to store slice of 2D hists
    TH1D *htemp = new TH1D("htemp","",
			   corrIt->second->GetNbinsY(),
			   corrIt->second->GetYaxis()->GetXmin(),
			   corrIt->second->GetYaxis()->GetXmax());
    htemp->Sumw2();

    const int nq = 2;
    double yq[2],xq[2];
    xq[0] = 0.5;
    xq[1] = 0.90;
    double yq_IQM[2],xq_IQM[2];
    xq_IQM[0] = 0.25;
    xq_IQM[1] = 0.75;
    // Loop over x bins
    for(int xBin = 1; xBin <= corrIt->second->GetNbinsX(); xBin++) {
      htemp->Reset();
      // Copy y bin content in this x bins to htemp
      for(int yBin = 1; yBin <= corrIt->second->GetNbinsY(); yBin++) {
	htemp->SetBinContent(yBin,corrIt->second->GetCellContent(xBin,yBin));
	htemp->SetBinError(yBin,corrIt->second->GetCellError(xBin,yBin));
      }
      htemp->SetEntries(htemp->GetEffectiveEntries());
      htemp->GetXaxis()->SetRange(0,-1);

      // Store distribution
      char name[100];
      sprintf(name,"%s_Distribution_%i",corrIt->second->GetName(),xBin);
      //      std::cout<< name << std::endl;
      htemp->SetName(name);
      htemp->SetTitle((config_->xBinTitle(xBin-1,min(),max())).c_str());
      htemp->SetXTitle((config_->yTitle()).c_str());
      htemp->SetYTitle("Number of jets");
      htemp->SetLineColor(config_->color(corrIt->first));
      if(htemp->GetEntries()>100){
	double mean = htemp->GetMean(); 
	double meanerror = htemp->GetMeanError();
	double width = htemp->GetRMS();
	if(width < 0.1) width = 0.1;
	if(htemp->GetSumOfWeights() <= 0) {
	  continue; 
	} else {
	  htemp->Fit("gaus","QNI","", mean - 3 * width,mean + 3 * width);
	  TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
	  mean = f->GetParameter(1);
	  meanerror = f->GetParError(1);
	  width = f->GetParameter(2);
	  if(width < 0.05) width = 0.05;
	  if( (htemp->Fit(f,"QI","goff",mean - 1.5 * width, mean + 1.5 * width) == 0) ) { //removed N option to store GaussFit with histogram
	    mean = f->GetParameter(1);
	    meanerror = f->GetParError(1);
	    width = f->GetParameter(2);
	    
	    hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->SetBinContent(xBin,mean);
	    hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->SetBinError(xBin,meanerror);
	    if( isResponse ) {
	      hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->SetBinContent(xBin,width/mean);
	      hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->SetBinError(xBin,f->GetParError(2)/mean);
	    } else {
	      hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->SetBinContent(xBin,width);
	      hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->SetBinError(xBin,f->GetParError(2));
	    }
	    hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfGaussFitMeans]->SetBinContent(xBin,(1+mean)/(1-mean));
	    hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfGaussFitMeans]->SetBinError(xBin,2 /((1-mean)*(1-mean))*meanerror);

	  }
	  hXProfile_[corrIt->first][ControlPlotsConfig::Chi2]->SetBinContent(xBin, f->GetChisquare() / f->GetNumberFreeParameters());
	  hXProfile_[corrIt->first][ControlPlotsConfig::Chi2]->SetBinError(xBin, 0.01);
	  hXProfile_[corrIt->first][ControlPlotsConfig::Probability]->SetBinContent(xBin, f->GetProb());
	  hXProfile_[corrIt->first][ControlPlotsConfig::Probability]->SetBinError(xBin, 0.01);
	  mean = htemp->GetMean();
	  meanerror = htemp->GetMeanError();
	  width = htemp->GetRMS();
	  //	  std::cout << "Mean: " << mean << " MeanError: "<< meanerror << std::endl;
	  hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->SetBinContent(xBin,mean);
	  hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->SetBinError(xBin,meanerror);
	  if( isResponse) {
	    hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->SetBinContent(xBin,width/mean); 
	    hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->SetBinError(xBin,htemp->GetRMSError()/mean);
	  } else {
	    hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->SetBinContent(xBin,width); 
	    hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->SetBinError(xBin,htemp->GetRMSError());
	  } 
	  hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfMeans]->SetBinContent(xBin,(1+mean)/(1-mean));
	  hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfMeans]->SetBinError(xBin,2 /((1-mean)*(1-mean))*meanerror);
	  htemp->ComputeIntegral(); //needed for TH1::GetQuantiles() according to documentation
	  htemp->GetQuantiles(nq,yq,xq);
	  hXProfile_[corrIt->first][ControlPlotsConfig::Median]->SetBinContent(xBin,yq[0]);
	  hXProfile_[corrIt->first][ControlPlotsConfig::Median]->SetBinError(xBin,0.0001);
	  hXProfile_[corrIt->first][ControlPlotsConfig::Quantiles]->SetBinContent(xBin,yq[1]/yq[0]-1);
	  hXProfile_[corrIt->first][ControlPlotsConfig::Quantiles]->SetBinError(xBin,0.0001);
	  htemp->GetQuantiles(nq,yq_IQM,xq_IQM);
	  Int_t IQ_low_bin_i = htemp->FindBin(yq_IQM[0]);
	  Int_t IQ_hig_bin_i = htemp->FindBin(yq_IQM[1]);
	  htemp->GetXaxis()->SetRange(IQ_low_bin_i,IQ_hig_bin_i);
	  mean = htemp->GetMean();
	  meanerror = htemp->GetMeanError();
	  //	  std::cout << "IQMean: " << mean << " IQMeanError: "<< meanerror << std::endl;
	  hXProfile_[corrIt->first][ControlPlotsConfig::IQMean]->SetBinContent(xBin,mean);
	  hXProfile_[corrIt->first][ControlPlotsConfig::IQMean]->SetBinError(xBin,meanerror);
	  hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfIQMeans]->SetBinContent(xBin,(1+mean)/(1-mean));
	  hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfIQMeans]->SetBinError(xBin,2 /((1-mean)*(1-mean))*meanerror);
	  //	  std::cout << "Mean, Error, GaussMean, GaussError, IQMean, IQError: " << std::endl;
	  //	  std::cout << hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->GetBinContent(xBin) << ", " << hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->GetBinError(xBin) << ", " 
	  //	       << hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->GetBinContent(xBin) << ", " << hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->GetBinError(xBin) << ", " 
	  //	       << hXProfile_[corrIt->first][ControlPlotsConfig::IQMean]->GetBinContent(xBin) << ", " << hXProfile_[corrIt->first][ControlPlotsConfig::IQMean]->GetBinError(xBin) << ", " 
	  //		    <<std::endl;
	  //	  std::cout << IQ_low_bin_i << " high bin i" << IQ_hig_bin_i << std::endl;

	  // <<  "Mean, Error, GaussMean, GaussError, IQMean, IQError: " << endl;
	  delete f;
	}
      } // do fits only if more than 100 entries
    TH1D *h = static_cast<TH1D*>(htemp->Clone(name));
    hYDistributions_[corrIt->first].at(xBin-1) = h;  
    } // End of loop over x bins
    delete htemp;
    hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitMean]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::GaussFitWidth]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfGaussFitMeans]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfGaussFitMeans]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfMeans]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfMeans]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::Chi2]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::Chi2]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::Probability]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::Probability]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::Mean]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::StandardDeviation]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::Median]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::Median]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::IQMean]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::IQMean]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfIQMeans]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::RatioOfIQMeans]->GetEffectiveEntries());
    hXProfile_[corrIt->first][ControlPlotsConfig::Quantiles]->SetEntries(hXProfile_[corrIt->first][ControlPlotsConfig::Quantiles]->GetEffectiveEntries());

  } // End of loop over CorrectionTypes

  return status;
}


// ----------------------------------------------------------------   
std::string ControlPlotsProfile::Bin::hist2DFileName(const ControlPlotsConfig::InputTag& tag) const {
  std::string name = config_->name();
  name += "_2D_"+config_->binName(idx_);
  name += "_"+config_->sampleName(tag.first);
  name += "_"+config_->correctionTypeName(tag.second);
  return name;
}



// ----------------------------------------------------------------   
std::string ControlPlotsProfile::Bin::profileFileName(ControlPlotsConfig::ProfileType type) const {
  std::string name = config_->name();
  name += "_"+config_->profileTypeName(type);
  name += "_"+config_->binName(idx_);

  return name;
}



// ----------------------------------------------------------------   
std::string ControlPlotsProfile::Bin::distributionFileName(int xBin, const ControlPlotsConfig::InputTag& tag) const {
  std::string name = config_->name();
  name += "_Dist_"+config_->binName(idx_);
  name += "_"+config_->sampleName(tag.first);
  name += "_"+config_->correctionTypeName(tag.second);
  name += "_"+config_->xBinName(xBin);

  return name;
}



// ----------------------------------------------------------------   
TLine *ControlPlotsProfile::Bin::createHorizontalLine() const {
  double y = 0.;
  if( (config_->yVariable()).find("Response") != std::string::npos )
    y = 1.;
  TLine *line = new TLine(config_->xMin(),y,config_->xMax(),y);
  line->SetLineStyle(2);
  line->SetLineColor(1);

  return line;
}



// ----------------------------------------------------------------   
TLegend *ControlPlotsProfile::Bin::createLegend() {
  size_t nEntries = hXProfile_.size();
  TLegend * leg = 0;
  if( nEntries ) {
    leg = new TLegend(0.4,0.85-nEntries*0.06,0.8,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);

    ControlPlotsConfig::ProfileTypeIt profTypeIt = config_->profileTypesBegin();
    for( ControlPlotsConfig::InputTagsIterator corrTypeIt = config_->inputTagsBegin() ; corrTypeIt != config_->inputTagsEnd(); corrTypeIt++) {
      leg->AddEntry(hXProfile(*corrTypeIt,*profTypeIt),
		    (config_->legendLabel(*corrTypeIt)).c_str(),
		    "PL");
    }
  }
  
  return leg;
}
