#include "OverlayPlainProfileAnd2DPlots.h"

OverlayPlainProfileAnd2DPlots::OverlayPlainProfileAnd2DPlots(TString titleLabel,TString yAxisLabel) : titleLabel_(titleLabel),yAxisLabel_(yAxisLabel){
  yminrange_= 0;
  ymaxrange_= -1;

}

OverlayPlainProfileAnd2DPlots::~OverlayPlainProfileAnd2DPlots(){
  std::cout << "calling destructor of OverlayPlainProfileAnd2DPlots" << std::endl;
  for(int n_i=0;n_i<BasePlotExtractors_.size();n_i++)delete BasePlotExtractors_.at(n_i);

}

void OverlayPlainProfileAnd2DPlots::addPlots(TString plotsnames,TString kalibriPlotsShortName, TString specificPlotName, TString specificSampleName, Int_t BinNumber, TString LegendEntry) {
  std::cout << "bla " << std::endl;
  BasePlotExtractors_.push_back(new BasePlotExtractor(plotsnames,kalibriPlotsShortName));
  BasePlotExtractors_.back()->init();
  SpecificPlotNames_.push_back(specificPlotName);
  SpecificSampleNames_.push_back(specificSampleName);
  BinNumbers_.push_back(BinNumber);
  chdir(BasePlotExtractors_.back()->outputPathROOT());
  LegendEntries_.push_back(LegendEntry);
}

int OverlayPlainProfileAnd2DPlots::IndexOfSpecificPlot(BasePlotExtractor* BPExtractor, TString specificPlotName) {
  for(int n_i=0;n_i<BPExtractor->names_.size();n_i++){
    if(specificPlotName==BPExtractor->names_.at(n_i))return n_i;
  }
  return -1;
}

int OverlayPlainProfileAnd2DPlots::IndexOfSpecificSample(BasePlotExtractor* BPExtractor, TString specificSampleName) {
  std::cout << SpecificPlotIndices_.back() << std::endl;
  ControlPlotsConfig *pConfig = BPExtractor->configs_.at(SpecificPlotIndices_.back());
  int idx = 0;
  for(ControlPlotsConfig::InputTagsIterator i = pConfig->inputTagsBegin();
      i != pConfig->inputTagsEnd(); ++i) {
    std::cout << pConfig->sampleName(i->first) << std::endl;
    if(pConfig->sampleName(i->first)==specificSampleName)return idx;
    idx++;
  }

  return -1;
}

void OverlayPlainProfileAnd2DPlots::Init() {
  //organizational stuff to select the correct plots (i.e. match specificPlotName to entries in list defined by plotsnames)
  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
    SpecificPlotIndices_.push_back(IndexOfSpecificPlot(BasePlotExtractors_.at(e_i),SpecificPlotNames_.at(e_i)));
    //    std::cout << SpecificPlotIndices_.at(e_i) << std::endl;
    assert (SpecificPlotIndices_.back()>-1);
    SpecificSampleIndices_.push_back(IndexOfSpecificSample(BasePlotExtractors_.at(e_i),SpecificSampleNames_.at(e_i)));
    std::cout << SpecificSampleIndices_.back() << std::endl;
    assert (SpecificSampleIndices_.back()>-1);
    //    BasePlotExtractors_.at(e_i)->configs_.at(SpecificPlotIndices_.at(e_i))

    collectOneDPlots_.push_back((TH1D*)BasePlotExtractors_.at(e_i)->AllPlots_.at(SpecificPlotIndices_.at(e_i)).at(BinNumbers_.at(e_i)).at(SpecificSampleIndices_.at(e_i))->Clone());
    collectTwoDPlots_.push_back((TH2D*)BasePlotExtractors_.at(e_i)->AllTwoDPlots_.at(SpecificPlotIndices_.at(e_i)).at(BinNumbers_.at(e_i)).at(SpecificSampleIndices_.at(e_i))->Clone());
    collectOneDPlots_.back()->UseCurrentStyle();
    collectTwoDPlots_.back()->UseCurrentStyle();
  }


//  std::cout << "reading in " << nExtrapols_ << " graphs with " << nXBins_ << " nXBins and " << nBins_ << " nBins_" <<  std::endl;
//  for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
//      VecOfTGErrvec_t allCollectedGraphs;
//      for(int bin_i=0;bin_i<nBins_;bin_i++){
//	int xBin_i=0;
//	
//	std::vector <TGraphErrors*> collectGraphs;
//	for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//	  
//	  //      std::cout <<BasePlotExtractors_.at(e_i)->outputPathROOT() <<std::endl;
//	  //      int bin_i=0;
//	  
//	  TString inName =   "ResolutionPlots_ExtrapolPtThree_MCDataRatio_"+BasePlotExtractors_.at(e_i)->plotsnames_+"_"+BasePlotExtractors_.at(e_i)->configs_.at(0)->binName(bin_i)+"_"+BasePlotExtractors_.at(e_i)->configs_.at(0)->xBinName(xBin_i);
//	  
//	  collectGraphs.push_back((TGraphErrors*) util::FileOps::readTGraphErrors(BasePlotExtractors_.at(e_i)->outputPathROOT()+"/Output"+BasePlotExtractors_.at(e_i)->kalibriPlotsShortName()+".root",BasePlotExtractors_.at(e_i)->configs_.at(0)->name()+"/"+inName));
//	  
//	}
//	
//	allCollectedGraphs.push_back(collectGraphs);
//      }
//      allCollectedGraphsXBins_.push_back(allCollectedGraphs);
//    }
}

void OverlayPlainProfileAnd2DPlots::checkConsistency() {
  std::cout << "checking consistency now..." << std::endl;
  nBins_=BasePlotExtractors_.at(0)->configs_.at(0)->nBins();
  nXBins_=BasePlotExtractors_.at(0)->configs_.at(0)->nXBins();
  nExtrapols_=BasePlotExtractors_.size();

  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
    assert(BasePlotExtractors_.at(e_i)->configs_.at(0)->nBins()==nBins_);
    assert(BasePlotExtractors_.at(e_i)->configs_.at(0)->nXBins()==nXBins_);
  }

//  for(int bin_i=0;bin_i<nBins_;bin_i++){
//    assert(allCollectedGraphs_.at(bin_i).size()==nExtrapols_);
//    //     assert(BasePlotExtractors_.at(e_i)->configs_.at(0)->nBins()==nBins_);
//  }
  std::cout << "checked consistency successfully..." << std::endl;


}


// use two subsequent plots to normalize (e.g. b response to inclusive response)
// in order to allow comparisons across MC generators
void OverlayPlainProfileAnd2DPlots::plotNormalizedComboOverlays() {

  util::StyleSettings::setStyleJMEPaperNoTitle();

  if(chdir("OverlayExtrapolPlots") != 0){ 
    mkdir("OverlayExtrapolPlots", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("OverlayExtrapolPlots"); 
  } 

  DefaultStyles style;

  TCanvas* c2 = new TCanvas("c2","",600,600);
  if( BasePlotExtractors_.at(0)->configs_.at(0)->logX() ) c2->SetLogx(1);
  if( BasePlotExtractors_.at(0)->configs_.at(0)->logY() ) c2->SetLogy(1);

  std::cout << "nExtrapols_ " << nExtrapols_ << " nExtrapols_%2: " << nExtrapols_%2 << " !nExtrapols_%2: " << !nExtrapols_%2<< std::endl;

  assert((!(nExtrapols_%2)));




//DEBUG snippet to plot all plots used for calculating the ratios
//  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++)    std::cout <<collectOneDPlots_.at(e_i)->GetName() << std::endl;
//  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//    collectOneDPlots_.at(e_i)->SetLineColor(style.getColor(e_i));
//    collectOneDPlots_.at(e_i)->SetMarkerColor(style.getColor(e_i));
//    collectOneDPlots_.at(e_i)->SetFillColor(style.getColor(e_i));
//    collectOneDPlots_.at(e_i)->SetFillStyle(1001);
//    collectOneDPlots_.at(e_i)->SetMarkerStyle(style.getMarker(e_i));
//    collectOneDPlots_.at(e_i)->GetListOfFunctions()->Delete();
//    collectOneDPlots_.at(e_i)->SetTitle("");
//  }
//
//  TLegend* leg2 = util::LabelFactory::createLegendWithOffset(nExtrapols_/2,util::LabelFactory::lineHeight());
//  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++)leg2->AddEntry(collectOneDPlots_.at(e_i),LegendEntries_.at(e_i),"PL");
//
//  collectOneDPlots_.at(0)->GetYaxis()->SetRangeUser(0.95,1.05);
//  collectOneDPlots_.at(0)->Draw();
//  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++)collectOneDPlots_.at(e_i)->Draw("same");
//
//
//  leg2->Draw();
//  BasePlotExtractors_.at(0)->drawCMSPrel();
//  TString outname = "OneDNormalizedComboRatiosDEBUGGINGPLOT_"+ util::cleanStringMore(titleLabel_)+".eps";
//  c2->RedrawAxis();
//  c2->SaveAs(outname);




  std::vector <TH1*> collectNormalizedComboRatios;
  for(unsigned int e_i=0;e_i<nExtrapols_/2;e_i++){
    if(DEBUG){
    std::cout << "Debug e_i: " << e_i << std::endl;
    std::cout <<collectOneDPlots_.at(2*e_i)->GetName() << std::endl;
    std::cout <<collectOneDPlots_.at(2*e_i+1)->GetName() << std::endl;
    }
    collectNormalizedComboRatios.push_back(util::HistOps::createRatioPlot(collectOneDPlots_.at(2*e_i+1),collectOneDPlots_.at(2*e_i),"Ratio to (i*2)+1th plot"));
//    for(int bin_i = 1; bin_i <= collectRatios_.at(e_i-1)->GetNbinsX(); bin_i++) {
//      collectRatios_.at(e_i-1)->SetBinContent(bin_i,(collectRatios_.at(e_i-1)->GetBinContent(bin_i)-1)*100);
//      collectRatios_.at(e_i-1)->SetBinError(bin_i,(collectRatios_.at(e_i-1)->GetBinError(bin_i))*100);
//    }

  }
 
 

  for(unsigned int e_i=0;e_i<nExtrapols_/2;e_i++){
    collectNormalizedComboRatios.at(e_i)->SetLineColor(style.getColor(e_i));
    collectNormalizedComboRatios.at(e_i)->SetMarkerColor(style.getColor(e_i));
    collectNormalizedComboRatios.at(e_i)->SetFillColor(style.getColor(e_i));
    collectNormalizedComboRatios.at(e_i)->SetFillStyle(1001);
    collectNormalizedComboRatios.at(e_i)->SetMarkerStyle(style.getMarker(e_i));
    collectNormalizedComboRatios.at(e_i)->GetListOfFunctions()->Delete();
    collectNormalizedComboRatios.at(e_i)->SetTitle("");
  }



  if(yAxisLabel_!="")collectNormalizedComboRatios.at(0)->GetYaxis()->SetTitle(yAxisLabel_);
  if(ymaxrange_!=-1)  collectNormalizedComboRatios.at(0)->GetYaxis()->SetRangeUser(yminrange_,ymaxrange_);
  else  collectNormalizedComboRatios.at(0)->GetYaxis()->SetRangeUser(BasePlotExtractors_.at(0)->yProfileMinMax().at(0),BasePlotExtractors_.at(0)->yProfileMinMax().at(1));
  collectNormalizedComboRatios.at(0)->Draw();
  
  for(unsigned int e_i=0;e_i<nExtrapols_/2;e_i++){
    std::cout << "plotting " << collectNormalizedComboRatios.at(e_i)->GetName() << std::endl;
    collectNormalizedComboRatios.at(e_i)->Draw("same");
  }
  TPaveText *label = util::LabelFactory::createPaveTextWithOffset(1,1.0,0.0);
  label->AddText(titleLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->Draw("same");

  TLegend* leg1 = util::LabelFactory::createLegendWithOffset(nExtrapols_/2,util::LabelFactory::lineHeight());
  for(unsigned int e_i=0;e_i<nExtrapols_/2;e_i++)leg1->AddEntry(collectNormalizedComboRatios.at(e_i),LegendEntries_.at(e_i*2),"PL");

  leg1->Draw();
  BasePlotExtractors_.at(0)->drawCMSPrel();
  TString outname = "OneDNormalizedComboRatios_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);


  chdir(".."); 
  

}



// use four subsequent plots to normalize (e.g. b response to inclusive response of Herwig to b response to inclusive response of PYTHIA)
//1: inclusive PYTHIA
//2: b PYTHIA
//3: inclusive Herwig
//4: b Herwig
//labelname is taken from the first entry (inclusive Herwig in this example)
void OverlayPlainProfileAnd2DPlots::plotNormalizedDoubleComboOverlays() {

  util::StyleSettings::setStyleJMEPaperNoTitle();

  if(chdir("OverlayExtrapolPlots") != 0){ 
    mkdir("OverlayExtrapolPlots", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("OverlayExtrapolPlots"); 
  } 

  DefaultStyles style;

  TCanvas* c2 = new TCanvas("c2","",600,600);
  if( BasePlotExtractors_.at(0)->configs_.at(0)->logX() ) c2->SetLogx(1);

  std::cout << "nExtrapols_ " << nExtrapols_ << " nExtrapols_%4: " << nExtrapols_%4 << " !nExtrapols_%4: " << !nExtrapols_%4<< std::endl;

  assert((!(nExtrapols_%4)));

  unsigned int nDoubleRatios = nExtrapols_/4;

  std::vector <TH1*> collectNormalizedDoubleComboRatios;
  for(unsigned int e_i=0;e_i<nDoubleRatios;e_i++){
    if(DEBUG){
    std::cout << "Debug e_i: " << e_i << std::endl;
    for(unsigned int c_i=0;c_i<4;c_i++)std::cout <<collectOneDPlots_.at(4*e_i+c_i)->GetName() << std::endl;
    }
    TH1* firstRatio = util::HistOps::createRatioPlot(collectOneDPlots_.at(4*e_i+1),collectOneDPlots_.at(4*e_i),"Ratio to (i*4)+1st plot");
    TH1* secondRatio = util::HistOps::createRatioPlot(collectOneDPlots_.at(4*e_i+3),collectOneDPlots_.at(4*e_i+2),"Ratio to (i*4)+1st plot");
    collectNormalizedDoubleComboRatios.push_back(util::HistOps::createRatioPlot(secondRatio,firstRatio,"Double ratio"));
    delete firstRatio; delete secondRatio;
  }
 
  for(unsigned int e_i=0;e_i<nDoubleRatios;e_i++){
    collectNormalizedDoubleComboRatios.at(e_i)->SetLineColor(style.getColor(e_i));
    collectNormalizedDoubleComboRatios.at(e_i)->SetMarkerColor(style.getColor(e_i));
    collectNormalizedDoubleComboRatios.at(e_i)->SetFillColor(style.getColor(e_i));
    collectNormalizedDoubleComboRatios.at(e_i)->SetFillStyle(1001);
    collectNormalizedDoubleComboRatios.at(e_i)->SetMarkerStyle(style.getMarker(e_i));
    collectNormalizedDoubleComboRatios.at(e_i)->GetListOfFunctions()->Delete();
    collectNormalizedDoubleComboRatios.at(e_i)->SetTitle("");
  }



  if(yAxisLabel_!="")collectNormalizedDoubleComboRatios.at(0)->GetYaxis()->SetTitle(yAxisLabel_);
  if(ymaxrange_!=-1)  collectNormalizedDoubleComboRatios.at(0)->GetYaxis()->SetRangeUser(yminrange_,ymaxrange_);
  else  collectNormalizedDoubleComboRatios.at(0)->GetYaxis()->SetRangeUser(BasePlotExtractors_.at(0)->yProfileMinMax().at(0),BasePlotExtractors_.at(0)->yProfileMinMax().at(1));
  collectNormalizedDoubleComboRatios.at(0)->Draw();
  
  for(unsigned int e_i=0;e_i<nDoubleRatios;e_i++){
    std::cout << "plotting " << collectNormalizedDoubleComboRatios.at(e_i)->GetName() << std::endl;
    collectNormalizedDoubleComboRatios.at(e_i)->Draw("same");
  }
  TPaveText *label = util::LabelFactory::createPaveTextWithOffset(1,1.0,0.0);
  label->AddText(titleLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->Draw("same");

  TLegend* leg1 = util::LabelFactory::createLegendWithOffset(nDoubleRatios,util::LabelFactory::lineHeight());
  for(unsigned int e_i=0;e_i<nDoubleRatios;e_i++)leg1->AddEntry(collectNormalizedDoubleComboRatios.at(e_i),LegendEntries_.at(e_i*4),"PL");

  leg1->Draw();
  BasePlotExtractors_.at(0)->drawCMSPrel();
  TString outname = "OneDNormalizedDoubleComboRatios_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);


  chdir(".."); 
  

}



//should be called by plotAll()
void OverlayPlainProfileAnd2DPlots::plotOverlays() {

  util::StyleSettings::setStyleJMEPaperNoTitle();

  if(chdir("OverlayExtrapolPlots") != 0){ 
    mkdir("OverlayExtrapolPlots", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("OverlayExtrapolPlots"); 
  } 
  
  DefaultStyles style;
  DefaultStyles PFFractionStyle;
  PFFractionStyle.setStyle("PFComp");
  TCanvas* c2 = new TCanvas("c2","",600,600);
  if( BasePlotExtractors_.at(0)->configs_.at(0)->logX() ) c2->SetLogx(1);
  if( BasePlotExtractors_.at(0)->configs_.at(0)->logY() ) c2->SetLogy(1);



  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
    collectOneDPlots_.at(e_i)->SetLineColor(style.getColor(e_i));
    collectOneDPlots_.at(e_i)->SetMarkerColor(style.getColor(e_i));
    collectOneDPlots_.at(e_i)->SetFillColor(style.getColor(e_i));
    collectOneDPlots_.at(e_i)->SetFillStyle(1001);
    collectOneDPlots_.at(e_i)->SetMarkerStyle(style.getMarker(e_i));
    collectOneDPlots_.at(e_i)->GetListOfFunctions()->Delete();
    collectOneDPlots_.at(e_i)->SetTitle("");

    collectTwoDPlots_.at(e_i)->SetLineColor(1);
    //    collectTwoDPlots_.at(e_i)->SetLineColor(PFFractionStyle.getColor(e_i));
    collectTwoDPlots_.at(e_i)->SetMarkerColor(PFFractionStyle.getColor(e_i));
    collectTwoDPlots_.at(e_i)->SetFillColor(PFFractionStyle.getColor(e_i));
    collectTwoDPlots_.at(e_i)->SetFillStyle(1001);
    collectTwoDPlots_.at(e_i)->SetMarkerStyle(PFFractionStyle.getMarker(e_i));
    collectTwoDPlots_.at(e_i)->GetListOfFunctions()->Delete();
    collectTwoDPlots_.at(e_i)->SetTitle("");

  }


  if(yAxisLabel_!="")collectOneDPlots_.at(0)->GetYaxis()->SetTitle(yAxisLabel_);
  if(ymaxrange_!=-1)  collectOneDPlots_.at(0)->GetYaxis()->SetRangeUser(yminrange_,ymaxrange_);
  else collectOneDPlots_.at(0)->GetYaxis()->SetRangeUser(BasePlotExtractors_.at(0)->yProfileMinMax().at(0),BasePlotExtractors_.at(0)->yProfileMinMax().at(1));
  collectOneDPlots_.at(0)->Draw();
  
  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
    std::cout << "plotting " << collectOneDPlots_.at(e_i)->GetName() << std::endl;
    collectOneDPlots_.at(e_i)->Draw("same");
  }
  TPaveText *label = util::LabelFactory::createPaveTextWithOffset(1,1.0,0.0);
  label->AddText(titleLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->Draw("same");

  TLegend* leg1 = util::LabelFactory::createLegendWithOffset(nExtrapols_,util::LabelFactory::lineHeight());
  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++)leg1->AddEntry(collectOneDPlots_.at(e_i),LegendEntries_.at(e_i),"PL");
  TLegend* leg1Fill = util::LabelFactory::createLegendWithOffset(nExtrapols_,util::LabelFactory::lineHeight());
  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++)leg1Fill->AddEntry(collectTwoDPlots_.at(e_i),LegendEntries_.at(e_i),"F");
  //    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(0),dataLabel()/*+" Data"*/,"P");
  leg1Fill->SetFillColor(0);
  leg1Fill->SetFillStyle(0);

  leg1->Draw();
  BasePlotExtractors_.at(0)->drawCMSPrel();
  TString outname = "OneD_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);

  std::cout << "test " << nExtrapols_ << std::endl;
  for(unsigned int e_i=1;e_i<nExtrapols_;e_i++){
    std::cout << "test" << std::endl;
    assert(nExtrapols_>=2);
    //    std::cout << collectRatios_.at(0)->GetNbinsX() << std::endl;
    assert(collectOneDPlots_.at(0)->GetNbinsX()==collectOneDPlots_.at(1)->GetNbinsX());
    std::cout << "test" << std::endl;
    collectRatios_.push_back(util::HistOps::createRatioPlot(collectOneDPlots_.at(e_i),collectOneDPlots_.at(0),"Deviation from "+LegendEntries_.at(0)+" [%]"));
    for(int bin_i = 1; bin_i <= collectRatios_.at(e_i-1)->GetNbinsX(); bin_i++) {
      collectRatios_.at(e_i-1)->SetBinContent(bin_i,(collectRatios_.at(e_i-1)->GetBinContent(bin_i)-1)*100);
      collectRatios_.at(e_i-1)->SetBinError(bin_i,(collectRatios_.at(e_i-1)->GetBinError(bin_i))*100);
    }

  }
 
  std::cout << "about to export ratios as well" <<std::endl;
  collectRatios_.at(0)->Draw("");
  collectRatios_.at(0)->GetYaxis()->SetRangeUser(-10,15);
  for(unsigned int i=0;i<collectRatios_.size();i++){
    collectRatios_.at(i)->Draw("same");
  }
  label->Draw("same");
  TLegend* leg2WoFirst = util::LabelFactory::createLegendWithOffset(nExtrapols_-1,util::LabelFactory::lineHeight());
  for(unsigned int e_i=1;e_i<nExtrapols_;e_i++)leg2WoFirst->AddEntry(collectOneDPlots_.at(e_i),LegendEntries_.at(e_i),"PL");

  leg2WoFirst->Draw();
  BasePlotExtractors_.at(0)->drawCMSPrel();
  outname = "OneDRatio_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);

  THStack *hs = new THStack("hs","");
  for(unsigned int i=0;i<collectTwoDPlots_.size();i++){
    hs->Add(collectTwoDPlots_.at(i)->ProjectionX());
  }
  THStack *hsNorm = new THStack("hsNorm","");
  for(unsigned int i=0;i<collectTwoDPlots_.size();i++){
    hsNorm->Add(collectTwoDPlots_.at(i)->ProjectionX());
  }

  hs->Draw();
  //  hsNorm->Draw();

  hs->GetHistogram()->GetYaxis()->SetTitle("Number of Entries");
  
  hs->GetHistogram()->GetXaxis()->SetTitle(collectOneDPlots_.at(0)->GetXaxis()->GetTitle());


  hs->Draw("hist");
  label->Draw("same");
  leg1Fill->Draw();
  //  collectTwoDPlots_.at(0)->Draw("colz");
  BasePlotExtractors_.at(0)->drawCMSPrel();
  outname = "TwoD_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);


  TIter next(hs->GetHists());
  TObject *obj;// =  next();
  while ((obj =  next())){
    TH1D* temp=(TH1D*)obj;
    temp->SetLineColor(temp->GetFillColor());
    temp->SetFillStyle(0);
  }

//  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//    collectTwoDPlots_.at(e_i)->SetLineColor(PFFractionStyle.getColor(e_i));
//    collectTwoDPlots_.at(e_i)->SetFillStyle(0);
//  }

  hs->Draw("nostack hist");
  label->Draw("same");
  leg1Fill->Draw();
  BasePlotExtractors_.at(0)->drawCMSPrel();
  //  collectTwoDPlots_.at(0)->Draw("colz");
  outname = "TwoDNoStack_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);




  util::HistOps::normalizeEachBinCompositionTHStackHistogram(hsNorm);

  hs->GetHistogram()->GetYaxis()->SetTitle("Relative Composition [%]");
  hs->GetHistogram()->GetXaxis()->SetTitle(collectOneDPlots_.at(0)->GetXaxis()->GetTitle());

//  hs->GetHistogram()->SetMinimum(1);
//  hs->GetHistogram()->SetMaximum(2000);
//  hs->GetHistogram()->GetYaxis()->SetRangeUser(1.,2000);

  TH1 * drawFrame = util::HistOps::createFrame(collectOneDPlots_.at(0), collectOneDPlots_.at(0)->GetXaxis()->GetTitle(), 1, 2000);

  drawFrame->Draw();
  c2->SetLogy(1);
  hs->Draw("nostack hist same");
  label->Draw("same");
  leg1Fill->Draw();
  //  collectTwoDPlots_.at(0)->Draw("colz");
  BasePlotExtractors_.at(0)->drawCMSPrel();
  outname = "TwoDEachBinNormalizedNoStack_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);

  next = hsNorm->GetHists();
  while ((obj =  next())){
    TH1D* temp=(TH1D*)obj;
    temp->SetLineColor(1);
    temp->SetFillStyle(1001);
  }
  c2->SetLogy(0);
  hs->GetHistogram()->GetYaxis()->SetRangeUser(0,-1);


  util::HistOps::normalizeEachBinCompositionTHStackHistogram(hsNorm);

  //  hs->Draw();
  hsNorm->Draw("hist");
  hsNorm->GetHistogram()->GetYaxis()->SetTitle("Relative Composition [%]");
  hsNorm->GetHistogram()->GetXaxis()->SetTitle(collectOneDPlots_.at(0)->GetXaxis()->GetTitle());

  //  hsNorm->Draw("hist same");
    label->Draw("same");
    leg1Fill->Draw();
  //  collectTwoDPlots_.at(0)->Draw("colz");
  BasePlotExtractors_.at(0)->drawCMSPrel();
  outname = "TwoDEachBinNormalized_"+ util::cleanStringMore(titleLabel_)+".eps";
  c2->RedrawAxis();
  c2->SaveAs(outname);




  chdir(".."); 
  


  
//  for(int bin_i=0;bin_i<nBins_;bin_i++){
//    //         int xBin_i=0;
//    for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
//      
//      std::pair <float,float> minMaxPair = BasePlotExtractors_.at(0).determineMinMax(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(0));
//      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//	std::pair <float,float> tempMinMaxPair = BasePlotExtractors_.at(e_i).determineMinMax(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
//	if(tempMinMaxPair.first<minMaxPair.first)minMaxPair.first=tempMinMaxPair.first;
//	if(tempMinMaxPair.second<minMaxPair.second)minMaxPair.second=tempMinMaxPair.second;
//      }
//      c2->DrawFrame(0,minMaxPair.first*0.85-0.05,BasePlotExtractors_.at(0).cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+BasePlotExtractors_.at(0).configs_.at(0)->cutAxisTitle()+";MC/Data Ratio"/*"#sqrt{2} #sigma"*/).c_str());
//      
//      //    collectGraphs.at(0)->Draw("P");
//      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//	BasePlotExtractors_.at(e_i).drawConfidenceIntervals(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
//      }
//      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->Draw("Psame");
//	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetLineColor(style.getColor(e_i));
//	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerColor(style.getColor(e_i));
//	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerStyle(style.getMarker(e_i));
//	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(e_i));
//	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
//	TF1* MCDataRatioTemp=(TF1*) allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->Clone();
//	MCDataRatioTemp->SetRange(0.1,1);
//	MCDataRatioTemp->SetLineStyle(1);
//	MCDataRatioTemp->Draw("same");
//      }
//      
//      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.65-nExtrapols_*util::LabelFactory::lineHeight());
//      label->AddText(BasePlotExtractors_.at(0).jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
//      label->AddText((BasePlotExtractors_.at(0).configs_.at(0)->xBinTitle(xBin_i,BasePlotExtractors_.at(0).configs_.at(0)->binEdges()->at(bin_i),BasePlotExtractors_.at(0).configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
//      label->Draw("same");
//      TLegend* leg1 = util::LabelFactory::createLegendWithOffset(nExtrapols_,0.75-nExtrapols_*util::LabelFactory::lineHeight());
//      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//	leg1->AddEntry(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i),BasePlotExtractors_.at(e_i).yProfileTitle()+" - " + BasePlotExtractors_.at(e_i).kalibriPlotsShortName(),"LP");
//      }
//      leg1->Draw();
//      BasePlotExtractors_.at(0).drawCMSPrel();
//      
//      TString outname = BasePlotExtractors_.at(0).kalibriPlotsShortName()+BasePlotExtractors_.at(0).plotsnames_;
//      for(unsigned int e_i=1;e_i<nExtrapols_;e_i++){
//	outname+="_"+BasePlotExtractors_.at(e_i).kalibriPlotsShortName()+BasePlotExtractors_.at(e_i).plotsnames_;
//      }
//      outname+="_ExtrapolPtThree_MCDataRatio_"+BasePlotExtractors_.at(0).configs_.at(0)->binName(bin_i)+"_"+BasePlotExtractors_.at(0).configs_.at(0)->xBinName(xBin_i);
//      outname+=".eps";
//      c2->SaveAs(outname);
//    }//end plotting in xbins
//    
//    
//    
//    if(nXBins_==1){//only plot when there are indeed more than one nXBins
//      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//        std::pair <float,float> minMaxPair = BasePlotExtractors_.at(e_i).determineMinMax(allCollectedGraphsXBins_.at(0).at(bin_i).at(e_i));
//        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
//  	std::pair <float,float> tempMinMaxPair = BasePlotExtractors_.at(e_i).determineMinMax(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
//  	if(tempMinMaxPair.first<minMaxPair.first)minMaxPair.first=tempMinMaxPair.first;
//  	if(tempMinMaxPair.second<minMaxPair.second)minMaxPair.second=tempMinMaxPair.second;
//        }
//	//        std::cout << "works here" << std::endl;
//        c2->DrawFrame(0,minMaxPair.first*0.85-0.05,BasePlotExtractors_.at(e_i).cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+BasePlotExtractors_.at(e_i).configs_.at(0)->cutAxisTitle()+";MC/Data Ratio"/*"#sqrt{2} #sigma"*/).c_str());
//        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
//  	BasePlotExtractors_.at(e_i).drawConfidenceIntervals(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
//        }
//        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
//  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->Draw("Psame");
//  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetLineColor(style.getColor(e_i));
//  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerColor(style.getColor(e_i));
//  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerStyle(style.getMarker(e_i));
//  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(e_i));
//  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
//  	TF1* MCDataRatioTemp=(TF1*) allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->Clone();
//  	MCDataRatioTemp->SetRange(0.1,1);
//  	MCDataRatioTemp->SetLineStyle(1);
//  	MCDataRatioTemp->Draw("same");
//        }
//        
//	//        std::cout << "works here2" << std::endl;
//        TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.6-(nXBins_/3)*util::LabelFactory::lineHeight());
//        label->AddText(BasePlotExtractors_.at(e_i).jetLabel_+"; "+ BasePlotExtractors_.at(e_i).yProfileTitle()+" - " + BasePlotExtractors_.at(e_i).kalibriPlotsShortName());//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
//        label->AddText((TString)"0/1 below #equiv 0#leq" + (BasePlotExtractors_.at(e_i).configs_.at(0)->xTitle()).c_str() + (TString) "<1" );
//        label->Draw("same");
//        TLegend* leg1 = util::LabelFactory::createLegendWithOffset(nExtrapols_,0.68-(nXBins_/3)*util::LabelFactory::lineHeight());
//        leg1->SetNColumns(3);
//      
//        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
//  	//	for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
//  	TString temp="";
//  	temp+=(BasePlotExtractors_.at(e_i).configs_.at(0)->xBinEdges()->at(xBin_i));
//  	temp+="/";
//  	temp+=(BasePlotExtractors_.at(e_i).configs_.at(0)->xBinEdges()->at(xBin_i+1));
//  	leg1->AddEntry(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i),temp ,"LP");
//        }
//        leg1->Draw();
//        BasePlotExtractors_.at(e_i).drawCMSPrel();
//        
//        TString outname = BasePlotExtractors_.at(e_i).kalibriPlotsShortName()+BasePlotExtractors_.at(e_i).plotsnames_;
//        outname+="_ExtrapolPtThree_XBins_MCDataRatio_"+BasePlotExtractors_.at(e_i).configs_.at(0)->binName(bin_i)+"_"+BasePlotExtractors_.at(e_i).configs_.at(0)->xBinName(0);
//        outname+=".eps";
//        c2->SaveAs(outname);
//      }//end plotting in BasePlotExtractors
//    }//end if    
//  }
//  chdir(".."); 
  
}
  
  
  void OverlayPlainProfileAnd2DPlots::plotAll() {
  checkConsistency();
  Init();
  plotOverlays();
}


  void OverlayPlainProfileAnd2DPlots::plotNormalizedComboRatiosAll() {
  checkConsistency();
  Init();
  plotNormalizedComboOverlays();
}

  void OverlayPlainProfileAnd2DPlots::plotNormalizedDoubleComboRatiosAll() {
  checkConsistency();
  Init();
  plotNormalizedDoubleComboOverlays();
}

void OverlayPlainProfileAnd2DPlots::ConfigureSetRangeUser(double ymin, double ymax){
  yminrange_= ymin;
  ymaxrange_= ymax;
}
