#include "OverlayExtrapolGraphs.h"

void OverlayExtrapolGraphs::addPlots(TString plotsnames,TString kalibriPlotsShortName) {
  std::cout << "bla " << std::endl;
  Extrapolations_.push_back(Extrapolation(plotsnames,kalibriPlotsShortName));
  chdir(Extrapolations_.back().outputPathROOT());

}


void OverlayExtrapolGraphs::readInGraphs() {
  std::cout << "reading in " << nExtrapols_ << " graphs with " << nXBins_ << " nXBins and " << nBins_ << " nBins_" <<  std::endl;
  for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
      VecOfTGErrvec_t allCollectedGraphs;
      for(int bin_i=0;bin_i<nBins_;bin_i++){
	int xBin_i=0;
	
	std::vector <TGraphErrors*> collectGraphs;
	for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
	  
	  //      std::cout <<Extrapolations_.at(e_i).outputPathROOT() <<std::endl;
	  //      int bin_i=0;
	  
	  TString inName =   "ResolutionPlots_ExtrapolPtThree_MCDataRatio_"+Extrapolations_.at(e_i).plotsnames_+"_"+Extrapolations_.at(e_i).configs_.at(0)->binName(bin_i)+"_"+Extrapolations_.at(e_i).configs_.at(0)->xBinName(xBin_i);
	  
	  collectGraphs.push_back((TGraphErrors*) util::FileOps::readTGraphErrors(Extrapolations_.at(e_i).outputPathROOT()+"/Output"+Extrapolations_.at(e_i).kalibriPlotsShortName()+".root",Extrapolations_.at(e_i).configs_.at(0)->name()+"/"+inName));
	  
	}
	
	allCollectedGraphs.push_back(collectGraphs);
      }
      allCollectedGraphsXBins_.push_back(allCollectedGraphs);
    }
}

void OverlayExtrapolGraphs::checkConsistency() {
  std::cout << "checking consistency now..." << std::endl;
  nBins_=Extrapolations_.at(0).configs_.at(0)->nBins();
  nXBins_=Extrapolations_.at(0).configs_.at(0)->nXBins();
  nExtrapols_=Extrapolations_.size();

  for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
    assert(Extrapolations_.at(e_i).configs_.at(0)->nBins()==nBins_);
    assert(Extrapolations_.at(e_i).configs_.at(0)->nXBins()==nXBins_);
  }

//  for(int bin_i=0;bin_i<nBins_;bin_i++){
//    assert(allCollectedGraphs_.at(bin_i).size()==nExtrapols_);
//    //     assert(Extrapolations_.at(e_i).configs_.at(0)->nBins()==nBins_);
//  }
  std::cout << "checked consistency successfully..." << std::endl;


}


//should be called by plotAll()
void OverlayExtrapolGraphs::plotOverlays() {
  if(chdir("OverlayExtrapolPlots") != 0){ 
    mkdir("OverlayExtrapolPlots", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("OverlayExtrapolPlots"); 
  } 
  
  DefaultStyles style;
  TCanvas* c2 = new TCanvas("c2","",600,600);
  
  for(int bin_i=0;bin_i<nBins_;bin_i++){
    //         int xBin_i=0;
    for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
      
      std::pair <float,float> minMaxPair = Extrapolations_.at(0).determineMinMax(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(0));
      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
	std::pair <float,float> tempMinMaxPair = Extrapolations_.at(e_i).determineMinMax(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
	if(tempMinMaxPair.first<minMaxPair.first)minMaxPair.first=tempMinMaxPair.first;
	if(tempMinMaxPair.second<minMaxPair.second)minMaxPair.second=tempMinMaxPair.second;
      }
      c2->DrawFrame(0,minMaxPair.first*0.85-0.05,Extrapolations_.at(0).cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+Extrapolations_.at(0).configs_.at(0)->cutAxisTitle()+";MC/Data Ratio"/*"#sqrt{2} #sigma"*/).c_str());
      
      //    collectGraphs.at(0)->Draw("P");
      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
	Extrapolations_.at(e_i).drawConfidenceIntervals(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
      }
      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->Draw("Psame");
	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetLineColor(style.getColor(e_i));
	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerColor(style.getColor(e_i));
	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerStyle(style.getMarker(e_i));
	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(e_i));
	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
	TF1* MCDataRatioTemp=(TF1*) allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->Clone();
	MCDataRatioTemp->SetRange(0.1,1);
	MCDataRatioTemp->SetLineStyle(1);
	MCDataRatioTemp->Draw("same");
      }
      
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.65-nExtrapols_*util::LabelFactory::lineHeight());
      label->AddText(Extrapolations_.at(0).jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
      label->AddText((Extrapolations_.at(0).configs_.at(0)->xBinTitle(xBin_i,Extrapolations_.at(0).configs_.at(0)->binEdges()->at(bin_i),Extrapolations_.at(0).configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
      label->Draw("same");
      TLegend* leg1 = util::LabelFactory::createLegendWithOffset(nExtrapols_,0.75-nExtrapols_*util::LabelFactory::lineHeight());
      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
	leg1->AddEntry(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i),Extrapolations_.at(e_i).yProfileTitle()+" - " + Extrapolations_.at(e_i).kalibriPlotsShortName(),"LP");
      }
      leg1->Draw();
      Extrapolations_.at(0).drawCMSPrel();
      
      TString outname = Extrapolations_.at(0).kalibriPlotsShortName()+Extrapolations_.at(0).plotsnames_;
      for(unsigned int e_i=1;e_i<nExtrapols_;e_i++){
	outname+="_"+Extrapolations_.at(e_i).kalibriPlotsShortName()+Extrapolations_.at(e_i).plotsnames_;
      }
      outname+="_ExtrapolPtThree_MCDataRatio_"+Extrapolations_.at(0).configs_.at(0)->binName(bin_i)+"_"+Extrapolations_.at(0).configs_.at(0)->xBinName(xBin_i);
      outname+=".eps";
      c2->SaveAs(outname);
    }//end plotting in xbins
    
    
    
    if(nXBins_==1){//only plot when there are indeed more than one nXBins
      for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
        std::pair <float,float> minMaxPair = Extrapolations_.at(e_i).determineMinMax(allCollectedGraphsXBins_.at(0).at(bin_i).at(e_i));
        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
  	std::pair <float,float> tempMinMaxPair = Extrapolations_.at(e_i).determineMinMax(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
  	if(tempMinMaxPair.first<minMaxPair.first)minMaxPair.first=tempMinMaxPair.first;
  	if(tempMinMaxPair.second<minMaxPair.second)minMaxPair.second=tempMinMaxPair.second;
        }
	//        std::cout << "works here" << std::endl;
        c2->DrawFrame(0,minMaxPair.first*0.85-0.05,Extrapolations_.at(e_i).cutNumbers_.back()+0.05,minMaxPair.second*1.1,(";cut on "+Extrapolations_.at(e_i).configs_.at(0)->cutAxisTitle()+";MC/Data Ratio"/*"#sqrt{2} #sigma"*/).c_str());
        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
  	Extrapolations_.at(e_i).drawConfidenceIntervals(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i));
        }
        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->Draw("Psame");
  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetLineColor(style.getColor(e_i));
  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerColor(style.getColor(e_i));
  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->SetMarkerStyle(style.getMarker(e_i));
  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(e_i));
  	allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
  	TF1* MCDataRatioTemp=(TF1*) allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i)->GetFunction("lin_extrapol")->Clone();
  	MCDataRatioTemp->SetRange(0.1,1);
  	MCDataRatioTemp->SetLineStyle(1);
  	MCDataRatioTemp->Draw("same");
        }
        
	//        std::cout << "works here2" << std::endl;
        TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.6-(nXBins_/3)*util::LabelFactory::lineHeight());
        label->AddText(Extrapolations_.at(e_i).jetLabel_+"; "+ Extrapolations_.at(e_i).yProfileTitle()+" - " + Extrapolations_.at(e_i).kalibriPlotsShortName());//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
        label->AddText((TString)"0/1 below #equiv 0#leq" + (Extrapolations_.at(e_i).configs_.at(0)->xTitle()).c_str() + (TString) "<1" );
        label->Draw("same");
        TLegend* leg1 = util::LabelFactory::createLegendWithOffset(nExtrapols_,0.68-(nXBins_/3)*util::LabelFactory::lineHeight());
        leg1->SetNColumns(3);
      
        for(Int_t xBin_i=0;xBin_i<nXBins_;xBin_i++){
  	//	for(unsigned int e_i=0;e_i<nExtrapols_;e_i++){
  	TString temp="";
  	temp+=(Extrapolations_.at(e_i).configs_.at(0)->xBinEdges()->at(xBin_i));
  	temp+="/";
  	temp+=(Extrapolations_.at(e_i).configs_.at(0)->xBinEdges()->at(xBin_i+1));
  	leg1->AddEntry(allCollectedGraphsXBins_.at(xBin_i).at(bin_i).at(e_i),temp ,"LP");
        }
        leg1->Draw();
        Extrapolations_.at(e_i).drawCMSPrel();
        
        TString outname = Extrapolations_.at(e_i).kalibriPlotsShortName()+Extrapolations_.at(e_i).plotsnames_;
        outname+="_ExtrapolPtThree_XBins_MCDataRatio_"+Extrapolations_.at(e_i).configs_.at(0)->binName(bin_i)+"_"+Extrapolations_.at(e_i).configs_.at(0)->xBinName(0);
        outname+=".eps";
        c2->SaveAs(outname);
      }//end plotting in Extrapolations
    }//end if    
  }
  chdir(".."); 
  
}
  
  
  void OverlayExtrapolGraphs::plotAll() {
  checkConsistency();
  readInGraphs();
  plotOverlays();
}
