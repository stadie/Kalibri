#include "Extrapolation.h"

//Default constructor
Extrapolation::Extrapolation(TString plotsnames,TString kalibriPlotsPath) : BasePlotExtractor(plotsnames,kalibriPlotsPath){
  //init config_file;

  //  init("GaussFitWidth");
  init(profileType_);
  //  if(plotsnames=="RelResponseVsPt")  init("RatioOfMeans");
  //  else init("StandardDeviation");
  extrapolInit();
}

void Extrapolation::extrapolInit() {
  
  cutNames_ = bag_of_string(ExternalConfig_.read<std::string>("TwoJetsPtBalanceEvent plots cut_list",""));
cutNumbers_ = bag_of<double>(ExternalConfig_.read<std::string>("TwoJetsPtBalanceEvent plots cut_no_list",""));

}



void Extrapolation::Plot() {
  MakeDateDir();
  if(chdir("Extrapol") != 0){ 
    mkdir("Extrapol", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("Extrapol"); 
  } 
  createPtRelExtrapol();


  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  c->SetLogx(configs_.at(0)->logX());
  DefaultStyles style;
  style.setStyle();

  //  createRatios();
  //  for(int bin_i=0;bin_i<pConfig->nBins();bin_i++){
  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
    int nEntries =2;//names_.size();
  TLegend *leg = new TLegend(0.2,0.9-nEntries*0.07,0.4,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  
    for(int conf_i=0;conf_i<configs_.size();conf_i++){
      for(int j=0;j<2;j++){
	//	if(j==0){
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetLineColor(style.getColor(0));
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerColor(style.getColor(0));
	  //	}
	  //	else {
	  //	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetLineColor(1);
	  //	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerColor(1);
	  //	  
	  //	}
	//	AllPlots_.at(conf_i).at(bin_i).at(j)->SetFillColor(style.getColor(conf_i));
	  //	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetFillStyle(3001);
	  AllPlots_.at(conf_i).at(bin_i).at(j)->SetMarkerStyle(style.getMarker(0));
      }
    }
    leg->AddEntry(AllPlots_.at(0).at(bin_i).at(0),(/*configs_.at(0)->yTitle()+*/" Data")/*.c_str()*/,"P");
      leg->AddEntry(AllPlots_.at(0).at(bin_i).at(1),(/*configs_.at(0)->yTitle()+*/" MC")/*.c_str()*/,"L");



  for(int i=0;i<configs_.size();i++){
    //Draw res in data and MC

    AllPlots_.at(i).at(bin_i).at(1)->Draw("hist");
    AllPlots_.at(i).at(bin_i).at(1)->GetYaxis()->SetRangeUser(yProfileMinMax().at(0),yProfileMinMax().at(1));
    AllPlots_.at(i).at(bin_i).at(1)->GetYaxis()->SetTitle(yProfileTitle());
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    AllPlots_.at(i).at(bin_i).at(1)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    AllPlots_.at(i).at(bin_i).at(1)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());
 
    AllPlots_.at(i).at(bin_i).at(0)->Draw("same");

    AllPlots_.at(i).at(bin_i).at(1)->Draw("hist same");
    leg->SetHeader((configs_.at(i)->binTitle(configs_.at(i)->binEdges()->at(bin_i),
					     configs_.at(i)->binEdges()->at(bin_i+1))).c_str());
    leg->Draw();
    cmsPrel();
    TString outname = "ResolutionPlots_"+plotsnames_+"_"+cutNames_.at(i)+"_"+configs_.at(0)->binName(bin_i);
    //  outname+=bin_i;
    outname+=".eps";
    c->SaveAs(outname);

    //Draw DataMC-ratios
    std::cout << AllRatiosDataMC_.size() << " and " << std::endl;
    AllRatiosDataMC_.at(i).at(bin_i)->Draw();
    drawConfidenceIntervals(AllRatiosDataMC_.at(i).at(bin_i));
    AllRatiosDataMC_.at(i).at(bin_i)->Draw("same");
    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetRangeUser(yRatioMinMax().at(0),yRatioMinMax().at(1));
    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetTitle(yRatioTitle());
    //    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetRangeUser(0.5,1.5);
    //    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetTitle("Resolution ratio");
    AllRatiosDataMC_.at(i).at(bin_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    AllRatiosDataMC_.at(i).at(bin_i)->GetXaxis()->SetRangeUser(configs_.at(0)->xMin(),configs_.at(0)->xMax());


  TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
  label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->AddText((configs_.at(i)->binTitle(configs_.at(i)->binEdges()->at(bin_i),
					   configs_.at(i)->binEdges()->at(bin_i+1))).c_str());//(configs_.at(0)->BinTitle(xBin_i,configs_.at(0)->binEdges()->at(bin_i),configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
  TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
  addFunctionLabelsToLegend( AllRatiosDataMC_.at(i).at(bin_i),leg1);
  //  leg1->AddEntry(MCExtrapols_.at(xBin_i),"Extrapolation (MC)","LP");
  //  leg1->AddEntry(DataExtrapols_.at(xBin_i),"Extrapolation (data)","LP");
  label->Draw("same");
  leg1->Draw();
  cmsPrel();





     cmsPrel();
    outname = "ResolutionPlots_"+plotsnames_+"_ratio"+cutNames_.at(i)+"_"+configs_.at(0)->binName(bin_i);
    //  outname+=bin_i;
    outname+=".eps";
    c->SaveAs(outname);

  }
  
  


  }

  for(int conf_i=0;conf_i<configs_.size();conf_i++){
  c->SetLogx(0);
  //  RatioVsBinVarHistos_.at(conf_i)
    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetRangeUser(0.7,1.3);
    //    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetTitle("Resolution ratio");
    TStyle *tdrStyle = (TStyle*)gROOT->FindObject("tdrStyle"); 
    RatioVsBinVarHistos_.at(conf_i)->GetYaxis()->SetTitleOffset(tdrStyle->GetTitleYOffset());
    RatioVsBinVarHistos_.at(conf_i)->Draw();
    cmsPrel();
    TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
    label->AddText(jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
    label->SetFillStyle(0);
    label->AddText(yProfileTitle()/*plotsnames_*/);
    label->Draw("same");
    TString outname = "ResolutionPlots_"+plotsnames_+"_RatioVsBinVar"+cutNames_.at(conf_i)+"_"+names_.at(conf_i);
    //  outname+=bin_i;
    outname+=".eps";
    c->SaveAs(outname);

  }


 chdir("../../."); 

}


void Extrapolation::createPtRelExtrapol() {


  VecOfTH1vec_t AllExtrapolatedRes;
  for(int bin_i=0;bin_i<configs_.at(0)->nBins();bin_i++){
    ExtrapolateBin ExtrapolationBin(this);
    
    for(int conf_i=0;conf_i<configs_.size();conf_i++){
      ExtrapolationBin.addMCHisto(AllPlots_.at(conf_i).at(bin_i).at(1));
      ExtrapolationBin.addDataHisto(AllPlots_.at(conf_i).at(bin_i).at(0));
    }
    for(Int_t xbin_i=0;xbin_i<configs_.at(0)->nXBins();xbin_i++){
      ExtrapolationBin.createTGraphErrors(xbin_i);
      ExtrapolationBin.plotExtrapol(xbin_i,bin_i);
    }
    ExtrapolationBin.produceExtrapolatedRes();
    TH1vec_t CollectExtrapolatedRes;
    CollectExtrapolatedRes.push_back(ExtrapolationBin.ExtrapolatedResData());
    CollectExtrapolatedRes.push_back(ExtrapolationBin.ExtrapolatedResMC());
    AllExtrapolatedRes.push_back(CollectExtrapolatedRes);

  }

  AllPlots_.push_back(AllExtrapolatedRes);
  cutNumbers_.push_back(0.0);
  cutNames_.push_back("00");
  configs_.push_back(configs_.at(0));
  functions_.push_back(functions_.at(0));
  profiles_.push_back(profiles_.at(0));
  names_.push_back(names_.at(0)+"_Extrapol");
  refreshRatiosDataMC();
  makeRatioVsBinVarHistos();
}


Extrapolation::ExtrapolateBin::ExtrapolateBin(Extrapolation* Outer) {
  Outer_=Outer;

}


void Extrapolation::ExtrapolateBin::addMCHisto(TH1D* MCHisto){
  MCHistos_.push_back(MCHisto);
}


void Extrapolation::ExtrapolateBin::addDataHisto(TH1D* DataHisto){
  DataHistos_.push_back(DataHisto);
}

void Extrapolation::ExtrapolateBin::createTGraphErrors(Int_t xBin_i){
  //  std::cout << Outer_->cutNumbers_.at(1) << std::endl;
  
  
  //  std::cout << Outer_->configs_.at(0)->xBinEdges()->at(xBin_i) << std::endl;
  std::vector<double> x,x_e,MCy,MCy_e,Datay,Datay_e;

  //  for()
  for (Int_t cut_i =0;cut_i<Outer_->cutNumbers_.size();cut_i++){
    x.push_back(Outer_->cutNumbers_.at(cut_i));
    x_e.push_back(0.);
    MCy.push_back(MCHistos_.at(cut_i)->GetBinContent(xBin_i));
    MCy_e.push_back(MCHistos_.at(cut_i)->GetBinError(xBin_i));
    Datay.push_back(DataHistos_.at(cut_i)->GetBinContent(xBin_i));
    Datay_e.push_back(DataHistos_.at(cut_i)->GetBinError(xBin_i));

       std::cout << Outer_->cutNumbers_.at(cut_i) << std::endl;
    
    
       std::cout<< MCy.back() << " and error: " << MCy_e.back() << std::endl;
  }

  
  TGraphErrors *extrapol_MC = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&MCy[0],&x_e[0],&MCy_e[0]);
  TGraphErrors *extrapol_Data = new TGraphErrors(Outer_->cutNumbers_.size(),&x[0],&Datay[0],&x_e[0],&Datay_e[0]);
  //  TGraphErrors *gr_res1 = new TGraphErrors(n,&x_ptthree_[0],&y_mean_ratio_res1_[0],&ex_ptthree_[0],&ey_mean_ratio_res1_[0]);
  TF1 *lin_extrapol = new TF1("lin_extrapol","[0]+[1]*x",0,0.5); //was used before...
  lin_extrapol->SetParameters(1,-0.1);
  lin_extrapol->SetParName(0,"ResZero");
  lin_extrapol->SetParName(1,"slope");
  
  extrapol_MC->Fit("lin_extrapol","Q","same",0,0.5);
  extrapol_Data->Fit("lin_extrapol","Q","same",0,0.5);

  MCExtrapols_.push_back(extrapol_MC);
  DataExtrapols_.push_back(extrapol_Data);


}

std::pair <float,float> Extrapolation::ExtrapolateBin::determineMinMax(TGraphErrors* graph){
  std::pair <float,float> minMaxPair(graph->GetY()[0],graph->GetY()[0]);
  for(Int_t i=0;i<graph->GetN();i++){
    if(graph->GetY()[i]<minMaxPair.first)minMaxPair.first=graph->GetY()[i];
    if(graph->GetY()[i]>minMaxPair.second)minMaxPair.second=graph->GetY()[i];
  }
  return minMaxPair;
}


void Extrapolation::ExtrapolateBin::plotExtrapol(Int_t xBin_i, Int_t bin_i){
  //  std::cout << Outer_->cutNumbers_.at(1) << std::endl;
  if(chdir("Resolution") != 0){ 
    mkdir("Resolution", S_IRWXU|S_IRWXG|S_IRWXO); 
    chdir("Resolution"); 
  } 
  
  setTDRStyle();
  DefaultStyles style;
  style.setStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  std::pair <float,float> minMaxPair = determineMinMax(MCExtrapols_.at(xBin_i));
  c->DrawFrame(0,minMaxPair.first*0.5-0.05,Outer_->cutNumbers_.back()+0.05,minMaxPair.second*1.2,(";"+Outer_->configs_.at(0)->cutAxisTitle()+";"+Outer_->yProfileTitle()/*"#sqrt{2} #sigma"*/)/*.c_str()*/);
  MCExtrapols_.at(xBin_i)->Draw("P");
  MCExtrapols_.at(xBin_i)->SetLineColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(0));
  MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(0));
  MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
  TF1* MCTemp=(TF1*) MCExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->Clone();
  MCTemp->SetRange(0.1,1);
  MCTemp->SetLineStyle(1);
  MCTemp->Draw("same");
  TH1F* DrawHist = (TH1F*) MCExtrapols_.at(xBin_i)->GetHistogram();
  DrawHist->GetXaxis()->SetRangeUser(0.0,Outer_->cutNumbers_.back()+0.05);
  DrawHist->GetYaxis()->SetRangeUser(0.0,0.15);
  DataExtrapols_.at(xBin_i)->Draw("Psame");
  DataExtrapols_.at(xBin_i)->SetLineColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->SetMarkerColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->SetMarkerStyle(style.getMarker(1));
  DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineColor(style.getColor(1));
  DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->SetLineStyle(2);
  TF1* DataTemp=(TF1*) DataExtrapols_.at(xBin_i)->GetFunction("lin_extrapol")->Clone();
  DataTemp->SetRange(0.1,1);
  DataTemp->SetLineStyle(1);
  DataTemp->Draw("same");
  int nEntries =2;


  TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
//  TPaveText *label = util::LabelFactory::createPaveText(1);
  label->AddText(Outer_->jetLabel_);//+  |#eta_{1,2}| > "+util::toTString(1.4)+",  L = "+util::StyleSettings::luminosity(4.6));
  label->AddText((Outer_->configs_.at(0)->xBinTitle(xBin_i,Outer_->configs_.at(0)->binEdges()->at(bin_i),Outer_->configs_.at(0)->binEdges()->at(bin_i+1),0)).c_str()/*+(TString)" and "*/);
//  else label->AddText(jetLabel+",  L = "+util::StyleSettings::luminosity(lumi));
  label->Draw("same");
  TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
  //  TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(4,0.5,1);

//  TLegend *leg = new TLegend(0.2,0.9-nEntries*0.07,0.4,0.9);
//  leg->SetBorderSize(0);
//  leg->SetFillColor(0);
//  leg->SetFillStyle(0);
//  leg->SetTextFont(42);
//  leg->SetTextSize(0.04);

  leg1->AddEntry(MCExtrapols_.at(xBin_i),"Extrapolation (MC)","LP");
  leg1->AddEntry(DataExtrapols_.at(xBin_i),"Extrapolation (data)","LP");

  leg1->Draw();
  cmsPrel();


  TString outname = "ResolutionPlots_ExtrapolPtThree_"+Outer_->plotsnames_+"_"+Outer_->configs_.at(0)->binName(bin_i)+"_"+Outer_->configs_.at(0)->xBinName(xBin_i);
  //  outname+=bin_i;
  outname+=".eps";
  c->SaveAs(outname);
  chdir(".."); 


}

void Extrapolation::ExtrapolateBin::produceExtrapolatedRes(){

  TH1D* tempExtrapolatedResMC = (TH1D*) MCHistos_.at(0)->Clone();
  TH1D* tempExtrapolatedResData = (TH1D*) DataHistos_.at(0)->Clone();

  for(Int_t xbin_i=0;xbin_i<Outer_->configs_.at(0)->nXBins();xbin_i++){
    tempExtrapolatedResMC->SetBinContent(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedResData->SetBinContent(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0) : 0.0);
    tempExtrapolatedResMC->SetBinError(xbin_i+1,MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? MCExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
    tempExtrapolatedResData->SetBinError(xbin_i+1,DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParameter(0)>0.01 ? DataExtrapols_.at(xbin_i)->GetFunction("lin_extrapol")->GetParError(0) : 0.0);
  }
  ExtrapolatedResMC_=(TH1D*) tempExtrapolatedResMC;
  ExtrapolatedResData_=(TH1D*) tempExtrapolatedResData;
    
}
