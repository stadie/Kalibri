// Adapted from Matthias Usercode
//  (computeTriggerEfficiencies.C,v 1.5 2012/01/24 10:16:38 mschrode Exp )

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TString.h"

#include "tdrstyle_mod.C"

#define UTILS_AS_HEADER_FILE
#include "/afs/naf.desy.de/user/k/kirschen/public/util/utils.h"
#include "/afs/naf.desy.de/user/k/kirschen/public/util/HistOps.h"
#include "/afs/naf.desy.de/user/k/kirschen/public/util/LabelFactory.h"
#include "/afs/naf.desy.de/user/k/kirschen/public/util/StyleSettings.h"


TString single_dijet;
TString single_dijet_ptVarLabel;
TString TriggerYearLabel;

//! Trigger efficiency plots from Kalibri ntuple
//!
//! Note: traditionally you would compute the efficiency of a
//! trigger A using a reference trigger B as
//!
//!   eff = N(passing A && passing B) / N(passing B) .
//!
//! However, due to the sometimes large pre-scales there might
//! only be very few events available for this method preventing
//! a reliable measurement. Instead, the efficiency is defined
//! as
//!
//!   eff = N(passing A) / N(passing B) .
//!
//! Of course, now 0 <= eff <= 1 is not true anymore but without
//! knowledge of the pre-scales the previous definition is not
//! exact either. The pt of the turn-on can still be correctly
//! determined as the pt where eff is 99% of the plateau value
//! of eff.
//!
//! TODO: Encapsulate information per trigger into type


//! Fit function for turn-on
// --------------------------------------------------
double eff(double *x, double *par) {
  return 0.5*par[2]*(erf(par[0]*(x[0]-par[1]))+1);
}



//! Create TChain from input root files. The root
//! files are expected to contain a TTree "DiJetTree".
//! There are two possible input options:
//!
//! 1) 'fileName' specifies a single root file; it ends
//!    with '.root';
//! 2) 'fileName' contains a list of root file names.
// --------------------------------------------------
TChain *createTChain(const TString &fileName) {
  TChain* chain = new TChain("DiJetTree"); 

  // Option 1: single root file
  if( fileName.EndsWith(".root") ) {
    chain->Add(util::absolutePath(fileName));
  }
  // Option 2: list of root files
  else {
    std::ifstream filelist;
    filelist.open(util::absolutePath(fileName));
    int nOpenedFiles = 0;
    if( filelist.is_open() ) {
      TString name = "";
      while( !filelist.eof() ) {
	filelist >> name;
	if( filelist.eof() ) break;
	chain->Add(name);
	nOpenedFiles++;
      }
    } else {
      std::cerr << "ERROR opening file '" << fileName << "'\n";
      exit(1);
    }
    filelist.close();
  }

  return chain;
}



// --------------------------------------------------
void getTriggerInfo(std::vector<TString> &tagTriggers, std::vector<TString> &probeTriggers, std::vector<double> &probeTriggerVals, std::vector<double> &tagTriggerVals) {
  std::vector<double> NominalPtOfTriggers;
  NominalPtOfTriggers.push_back(30);
  NominalPtOfTriggers.push_back(60);
  NominalPtOfTriggers.push_back(80);
  NominalPtOfTriggers.push_back(110);
  NominalPtOfTriggers.push_back(150);
  NominalPtOfTriggers.push_back(190);
  NominalPtOfTriggers.push_back(240);
  NominalPtOfTriggers.push_back(300);
  NominalPtOfTriggers.push_back(370);
  std::vector<TString> triggers;
  if(TriggerYearLabel=="2011"){
    if(single_dijet=="DiJetAve"||single_dijet=="Jet"){

      for(size_t i =0; i<NominalPtOfTriggers.size();i++){
	triggers.push_back("Hlt"+single_dijet+util::toTString(NominalPtOfTriggers.at(i))); 
      }
    }
    else{exit;}
    
    for(size_t i = 0; i < triggers.size()-1; ++i) {
      tagTriggers.push_back(triggers.at(i));
      tagTriggerVals.push_back(NominalPtOfTriggers.at(i));
      probeTriggers.push_back(triggers.at(i+1));
      probeTriggerVals.push_back(NominalPtOfTriggers.at(i+1));
    }
    
    assert( probeTriggerVals.size() == probeTriggers.size() );
  }
}


//TH1D

// --------------------------------------------------
void computeTriggerEfficiences(const TString &fileNames, double lumi, int nEvts = -1, double etaMin = 0.) {

  std::cout << "Preparing input...\n";

  // Style
  //  util::StyleSettings::setStyleNote();
  util::StyleSettings::setStyleJMEPaper();
  TString jetType = util::LabelFactory::jetAlgo(fileNames)+util::LabelFactory::jetType(fileNames);
  TString jetLabel = util::LabelFactory::labelJet(fileNames);
  int nBins = 250;
  if( etaMin > 0. ) nBins = 100;
  double histMin = 0.;
  double histMax = 500.;

  // The trigger names
  std::vector<TString> tagTriggers;
  std::vector<TString> probeTriggers;
  std::vector<double> probeTriggerVals;
  std::vector<double> tagTriggerVals;
  getTriggerInfo(tagTriggers,probeTriggers,probeTriggerVals,tagTriggerVals);

  // Fit ranges
  std::vector<double> ptMin;
  std::vector<double> ptMax;
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
        ptMin.push_back(0.7*probeTriggerVals.at(i));
    //    ptMax.push_back(std::min(1.5*probeTriggerVals.at(i),histMax));
    //    ptMin.push_back(0.5*probeTriggerVals.at(i));
    ptMax.push_back(std::min(2*probeTriggerVals.at(i),histMax));
  }

  // Event counters
  std::vector<TH1*> hNTag;
  std::vector<TH1*> hNProbe;
  std::vector<TH1*> hEff;
  std::vector<TF1*> fEff;
  std::vector<TLine*> lPtVar99;

  //  single_dijet_ptVarLabel

  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    TH1* h = util::HistOps::createTH1D("hNTag_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin),nBins,histMin,histMax,""+single_dijet_ptVarLabel+"","GeV","events");
    h->Sumw2();
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.1);
    hNTag.push_back(h);

    h = static_cast<TH1*>(hNTag.back()->Clone("hNProbe_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin)));
    hNProbe.push_back(h);

    h = static_cast<TH1*>(hNTag.back()->Clone("hEff_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin)));
    hEff.push_back(h);

    TF1* f = new TF1("efficiency_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin),eff,ptMin.at(i),ptMax.at(i),3);
    f->SetLineColor(kBlue);
    f->SetLineWidth(2);
    fEff.push_back(f);

    TLine* l = new TLine(1.,0.,1.,1.);
    l->SetLineWidth(fEff.back()->GetLineWidth());
    l->SetLineColor(fEff.back()->GetLineColor());
    lPtVar99.push_back(l);
  }
  std::vector<double> plateauVals;
  std::vector<double> ptVar99;


  std::cout << "Preparing chain...\n";
  // Set up chain
  TChain *chain = createTChain(fileNames);

  // Set branch addresses
  const int maxNJet = 50;
  
  std::vector<double> corrJetPt(20);

  int nObjJet = 0;
  float jetPt[maxNJet];
  float jetEta[maxNJet];
  float jetCorrL1[maxNJet];
  float jetCorrL2L3[maxNJet];
  bool hltphys = false;
  bool hlt30 = false; 
  bool hlt60 = false; 
  bool hlt80 = false; 
  bool hlt110 = false;
  bool hlt150 = false;
  bool hlt190 = false;
  bool hlt240 = false;
  bool hlt300 = false;
  bool hlt370 = false;
  chain->SetBranchAddress("NobjJet",&nObjJet);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetCorrL1",jetCorrL1);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  chain->SetBranchAddress("HltJet30",&hltphys); 
  chain->SetBranchAddress("Hlt"+single_dijet+"30",&hlt30); 
  chain->SetBranchAddress("Hlt"+single_dijet+"60",&hlt60); 
  chain->SetBranchAddress("Hlt"+single_dijet+"80",&hlt80); 
  chain->SetBranchAddress("Hlt"+single_dijet+"110",&hlt110);
  chain->SetBranchAddress("Hlt"+single_dijet+"150",&hlt150);
  chain->SetBranchAddress("Hlt"+single_dijet+"190",&hlt190);
  chain->SetBranchAddress("Hlt"+single_dijet+"240",&hlt240);
  chain->SetBranchAddress("Hlt"+single_dijet+"300",&hlt300);
  chain->SetBranchAddress("Hlt"+single_dijet+"370",&hlt370);

  // Association of trigger names and booleans
  std::map<TString,bool> triggerDecisions;
  triggerDecisions["HltJet30"] = false;
  triggerDecisions["Hlt"+single_dijet+"30"] = false;
  triggerDecisions["Hlt"+single_dijet+"60"] = false;
  triggerDecisions["Hlt"+single_dijet+"80"] = false;
  triggerDecisions["Hlt"+single_dijet+"110"] = false;
  triggerDecisions["Hlt"+single_dijet+"150"] = false;
  triggerDecisions["Hlt"+single_dijet+"190"] = false;
  triggerDecisions["Hlt"+single_dijet+"240"] = false;
  triggerDecisions["Hlt"+single_dijet+"300"] = false;
  triggerDecisions["Hlt"+single_dijet+"370"] = false;

  // Loop over tree entries and fill histograms
  if( nEvts < 0 || nEvts > chain->GetEntries() ) nEvts = chain->GetEntries();
  std::cout << "Reading " << nEvts << " entries...\n";
  for(int n = 0; n < nEvts; ++n) {
    if( n%50000 == 0 ) std::cout << "  " << n << std::endl;
    chain->GetEntry(n);

    if( nObjJet > maxNJet ) {
      std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet << ". Skipping event.\n";
      continue;
    }

    if( nObjJet > 1 ) {
      // Sort jets by corrected pt
      for(int i = 0; i < static_cast<int>(corrJetPt.size()); ++i) {
	if( i < nObjJet ) corrJetPt[i] = jetCorrL1[i]*jetCorrL2L3[i]*jetPt[i];
	else corrJetPt[i] = 0.;      
      }
      std::sort(corrJetPt.begin(),corrJetPt.end());
 
      //      std::cout << corrJetPt[corrJetPt.size()-1] << std::endl;      
      //      std::cout << corrJetPt[corrJetPt.size()-2] << std::endl;      

     
      // Compute average corrected pt
      double ptVar=0;

      if(single_dijet=="DiJetAve"){
      ptVar = 0.5*(corrJetPt[corrJetPt.size()-1]+corrJetPt[corrJetPt.size()-2]);
      }
      else if(single_dijet=="Jet"){
		ptVar = corrJetPt[corrJetPt.size()-1];
	//ptVar = jetPt[0];
      }
      //      corrJetPt.Print();

      if(std::abs(jetEta[0]) < etaMin && std::abs(jetEta[1]) < etaMin ) continue;
      //      if( std::abs(jetEta[0]) > 1.3 || std::abs(jetEta[1]) > 1.3) continue;
	
      triggerDecisions["HltJet30"] = hltphys; 
      triggerDecisions["Hlt"+single_dijet+"30"] = hlt30; 
      triggerDecisions["Hlt"+single_dijet+"60"] = hlt60; 
      triggerDecisions["Hlt"+single_dijet+"80"] = hlt80; 
      triggerDecisions["Hlt"+single_dijet+"110"] = hlt110;
      triggerDecisions["Hlt"+single_dijet+"150"] = hlt150;
      triggerDecisions["Hlt"+single_dijet+"190"] = hlt190;
      triggerDecisions["Hlt"+single_dijet+"240"] = hlt240;
      triggerDecisions["Hlt"+single_dijet+"300"] = hlt300;
      triggerDecisions["Hlt"+single_dijet+"370"] = hlt370;

      // Loop over triggers
      for(size_t i = 0; i < probeTriggers.size(); ++i) {
	if( ptVar > ptMin.at(i) && ptVar < ptMax.at(i) ) {
	  // Count events passing reference (tag) trigger
	  if( triggerDecisions[tagTriggers.at(i)] ) hNTag.at(i)->Fill(ptVar);
	  // Count events passing probe trigger
	  if( triggerDecisions[probeTriggers.at(i)] ) hNProbe.at(i)->Fill(ptVar);
	}
      }
    } // End if( nObjJet > 1 )
  } // End loop over entries


  // Compute and fit efficiencies,
  // determine 99% threshold
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    hEff.at(i)->Divide(hNProbe.at(i),hNTag.at(i),1,1);
    fEff.at(i)->SetParameter(0,1.);
    fEff.at(i)->SetParameter(1,0.5*(ptMin.at(i)+ptMax.at(i)));
    fEff.at(i)->SetParameter(2,1.);
    hEff.at(i)->Fit(fEff.at(i),"NQRI","",ptMin.at(i),ptMax.at(i));

    plateauVals.push_back(fEff.at(i)->GetParameter(2));
    ptVar99.push_back(fEff.at(i)->GetX(0.99*plateauVals.at(i),ptMin.at(i),ptMax.at(i)));
    lPtVar99.at(i)->SetX1(ptVar99.back());
    lPtVar99.at(i)->SetX2(ptVar99.back());
    lPtVar99.at(i)->SetY2(0.99*plateauVals.at(i));
  }


  // Plot efficiencies
  TString outNamePrefix = "TurnOn_Full2011data_"+jetType;
  if( etaMin > 0. ) outNamePrefix += "_MinEta"+util::toTStringNoPoint(etaMin,1);

  for(size_t i = 0; i < probeTriggers.size(); ++i) {

    // Labels
    TPaveText *label = util::LabelFactory::createPaveText(3);
    if( etaMin > 0. ) label->AddText(jetLabel+",  |#eta_{1,2}| > "+util::toTString(etaMin)+",  L = "+util::StyleSettings::luminosity(lumi));
    else label->AddText(jetLabel+",  L = "+util::StyleSettings::luminosity(lumi));
    label->AddText("Reference: "+tagTriggers.at(i));
    label->AddText(""+single_dijet_ptVarLabel+"(#epsilon > 99%) = "+util::toTString(ptVar99.at(i))+" GeV");
    label->GetLine(2)->SetTextColor(lPtVar99.at(i)->GetLineColor());
  
    TH1 *hFrame = util::HistOps::createRatioFrame(ptMin.at(i),ptMax.at(i),0.,1.6*plateauVals.at(i),
						  ""+single_dijet_ptVarLabel+" (GeV)",probeTriggers.at(i)+" Efficiency");
    hFrame->SetLineColor(fEff.at(i)->GetLineColor());
    for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
      hFrame->SetBinContent(bin,plateauVals.at(i));
    }
    
    TCanvas *canEff = new TCanvas("canEff_"+probeTriggers.at(i)+"_Eta"+util::toTString(etaMin),probeTriggers.at(i)+" efficiency",600,600);
    canEff->cd();
    hFrame->Draw();
    hEff.at(i)->Draw("PE1same");
    fEff.at(i)->Draw("same");
    lPtVar99.at(i)->Draw("same");
    label->Draw("same");
    cmsPrel();
    canEff->SaveAs(outNamePrefix+"_"+probeTriggers.at(i)+".eps","eps");
  }

  // All in one plot
  TH1 *hFrame = util::HistOps::createRatioFrame(histMin,histMax,0.,2.45,""+single_dijet_ptVarLabel+" (GeV)","Efficiency");
  for(int bin = 1; bin <= hFrame->GetNbinsX(); ++bin) {
    hFrame->SetBinContent(bin,1.);
  }
  TCanvas *canEff = new TCanvas("canEffAll_Eta"+util::toTString(etaMin),"efficiency",600,600);
  canEff->cd();
  hFrame->Draw();
  for(int i = static_cast<int>(probeTriggers.size()-1); i >= 0; --i) {
    if( plateauVals.at(i) > 0. ) hEff.at(i)->Scale(1./plateauVals.at(i));
    fEff.at(i)->SetParameter(2,1.);

    hEff.at(i)->SetMarkerStyle(20+(i%4));
    hEff.at(i)->SetMarkerColor(util::StyleSettings::color(i));
    hEff.at(i)->SetLineColor(hEff.at(i)->GetMarkerColor());
    fEff.at(i)->SetLineColor(hEff.at(i)->GetMarkerColor());

    hEff.at(i)->Draw("PE1same");
    fEff.at(i)->Draw("same");
  }
  TPaveText *label = util::LabelFactory::createPaveText(1);
  if( etaMin > 0. ) label->AddText(jetLabel+",  |#eta_{1,2}| > "+util::toTString(etaMin)+",  L = "+util::StyleSettings::luminosity(lumi));
  else label->AddText(jetLabel+",  L = "+util::StyleSettings::luminosity(lumi));
  label->Draw("same");
  TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(4,-0.5,1);
  TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(4,0.5,1);
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    if( i < 4 ) {
      leg1->AddEntry(hEff.at(i),probeTriggers.at(i),"LP");
    } else {
      leg2->AddEntry(hEff.at(i),probeTriggers.at(i),"LP");
    }
  }
  leg1->Draw("same");
  leg2->Draw("same");
  cmsPrel();
  canEff->SaveAs(outNamePrefix+".eps","eps");


  //Trigger threshold linear extrapolation plot
  TH1 *hExtrapolFrame = util::HistOps::createFrame(histMin,histMax,histMin,1.5*histMax," Nominal trigger threshold (GeV)","99% efficiency threshold (GeV)");
  TCanvas *canExtrapol = new TCanvas("canExtrapolAll_Eta"+util::toTString(etaMin),"efficiency",600,600);
  canExtrapol->cd();
  hExtrapolFrame->Draw();

  TGraphErrors* extrapol = new TGraphErrors(probeTriggerVals.size(),&probeTriggerVals[0],&ptVar99[0]);
  extrapol->SetMarkerStyle(20);
  extrapol->SetMarkerSize(1.1);
  extrapol->Fit("pol1");

  label = util::LabelFactory::createPaveText(3);
  if( etaMin > 0. ) label->AddText(jetLabel+",  |#eta_{1,2}| > "+util::toTString(etaMin)+",  L = "+util::StyleSettings::luminosity(lumi));
  else label->AddText(jetLabel+",  L = "+util::StyleSettings::luminosity(lumi));
  label->AddText("Lowest trigger: "+tagTriggers.at(0));
  label->AddText("extrapolated "+single_dijet_ptVarLabel+"(#epsilon > 99%) = "+util::toTString(extrapol->GetFunction("pol1")->Eval(tagTriggerVals.at(0)),1)+" GeV");
  
  extrapol->Draw("P");
  label->Draw("same");
  cmsPrel();
  canExtrapol->SaveAs(outNamePrefix+"Extrapol.eps","eps");


  // Print efficiencies
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    std::cout << "\\texttt{" << probeTriggers.at(i) << "} & " << util::toTString(ptVar99.at(i),1) << " \\\\"  << std::endl;
  }

  std::cout << std::endl << std::endl;
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
    std::cout << probeTriggers.at(i) << " : " << util::toTString(ptVar99.at(i),1) << std::endl;
  }


  TFile *outf = new TFile("TriggerEfficiencies_"+outNamePrefix+"_"+single_dijet+".root","RECREATE");
  for(size_t i = 0; i < probeTriggers.size(); ++i) {
  hNTag.at(i)->Write();
  hNProbe.at(i)->Write();
  hEff.at(i)->Write();
  fEff.at(i)->Write();
  }
  outf->Close();

}


// --------------------------------------------------
void run(int nEvts = -1) {
  TString input = "/afs/naf.desy.de/user/k/kirschen/scratch/Kalibri2/L2andJERScripts/filelist_Full2011_L1FastJet_AK5PF";

  double lumi = 4965;
  
  single_dijet="DiJetAve";
  //  single_dijet="Jet";
  TriggerYearLabel="2011";

  if(single_dijet=="DiJetAve"){
    single_dijet_ptVarLabel = "p^{ave}_{T}";
  }
  else if(single_dijet=="Jet"){
    single_dijet_ptVarLabel = "p^{lead}_{T}";
  }
  
  computeTriggerEfficiences(input,lumi,nEvts);
//   computeTriggerEfficiences(input,lumi,nEvts,1.7);
//   computeTriggerEfficiences(input,lumi,nEvts,2.3);


}
