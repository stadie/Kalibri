// $Id: ControlPlotsJetSmearing.cc,v 1.3 2009/07/16 14:49:10 mschrode Exp $

#include "ControlPlotsJetSmearing.h"

#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TVector2.h"

#include "SmearData.h"
#include "SmearDiJet.h"
#include "SmearPhotonJet.h"




// --------------------------------------------------
ControlPlotsJetSmearing::ControlPlotsJetSmearing(const std::string& configfile, const std::vector<TData*> * data, TParameters * param)
  : data_(data), config_(new ConfigFile(configfile.c_str())), param_(param),
    respNBins_(150), respMin_(0.), respMax_(6.),
    dir_("./controlplots")
{
  setGStyle();
}



//!  \brief Draw response control plots for events
//!         of type \p SmearData
// --------------------------------------------------
void ControlPlotsJetSmearing::plotResponse() const
{
  std::cout << "Creating response control plots..." << std::endl;

  // Create histograms of response
  TH1D* hRespMeas = new TH1D("hRespMeas",";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",respNBins_,respMin_,respMax_);
  hRespMeas->Sumw2();

  TH1D* hRespFit = new TH1D("hRespFit",";R = p^{jet}_{T} / p^{true}_{T};1/(Nw)  dN / dR",5*respNBins_,respMin_,respMax_);
  hRespFit->SetLineColor(2);
  hRespFit->SetLineWidth(2);
  hRespFit->Sumw2();

  TH1D* hRespFitStep = new TH1D("hRespFitStep",";R = p^{jet}_{T} / p^{true}_{T};1 / (Nw)  dN / dR",
				config_->read<int>("Response pdf nsteps",10),
				config_->read<double>("Response pdf min",0.),
				config_->read<double>("Response pdf max",1.8));
  hRespFitStep->Sumw2();
  hRespFitStep->SetLineColor(9);
  hRespFitStep->SetLineWidth(2);

  TH1D* hRespFitGaus = new TH1D("hRespFitGaus",";R = p^{jet}_{T} / p^{true}_{T};1 / (Nw)  dN / dR",5*respNBins_,respMin_,respMax_);
  hRespFitGaus->Sumw2();
  hRespFitGaus->SetLineColor(8);
  hRespFitGaus->SetLineWidth(2);

  TH1D* hRespFitSum = new TH1D("hRespFitSum",";R = p^{jet}_{T} / p^{true}_{T};1 / (Nw)  dN / dR",5*respNBins_,respMin_,respMax_);
  hRespFitSum->Sumw2();
  hRespFitSum->SetLineColor(1);
  hRespFitSum->SetLineWidth(2);


  // Create histograms of genjet pdf
  double mingpt = config_->read<double>("DiJet integration min truth",0);
  double maxgpt = config_->read<double>("DiJet integration max truth",1);

  TH1D * hDijetGenJetPt     = new TH1D("hDijetGenJetPt",
				       ";p^{true}_{T} (GeV);1 / (Nw)  dN / dp^{true}_{T}",
				       25,0.9*mingpt,1.1*maxgpt);
  hDijetGenJetPt->GetXaxis()->SetNdivisions(505);
  hDijetGenJetPt->Sumw2();

  TH1D * hDijetTruthPDF     = new TH1D("hDijetTruthPDF",
				       ";p^{true}_{T} (GeV);1 / (Nw)  dN / dp^{true}_{T}",
				       5*respNBins_,0.9*mingpt,1.1*maxgpt);
  hDijetTruthPDF->SetLineColor(2);
  hDijetTruthPDF->SetLineWidth(2);

  // Fill histogram of measured response
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
      for(int i = 0; i < 2; i++) {        // Loop over both jets
	TJet * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 ) jet = static_cast<TJet*>(dijet->GetSecondMess());
	hDijetGenJetPt->Fill(jet->genPt);
	hRespMeas->Fill(jet->pt / jet->genPt);
      }
    }

    // Select PhotonJet events
    if( (*datait)->GetType() == TypeSmearPhotonJet ) {
      SmearPhotonJet * pjet = static_cast<SmearPhotonJet*>(*datait);  
      double photonpt = pjet->GetTruth();
      double jetpt    = pjet->GetMess()->pt;
      hRespMeas->Fill( jetpt / photonpt );		    
    }
  } // End of loop over data

  if( hRespMeas->Integral("width") ) hRespMeas->Scale(1./hRespMeas->Integral("width"));
  if( hDijetGenJetPt->Integral("width") ) hDijetGenJetPt->Scale(1./hDijetGenJetPt->Integral("width"));


  // Fill histogram of fitted response
  SmearData * smeardata = dynamic_cast<SmearData*>(data_->front());
  if( smeardata ) {
    // Interpolated fit function
    for(int bin = 1; bin <= hRespFit->GetNbinsX(); bin++) {
      double r = hRespFit->GetBinCenter(bin);
      hRespFit->SetBinContent(bin,smeardata->RespPDF(r));
    }

    // In case of step + gauss parametrizations
    std::string param = config_->read<std::string>("Parametrization Class","");
    if( param == "SmearParametrizationStepGauss" || param == "SmearParametrizationStepGaussInter" ) {
      std::vector<double> scale = bag_of<double>(config_->read<string>("Jet parameter scales",""));

      // Step part of fit function
      for(int bin = 1; bin <= hRespFitStep->GetNbinsX(); bin++) {
	double val  = scale.at(bin+2)*(smeardata->GetRespPar()[bin+2]);
	hRespFitStep->SetBinContent(bin,val);
      }
      hRespFitStep->Scale(1./hRespFitStep->Integral("width"));
      hRespFitStep->Scale(1. - scale.at(0)*(smeardata->GetRespPar()[0]));
      
      // Gauss part of fit function
      for(int bin = 1; bin <= hRespFitGaus->GetNbinsX(); bin++) {
	double c     = scale.at(0)*(smeardata->GetRespPar()[0]);
	double mu    = scale.at(1)*(smeardata->GetRespPar()[1]);
	double sigma = scale.at(2)*(smeardata->GetRespPar()[2]);
	double r     = hRespFitGaus->GetBinCenter(bin);
	double val   = c * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
	hRespFitGaus->SetBinContent(bin,val);
      }
      
      // Sum
      for(int binGaus = 1; binGaus <= hRespFitGaus->GetNbinsX(); binGaus++) {
	int    binStep = hRespFitStep->FindBin(hRespFitGaus->GetBinCenter(binGaus));
	double val     = hRespFitStep->GetBinContent(binStep) + hRespFitGaus->GetBinContent(binGaus);
	hRespFitSum->SetBinContent(binGaus,val);
      }
    }

    // In case of two gauss parametrizations
    else if( param == "SmearParametrizationTwoGauss" ) {
    
      std::vector<double> scale = bag_of<double>(config_->read<string>("Jet parameter scales",""));
      
      // Central Gauss of fit function
      for(int bin = 1; bin <= hRespFitGaus->GetNbinsX(); bin++) {
	double c     = scale.at(0)*(smeardata->GetRespPar()[0]);
	double mu    = scale.at(1)*(smeardata->GetRespPar()[1]);
	double sigma = scale.at(2)*(smeardata->GetRespPar()[2]);
	double r     = hRespFitGaus->GetBinCenter(bin);
	double val   = c * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
	hRespFitGaus->SetBinContent(bin,val);
      }
      
      // Tail Gauss of fit function
      for(int bin = 1; bin <= hRespFitSum->GetNbinsX(); bin++) {
	double c     = scale.at(0)*(smeardata->GetRespPar()[0]);
	double mu    = scale.at(3)*(smeardata->GetRespPar()[3]);
	double sigma = scale.at(4)*(smeardata->GetRespPar()[4]);
	double r     = hRespFitGaus->GetBinCenter(bin);
	double val   = (1.-c) * exp( -pow((mu-r)/sigma,2) / 2. ) / sqrt(2.*M_PI) / sigma;
	hRespFitSum->SetBinContent(bin,val);
      }
    }
  }

  // Fill histogram of assumed dijet truth pdf
  // and fit with 1/x^n function
  double n      = 6.5;
  DataIt datait = data_->begin();
  while( (*datait)->GetType() != TypeSmearDiJet  &&  datait != data_->end() ) datait++;
  if( datait != data_->end() ) {
    SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  
    for(int bin = 1; bin <= hDijetTruthPDF->GetNbinsX(); bin++) {
      double t = hDijetTruthPDF->GetBinCenter(bin);
      hDijetTruthPDF->SetBinContent(bin,dijet->TruthPDF(t));
    }

    n            = (dijet->GetTruthPar()[0]);
    double k     = n - 3.;
    double m     = n - 1.;
    double norm1 = k / ( pow(mingpt,-k) - pow(maxgpt,-k) );
    double norm2 = m / ( pow(mingpt,-m) - pow(maxgpt,-m) );
    hDijetTruthPDF->Scale(norm2/norm1);

    //   double norm = ( pow(maxgpt,3) + pow(mingpt,3) ) / 3.;
    //   hDijetTruthPDF->Scale(norm/(maxgpt-mingpt));
  }

  // Find populated x-axis range
  int maxBin = 0;
  int minBin = 1;
  for(int bin = 1; bin <= hRespMeas->GetNbinsX(); bin++) {
    if( hRespMeas->GetBinContent(bin) > 0  ) maxBin = bin;
    if( minBin > maxBin ) minBin = bin;
  }
  if( maxBin < hRespMeas->GetNbinsX() ) maxBin++;
  if( minBin > 1 )                      minBin--;
  hRespMeas->GetXaxis()->SetRange(minBin,maxBin);


  // Write histos to eps file
  TCanvas *c1 = new TCanvas("c1","Jet Response",0,0,600,600);
  c1->cd();
  hRespMeas->Draw();
  hRespFit->Draw("Lsame");
  c1->SetLogy();
  c1->SaveAs((dir_+"/JetResponse.eps").c_str());

  c1->cd();
  hRespMeas->Draw();
  std::string param = config_->read<std::string>("Parametrization Class","");
  if( param == "SmearParametrizationStepGauss" || param == "SmearParametrizationStepGaussInter" ) {
    hRespFitStep->Draw("same");
    hRespFitGaus->Draw("same");
    hRespFitSum->Draw("same");
  } else if( param == "SmearParametrizationTwoGauss" ) {
    hRespFitGaus->Draw("same");
    hRespFitSum->Draw("same");
  }
  hRespFit->Draw("Lsame");
  c1->SetLogy();
  c1->SaveAs((dir_+"/JetResponseDetail.eps").c_str());

  c1->cd();
  hDijetGenJetPt->Draw();
  hDijetTruthPDF->Draw("SAME");

  TLegend *leg = new TLegend(0.4,0.75,0.89,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  char entry[50];
  sprintf(entry,"Fit #propto 1 / (p^{true}_{T})^{%.1f}",n);
  leg->AddEntry(hDijetTruthPDF,entry,"L");
  leg->Draw("same");

  c1->SetLogy();
  c1->SaveAs((dir_+"/DiJetTruth.eps").c_str());

  // Write histos to root file
  TFile rootfile((dir_+"/JetResponse.root").c_str(),"RECREATE");
  rootfile.WriteTObject(hRespMeas);
  rootfile.WriteTObject(hRespFit);
  rootfile.WriteTObject(hRespFitStep);
  rootfile.WriteTObject(hRespFitGaus);
  rootfile.WriteTObject(hRespFitSum);
  rootfile.WriteTObject(hDijetGenJetPt);
  rootfile.WriteTObject(hDijetTruthPDF);
  rootfile.Close();

  // Clean up
  delete c1;
  delete hRespMeas;
  delete hRespFit;
  delete hRespFitStep;
  delete hRespFitGaus;
  delete hRespFitSum;
  delete hDijetGenJetPt;
  delete hDijetTruthPDF;
  delete leg;

  std::cout << "Done" << std::endl;
}



// --------------------------------------------------
void ControlPlotsJetSmearing::plotParameterScan(const std::vector<unsigned int>& pars) const {
  std::cout << "Creating parameter scan control plots..." << std::endl;

  // Create one histogram of likelihood per parameter
  std::vector<TH1D*> hLikelihood;
  for(int i = 0; i < param_->GetNumberOfParameters(); i++) {
    char name[50];
    sprintf(name,"hNLogL%i",i);
    char title[50];
    sprintf(title,";Parameter %i;-ln(L)",i);
    double min = param_->GetPars()[i] - 0.5;
    if( min < 0. ) min = 0.;
    double max = param_->GetPars()[i] + 0.5;
    TH1D * h = new TH1D(name,title,10,min,max);
    h->SetMarkerStyle(24);
    hLikelihood.push_back(h);
  }

  // Do parameter scan
  for(unsigned int i = 0; i < pars.size(); i++) {
    int idx = pars.at(i);
    if( idx < param_->GetNumberOfParameters() ) {
      std::cout << "  Par " << idx << "... " << std::flush;
      double oldpar = param_->GetPars()[idx];
      TH1D * h = hLikelihood.at(idx);
      for(int bin = 1; bin <= h->GetNbinsX(); bin++) {
	param_->GetPars()[idx] = h->GetBinCenter(bin);
	double fsum = 0.;
	for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
	  fsum += (*datait)->chi2();
	}
	h->SetBinContent(bin,fsum);
      }
      param_->GetPars()[idx] = oldpar;
      std::cout << "ok" << std::endl;
    }
  }

  // Plot histograms
  TPostScript * const ps = new TPostScript((dir_+"/ParameterScan.ps").c_str(),111);
  TCanvas *c1 = new TCanvas("c1","Parameter Scan",0,0,600,600);
  ps->NewPage();
  c1->cd();
  for(unsigned int i = 0; i < pars.size(); i++) {
    int idx = pars.at(i);
    if( idx < param_->GetNumberOfParameters() ) {
      TH1D * h = hLikelihood.at(idx);
      h->Draw("P");
      TLine * line = new TLine(param_->GetPars()[idx],h->GetMinimum(),
			       param_->GetPars()[idx],h->GetMaximum());
      line->SetLineStyle(2);
      line->SetLineWidth(2);
      line->Draw("same");
      c1->Draw();
      ps->NewPage();
      delete line;
    }
  }
  ps->Close();

  // Clean up
  for(std::vector<TH1D*>::iterator it = hLikelihood.begin(); it != hLikelihood.end(); it++) {
    delete *it;
  }
  hLikelihood.clear();
  delete c1;
  delete ps;

  std::cout << "Done" << std::endl;
}



//!  \brief Draw control plots for events
//!         of type \p TypeSmearDiJet
// --------------------------------------------------
void ControlPlotsJetSmearing::plotDijets() const
{
  std::cout << "Creating dijet control plots..." << std::endl;

  // --- Create histograms --------------------------

  // Find pt ranges
  double minGenJetPt = 10000.;
  double maxGenJetPt = 0.;
  double minCalJetPt = 10000.;
  double maxCalJetPt = 0.;
  double minDijetPt  = 10000.;
  double maxDijetPt  = 0.;
  double min3rdJetPt = 10000.;
  double max3rdJetPt = 0.;
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    // Select DiJet events
    if( (*datait)->GetType() == TypeSmearDiJet )  {
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      // Loop over both jets
      for(int i = 0; i < 2; i++) {        
	TJet         * jet = static_cast<TJet*>(dijet->GetMess());
	if( i == 1 )   jet = static_cast<TJet*>(dijet->GetSecondMess());

	if( jet->genPt < minGenJetPt ) minGenJetPt = jet->genPt;
	if( jet->genPt > maxGenJetPt ) maxGenJetPt = jet->genPt;

	if( jet->pt < minCalJetPt ) minCalJetPt = jet->pt;
	if( jet->pt > maxCalJetPt ) maxCalJetPt = jet->pt;
      }

      if( dijet->dijetPt() < minDijetPt ) minDijetPt = dijet->dijetPt();
      if( dijet->dijetPt() > maxDijetPt ) maxDijetPt = dijet->dijetPt();

      TJet * jet3 = static_cast<TJet*>(dijet->GetThirdMess());
      if( jet3->pt < min3rdJetPt ) min3rdJetPt = jet3->pt;
      if( jet3->pt > max3rdJetPt ) max3rdJetPt = jet3->pt;
    }
  }

  // Pt distributions
  std::vector<TH1F*> hGenJetPt;
  std::vector<TH1F*> hCalJetPt;
  std::string genJetNames[3] = { "hGenJetPtBothJets",
				 "hGenJetPtJet1", 
				 "hGenJetPtJet2" };
  std::string calJetNames[3] = { "hCalJetPtBothJets",
				 "hCalJetPtJet1", 
				 "hCalJetPtJet2" };
  int color[3] = { 1, 2, 4 };
  for(int i = 0; i < 3; i++) {
    TH1F * h = 0;

    h = new TH1F(genJetNames[i].c_str(),";p^{gen}_{T} (GeV);1 / N  dN / dp^{gen}_{T}  1 / (GeV)",
		 50,0.9*minGenJetPt,1.1*maxGenJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
    hGenJetPt.push_back(h);

    h = new TH1F(calJetNames[i].c_str(),";p^{jet}_{T} (GeV);1 / N  dN / dp^{jet}_{T}  1 / (GeV)",
		 50,0.9*minCalJetPt,1.1*maxCalJetPt);
    h->Sumw2();
    h->SetLineWidth(2);
    h->SetLineColor(color[i]);
    hCalJetPt.push_back(h);
  }

  TH2F * hCalJet2vsCalJet1Pt = new TH2F("hCalJet2vsCalJet1Pt",
					";p^{jet1}_{T} (GeV);p^{jet2}_{T} (GeV)",
					50,0.9*minCalJetPt,1.1*maxCalJetPt,
					50,0.9*minCalJetPt,1.1*maxCalJetPt);

  TH1F * hDijetPt = new TH1F("hDijetPt",";p^{dijet}_{T} (GeV);1 / N  dN / dp^{dijet}_{T}  1 / (GeV)",
			     50,0.9*minDijetPt, 1.1*maxDijetPt);
  hDijetPt->SetLineWidth(2);
  hDijetPt->Sumw2();

  TH1F * h3rdJetPt = new TH1F("h3rdJetPt",";p^{jet3}_{T} (GeV);1 / N  dN / dp^{jet3}_{T}  1 / (GeV)",
			     50,0.9*min3rdJetPt, 1.1*max3rdJetPt);
  h3rdJetPt->SetLineWidth(2);
  h3rdJetPt->Sumw2();

  TH2F * h3rdJetvsDijetPt = new TH2F("h3rdJetvsDijetPt",
				     ";p^{dijet}_{T} (GeV);p^{jet3}_{T} (GeV)",
				     50,0.9*minDijetPt, 1.1*maxDijetPt,
				     50,0.9*min3rdJetPt, 1.1*max3rdJetPt);

  TH1F * hRel3rdJetPt = new TH1F("hRel3rdJetPt",
				 ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};1 / N  dN / dp^{jet3}_{T,rel}",
			     50,0,1.4);
  hRel3rdJetPt->SetLineWidth(2);
  hRel3rdJetPt->Sumw2();

  TH1F * hDeltaPhi = new TH1F("hDeltaPhi",";#Delta#Phi;1 / N  dN / d#Delta#Phi",
			      25,1.5,M_PI);
  hDeltaPhi->SetLineWidth(2);
  hDeltaPhi->Sumw2();


  // Response correlations
  TH2F * hRvsDeltaPhi = new TH2F("hRvsDeltaPhi",";#Delta#Phi;p^{jet}_{T} / p^{gen}_{T}",
				 25,1.5,M_PI,25,0,1.4);
  TH2F * hRvsRel3rdJetPt = new TH2F("hRvsRel3rdJetPt",
				    ";p^{jet3}_{T,rel} = p^{jet3}_{T} / p^{dijet}_{T};p^{jet}_{T} / p^{gen}_{T}",
				    50,0,1.4,25,0,1.4);
  TH2F * hRvsEMF = new TH2F("hRvsEMF",";EMF;p^{jet}_{T} / p^{gen}_{T}",
			    25,0,1,25,0,1.4);



  // --- Fill histograms ----------------------------
  for(DataIt datait = data_->begin(); datait != data_->end(); datait++) {
    if( (*datait)->GetType() == TypeSmearDiJet )  { // Select DiJet events
      SmearDiJet * dijet = static_cast<SmearDiJet*>(*datait);  

      TJet * jet1 = static_cast<TJet*>(dijet->GetMess());
      TJet * jet2 = static_cast<TJet*>(dijet->GetSecondMess());
      TJet * jet3 = static_cast<TJet*>(dijet->GetThirdMess());

      double dPhi = std::abs(TVector2::Phi_mpi_pi( jet1->phi - jet2->phi ));
      hDeltaPhi->Fill(dPhi);

      hRvsDeltaPhi->Fill( dPhi, jet1->pt / jet1->genPt );
      hRvsDeltaPhi->Fill( dPhi, jet2->pt / jet2->genPt );

      hRvsEMF->Fill( jet1->EMF / jet1->pt, jet1->pt / jet1->genPt);
      hRvsEMF->Fill( jet2->EMF / jet2->pt, jet2->pt / jet2->genPt);

      hRvsRel3rdJetPt->Fill( jet3->pt / dijet->dijetPt(), jet1->pt / jet1->genPt );
      hRvsRel3rdJetPt->Fill( jet3->pt / dijet->dijetPt(), jet2->pt / jet2->genPt );

      hGenJetPt.at(0)->Fill( jet1->genPt );
      hGenJetPt.at(0)->Fill( jet2->genPt );
      hGenJetPt.at(1)->Fill( jet1->genPt );
      hGenJetPt.at(2)->Fill( jet2->genPt );

      hCalJetPt.at(0)->Fill( jet1->pt );
      hCalJetPt.at(0)->Fill( jet2->pt );
      hCalJetPt.at(1)->Fill( jet1->pt );
      hCalJetPt.at(2)->Fill( jet2->pt );

      hCalJet2vsCalJet1Pt->Fill( jet1->pt, jet2->pt );
      hDijetPt->Fill( dijet->dijetPt() );

      h3rdJetPt->Fill( jet3->pt );
      h3rdJetvsDijetPt->Fill( dijet->dijetPt(), jet3->pt);
      hRel3rdJetPt->Fill( jet3->pt / dijet->dijetPt() );
    }
  }


  // Normalizing histograms
  for(size_t i = 0; i < hGenJetPt.size(); i++) {
    normHist( hGenJetPt.at(i) );
    normHist( hCalJetPt.at(i) );
  }
  normHist( hDijetPt );
  normHist( h3rdJetPt );
  normHist( hRel3rdJetPt );
  normHist( hDeltaPhi );



  // --- Plot histograms ----------------------------
  TPostScript * const ps = new TPostScript((dir_+"/Dijets.ps").c_str(),111);
  TCanvas           * c1 = new TCanvas("c1","Dijets",0,0,600,600);

  TLegend * leg = new TLegend(0.25,0.75,0.45,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(hGenJetPt.at(1),"Jet 1","L");
  leg->AddEntry(hGenJetPt.at(2),"Jet 2","L");

  std::vector<TObject*> objs;
  for(size_t i = 1; i < hGenJetPt.size(); i++) {
    objs.push_back(hGenJetPt.at(i));
  }
  objs.push_back(leg);
  drawPSPage(ps,c1,hGenJetPt.at(0),"",true);
  drawPSPage(ps,c1,objs,"",true);

  objs.clear();
  for(size_t i = 1; i < hCalJetPt.size(); i++) {
    objs.push_back(hCalJetPt.at(i));
  }
  objs.push_back(leg);
  drawPSPage(ps,c1,hCalJetPt.at(0),"",true);
  drawPSPage(ps,c1,objs,"",true);

  drawPSPage(ps,c1,hCalJet2vsCalJet1Pt,"BOX",false);
  drawPSPage(ps,c1,hDijetPt,"",true);
  drawPSPage(ps,c1,h3rdJetPt,"",true);
  drawPSPage(ps,c1,h3rdJetvsDijetPt,"BOX",false);
  drawPSPage(ps,c1,hRel3rdJetPt,"",true);
  drawPSPage(ps,c1,hDeltaPhi,"",true);
  drawPSPage(ps,c1,hRvsRel3rdJetPt,"BOX",false);
  drawPSPage(ps,c1,hRvsDeltaPhi,"BOX",false);
  drawPSPage(ps,c1,hRvsEMF,"BOX",false);

  ps->Close();

  // Clean up
  for(size_t i = 0; i < hGenJetPt.size(); i++) {
    delete hGenJetPt.at(i);
    delete hCalJetPt.at(i);
  }
  delete leg;
  delete hCalJet2vsCalJet1Pt;
  delete hDijetPt;
  delete h3rdJetPt;
  delete h3rdJetvsDijetPt;
  delete hRel3rdJetPt;
  delete hDeltaPhi;
  delete hRvsRel3rdJetPt;
  delete hRvsDeltaPhi;
  delete hRvsEMF;
  delete c1;
  delete ps;

  std::cout << "Done" << std::endl;
}



//!  \brief Set default gStyle options
// --------------------------------------------------
void ControlPlotsJetSmearing::setGStyle() const
{
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  // For the frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the histo:
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);

  // For the statistics box:
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("neMR");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.92);              
  gStyle->SetStatY(0.86);              
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.22);

  // For the legend
  gStyle->SetLegendBorderSize(1);

  //  Margins
  // -------------------------------------------
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadRightMargin(0.04);

  // For the Global title:
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.58);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleXOffset(0);
  gStyle->SetTitleYOffset(0);
  gStyle->SetTitleBorderSize(0);

  // For the axis labels:
  //  For the axis labels and titles
  // -------------------------------------------
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetLabelColor(1,"XYZ");
  // For the axis labels:
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  
  // For the axis titles:
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.5);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}



//!  \brief Draw TObject on one page of a ps file
//!
//!  Opens a new page in the PostScript file \p ps
//!  and draws the TObject \p obj on it.
//!
//!  \param ps     A new page is added to this file
//!                containing the TObject \p obj
//!  \param can    The \p obj is drawn on this canvas
//!  \param obj    The TObject to be drawn
//!  \param option Draw options
//!  \param logy   Sets log-scale on y-axis if true
// --------------------------------------------------
void ControlPlotsJetSmearing::drawPSPage(TPostScript * ps, TCanvas * can, TObject * obj, std::string option, bool logy) const {
  std::vector<TObject*> objs;
  objs.push_back(obj);
  drawPSPage(ps,can,objs,option,logy);
}



//!  \brief Draw TObjects on one page of a ps file
//!
//!  Opens a new page in the PostScript file \p ps
//!  and draws all TObjects in \p obs on it. If \p objs
//!  contains more than one TObject, the draw option
//!  "same" is automatically used.
//!
//!  \param ps     A new page is added to this file
//!                containing all TObjects in \p objs
//!  \param can    The \p objs are drawn on this canvas
//!  \param objs   Contains the objects that are drawn
//!  \param option Draw options (except for "same",
//!                see above)
//!  \param logy   Sets log-scale on y-axis if true
// --------------------------------------------------
void ControlPlotsJetSmearing::drawPSPage(TPostScript * ps, TCanvas * can, std::vector<TObject*> objs, std::string option, bool logy) const {
  ps->NewPage();
  can->cd();
  for( size_t i = 0; i < objs.size(); i++ ) {
    if( i == 0 ) objs.at(i)->Draw(option.c_str());
    else         objs.at(i)->Draw((option.append("same")).c_str());
  }
  if( logy ) can->SetLogy(1);
  else       can->SetLogy(0);
  can->Draw();
}
