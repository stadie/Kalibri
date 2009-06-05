//!
//!  Compare L2 and L3 correction factors from
//!  Kalibri with official JetMET factors.
//!
//!  \author Matthias Schroeder
//!  \date   Thu Jun  4 14:35:41 CEST 2009
//!
//!  $Id: $
// -------------------------------------------------------------

#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TPaveText.h"
#include "TPostScript.h"
#include "TROOT.h"
#include "TStyle.h"


typedef std::map<unsigned int, std::vector<double> > parmap;
typedef std::map<unsigned int, std::vector<double> >::const_iterator parmapIt;


//!  Read parameters from an ascii file with the
//!  following syntax:
//!
//!   EtaMin EtaMax NumPar Par1 Par2 ... ParNumPar
parmap readParameters(const std::string& filename) {
  parmap corr;
  std::ifstream file;
  file.open(filename.c_str());
  unsigned int line = 0;
  double       val  = -1.;
  while( !file.eof() ) {
    std::vector<double> par;        // Store parameters
    for(int i = 0; i < 2; i++) {    // Loop over eta edges of bin
      file >> val;
      par.push_back(val);
    }
    int n = 0;
    file >> n;                      // Number of parameters
    for(int i = 0; i < n; i++) {    // Loop over parameters
      file >> val;
      par.push_back(val);
    }
    if( n != 0 ) corr[line] = par;  // Some weird last line is read otherwise
    line++;
  }
  file.close();

  return corr;
}


double corrL2(double * x, double * par) {
  double pt    = (x[0] < 4.0) ? 4.0 : (x[0] > 2000.0) ? 2000.0 : x[0]; 
  double logpt = log10(pt);
  double c1    = par[0] + logpt * ( par[1] + logpt * par[2] );

  return c1;
}


double corrL3(double * x, double * par) {
  double pt    = (x[0] < 4.0) ? 4.0 : (x[0] > 2000.0) ? 2000.0 : x[0];
  double logpt = log10(pt);
  double c2    = par[0] + par[1]/(pow(logpt,par[2]) + par[3]);

  return  c2;
}


int compareL2L3Correction(const std::string& kalibri, const std::string& jetMETL2, const std::string& jetMETL3) {
//int compareL2L3Correction(const std::string& jetMETL2, const std::string& jetMETL3) {
  // -- Read calibration constants ---------------
  std::cout << "Reading calibration constants...\n";
  parmap corrKalibri  = readParameters(kalibri);
  parmap corrJetMETL2 = readParameters(jetMETL2);
  parmap corrJetMETL3 = readParameters(jetMETL3);
  if( (1 + corrKalibri.size() ) != ( corrJetMETL2.size() + corrJetMETL3.size() ) ) {
    std::cerr << "Linker Kompaktierer defekt. Kontaktieren Sie das Personal (68/111)." << std::endl;
    return 1;
  }


  // -- Global variables -------------------------
  std::vector<TF1*> fCorrKalibri(1+corrKalibri.size());
  std::vector<TF1*> fCorrJetMET(1+corrKalibri.size());
  TF1 * f = 0;
  char name[50];
  std::string outFileName = "comparisonL2L3_Summer08_IC5Calo.ps";
  bool printPar = false;


  // -- Create correction plots for JetMET -------
  std::cout << "Creating JetMET plots...\n";
  // L3 correction
  if( printPar ) std::cout << "  JetMETL3: " << std::flush;
  std::vector<double>& par = corrJetMETL3[0];
  f = new TF1("fJetMETL3",corrL3,par.at(2),par.at(3),4);
  for(int i = 0; i < 4; i++) {
    f->SetParameter(i,par.at(4+i));
    if( printPar ) std::cout << par.at(4+i) << " ";
  }
  if( printPar ) std::cout << "\n";
  fCorrJetMET.at(0) = f;

  // L2 correction
  for(unsigned int bin = 0; bin < corrJetMETL2.size(); bin++) {
    if( printPar ) std::cout << "  JetMETL2 " << bin << ": " << std::flush;
    par = corrJetMETL2[bin];
    sprintf(name,"fJetMETL2_%i",bin);
    f = new TF1(name,corrL2,par.at(2),par.at(3),3);
    sprintf(name,"Bin %i: %.3f <  #eta < %.3f",bin,par.at(0),par.at(1));
    f->SetTitle(name);
    for(int i = 0; i < 3; i++) {
      f->SetParameter(i,par.at(4+i));
      if( printPar ) std::cout << par.at(4+i) << " ";
    }
    if( printPar ) std::cout << "\n";
    fCorrJetMET.at(1+bin) = f;
  }

  // Nice plots
  for(unsigned i = 0; i < fCorrJetMET.size(); i++) {
    fCorrJetMET.at(i)->SetLineWidth(2);
    fCorrJetMET.at(i)->SetLineColor(4);
  }


  // -- Create correction plots for Kalibri ------
  std::cout << "Creating Kalibri plots...\n";
  int nTowerPar     = 0;
  int nJetPar       = 3;
  int nTrackPar     = 0;
  int nGlobalJetPar = 4;
  // Global correction
  if( printPar ) std::cout << "  Kalibri: " << std::flush;
  par = corrKalibri[0];
  f = new TF1("fKalibri_global",corrL3,par.at(2),par.at(3),4);
  for(int i = 0; i < nGlobalJetPar; i++) {
    f->SetParameter(i,par.at(4+nTowerPar+nJetPar+nTrackPar+i));
    if( printPar ) std::cout << par.at(4+nTowerPar+nJetPar+nTrackPar+i) << " ";
  }
//   f = new TF1("fKalibri_global",corrL3,4,2000,4);
//   f->SetParameter(0,0.998293);
//   f->SetParameter(1,5.43056);
//   f->SetParameter(2,3.3444);
//   f->SetParameter(3,2.39809);
  if( printPar ) std::cout << "\n";
  fCorrKalibri.at(0) = f;

  // Local correction
  double scale[3] = { 1., 0.1, 0.01 };
  for(unsigned int bin = 0; bin < corrKalibri.size(); bin++) {
    if( printPar ) std::cout << "  Kalibri " << bin << ": " << std::flush;
    par = corrKalibri[bin];
    sprintf(name,"fKalibri_%i",bin);
    f = new TF1(name,corrL2,par.at(2),par.at(3),3);
    sprintf(name,"Bin %i: %.3f <  #eta < %.3f",bin,par.at(0),par.at(1));
    f->SetTitle(name);
    for(int i = 0; i < nJetPar; i++) {
      f->SetParameter(i,scale[i]*par.at(4+nTowerPar+i));
      if( printPar ) std::cout << scale[i]*par.at(4+nTowerPar+i) << " ";
    }
    if( printPar ) std::cout << "\nDone";
    fCorrKalibri.at(1+bin) = f;
  } 

  // Nice plots
  for(unsigned i = 0; i < fCorrKalibri.size(); i++) {
    fCorrKalibri.at(i)->SetLineWidth(2);
    fCorrKalibri.at(i)->SetLineColor(2);
  }


  // -- Writing plots to file --------------------
  std::cout << "Writing plots to file...\n";
  gROOT->ProcessLine(".x ./gStyleSettings.C");
  gStyle->SetOptStat(0);
  TPostScript * const ps = new TPostScript(outFileName.c_str(),111);
  TCanvas     * const c1 = new TCanvas("c1","",500,500);
  TH1F          * hFrame = new TH1F("hFrame",";p_{T} (GeV)",100,0,3000);
  c1->cd();
  c1->SetLogx(1);
  for(unsigned i = 0; i < fCorrKalibri.size(); i++) {
    if( i == 0 ) {
      hFrame->GetYaxis()->SetRangeUser(1.0,5.0);
      hFrame->GetYaxis()->SetTitle("L3 correction factor");
    }
    else {
      hFrame->GetYaxis()->SetRangeUser(0.5,2.0);
      hFrame->GetYaxis()->SetTitle("L2 correction factor");
    }
    hFrame->SetTitle(fCorrKalibri.at(i)->GetTitle());
    hFrame->Draw();
    fCorrKalibri.at(i)->Draw("same");
    fCorrJetMET.at(i)->Draw("same");

    TPaveText * t[3];
    double maxy = 0.86;
    double miny = maxy - 0.05*(1+fCorrKalibri.at(i)->GetNpar());
    t[0] = new TPaveText(0.3,miny,0.38,maxy,"NDC");
    sprintf(name,"   ");
    t[0]->AddText(name);
    t[1] = new TPaveText(0.38,miny,0.58,maxy,"NDC");
    sprintf(name,"#color[2]{Kalibri}");
    t[1]->AddText(name);
    t[2] = new TPaveText(0.58,miny,0.78,maxy,"NDC");
    sprintf(name,"#color[4]{JetMET}");
    t[2]->AddText(name);
    for(int p = 0; p < fCorrKalibri.at(i)->GetNpar(); p++) {
      sprintf(name,"b_{%i}:",p);
      t[0]->AddText(name);
      sprintf(name,"%.4f",fCorrKalibri.at(i)->GetParameter(p));
      t[1]->AddText(name);
      sprintf(name,"%.4f",fCorrJetMET.at(i)->GetParameter(p));
      t[2]->AddText(name);
    }
    for(int p = 0; p < 3; p++) {
      t[p]->SetBorderSize(0);
      t[p]->SetFillColor(0);
      t[p]->SetTextFont(42);
      t[p]->SetTextAlign(11);
      t[p]->Draw("same");
    }

    c1->Draw();
    ps->NewPage();
    for(int p = 0; p < 3; p++) {
      delete t[p];
    }
  }
  ps->Close();
  

  // -- Cleaning up ------------------------------
  for(unsigned int i = 0; i < fCorrKalibri.size(); i++) {
    delete fCorrKalibri.at(i);
    delete fCorrJetMET.at(i);
  }
  delete hFrame;
  delete c1;
  delete ps;

  std::cout << "Wrote file " << outFileName << ".\n";

  return 0;
}
