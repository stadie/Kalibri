//!
//!  \brief Functions to quickly plot control plots
//!  from root-file
//!
//!  \author Matthias Schroeder
//!  \date   Wed Apr  1 18:28:02 CEST 2009
//!  $Id$
//!

#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"
#include "TStyle.h"


namespace makePlots
{
  // - Function declarations ----------------------------------------------------
  void AddFileName(const std::string& name);
  void Resolution(const std::string& type = "GammaJet");
  void Response(const std::string& type = "GammaJet");
  void ResponseDistribution(double min, double max, const std::string& type = "GammaJet");
  void ResponseParametrization(const std::string& className, const std::vector<double>& par, const std::vector<double>& parGlobal);
  void CorrectionParametrization(const std::string& className, const std::vector<double>& par, const std::vector<double>& parGlobal);



  // - Global variables ---------------------------------------------------------

  //!  List of files to be included
  std::vector<std::string> fileNames;




  // - Function implementation --------------------------------------------------

  //!  \brief Add a file name to fileNames
  //!  
  //!  \param name File name
  void AddFileName(const std::string& name)
  {
    fileNames.push_back(name);
  }


  //!  \brief Plot resolution
  //!
  //!  Plot relative resolution of \f$ E^{meas}_{T} / E^{gen}_{T} \f$
  //!  vs \f$ E^{gen}_{T} \f$ for all files in fileNames .
  //!  The relative resolution is calculated once as RMS and once
  //!  as sigma of a Gaussian fit of the
  //!  \f$ E^{meas}_{T} / E^{gen}_{T} \f$ distribution.
  //!
  //!  \param type Name of branch, default is "GammaJet"
  void Resolution(const std::string& type)
  {
    std::vector<TH1F*> histsUncorr;
    std::vector<TH1F*> histsCorr;
    std::vector<TH1F*> histsGaussUncorr;
    std::vector<TH1F*> histsGaussCorr;

    // Get histograms from file
    for(std::vector<std::string>::const_iterator it = fileNames.begin();
	it != fileNames.end(); it++)
      {
	TH1F *hUncorr       = 0;
	TH1F *hCorr         = 0;
	TH1F *hGaussUncorr  = 0;
	TH1F *hGaussCorr    = 0;

	TFile *file = new TFile((*it).c_str(),"READ");
	std::string histName = type;
	histName.append("/hpt0_result1");
	file->GetObject(histName.c_str(),hUncorr);
	histName = type;
	histName.append("/hpt1_result1");
	file->GetObject(histName.c_str(),hCorr);
	histName = type;
	histName.append("/hpt0_result3");
	file->GetObject(histName.c_str(),hGaussUncorr);
	histName = type;
	histName.append("/hpt1_result3");
	file->GetObject(histName.c_str(),hGaussCorr);

	if(  hUncorr && hCorr && hGaussUncorr && hGaussCorr )
	  {
	    hUncorr->SetDirectory(0);
	    histsUncorr.push_back(hUncorr);

	    hCorr->SetDirectory(0);
	    histsCorr.push_back(hCorr);

	    hGaussUncorr->SetDirectory(0);
	    histsGaussUncorr.push_back(hGaussUncorr);

	    hGaussCorr->SetDirectory(0);
	    histsGaussCorr.push_back(hGaussCorr);
	  }
	else
	  {
	    std::cerr << "ERROR getting histogram from file " << *it << std::endl;
	  }

	delete file;
      }


    // Set nice style
    for(unsigned int i = 0; i < histsUncorr.size(); i++)
      {
	histsUncorr.at(i)->GetYaxis()->SetRangeUser(0,0.5);
	histsUncorr.at(i)->SetTitle("Standard deviation");
	histsUncorr.at(i)->GetYaxis()->SetTitle("Relative resolution  < #sigma / E^{gen}_{T} >");
	histsUncorr.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");

	histsCorr.at(i)->GetYaxis()->SetRangeUser(0,0.5);
	histsCorr.at(i)->SetTitle("Standard deviation");
	histsCorr.at(i)->GetYaxis()->SetTitle("Relative resolution  < #sigma / E^{gen}_{T} >");
	histsCorr.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");
	histsCorr.at(i)->SetMarkerColor(2+i);
	histsCorr.at(i)->SetLineColor(2+i);

	histsGaussUncorr.at(i)->GetYaxis()->SetRangeUser(0,0.5);
	histsGaussUncorr.at(i)->SetTitle("Width of Gaussian fit");
	histsGaussUncorr.at(i)->GetYaxis()->SetTitle("Relative resolution  < #sigma / E^{gen}_{T} >");
	histsGaussUncorr.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");

	histsGaussCorr.at(i)->GetYaxis()->SetRangeUser(0,0.5);
	histsGaussCorr.at(i)->SetTitle("Width of Gaussian fit");
	histsGaussCorr.at(i)->GetYaxis()->SetTitle("Relative resolution  < #sigma / E^{gen}_{T} >");
	histsGaussCorr.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");
	histsGaussCorr.at(i)->SetMarkerColor(2+i);
	histsGaussCorr.at(i)->SetLineColor(2+i);
      }


    // Plot histograms
    TCanvas *can = new TCanvas("canResolution","Resolution",1000,500);
    can->Divide(2,1);

    can->cd(1);
    histsUncorr.at(0)->Draw();
    for(std::vector<TH1F*>::iterator it = histsCorr.begin();
	it != histsCorr.end(); it++)
      {
	(*it)->Draw("same");
      }
    gPad->SetGrid();

    can->cd(2);
    histsGaussUncorr.at(0)->Draw();
    for(std::vector<TH1F*>::iterator it = histsGaussCorr.begin();
	it != histsGaussCorr.end(); it++)
      {
	(*it)->Draw("same");
      }
    gPad->SetGrid();

  }



  //!  \brief Plot response
  //!
  //!  Plot response \f$ E^{meas}_{T} / E^{gen}_{T} \f$
  //!  vs \f$ E^{gen}_{T} \f$ and vs \f$ \eta_{T} \f$
  //!  for all files in fileNames .
  //!
  //!  \param type Name of branch, default is "GammaJet"
  void Response(const std::string& type)
  {
    // Response vs pt
    std::vector<TH1F*> histsUncorrPt;
    std::vector<TH1F*> histsCorrPt;
    std::vector<TH1F*> histsGaussUncorrPt;
    std::vector<TH1F*> histsGaussCorrPt;

    // Response vs eta
    std::vector<TH1F*> histsUncorrEta;
    std::vector<TH1F*> histsCorrEta;
    std::vector<TH1F*> histsGaussUncorrEta;
    std::vector<TH1F*> histsGaussCorrEta;

    // Get histograms from file
    for(std::vector<std::string>::const_iterator it = fileNames.begin();
	it != fileNames.end(); it++)
      {
	TH1F *hUncorrPt        = 0;
	TH1F *hCorrPt          = 0;
	TH1F *hGaussUncorrPt   = 0;
	TH1F *hGaussCorrPt     = 0;

	TH1F *hUncorrEta       = 0;
	TH1F *hCorrEta         = 0;
	TH1F *hGaussUncorrEta  = 0;
	TH1F *hGaussCorrEta    = 0;


	TFile *file = new TFile((*it).c_str(),"READ");

	std::string histName = type;
	histName.append("/hpt0_result0");
	file->GetObject(histName.c_str(),hUncorrPt);
	histName = type;
	histName.append("/hpt1_result0");
	file->GetObject(histName.c_str(),hCorrPt);
	histName = type;
	histName.append("/hpt0_result2");
	file->GetObject(histName.c_str(),hGaussUncorrPt);
	histName = type;
	histName.append("/hpt1_result2");
	file->GetObject(histName.c_str(),hGaussCorrPt);

	histName = type;
	histName.append("/heta0_result0");
	file->GetObject(histName.c_str(),hUncorrEta);
	histName = type;
	histName.append("/heta1_result0");
	file->GetObject(histName.c_str(),hCorrEta);
	histName = type;
	histName.append("/heta0_result2");
	file->GetObject(histName.c_str(),hGaussUncorrEta);
	histName = type;
	histName.append("/heta1_result2");
	file->GetObject(histName.c_str(),hGaussCorrEta);


	if(  hUncorrPt && hCorrPt && hGaussUncorrPt && hGaussCorrPt
	     &&  hUncorrPt && hCorrPt && hGaussUncorrPt && hGaussCorrPt  )
	  {
	    hUncorrPt->SetDirectory(0);
	    histsUncorrPt.push_back(hUncorrPt);

	    hCorrPt->SetDirectory(0);
	    histsCorrPt.push_back(hCorrPt);

	    hGaussUncorrPt->SetDirectory(0);
	    histsGaussUncorrPt.push_back(hGaussUncorrPt);

	    hGaussCorrPt->SetDirectory(0);
	    histsGaussCorrPt.push_back(hGaussCorrPt);

	    hUncorrEta->SetDirectory(0);
	    histsUncorrEta.push_back(hUncorrEta);

	    hCorrEta->SetDirectory(0);
	    histsCorrEta.push_back(hCorrEta);

	    hGaussUncorrEta->SetDirectory(0);
	    histsGaussUncorrEta.push_back(hGaussUncorrEta);

	    hGaussCorrEta->SetDirectory(0);
	    histsGaussCorrEta.push_back(hGaussCorrEta);
	  }
	else
	  {
	    std::cerr << "ERROR getting histogram from file " << *it << std::endl;
	  }

	delete file;
      }


    // Set nice style
    for(unsigned int i = 0; i < histsUncorrPt.size(); i++)
      {
	histsUncorrPt.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsUncorrPt.at(i)->SetTitle("Mean");
	histsUncorrPt.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsUncorrPt.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");

	histsCorrPt.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsCorrPt.at(i)->SetTitle("Mean");
	histsCorrPt.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsCorrPt.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");
	histsCorrPt.at(i)->SetMarkerColor(2+i);
	histsCorrPt.at(i)->SetLineColor(2+i);

	histsGaussUncorrPt.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsGaussUncorrPt.at(i)->SetTitle("Mean of Gaussian fit");
	histsGaussUncorrPt.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsGaussUncorrPt.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");

	histsGaussCorrPt.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsGaussCorrPt.at(i)->SetTitle("Mean of Gaussian fit");
	histsGaussCorrPt.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsGaussCorrPt.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");
	histsGaussCorrPt.at(i)->SetMarkerColor(2+i);
	histsGaussCorrPt.at(i)->SetLineColor(2+i);

	histsUncorrEta.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsUncorrEta.at(i)->SetTitle("Mean");
	histsUncorrEta.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsUncorrEta.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");

	histsCorrEta.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsCorrEta.at(i)->SetTitle("Mean");
	histsCorrEta.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsCorrEta.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");
	histsCorrEta.at(i)->SetMarkerColor(2+i);
	histsCorrEta.at(i)->SetLineColor(2+i);

	histsGaussUncorrEta.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsGaussUncorrEta.at(i)->SetTitle("Mean of Gaussian fit");
	histsGaussUncorrEta.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsGaussUncorrEta.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");

	histsGaussCorrEta.at(i)->GetYaxis()->SetRangeUser(0,1.5);
	histsGaussCorrEta.at(i)->SetTitle("Mean of Gaussian fit");
	histsGaussCorrEta.at(i)->GetYaxis()->SetTitle("Response  < E_{T} / E^{gen}_{T} >");
	histsGaussCorrEta.at(i)->GetXaxis()->SetTitle("E^{gen}_{T}  (GeV)");
	histsGaussCorrEta.at(i)->SetMarkerColor(2+i);
	histsGaussCorrEta.at(i)->SetLineColor(2+i);
      }


    // Plot histograms
    TCanvas *canPt = new TCanvas("canResponsePt","Response Pt",1000,500);
    canPt->Divide(2,1);

    canPt->cd(1);
    histsUncorrPt.at(0)->DrawClone();
    for(std::vector<TH1F*>::iterator it = histsCorrPt.begin();
	it != histsCorrPt.end(); it++)
      {
	(*it)->DrawClone("same");
      }
    gPad->SetGrid();

    canPt->cd(2);
    histsGaussUncorrPt.at(0)->DrawClone();
    for(std::vector<TH1F*>::iterator it = histsGaussCorrPt.begin();
	it != histsGaussCorrPt.end(); it++)
      {
	(*it)->DrawClone("same");
      }
    gPad->SetGrid();


    TCanvas *canPtZ = new TCanvas("canResponsePtZoom","Response Pt Zoom",1000,500);
    canPtZ->Divide(2,1);

    canPtZ->cd(1);
    histsUncorrPt.at(0)->GetYaxis()->SetRangeUser(0.8,1.2);
    histsUncorrPt.at(0)->Draw();
    for(std::vector<TH1F*>::iterator it = histsCorrPt.begin();
	it != histsCorrPt.end(); it++)
      {
	(*it)->Draw("same");
      }
    gPad->SetGrid();

    canPtZ->cd(2);
    histsGaussUncorrPt.at(0)->GetYaxis()->SetRangeUser(0.8,1.2);
    histsGaussUncorrPt.at(0)->Draw();
    for(std::vector<TH1F*>::iterator it = histsGaussCorrPt.begin();
	it != histsGaussCorrPt.end(); it++)
      {
	(*it)->Draw("same");
      }
    gPad->SetGrid();


    TCanvas *canEta = new TCanvas("canResponseEta","Response Eta",1000,500);
    canEta->Divide(2,1);

    canEta->cd(1);
    histsUncorrEta.at(0)->DrawClone();
    for(std::vector<TH1F*>::iterator it = histsCorrEta.begin();
	it != histsCorrEta.end(); it++)
      {
	(*it)->DrawClone("same");
      }
    gPad->SetGrid();

    canEta->cd(2);
    histsGaussUncorrEta.at(0)->DrawClone();
    for(std::vector<TH1F*>::iterator it = histsGaussCorrEta.begin();
	it != histsGaussCorrEta.end(); it++)
      {
	(*it)->DrawClone("same");
      }
    gPad->SetGrid();


    TCanvas *canEtaZ = new TCanvas("canResponseEtaZoom","Response Eta Zoom",1000,500);
    canEtaZ->Divide(2,1);

    canEtaZ->cd(1);
    histsUncorrEta.at(0)->GetYaxis()->SetRangeUser(0.8,1.2);
    histsUncorrEta.at(0)->Draw();
    for(std::vector<TH1F*>::iterator it = histsCorrEta.begin();
	it != histsCorrEta.end(); it++)
      {
	(*it)->Draw("same");
      }
    gPad->SetGrid();

    canEtaZ->cd(2);
    histsGaussUncorrEta.at(0)->GetYaxis()->SetRangeUser(0.8,1.2);
    histsGaussUncorrEta.at(0)->Draw();
    for(std::vector<TH1F*>::iterator it = histsGaussCorrEta.begin();
	it != histsGaussCorrEta.end(); it++)
      {
	(*it)->Draw("same");
      }
    gPad->SetGrid();
  }



  void ResponseDistribution(double min, double max, const std::string& type)
  {
    // Store the projections of corrected / uncorrected
    // response distributions
    std::vector<TH1F*> histsUncorr;
    std::vector<TH1F*> histsCorr;

    // Get 2D histograms resp vs Et gen from file
    // and extract projection
    for(std::vector<std::string>::const_iterator it = fileNames.begin();
	it != fileNames.end(); it++)
      {
	// Get histograms from file
	TH2F *h2Uncorr = 0;
	TH2F *h2Corr   = 0;
	TFile *file = new TFile((*it).c_str(),"READ");
	std::string histName = type;
	histName.append("/hpt0");
	file->GetObject(histName.c_str(),h2Uncorr);
	histName = type;
	histName.append("/hpt1");
	file->GetObject(histName.c_str(),h2Corr);

	// Extract projection
	if(  h2Uncorr && h2Corr  )
	  {
	    // Declare histograms holding projections
	    TH1F *hUncorr = new TH1F("ProjUncorr","",
				     h2Uncorr->GetNbinsY(),
				     h2Uncorr->GetYaxis()->GetBinLowEdge(1),
				     h2Uncorr->GetYaxis()->GetBinUpEdge(h2Uncorr->GetNbinsY()));
	    TH1F *hCorr = static_cast<TH1F*>(hUncorr->Clone("ProjCorr"));

	    // Find x-bins of min and max
	    int minXbin = h2Uncorr->GetXaxis()->FindBin(min);
	    int maxXbin = h2Uncorr->GetXaxis()->FindBin(max);

	    // Loop over y-bins
	    for(int ybin = 1; ybin <= h2Uncorr->GetNbinsY(); ybin++)
	      {
		// Store content of bins (xbin,ybin)
		// with xbin between min, max
		double contentUncorr = 0;
		double contentCorr   = 0;

		// Loop over x-bins between min and max
		for(int xbin = minXbin; xbin <= maxXbin; xbin++)
		  {
		    // Add up content
		    contentUncorr += h2Uncorr -> GetBinContent( h2Uncorr -> GetBin(xbin,ybin) );
		    contentCorr   += h2Corr   -> GetBinContent( h2Corr   -> GetBin(xbin,ybin) );
		  }

		// Fill content into projections
		hUncorr -> SetBinContent( ybin, contentUncorr );
		hCorr   -> SetBinContent( ybin, contentCorr   );
	      }

	    // Fill projections into vectors
	    hUncorr->SetDirectory(0);
	    histsUncorr.push_back(hUncorr);

	    hCorr->SetDirectory(0);
	    histsCorr.push_back(hCorr);
	  }
	else
	  {
	    std::cerr << "ERROR getting histogram from file " << *it << std::endl;
	  }

	delete file;
      }


    // Set nice style
    // Find maximum
    double ymax = 0;
    for(unsigned int i = 0; i < histsUncorr.size(); i++)
      {
	if( histsUncorr.at(i) -> GetMaximum() > ymax ) ymax = histsUncorr.at(i) -> GetMaximum();
	if( histsCorr.at(i)   -> GetMaximum() > ymax ) ymax = histsCorr.at(i)   -> GetMaximum();
      }
    for(unsigned int i = 0; i < histsUncorr.size(); i++)
      {
	histsUncorr.at(i)->GetYaxis()->SetRangeUser(0,1.1*ymax);
 	histsUncorr.at(i)->GetXaxis()->SetTitle("E^{meas}_{T} / E^{gen}_{T}");
	char title[50];
	sprintf(title,"%.1f < E^{gen}_{T} < %.1f GeV",min,max);
	histsUncorr.at(i)->SetTitle(title);

	histsCorr.at(i)->GetYaxis()->SetRangeUser(0,1.1*ymax);
 	histsCorr.at(i)->GetXaxis()->SetTitle("E^{meas}_{T} / E^{gen}_{T}");
	histsCorr.at(i)->SetTitle(title);
	histsCorr.at(i)->SetMarkerColor(2+i);
	histsCorr.at(i)->SetLineColor(2+i);
      }


    // Plot histograms
    gStyle->SetOptStat(0);
    TCanvas *can = new TCanvas("canProjection","Projection",500,500);
    can->cd();
    histsUncorr.at(0)->Draw();
    for(std::vector<TH1F*>::iterator it = histsCorr.begin();
	it != histsCorr.end(); it++)
      {
	(*it)->Draw("same");
      }
    gPad->SetGrid();
  }



  //!  \brief Plot response functions of a certain parametrization class
  //!
  //!  The supported parametrization classes are
  //!  - L2L3: JetMET L2 and L3 parametrization
  //!  If no parameters are given (empty or wrongly sized vectors),
  //!  default values are used.
  //!
  //!  \param className  Name of parametrization class
  //!  \param par        Parameter values for tower parametrization
  //!  \param parGlobal  Parameter values for global parametrization
  void ResponseParametrization(const std::string& className, const std::vector<double>& par, const std::vector<double>& parGlobal)
  {
    if( className.compare("L2L3") == 0  )
      {
	// Check if parameter vectors have correct size
	std::vector<double> tmpParGlobal;
	if( parGlobal.size() >= 5 )
	  {
	    tmpParGlobal = parGlobal;
	  }
	// if not, use default values
	else
	  {
	    tmpParGlobal.push_back(0.99);
	    tmpParGlobal.push_back(7.24);
	    tmpParGlobal.push_back(3.50);
	    tmpParGlobal.push_back(6.77);
	    tmpParGlobal.push_back(2.64);
	    std::cout << "makePlots::ResponseParametrization: Using default parameter values" << std::endl;
	  }

	// Create global response function
	TF1 *fGlobal = new TF1("fGlobalResponseL2L3","[0] - [1]/(pow(log10(x),[2]) + [3]) + [4]/x",5,1000);
	for(unsigned int i = 0; i < tmpParGlobal.size(); i++)
	  {
	    fGlobal->SetParameter(i,tmpParGlobal.at(i));
	  }
	fGlobal->SetTitle("L3 response function");
	fGlobal->GetXaxis()->SetTitle("p^{gen}_{T,had} (GeV)");

	// Plot response function
	TCanvas *can = new TCanvas("canGlobalParametrizationL2L3","L2L3 Response",500,500);
	can->cd();
	fGlobal->Draw();
	gPad->SetLogx();
      }
    else
      {
	std::cerr << "ERROR makePlots::ResponseParametrization: " << className << " is not a known parametrization class." << std::endl;
      }
  }


  //!  \brief Plot correction functions of a certain parametrization class
  //!
  //!  The supported parametrization classes are
  //!  - L2L3: JetMET L2 and L3 parametrization
  //!  If no parameters are given (empty or wrongly sized vectors),
  //!  default values are used.
  //!
  //!  \param className  Name of parametrization class
  //!  \param par        Parameter values for tower parametrization
  //!  \param parGlobal  Parameter values for global parametrization
  void CorrectionParametrization(const std::string& className, const std::vector<double>& par, const std::vector<double>& parGlobal)
  {
    if( className.compare("L2L3") == 0  )
      {
	// Check if parameter vectors have correct size
	std::vector<double> tmpParGlobal;
	if( parGlobal.size() >= 5 )
	  {
	    tmpParGlobal = parGlobal;
	  }
	// if not, use default values
	else
	  {
	    tmpParGlobal.push_back(1.00);
	    tmpParGlobal.push_back(5.00);
	    tmpParGlobal.push_back(3.08);
	    tmpParGlobal.push_back(2.05);
	    tmpParGlobal.push_back(8.00);
	    std::cout << "makePlots::CorrectionParametrization: Using default parameter values" << std::endl;
	  }

	// Create global response function
	TF1 *fGlobal = new TF1("fGlobalCorrectionL2L3","[0] + [1]/(pow(log10(x),[2]) + [3]) - [4]/x",7,1000);
	for(unsigned int i = 0; i < tmpParGlobal.size(); i++)
	  {
	    fGlobal->SetParameter(i,tmpParGlobal.at(i));
	  }
	fGlobal->SetTitle("L3 correction function");
	fGlobal->GetXaxis()->SetTitle("p^{gen}_{T,had} (GeV)");

	// Plot response function
	TCanvas *can = new TCanvas("canGlobalCorrectionL2L3","L2L3 Correction",500,500);
	can->cd();
	fGlobal->Draw();
	gPad->SetLogx();
      }
    else
      {
	std::cerr << "ERROR makePlots::CorrectionParametrization: " << className << " is not a known parametrization class." << std::endl;
      }
  }
}


