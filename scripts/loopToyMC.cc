#include <vector>
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TStyle.h"


// Always COMPILE script before execution!

// Generate a toy mc sample, run calibration,
// and store fitted parameters. Loop nLoops
// times
void loopToyMC(int nLoops = 1)
{
  // Two dimensional vector (n,m) storing the
  // parameters with index n of loop m
  std::vector< std::vector<double> > paramLoops;
  int nPar = -1;
  for(int n = 0; n < nLoops; n++)
    {
      cout << endl << "Processing loop " << n+1 << " of " << nLoops << ":" << endl;

      gROOT->ProcessLine(".! ./../toy ./../config/test_tower_bias.cfg");
      gROOT->ProcessLine(".! ./../junk ./../config/test_tower_bias.cfg");

      // Read fitted parameters
      // Needs txt-output from caliber
      // TODO: extend to read different eta bins
      TString fileName = "./../CalibMaker.txt";
      std::ifstream read;
      read.open(fileName);
      double val = -1.;
      if( read >> val )		// Eta low
	{
	  read >> val;		// Eta up
	  read >> nPar;		// number parameters
	  if( n==0 ) paramLoops = std::vector< std::vector<double> >(nPar,std::vector<double>(nLoops));

	  for(int i = 0; i < nPar; i++)
	    {
	      read >> val;
	      paramLoops.at(i).at(n) = val;
	    }
	}
      read.close();
      gROOT->ProcessLine(".! rm "+fileName);
    }


//   for(int p = 0; p < nPar; p++)
//     {
//       cout << "Par " << p << ":  " << flush;
//       for(int n = 0; n < nLoops; n++)
// 	{
// 	  cout << paramLoops.at(p).at(n) << "  " << flush;
// 	}
//       cout << endl;
//     }


  // One histogram per parameter
  cout << "Plotting parameter values.... " << flush;

  std::vector<TH1F*> hPar;
  std::vector<TCanvas*> can;
  for(int p = 0; p < nPar; p++)
    {
      std::vector<double> params = paramLoops.at(p);
      std::vector<double>::const_iterator minParam = std::min_element(params.begin(),params.end());
      std::vector<double>::const_iterator maxParam = std::max_element(params.begin(),params.end());

      TString name = "hPar";
      name += p;
      TString title = "Parameter ";
      title += p;
      title += ";p_{";
      title += p;
      title += "}";
      hPar.push_back(new TH1F(name,title,50,(*minParam)-0.01,(*maxParam)+0.01));
      hPar.back()->SetNdivisions(5);
      for(int n = 0; n < nLoops; n++)
	{
	  hPar.back()->Fill(params.at(n));
	}

      if( p%4 == 0 )		// Put 4 histograms on one canvas
	{
	  name = "can";
	  name += p/4;
	  can.push_back(new TCanvas(name,name,20*p,10*p,800,800));
	  can.back()->Divide(2,2);
	}
      can.back()->cd(1 + p%4);
      hPar.back()->Draw();
    }

  cout << "ok" << endl;


  // Mean values of parameters
  cout << "Plotting mean parameter values.... " << flush;

  gStyle->SetErrorX(0);

  TH1F *hMeanPar = new TH1F("hMeanPar","Parameter mean values;Parameter index i;#bar{p}_{i}",nPar,-0.5,nPar-0.5);
  hMeanPar->SetNdivisions(nPar);
  hMeanPar->SetMarkerStyle(24);

  TH1F *hMeanGaussPar = new TH1F("hMeanGaussPar","Parameter mean values (Gauss fit);Parameter index i;#bar{p}_{i}",nPar,-0.5,nPar-0.5);
  hMeanGaussPar->SetNdivisions(nPar);
  hMeanGaussPar->SetMarkerStyle(24);

  for(int p = 0; p < nPar; p++)
    {
      hMeanPar->SetBinContent(1+p,hPar.at(p)->GetMean());
      hMeanPar->SetBinError(1+p,hPar.at(p)->GetRMS());

      hPar.at(p)->Fit("gaus","0QI");
      TF1 *fit = hPar.at(p)->GetFunction("gaus");
      hMeanGaussPar->SetBinContent(1+p,fit->GetParameter(1));
      hMeanGaussPar->SetBinError(1+p,fit->GetParameter(2));
    }

  TCanvas *canMean = new TCanvas("canMean","Mean Par",500,500);
  canMean->cd();
  hMeanPar->Draw("PE1");

  TCanvas *canMeanGauss = new TCanvas("canMeanGauss","MeanGauss Par",500,500);
  canMeanGauss->cd();
  hMeanGaussPar->Draw("PE1");

  cout << "ok" << endl;
}
