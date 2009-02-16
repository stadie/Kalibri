#ifndef TControlPlots_NEW_h
#define TControlPlots_NEW_h

#include <set>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObject.h"
#include "TStyle.h"

class TData;
class TParameters;

class TControlPlots
{
public:
  TControlPlots(const std::string& configfile, const std::vector<TData*> *data, TParameters *par);
  ~TControlPlots();

  bool OutputFormatRoot() const { return mOutputROOT; }

  void MakePlots();
  void MakeControlPlotsDiJet();
  void MakeControlPlotsGammaJet(const std::set<std::string>& plottedQuant);
  void MakeControlPlotsGammaJetPerJetBin();
  void MakeControlPlotsGammaJetPerTowerBin();
  void MakeControlPlotsGammaJetSigmas();
  void MakeControlPlotsParameterScan();
  void MakeControlPlotsTop();
  void MakeControlPlotsTowers();
 private:  
  void Fit2D(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] ) const;
  void SetGStyle() const;
  void WriteToRootFile(std::vector<TObject*> obj, std::string dir);

  const std::vector<TData*> *mData; 
  TParameters *mPar;
  TFile *mOutFile;
  TString mPtRatioName[3];	  // For histo titles etc
  TString mControlQuantityName[8]; // For histo titles etc
  bool mOutputROOT; 
  bool makeControlPlotsTowers;
  bool makeControlPlotsGammaJet;
  bool makeControlPlotsGammaJet2;
  bool makeControlPlotsDiJet;
  bool makeControlPlotsTop;
  bool makeControlPlotsParScan;

  std::set<std::string> mPlottedQuant;
};
#endif
