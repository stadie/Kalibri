#ifndef TControlPlots_NEW_h
#define TControlPlots_NEW_h

#include <string>
#include <vector>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TObject.h>
#include <TStyle.h>

class TData;
class TParameters;

class TControlPlots
{
public:
  TControlPlots(const std::vector<TData*> *data, TParameters *par, int outputFormat = 0);
  ~TControlPlots();

  void MakeControlPlotsDiJet();
  void MakeControlPlotsGammaJet();
  void MakeControlPlotsGammaJetPerJetBin();
  void MakeControlPlotsGammaJetPerTowerBin();
  void MakeControlPlotsGammaJetSigmas();
  void MakeControlPlotsParameterScan();
  void MakeControlPlotsTowers();
  bool OutputFormatRoot() const { return mOutputROOT; }
  
private:
  void Fit2D(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] ) const;
  void SetGStyle() const;
  void WriteToRootFile(std::vector<TObject*> obj, std::string dir);

  const std::vector<TData*> *mData; 
  TParameters *mPar;
  TFile * const mOutFile;
  TString mPtRatioName[3];	  // For histo titles etc
  TString mControlQuantityName[8]; // For histo titles etc
  bool mOutputROOT;
};
#endif
