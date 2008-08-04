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
  bool OutputFormatRoot() const { return _outputROOT; }
  
private:
  void Fit2D(TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] );
  void SetGStyle();
  void WriteToRootFile(std::vector<TObject*> obj, std::string dir);

  const std::vector<TData*> *_data; 
  TParameters *_par;
  TFile * const _outFile;
  TString ptRatioName[3];	  // For histo titles etc
  TString controlQuantityName[8]; // For histo titles etc
  bool _outputROOT;
};
#endif
