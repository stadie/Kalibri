#ifndef TControlPlots_h
#define TControlPlots_h

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


//!  \brief Create control plots
//!
//!  Objects of 'TControlPlots' can create control plots via
//!  the 'MakeControlPlots()' method from several TData
//!  objects. The output is in .ps or both .ps and .root format.
//!  The kind of the control plots and the output format is
//!  specified via the config file.
//!
//!  \author Christian Autermann
//!  \date Fri Jan 18 13:55:15 2008 UTC
//!  $Id$
// -------------------------------------------------------------
class TControlPlots
{
public:
  TControlPlots(const std::string& configfile, const std::vector<TData*> *data, TParameters *par);
  ~TControlPlots();

  bool OutputFormatRoot() const { return mOutputROOT; }

  void MakePlots();

 private:  
  void MakeControlPlotsBinnedResponse();
  void MakeControlPlotsChi2();
  void MakeControlPlotsDiJet();
  void MakeControlPlotsGammaJet(const std::set<std::string>& plottedQuant);
  void MakeControlPlotsGammaJetPerJetBin();
  void MakeControlPlotsGammaJetPerTowerBin();
  void MakeControlPlotsGammaJetSigmas();
  void MakeControlPlotsParameterScan();
  void MakeControlPlotsTop();
  void MakeControlPlotsTowers();

  void Fit2D(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] ) const;
  void Fit2DRes(const TH2F* hist, TH1F* hresults[8], TH1F* gaussplots[4], TF1* gf[4] ) const;
  void SetGStyle() const;
  void WriteToRootFile(std::vector<TObject*> obj, std::string dir);

  const        std::vector<TData*> *mData;      //!< Pointer to data
  TParameters *mPar;                            //!< Pointer to parameter values
  TFile       *mOutFile;                        //!< Pointer to root output file
  TString      mPtRatioName[3];	                //!< For histo titles etc
  TString      mControlQuantityName[8];         //!< For histo titles etc
  bool         mOutputROOT;                     //!< If true, histograms are written to ROOT file
  bool         makeControlPlotsBinnedResponse;
  bool         makeControlPlotsChi2;
  bool         makeControlPlotsTowers;
  bool         makeControlPlotsGammaJet;
  bool         makeControlPlotsGammaJet2;
  bool         makeControlPlotsDiJet;
  bool         makeControlPlotsTop;
  bool         makeControlPlotsParScan;

  std::set<std::string> mPlottedQuant;



  //!  \brief A two-dimensional grid
  //!
  //!  The two dimensions of the grid are named 'x' and 'y'.
  //!  They are divided into bins whose borders can be specified
  //!  i.e. the binsize can vary. The index of the bins in x
  //!  direction is named 'ix' and ranges from 0 to NBinsX() - 1
  //!  and likewise for the y direction. Additionally, there is
  //!  a global bin numbering scheme 'bin' from 0 to NBins() - 1,
  //!  counting the bins first in x and then in y direction i.e.
  //!  
  //!     bin    ix  iy
  //!  -------------------
  //!      0      0   0
  //!      1      1   0
  //!     ...
  //!   NBinsX()  0   1
  //!  
  //!  \author Matthias Schroeder
  //!  \date Thu Apr 23 13:05:54 CEST 2009
  // -------------------------------------------------------------
  class Binning
    {
    public:
      //!  \brief Creates a binning in (x,y)
      //!
      //!  The bins are created from the given bin edges.
      //!  \note The vectors of bin edges must have at least two entries
      //!        and must be ordered from lowest to highest value.
      //!  \param binEdgesX Bin edges in x direction (binEdgesX.size() == NBinsX()+1)
      //!  \param binEdgesY Bin edges in y direction (binEdgesY.size() == NBinsY()+1)
      // -------------------------------------------------------------
      Binning(const std::vector<double>& binEdgesX, const std::vector<double>& binEdgesY);
      ~Binning() {};

      //!  \brief Lower bin edge in x direction
      //!  \param bin Global bin index
      //!  \return Lower bin edge in x direction
      // -------------------------------------------------------------
      double XLow(int bin) const { return mEdgesX.at(IX(bin)); }

      //!  \brief Upper bin edge in x direction
      //!  \param bin Global bin index
      //!  \return Upper bin edge in x direction
      // -------------------------------------------------------------
      double XUp(int bin) const { return mEdgesX.at(IX(bin)+1); }

      //!  \brief Lower bin edge in y direction
      //!  \param bin Global bin index
      //!  \return Lower bin edge in y direction
      // -------------------------------------------------------------
      double YLow(int bin) const { return mEdgesY.at(IY(bin)); }

      //!  \brief Upper bin edge in y direction
      //!  \param bin Global bin index
      //!  \return Upper bin edge in y direction
      // -------------------------------------------------------------
      double YUp(int bin) const { return mEdgesY.at(IY(bin)+1); }

      //!  \brief X bin index ix of global bin
      //!
      //!  Finds the index ix of the bin in x direction corresponding
      //!  to a global bin.
      //!  \param bin Global bin index
      //!  \return ix of global bin
      // -------------------------------------------------------------
      int IX(int bin) const { return bin % NBinsX(); }

      //!  \brief X bin index ix of value x
      //!
      //!  Finds the index ix of the bin in x direction that contains
      //!  the value x.
      //!  \param bin Global bin index
      //!  \return ix of the bin containing x
      //!          ( ix == -1 for x < XLow(0): Underflow, 
      //!            ix == NBinsX() for x > XUp(NBins()-1): Overflow )
      // -------------------------------------------------------------
      int IX(double x) const;

      //!  \brief Y bin index iy of global bin
      //!
      //!  Finds the index iy of the bin in y direction corresponding
      //!  to a global bin.
      //!  \param bin Global bin index
      //!  \return iy of global bin
      // -------------------------------------------------------------
      int IY(int bin) const { return bin / NBinsX(); }

      //!  \brief Y bin index iy of value y
      //!
      //!  Finds the index iy of the bin in y direction that contains
      //!  the value y.
      //!  \param bin Global bin index
      //!  \return iy of the bin containing y
      //!          ( iy == -1 for y < YLow(0): Underflow, 
      //!            iy == NBinsyY() for y > YUp(NBins()-1): Overflow )
      // -------------------------------------------------------------
      int IY(double y) const;

      //!  \brief Global bin index of bin with x and y indices ix and iy
      //!  \param ix Bin index in x direction
      //!  \param iy Bin index in y direction
      //!  \return Global bin index
      // -------------------------------------------------------------
      int Bin(int ix, int iy) const { return ix + iy*NBinsX(); }

      //!  \brief Global bin index of bin containing values x and y
      //!  \param x x value
      //!  \param y y value
      //!  \return Global bin index
      // -------------------------------------------------------------
      int Bin(double x, double y) const { return Bin(IX(x),IY(y)); }

      //!  \brief Number of global bins
      //!
      //!  NBins() == NBinsX() * NBinsY()
      //!  \return Number of global bins
      // -------------------------------------------------------------
      int NBins() const { return  NBinsX()*NBinsY(); }

      //!  \brief Number of x bins
      //!  \return Number of x bins
      // -------------------------------------------------------------
      int NBinsX() const { return static_cast<int>(mEdgesX.size()) - 1; }

      //!  \brief Number of y bins
      //!  \return Number of y bins
      // -------------------------------------------------------------
      int NBinsY() const { return static_cast<int>(mEdgesY.size()) - 1; }

      //!  \brief Print the binning to std
      // -------------------------------------------------------------
      void Print() const;

    private:
      std::vector<double> mEdgesX;   //!< Bin edges in x direction
      std::vector<double> mEdgesY;   //!< Bin edges in y direction
    };
};
#endif
