//
// $Id: CalibData.h,v 1.74 2009/11/20 12:25:41 stadie Exp $
//
#ifndef ControlPlotsComparison_h
#define ControlPlotsComparison_h

#include <string>
#include <vector>

#include <TFile.h>

class TControlPlotsComparison
{
public:
  TControlPlotsComparison();
  ~TControlPlotsComparison() {};

  void CompareControlPlots(const std::vector<std::string> &fileName, std::vector<std::string> descr = 0) const;

private:
  void CompareControlPlotsDiJet(const std::vector<TFile*> &file, const std::vector<std::string> &descr = 0) const;
  void CompareControlPlotsGammaJet(const std::vector<TFile*> &file, const std::vector<std::string> &descr = 0) const;
  void SetGStyle() const;
};
#endif
