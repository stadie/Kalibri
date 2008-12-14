//
//    Possible base class for all event processors
//    
//    Right now we only have the spectrum correctors.
//    Thus they are implemented directly in this class
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventReader.h,v 1.1 2008/12/12 13:43:15 stadie Exp $
//   
#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

class TData;
class TParameters;
class TH1;
class TF1;

#include <vector>
#include <string>

class EventProcessor
{
 public:
  EventProcessor(const std::string& configfile, TParameters* param);
  ~EventProcessor();
  int process(std::vector<TData*>& data);
 private:
  void FlattenSpectra(std::vector<TData*>& data);
  void BalanceSpectra(std::vector<TData*>& data);
  int GetSpectraBin(double m1, double m2, double m3);
  static void FitWithoutBottom(TH1 * hist, TF1 * func, double bottom=0.33);
  static double gauss_step(double *x, double *par);

  TParameters* p;
  double Et_cut_on_gamma, Et_cut_on_jet;
  bool flatten_spectra;
  double RelWeight[7];//@@ Replace 7 by something meaningful
  std::vector<int> _residualScalingScheme;          // Iteration scheme of scaling of residuals

  typedef std::vector<TData*>::iterator DataIter;
  typedef std::vector<TData*>::const_iterator DataConstIter;
  //"Not-Balanced" Rejection: Make average-fitting equal to peak-fitting
  struct NotBalancedRejection {
  NotBalancedRejection(double *cut, double min, double max):
    _cut(cut),_min(min),_max(max){};
    bool operator()(TData *d);
    double *_cut;
    double _min, _max;
  };


};


#endif
