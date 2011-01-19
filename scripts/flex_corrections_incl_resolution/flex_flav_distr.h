#ifndef flex_flav_distr_h
#define flex_flav_distr_h

#include "THelpers.h"

class flex_flav_distr : public base_corr{

public :

  //std::vector < std::vector < TH1D* > > tlj_Corr_Vars_counts_all_;

  std::vector < std::vector < std::vector <TH1D*> > > tlj_Corr_Vars_counts_all_pt_all_PDG_;
  std::vector < std::vector < TH1D* > > tlj_L2L3_Response_all_PDG_;

  //  std::vector < std::vector  <THStack*>  > tlj_Corr_Vars_counts_all_pt_all_Corr_stacked_;


  virtual void     Loop();
 virtual void Book_Histos();
 virtual void Write_Histos();

};

#endif
