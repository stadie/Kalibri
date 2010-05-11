#ifndef flex_res_h
#define flex_res_h

#include "THelpers.h"

class flex_res : public base_corr{

public :
  std::vector < std::vector < TH1D* > > tlj_X_counts_all_;

  std::vector < std::vector < TH2D* > > tlj_X_response_2D_all_;
  std::vector < std::vector < TProfile* > > tlj_X_response_prof_all_;
  std::vector < std::vector < TH1D* > > tlj_X_response_GMP_mean_all_;
  std::vector < std::vector < TH1D* > > tlj_X_response_all_corrections_all_;


  ////RESOLUTION PART

  std::vector < std::vector < TH1D* > > tlj_corrected_response_barrel_all_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_mean_all_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_sigma_all_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_resol_all_;

  std::vector < TGraphErrors* > tlj_Response_Graphs_mean_selec_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_sigma_selec_;
  std::vector < TGraphErrors* > tlj_Response_Graphs_rel_resol_selec_;

  TString root_resol_name_binning;


  virtual void     Loop();
  virtual void Import_Histos();
  virtual void Book_Histos();
  virtual void Write_Histos();
  virtual void Write_TGraphErrors();
   TGraphErrors* flex_res::make_graph(std::vector < std::vector < Double_t > >  Double_gauss_, Int_t y_para, TString title, Int_t X_par=-1);
  virtual void flex_res::draw_graphs(std::vector <TGraphErrors*> graphs_, Double_t ylow, Double_t yhigh, TLegend *legend, TString PDF_PNG_name);


};

#endif
