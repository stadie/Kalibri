#define flex_res_cxx
#include "THelpers.h"
#include "base_corr.C"
#include "flex_res.h"

void flex_res::Import_Histos()
{

  TString root_name_binning= "histos_" + X_labels_[bin_choice] + ".root";

  cout << "DEBUG: " << root_name_binning << endl;

  TFile* _file=new TFile(root_name_binning,"OPEN");

  //    TFile f("file.root");
    if (_file->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;



  //  std::vector < std::vector < TH1D* > > tlj_X_response_all_corrections_all_;


  for (Int_t X_i=0;X_i<no_X_labels_;X_i++)
    {
      std::vector < TH1D* >  tlj_X_counts_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   //	   cout << define_pt_histo_name("tlj_counts_"+(X_labels_[X_i]+"_"), pt_bins_, pt_i).Data() << endl;
      	   tlj_X_counts_.push_back((TH1D*)_file->Get(define_pt_histo_name("tlj_counts_"+(X_labels_[X_i]+"_"), pt_bins_, pt_i).Data()));
	 }
       tlj_X_counts_all_.push_back(tlj_X_counts_);
    }

  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < TH2D* >  tlj_X_response_2D_;
      std::vector < TProfile* >  tlj_X_response_prof_;
      std::vector < TH1D* >  tlj_X_response_GMP_mean_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_X_response_2D_.push_back((TH2D*)_file->Get(define_pt_histo_name( "tlj_X_response_2D_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i)));
	   tlj_X_response_prof_.push_back((TProfile*)_file->Get(define_pt_histo_name( "tlj_X_response_prof_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i)));
	   TString temp =  (TString) tlj_X_response_2D_.back()->GetName() + "_1";
	   tlj_X_response_GMP_mean_.push_back((TH1D*)_file->Get(temp.Data()));

	 }

              tlj_X_response_2D_all_.push_back(tlj_X_response_2D_);
              tlj_X_response_prof_all_.push_back(tlj_X_response_prof_);
              tlj_X_response_GMP_mean_all_.push_back(tlj_X_response_GMP_mean_);

    }

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* > tlj_X_response_all_corrections_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   if(Corr_i<no_Corr_labels_) tlj_X_response_all_corrections_.push_back((TH1D*)tlj_X_response_prof_all_[Corr_i].at(pt_i));
	   else tlj_X_response_all_corrections_.push_back(tlj_X_response_GMP_mean_all_[Corr_i-no_Corr_labels_].at(pt_i));
	   //	   cout << "ACHTUNG!!!!!!!!!!!!!!" << tlj_X_response_all_corrections_.back()->GetName() << endl;

	 }
       tlj_X_response_all_corrections_all_.push_back(tlj_X_response_all_corrections_);
    }


  //  tlj_X_response_all_corrections_all_[9].at(0)->Draw();

}


void flex_res::Book_Histos()
{

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* >  tlj_corrected_response_barrel_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_corrected_response_barrel_.push_back(define_pt_histo( "tlj_corrected_response_barrel_"+Corr_labels_[Corr_i].first+"_", ";corrected Response", pt_bins_, pt_i, response_bins, response_low, response_high));
	   //	   cout << tlj_corrected_response_barrel_.back()->GetName() << endl;
	   //	   cout << "Corr_i: " <<Corr_i <<" pt_i: " <<pt_i << endl;
	 }
       tlj_corrected_response_barrel_all_.push_back(tlj_corrected_response_barrel_);
    }


 for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
   {
  std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_counts_;
  std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_counts_helper_;
  std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_response_;
  std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_response2_;
  std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_mean_response_;
   for(unsigned int Sel_Corr_j=Sel_Corr_i+1;Sel_Corr_j<Corr_selected_labels_.size();Sel_Corr_j++)
     {
	   std::vector < TH2D* >  tlj_Sel_Correlations_2D_counts_;
	   std::vector < TH2D* >  tlj_Sel_Correlations_2D_counts_helper_;
	   std::vector < TH2D* >  tlj_Sel_Correlations_2D_response_;
	   std::vector < TH2D* >  tlj_Sel_Correlations_2D_response2_;
	   std::vector < TH2D* >  tlj_Sel_Correlations_2D_mean_response_;
//         if(Sel_Corr_i!=Sel_Corr_j)
//  	 {
	   for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	     {
	       tlj_Sel_Correlations_2D_counts_.push_back(define_pt_histo( "tlj_Sel_Correlations_2D_counts_"+
	       Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+"_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first + 
               "_", ";"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+";"+  
               Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].second,
                s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].second));
	       tlj_Sel_Correlations_2D_counts_.back()->Sumw2();

	       tlj_Sel_Correlations_2D_counts_helper_.push_back(define_pt_histo( "tlj_Sel_Correlations_2D_counts_helper_"+
	       Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+"_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first + 
               "_", ";"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+";"+  
               Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].second,
                s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].second));
	       tlj_Sel_Correlations_2D_counts_helper_.back()->Sumw2();

	       tlj_Sel_Correlations_2D_response_.push_back(define_pt_histo( "tlj_Sel_Correlations_2D_response_"+
	       Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+"_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first + 
               "_", ";"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+";"+  
               Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].second,
                s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].second));
	       tlj_Sel_Correlations_2D_response_.back()->Sumw2();

	       tlj_Sel_Correlations_2D_response2_.push_back(define_pt_histo( "tlj_Sel_Correlations_2D_response2_"+
	       Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+"_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first + 
               "_", ";"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+";"+  
               Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].second,
                s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].second));
	       tlj_Sel_Correlations_2D_response2_.back()->Sumw2();

	       tlj_Sel_Correlations_2D_mean_response_.push_back(define_pt_histo( "tlj_Sel_Correlations_2D_mean_response_"+
	       Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+"_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first + 
               "_", ";"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+";"+  
               Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].second,
                s_phi_bins_x,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first,corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].second));
	       tlj_Sel_Correlations_2D_mean_response_.back()->Sumw2();

	     }
	   tlj_Sel_Correlations_2D_one_corr_counts_helper_.push_back(tlj_Sel_Correlations_2D_counts_helper_);
	   tlj_Sel_Correlations_2D_one_corr_counts_.push_back(tlj_Sel_Correlations_2D_counts_);
	   tlj_Sel_Correlations_2D_one_corr_response_.push_back(tlj_Sel_Correlations_2D_response_);
	   tlj_Sel_Correlations_2D_one_corr_response2_.push_back(tlj_Sel_Correlations_2D_response2_);
	   tlj_Sel_Correlations_2D_one_corr_mean_response_.push_back(tlj_Sel_Correlations_2D_mean_response_);
	   //	 }
     }
   tlj_Sel_Correlations_2D_counts_helper_all_.push_back(tlj_Sel_Correlations_2D_one_corr_counts_helper_);
   tlj_Sel_Correlations_2D_counts_all_.push_back(tlj_Sel_Correlations_2D_one_corr_counts_);
   tlj_Sel_Correlations_2D_response_all_.push_back(tlj_Sel_Correlations_2D_one_corr_response_);
   tlj_Sel_Correlations_2D_response2_all_.push_back(tlj_Sel_Correlations_2D_one_corr_response2_);
   tlj_Sel_Correlations_2D_mean_response_all_.push_back(tlj_Sel_Correlations_2D_one_corr_mean_response_);
   }


}

void flex_res::Write_Histos()
{


 root_resol_name_binning= "histos_resolution_" + X_labels_[bin_choice] + ".root";

   TFile *outf = new TFile(root_resol_name_binning,"RECREATE");

     
    TString saveFolder = "Resolution_Plots";
    if(chdir(saveFolder) != 0){
      mkdir(saveFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(saveFolder);
    }

    cout <<"created folder, writing corrected response"<< endl;

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->Write();
	  cout << tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetName() << endl;
	  //	  tlj_X_response_all_ tlj_X_response_all_corrections_.back()->Draw();
	}
      draw_TH1D_save_with_Gauss_Fit_PS("PS_Corrected_Response",img_choice, tlj_corrected_response_barrel_all_[Corr_i], X_labels_[bin_choice] +tlj_corrected_response_barrel_all_[Corr_i].at(0)->GetName(),"nice");
      
    }
   
    cout <<" writing correction functions"<< endl;


  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      draw_TH1D_save_PS("PS_Correction_Functions",img_choice,tlj_X_response_all_corrections_all_[Corr_i], X_labels_[bin_choice] +"tlj_X_response_all_corrections_all_"+Corr_labels_[Corr_i].first,"nice");  
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->Write();
	  //	  tlj_X_response_all_ tlj_X_response_all_corrections_.back()->Draw();
	}


    }

  for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
    {
      //			 std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_;
      for(unsigned int Sel_Corr_j=Sel_Corr_i+1;Sel_Corr_j<Corr_selected_labels_.size();Sel_Corr_j++)
	{
	  //			     std::vector < TH2D* >  tlj_Sel_Correlations_2D_;

	  for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	    {
	      //	      tlj_Sel_Correlations_2D_counts_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Rebin2D();

	      tlj_Sel_Correlations_2D_counts_helper_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Rebin2D(4);
	      tlj_Sel_Correlations_2D_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Rebin2D(4);
	      tlj_Sel_Correlations_2D_response2_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Rebin2D(4);
	      tlj_Sel_Correlations_2D_mean_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Rebin2D(4);

	      tlj_Sel_Correlations_2D_counts_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Write();
	      tlj_Sel_Correlations_2D_counts_helper_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Write();
	      tlj_Sel_Correlations_2D_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Write();
	      tlj_Sel_Correlations_2D_response2_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Write();

	      tlj_Sel_Correlations_2D_mean_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Divide(tlj_Sel_Correlations_2D_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i),tlj_Sel_Correlations_2D_counts_helper_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i));
	      tlj_Sel_Correlations_2D_mean_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Write();

	      TMatrixDSym cov(1, n);
	      
	      for(unsigned int i = 1; i <= 2; i++)
		for(unsigned int j = 1; j <= 2; j++)
		  cov(i, j) = tlj_Sel_Correlations_2D_counts_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->GetCovariance(i,j);
	      std:: cout << tlj_Sel_Correlations_2D_counts_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->GetName() << std::endl;
	      std::cout << "Covariance matrix:" << std::endl;
	      cov.Print();
	      TDecompSVD svd(cov);
	      TMatrixD rotation = svd.GetU();
	      std::cout << "Matrix to diagonalize covariance:" << std::endl;
	      rotation.Print();

	    }

	      draw_TH2D_save_PS("PS_Correlation_Plots_counts",true,img_choice,tlj_Sel_Correlations_2D_counts_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1), X_labels_[bin_choice] +" tlj_Sel_Correlations_counts_"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first + "_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first ,"nice","colz");

	      gStyle->SetOptStat(0);
	      // gRoot->SetStatBox(0);
	      draw_TH2D_save_PS("PS_Correlation_Plots_mean_response",false,img_choice,tlj_Sel_Correlations_2D_mean_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1), X_labels_[bin_choice] +" tlj_Sel_Correlations_mean_response_"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first + "_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first ,"nice","colz","",0,-1,0,-1,0.8,1.2);
	      gStyle->SetOptStat(1);
	      //	      gRoot->SetStatBox(1);
	      //void draw_TH2D_save_PS(TString SaveTo,Bool_t with_correlation_factors, TString img_exp, std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {
	}
    }

  chdir("..");


   outf->Close();


}

void flex_res::Write_TGraphErrors()
{

   TFile *outf = new TFile(root_resol_name_binning,"UPDATE");

    TString saveFolder = "Resolution_Plots";
    if(chdir(saveFolder) != 0){
      mkdir(saveFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(saveFolder);
    }

      draw_TGraphErrors_save_PS("PS_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_mean_all_,X_labels_[bin_choice] +tlj_Response_Graphs_mean_all_.at(0)->GetName(),"nice","ALP","no","x1");
      draw_TGraphErrors_save_PS("PS_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_mean_all_,X_labels_[bin_choice] +tlj_Response_Graphs_mean_all_.at(0)->GetName(),"nice");
      draw_TGraphErrors_save_PS("PS_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_sigma_all_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_sigma_all_.at(0)->GetName(),"nice","ALP", "yes","x1",0,-1,0,0.5);
      draw_TGraphErrors_save_PS("PS_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_sigma_all_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_sigma_all_.at(0)->GetName(),"nice","ALP", "yes","x0",0,-1,0,0.5);
      draw_TGraphErrors_save_PS("PS_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_resol_all_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_resol_all_.at(0)->GetName(),"nice","ALP", "no","x1");
      draw_TGraphErrors_save_PS("PS_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_resol_all_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_resol_all_.at(0)->GetName(),"nice","ALP", "no","x0");

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      tlj_Response_Graphs_mean_all_.at(Corr_i)->Write();
      tlj_Response_Graphs_rel_sigma_all_.at(Corr_i)->Write();
      tlj_Response_Graphs_rel_resol_all_.at(Corr_i)->Write();
    }


  if(Corr_selected_labels_.size()>0)
{
      draw_TGraphErrors_save_PS("PS_Selected_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_mean_selec_,X_labels_[bin_choice] +tlj_Response_Graphs_mean_selec_.at(0)->GetName(),"nice","ALP","no","x1");
      draw_TGraphErrors_save_PS("PS_Selected_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_mean_selec_,X_labels_[bin_choice] +tlj_Response_Graphs_mean_selec_.at(0)->GetName(),"nice","ALP","no","x0");
      draw_TGraphErrors_save_PS("PS_Selected_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_sigma_selec_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_sigma_selec_.at(0)->GetName(),"nice","ALP","yes","x0",0,-1,0,0.5);
      draw_TGraphErrors_save_PS("PS_Selected_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_sigma_selec_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_sigma_selec_.at(0)->GetName(),"nice","ALP","yes","x1",0,-1,0,0.5);
      draw_TGraphErrors_save_PS("PS_Selected_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_resol_selec_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_resol_selec_.at(0)->GetName(),"nice","ALP","no","x0");
      draw_TGraphErrors_save_PS("PS_Selected_Mean_RELSIGMA_RELRES",img_choice,tlj_Response_Graphs_rel_resol_selec_,X_labels_[bin_choice] +tlj_Response_Graphs_rel_resol_selec_.at(0)->GetName(),"nice","ALP","no","x1");
}


   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
     {
       tlj_Response_Graphs_mean_selec_.at(Sel_Corr_i)->Write();
       tlj_Response_Graphs_rel_sigma_selec_.at(Sel_Corr_i)->Write();
       tlj_Response_Graphs_rel_resol_selec_.at(Sel_Corr_i)->Write();
     }



 TLegend *leg_all;
 leg_all = new TLegend(0.1,0.1,0.48,0.3);
 //   leg->SetHeader("Legende");
 for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
   {
   leg_all->AddEntry(tlj_Response_Graphs_mean_all_[Corr_i],Corr_labels_[Corr_i].first,"lep");
   }

 TLegend *leg_selected;
 leg_selected = new TLegend(0.1,0.1,0.48,0.3);
 //   leg->SetHeader("Legende");
 for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
   {
   leg_selected->AddEntry(tlj_Response_Graphs_mean_selec_[Sel_Corr_i],Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first,"lep");
   }



 draw_graphs( tlj_Response_Graphs_mean_all_, 0., 1.3, leg_all,  X_labels_[bin_choice] +"RESPONSE_res_All_response_graph_mean_w_mean_e");
 draw_graphs( tlj_Response_Graphs_rel_sigma_all_, -0.1, 0.5, leg_all,  X_labels_[bin_choice] +"RESPONSE_res_All_response_graph_rel_sigma");
 draw_graphs( tlj_Response_Graphs_rel_resol_all_, 0.6, 1.1, leg_all, X_labels_[bin_choice] + "RESPONSE_res_All_response_graph_rel_sigma_relative_to_L2L3");


 TString selection;
   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
     {
       selection+=Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first +"-";
     }
   cout << selection << endl;


  if(Corr_selected_labels_.size()>0)
{
 draw_graphs( tlj_Response_Graphs_mean_selec_, 0.85, 1.15, leg_selected,  X_labels_[bin_choice] +"RESPONSE_res_selected_" + selection + "_response_graph_mean");
 draw_graphs( tlj_Response_Graphs_rel_sigma_selec_, -0.1, 0.5, leg_selected,  X_labels_[bin_choice] +"RESPONSE_res_selected_" + selection + "_response_graph_rel_sigma");
 draw_graphs( tlj_Response_Graphs_rel_resol_selec_, 0.6, 1.1, leg_selected,  X_labels_[bin_choice] +"RESPONSE_res_selected_" + selection + "_response_graph_rel_sigma_relative_to_L2L3");
}


  chdir("..");

   outf->Close();

}

void flex_res::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L flex_res.C
//      Root > flex_res t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;




   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   Declare_Labels();

   cout << "Bitte Binningvariable wählen (ACHTUNG: I.A. sollte die Binningvariable hier die gleiche sein wie beim Rausschreiben der Korrekturplots (bisher Genpt)!!!):" <<endl;
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>bin_choice;
   fflush(stdin);
   cout << "Es wurde " << bin_choice << ", also " << X_labels_[bin_choice]<< ", ausgewaehlt. Danke..." << endl;

   

   cout << "Bitte Parametrisierung wählen:" <<endl;
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>X_choice;
   fflush(stdin);
   cout << "Es wurde " << X_choice << ", also " << X_labels_[X_choice]<< ", ausgewaehlt. Danke..." << endl;



   cout << "Bitte Korrekturen zum Vergleich auswaehlen:" <<endl;
   //   cout << "Test: " << (char)66 << endl;

   //   cout << Corr_labels_[2].first << endl;

   for(Int_t i=0;i<no_all_Corr_labels_;i++)  
     {
       cout << (char)(i+65) << ": " << Corr_labels_[i].first << endl;
     }
   cin >> Corr_choice;
   fflush(stdin);
   cout << "Es wurde " << Corr_choice << " ausgewaehlt. Danke... Gewaehlt wurde also: " << endl;


   for(Int_t i=0;i<no_all_Corr_labels_;i++)  
     {
       if(Corr_choice.Contains((char)(i+65)))
	 {
	   cout << (char)(i+65) << ": " << Corr_labels_[i].first << endl;
	   Corr_selected_labels_.push_back(i);
	 }
     }
   if(Corr_selected_labels_.size()==0)cout << "Es wurde nichts ausgewaehlt!!!" << endl;

   cout << "Bitte Exportformat für Bilder wählen ('A' für keinen Export)" <<endl;
   cin >>img_choice;
   fflush(stdin);
   cout << "Es wurde " << img_choice << " ausgewaehlt. Danke..." << endl;

   Import_Histos();
   Book_Histos();





   
   Long64_t nentries = fChain->GetEntriesFast();


   if(entries_to_run>0)nentries=entries_to_run;
   //           nentries=200000;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if(jentry % 1000 == 0) printf("event %d\n", jentry); //Alle tausend Events was sagen
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //Jet quality variables...
      Bool_t towers_gt_1=0;
      Double_t genjet_jet_deltar=0;


      //loop over two leading jets and 
      //1. calculate jet quality variables
      //2. fill histos with correction variables of good jets in (pt)-bins

          for(Int_t genjet_i =0; genjet_i<2;genjet_i++)//über erste zwei colgenjet
	{
	  //DeltaR ausrechnen (Gen -Jet Nummerierung und passender CaloJet)
	  genjet_jet_deltar= delta_r(GenJetColEta[genjet_i],JetEta[GenJetColJetIdx[genjet_i]],GenJetColPhi[genjet_i],JetPhi[GenJetColJetIdx[genjet_i]]); 

	  Bool_t isolated=1;//Isolation ueberpruefen
	  for(Int_t caljet_i=0;caljet_i<NobjJet;caljet_i++)
	    {
		  if(caljet_i!=GenJetColJetIdx[genjet_i]&&delta_r(GenJetColEta[genjet_i],JetEta[caljet_i],GenJetColPhi[genjet_i],JetPhi[caljet_i])<iso_max)
		    {
		      isolated=0;
		    }
	    }

	

	  
	  towers_gt_1=0;
	  Int_t no_towers=0;

	  for(Int_t a_i=0;a_i<NobjTow;a_i++) // Ueber alle Tower loopen und alle EM und Had-Tower zum jeweiligen GenJet gehoeren zusammenzaehlen
	    {
	      if(Tow_jetidx[a_i]==GenJetColJetIdx[genjet_i])
		{
		  no_towers++;
		}
	    }
	  if(no_towers>1)towers_gt_1=1;// ueberpruefe explizit, dass Jet mehr als einen Tower hat...


	  if(std::abs(GenJetColEta[genjet_i])<1.1&&genjet_jet_deltar<0.25&&isolated&&towers_gt_1)//demand Eta-selection, deltar (genjet,jet) <.25, isolated and more than one tower per jet
	    {
	      Int_t match=GenJetColJetIdx[genjet_i];


	      //Get variables needed for filling...
       	      base_corr::Fill_obvious_vars(genjet_i, match);
	      std::vector < Double_t > X_Var_=get_X_Vars_(genjet_i,match);
	      std::vector < Double_t > Correction_Var_=get_Correction_Vars_(genjet_i,match);

	     for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	       {//Jeder der zwei l. Jet wird dem pt-bin zugeordnet, der zu ihm passt
		 if(X_Var_[bin_choice]>pt_bins_[pt_i].first&&X_Var_[bin_choice]<pt_bins_[pt_i].second)
		   {
		     std::vector < Double_t > Corrected_Response_;
		     Corrected_Response_.push_back(_JetResponse);				
		     Corrected_Response_.push_back(_L2L3JetResponse);
		     for(Int_t Corr_i=2;Corr_i<no_Corr_labels_;Corr_i++)
		       {
			 Corrected_Response_.push_back(_L2L3JetResponse/tlj_X_response_prof_all_[Corr_i].at(pt_i)->GetBinContent(tlj_X_response_prof_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(Correction_Var_[Corr_i])));				
		       }
		     Corrected_Response_.push_back(_JetResponse);	
		     Corrected_Response_.push_back(_L2L3JetResponse);
		     for(Int_t Corr_i=2+no_Corr_labels_;Corr_i<no_all_Corr_labels_;Corr_i++)
		       {
			 Corrected_Response_.push_back(_L2L3JetResponse/tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->GetBinContent(tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(Correction_Var_[Corr_i-no_Corr_labels_])));				
		       }



		     //  cout << "test" << endl;
		     for(Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
		       {
			 tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->Fill(Corrected_Response_[Corr_i]);
			 //	   tlj_X_response_2D_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],L2L3JetResponse_.back());
			 //			 tlj_X_response_prof_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],L2L3JetResponse_.back());
			 //  cout << "test" << endl;
		       }


		     for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
		       {
			 //			 std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_;
			 for(unsigned int Sel_Corr_j=Sel_Corr_i+1;Sel_Corr_j<Corr_selected_labels_.size();Sel_Corr_j++)
			   {
			     //			     std::vector < TH2D* >  tlj_Sel_Correlations_2D_;
			     if(Corr_selected_labels_[Sel_Corr_i]<no_Corr_labels_&&Corr_selected_labels_[Sel_Corr_j]<no_Corr_labels_)
			       {
				 tlj_Sel_Correlations_2D_counts_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],Correction_Var_[Corr_selected_labels_[Sel_Corr_j]]);


				 tlj_Sel_Correlations_2D_counts_helper_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],Correction_Var_[Corr_selected_labels_[Sel_Corr_j]]);
				 tlj_Sel_Correlations_2D_counts_helper_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first+0.01,Correction_Var_[Corr_selected_labels_[Sel_Corr_j]]);
				 tlj_Sel_Correlations_2D_counts_helper_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first+0.01);


				 tlj_Sel_Correlations_2D_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],Correction_Var_[Corr_selected_labels_[Sel_Corr_j]],_L2L3JetResponse);
				 tlj_Sel_Correlations_2D_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first+0.01,Correction_Var_[Corr_selected_labels_[Sel_Corr_j]],_L2L3JetResponse);
				 tlj_Sel_Correlations_2D_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first+0.01,_L2L3JetResponse);



				 tlj_Sel_Correlations_2D_response2_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],Correction_Var_[Corr_selected_labels_[Sel_Corr_j]],_L2L3JetResponse*_L2L3JetResponse);
				 tlj_Sel_Correlations_2D_response2_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_i]].first+0.01,Correction_Var_[Corr_selected_labels_[Sel_Corr_j]],_L2L3JetResponse*_L2L3JetResponse);
				 tlj_Sel_Correlations_2D_response2_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],corr_var_x_edges_[Corr_selected_labels_[Sel_Corr_j]].first+0.01,_L2L3JetResponse*_L2L3JetResponse);
			       }
			     else
			       {
				 tlj_Sel_Correlations_2D_counts_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(0.5,0.5);	 
				 tlj_Sel_Correlations_2D_counts_helper_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(0.5,0.5);	 
				 tlj_Sel_Correlations_2D_response_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(0.5,0.5);
				 tlj_Sel_Correlations_2D_response2_all_[Sel_Corr_i].at(Sel_Corr_j-Sel_Corr_i-1).at(pt_i)->Fill(0.5,0.5);
			       }
			   }
		       }


		   }
	       }



	    }

	}



   }


   std::cout << "Funzt...auch hier :)" << std::endl;

   gStyle->SetOptStat(1);
   gStyle->SetOptFit(1);
   gStyle->SetPalette(1);

   test->Range(-1.5,-0.625,13.5,5.625);
   test->SetBorderSize(2);
   test->SetFrameFillColor(0);
//   test->SetGridy();
   test->SetSelected(test);





   std::vector < std::vector < std::vector < Double_t > > > Double_gauss_all_;
   std::vector < std::vector < std::vector < Double_t > > > Rel_resol_all_;

   for(Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
     {//Jeder der zwei l. Jet wird dem pt-bin zugeordnet, der zu ihm passt

       std::vector < std::vector < Double_t > > Double_gauss_;
       std::vector < std::vector < Double_t > > Rel_resol_;

        for(Int_t pt_i=no_small_pt_bins_;pt_i<no_pt_bins_;pt_i++)
 	 {
 	   if(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetEntries()>150)doublefit_gaus(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i));
 	 }
       for(Int_t pt_i=0;pt_i<no_small_pt_bins_;pt_i++)
	 {
	   cout << "Name: " << tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetName()  << " Anzahl der Einträge: " << tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetEntries() << " Corr_Nummer: " << Corr_i << "  Pt_i " << pt_i << endl;

	   //	   tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->Draw();
	   if(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetEntries()>150)
	     {
	   cout << "Name: " << tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetName()  << "Anzahl der Einträge: " << tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetEntries() << " Corr_Nummer: " << Corr_i << "  Pt_i " << pt_i << endl;
	       Double_gauss_.push_back(doublefit_gaus(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)));
	       std::vector < Double_t > rel_resol_pair_;
	       rel_resol_pair_.push_back(doublefit_gaus(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)).at(2)/doublefit_gaus(tlj_corrected_response_barrel_all_[1].at(pt_i)).at(2));
	       rel_resol_pair_.push_back(doublefit_gaus(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)).at(3)/doublefit_gaus(tlj_corrected_response_barrel_all_[1].at(pt_i)).at(2));

	       Rel_resol_.push_back(rel_resol_pair_);
	     }
	   else break;
	 }
       Double_gauss_all_.push_back(Double_gauss_);
       Rel_resol_all_.push_back(Rel_resol_);
     }



   for(Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
     {
       //Mean, default X
       tlj_Response_Graphs_mean_all_.push_back(make_graph(Double_gauss_all_[Corr_i], 0, X_labels_[X_choice] +"vs.Mean_" + Corr_labels_[Corr_i].first, X_choice));
       //Rel. Sigma, default X
       tlj_Response_Graphs_rel_sigma_all_.push_back(make_graph(Double_gauss_all_[Corr_i], 2, X_labels_[X_choice] +"vs.#sigma_{Rel}_" + Corr_labels_[Corr_i].first, X_choice));
       //Rel. Sigma relative to L2L3, default X
       tlj_Response_Graphs_rel_resol_all_.push_back(make_graph(Rel_resol_all_[Corr_i], 0, X_labels_[X_choice] +"vs.#sigma_{Rel}_over_#sigma_{Rel,L2L3}_" + Corr_labels_[Corr_i].first, X_choice));
       //       tlj_Response_Graphs_rel_resol_all_[Corr_i]->GetYaxis()->SetTitle("BLFH");
       //       tlj_Response_Graphs_rel_resol_all_.back()->GetYaxis()->SetTitle("#frac{#sigma}{Mean} / #frac{#sigma_{L2L3}}{Mean_{L2L3}}");

       if(Corr_i>=no_Corr_labels_)
	 {
	   tlj_Response_Graphs_mean_all_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_sigma_all_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_resol_all_.back()->SetLineStyle(2);
	 }
     }
   
   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
     {
       //Personal selection Mean, default X
       tlj_Response_Graphs_mean_selec_.push_back(make_graph(Double_gauss_all_[Corr_selected_labels_[Sel_Corr_i]], 0, "Selected_" + X_labels_[X_choice] +"vs.Mean_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first, X_choice));
       //Rel. Sigma, default X
       tlj_Response_Graphs_rel_sigma_selec_.push_back(make_graph(Double_gauss_all_[Corr_selected_labels_[Sel_Corr_i]], 2, "Selected_" + X_labels_[X_choice] +"vs.#sigma_{Rel}_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first, X_choice));
       //Rel. Sigma relative to L2L3, default X
       tlj_Response_Graphs_rel_resol_selec_.push_back(make_graph(Rel_resol_all_[Corr_selected_labels_[Sel_Corr_i]], 0, "Selected" + X_labels_[X_choice] +"vs.#sigma_{Rel}_over_#sigma_{Rel,L2L3}_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first, X_choice));

       //       tlj_Response_Graphs_rel_resol_selec_[Sel_Corr_i]->GetYaxis()->SetTitle("BLAH");
       //   gr->GetYaxis()->SetTitle("Y title");
       if(Corr_selected_labels_[Sel_Corr_i]>=no_Corr_labels_)
	 {
	   tlj_Response_Graphs_mean_selec_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_sigma_selec_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_resol_selec_.back()->SetLineStyle(2);
	 }

     }
   



   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!
         Write_Histos();
         Write_TGraphErrors();







}
