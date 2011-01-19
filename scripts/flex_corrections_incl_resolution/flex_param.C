#define flex_param_cxx
#include "THelpers.h"
#include "base_corr.C"
#include "flex_param.h"

void flex_param::Import_Histos()
{
  TString root_name_binning= "histos_" + eta_region_labels_[eta_choice] + X_labels_[bin_choice] + ".root";

  cout << "DEBUG: " << root_name_binning << endl;

  TFile* _file=new TFile(root_name_binning,"OPEN");

  //    TFile f("file.root");
    if (_file->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;


  for (Int_t X_i=0;X_i<no_X_labels_;X_i++)
    {
      std::vector < TH1D* >  tlj_X_counts_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   //	   cout << define_pt_histo_name("tlj_counts_"+(X_labels_[X_i]+"_"), pt_bins_, pt_i).Data() << endl;
      	   tlj_X_counts_.push_back((TH1D*)_file->Get(define_pt_histo_name("tlj_counts_"+(X_labels_[X_i]+"_"), pt_bins_, pt_i).Data()));
	 }
       tlj_X_counts_all_.push_back(tlj_X_counts_);
	   if(titleprint==0){
	     tlj_X_counts_.back()->SetTitle("");
	   }
    }

  cout << "loaded counts" << endl;

  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < TH2D* >  tlj_X_response_2D_;
      std::vector < TProfile* >  tlj_X_response_prof_;
      std::vector < TH1D* >  tlj_X_response_GMP_mean_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   //  cout << "trying to load response0 " << endl;
  cout <<define_pt_histo_name( "tlj_X_response_2D_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i) << endl;
	   tlj_X_response_2D_.push_back((TH2D*)_file->Get(define_pt_histo_name( "tlj_X_response_2D_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i)));
	   //	   cout <<tlj_X_response_2D_.back()->GetName() << endl;
	   //  cout << "trying to load response1 " << endl;
	   tlj_X_response_prof_.push_back((TProfile*)_file->Get(define_pt_histo_name( "tlj_X_response_prof_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i)));
	   TString temp =  (TString) tlj_X_response_2D_.back()->GetName() + "_1";
	   //  cout << "trying to load response2 " << endl;
	   tlj_X_response_GMP_mean_.push_back((TH1D*)_file->Get(temp.Data()));

	   if(titleprint==0){
	     tlj_X_response_GMP_mean_.back()->SetTitle("");
	     tlj_X_response_2D_.back()->SetTitle("");
	     tlj_X_response_prof_.back()->SetTitle("");
	   }
	 }

              tlj_X_response_2D_all_.push_back(tlj_X_response_2D_);
              tlj_X_response_prof_all_.push_back(tlj_X_response_prof_);
              tlj_X_response_GMP_mean_all_.push_back(tlj_X_response_GMP_mean_);

    }

  cout << "loaded response " << endl;



  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* > tlj_X_response_all_corrections_;
      std::vector < TH1D* > tlj_Corr_Vars_counts_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   if(Corr_i<no_Corr_labels_) tlj_X_response_all_corrections_.push_back((TH1D*)tlj_X_response_prof_all_[Corr_i].at(pt_i));
	   else tlj_X_response_all_corrections_.push_back(tlj_X_response_GMP_mean_all_[Corr_i-no_Corr_labels_].at(pt_i));
	   
	   //	   cout << "ACHTUNG!!!!!!!!!!!!!!" << tlj_X_response_all_corrections_.back()->GetName() << endl;
	   if(Corr_i<no_Corr_labels_)tlj_Corr_Vars_counts_.push_back((TH1D*)_file->Get(define_pt_histo_name( "tlj_Corr_Vars_counts_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i)));
	   else tlj_Corr_Vars_counts_.push_back((TH1D*)_file->Get(define_pt_histo_name( "tlj_Corr_Vars_counts_"+Corr_labels_[Corr_i-no_Corr_labels_].first+"_", pt_bins_, pt_i)));

	   if(titleprint==0){
	     tlj_Corr_Vars_counts_.back()->SetTitle("");
	     tlj_X_response_all_corrections_.back()->SetTitle("");
	   }

	 }
       tlj_X_response_all_corrections_all_.push_back(tlj_X_response_all_corrections_);
       tlj_Corr_Vars_counts_all_.push_back(tlj_Corr_Vars_counts_);
    }


  //  tlj_X_response_prof_all_[2].at(5)->Draw();
  //  tlj_X_response_prof_all_[2].at(5)->Fit("[0]*(x-[1]) + [2]","EM","same",0,1);
  //  tlj_X_response_prof_all_[2].at(5)->Fit("pol2","EM","same",0,1);


  /*
  TFile* _file2=new TFile("Fitted_paras.root","OPEN");

  //    TFile f("file.root");
    if (_file2->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;
  */
}


void flex_param::Book_Histos()
{



}

void flex_param::Fit_Histos()
{

  Double_t first_bin_pos = 0.; 
  Double_t last_bin_pos = 1.;
  const Int_t nq = 2;

  Double_t xq[nq];  // position where to compute the quantiles in [0,1]
  Double_t yq[nq];  // array to contain the quantiles
	  

	  xq[0] = 0.005; xq[1]=0.995;


  for(Int_t pt_i=0;pt_i<no_small_pt_bins_;pt_i++)
    {
      if(tlj_Corr_Vars_counts_all_[0].at(pt_i)->GetEntries()>150)no_of_populated_bins=pt_i;
    }


  for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
    {
      std::vector < TF1* > fit_functions_;
      std::vector < TF1* > fit_functions_GMP_;
      std::vector < std::vector < Double_t > > Corr_coeffs_;
      //      cout << "test" << endl;


	  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
	    {
	      std::vector < Double_t > corr_value_and_error_;
	      tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetQuantiles(nq,yq,xq);
	      cout << yq[0] << " and " << yq[1] << endl;
	      first_bin_pos=yq[0];
	      last_bin_pos=yq[1];
	      
	      
	      fit_functions_.push_back(  new TF1(eta_region_labels_[eta_choice] + X_labels_[bin_choice] + "Fit_" + tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->GetName() ,
						 "1/([0]+(1-[3]*[4]*[4])+[2]*(x-[1])+[3]*pow((x-[1]),2))",corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second)) ;

	   if(titleprint==0){
	     fit_functions_.back()->SetTitle("");
	   }

	      //	fit_functions_.push_back(  new TF1(eta_region_labels_[eta_choice] + X_labels_[bin_choice] + "Fit_" + tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->GetName() ,"[0]+(1-[3]*[4]*[4])+[2]*(x-[1])+[3]*pow((x-[1]),2)",first_bin_pos,last_bin_pos)) ;
	      fit_functions_.back()->SetParameters(1,2,3,4,5);
	      fit_functions_.back()->SetParNames("Const","X0","B*(X-X0)","C*(X-X_0)^2","Sigma"); 
	      //	fit_functions_.back()->SetParNames(param_fit_labels_[0],param_fit_labels_[1],param_fit_labels_[2],param_fit_labels_[3]);
	      fit_functions_.back()->FixParameter(0,0.);
	      fit_functions_.back()->FixParameter(1,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetMean());
	      fit_functions_.back()->SetParError(1,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetMeanError());
	      cout <<"MEANERROR:" << tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetMeanError() << endl;
	      fit_functions_.back()->FixParameter(4,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetRMS());
	      fit_functions_.back()->SetParError(4,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetRMSError());
	      cout <<"RMSERROR:" << tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetRMSError() << endl;
	      fit_functions_.back()->SetLineColor(kRed);
	      cout <<	fit_functions_.back()->GetName() << endl;
	      
	      TFitResultPtr r;
	      r=tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->Fit(fit_functions_.back(),"ES","same",first_bin_pos,last_bin_pos);
	      if(r==0)
		{
		  cout << fit_functions_.back()->GetName() << endl;
		  for (Int_t i=0;i<4;i++)cout << fit_functions_.back()->GetParameter(i) << endl;
		  TMatrixDSym correlation = r->GetCorrelationMatrix();
		  cout << "Correlation matrix: " <<correlation(2,3)<< endl;
		  r->GetCorrelationMatrix().Print();
		  

		  TMatrixDSym cov = r->GetCovarianceMatrix();	
		  cout << "Covariance matrix:"  << endl;
		  r->GetCovarianceMatrix().Print();
		  
		  TDecompSVD svd(cov);
		  TMatrixD rotation = svd.GetU();
		  std::cout << "Matrix to diagonalize covariance:" << std::endl;
		  rotation.Print();
		  
		  corr_value_and_error_.push_back(correlation(2,3));
		  corr_value_and_error_.push_back(0.);
		  Corr_coeffs_.push_back(corr_value_and_error_);
		  //	    cov(2, 3)
		  
		}
	      else
		{	    corr_value_and_error_.push_back( 0.);
		corr_value_and_error_.push_back(0.);
		Corr_coeffs_.push_back(corr_value_and_error_);
		}
	      //Same for GMP
	      fit_functions_GMP_.push_back(  new TF1(eta_region_labels_[eta_choice] + X_labels_[bin_choice] + "Fit_GMP_" + tlj_X_response_all_corrections_all_[Corr_i+no_Corr_labels_].at(pt_i)->GetName() ,
						     "1/([0]+(1-[3]*[4]*[4])+[2]*(x-[1])+[3]*pow((x-[1]),2))",corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second)) ;
	      //	fit_functions_GMP_.push_back(  new TF1(eta_region_labels_[eta_choice] + X_labels_[bin_choice] + "Fit_GMP_" + tlj_X_response_all_corrections_all_[Corr_i+no_Corr_labels_].at(pt_i)->GetName() ,"[0]+(1-[3]*[4]*[4])+[2]*(x-[1])+[3]*pow((x-[1]),2)",first_bin_pos,last_bin_pos)) ;
	   if(titleprint==0){
	     fit_functions_.back()->SetTitle("");
	   }

	      fit_functions_GMP_.back()->SetParameters(1,2,3,4,5);
	      fit_functions_GMP_.back()->SetParNames("Const","X0","B*(X-X0)","C*(X-X_0)^2","Sigma");
	      fit_functions_GMP_.back()->FixParameter(0,0.);
	      fit_functions_GMP_.back()->FixParameter(1,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetMean());
	      fit_functions_GMP_.back()->SetParError(1,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetMeanError());
	      fit_functions_GMP_.back()->FixParameter(4,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetRMS());
	      fit_functions_GMP_.back()->SetParError(4,tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetRMSError());
	      fit_functions_GMP_.back()->SetLineColor(kRed);
	      cout <<	fit_functions_GMP_.back()->GetName() << endl;
	      
	      tlj_X_response_all_corrections_all_[Corr_i+no_Corr_labels_].at(pt_i)->Fit(fit_functions_GMP_.back(),"E","same",first_bin_pos,last_bin_pos);


	      //Eventuell "EMV" als Fit-OPtion nehmen...
	    }

	  Corr_coeffs_all_.push_back(Corr_coeffs_);
	  
	  for(unsigned int GMP_i=0;GMP_i<fit_functions_GMP_.size();GMP_i++)fit_functions_.push_back(fit_functions_GMP_[GMP_i]);

   
      fit_functions_all_.push_back(fit_functions_);
    }



}

void flex_param::Write_Histos()
{
TString root_resol_name_binning= "histos_parametrisations_" + eta_region_labels_[eta_choice] + X_labels_[bin_choice] + X_labels_[X_choice]+ ".root";

   TFile *outf = new TFile(root_resol_name_binning,"RECREATE");

     
    TString saveFolder = "Parametrization_Plots";
    if(chdir(saveFolder) != 0){
      mkdir(saveFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(saveFolder);
    }

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      cout << Corr_i+1 << " of " << no_all_Corr_labels_ << endl;
      draw_TH1D_save_PS("PS_Fitted_Corr_Vars_vs_Response",img_choice, tlj_X_response_all_corrections_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_X_response_all_corrections_all_[Corr_i].at(0)->GetName(),"nice");
      cout << "writing PS done... trying to save to ROOT " << endl;
      
      for(Int_t pt_i=0;pt_i<no_small_pt_bins_;pt_i++)
	{
	  tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->Write();
	}
      
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
      	  fit_functions_all_[pt_i].at(Corr_i)->Write();
	}
      


    }
   

  chdir("..");


   outf->Close();


}




void flex_param::Create_TGraphErrors()
{


  cout << "test....................." <<fit_functions_all_[4].at(5)->GetParameter(0) <<fit_functions_all_[4].at(5)->GetParError(0)<< endl;

  std::vector < std::vector < std::vector < Double_t > > > values_fit_params_in_X_all_;
  //  std::vector < std::vector < std::vector < Double_t > > > values_fit_params_in_X_GMP_all_;


  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      std::vector < std::vector < Double_t > > values_fit_params_in_X_;
      //      std::vector < std::vector < Double_t > > values_fit_params_in_X_GMP_;

      //      cout << "test" << endl;
      for(Int_t pt_i=0;pt_i<no_of_populated_bins;pt_i++)
	{
	  offset_for_X = 0;
	  if(tlj_X_counts_all_[X_choice].at(pt_i)->GetEntries()>150){
	  std::vector < Double_t >  values_fit_params_in_bins_in_X_;

	  for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
	    {
	      values_fit_params_in_bins_in_X_.push_back(fit_functions_all_[pt_i].at(Corr_i)->GetParameter(par_i));

	      if(par_i==1)
		{values_fit_params_in_bins_in_X_.push_back(tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetMeanError());}
	      else if(par_i==4)
		{values_fit_params_in_bins_in_X_.push_back(tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->GetRMSError());}
	      else values_fit_params_in_bins_in_X_.push_back(fit_functions_all_[pt_i].at(Corr_i)->GetParError(par_i));


	      cout << "test......." +(fit_functions_all_[pt_i].at(Corr_i)->GetName() + param_fit_labels_[par_i]) +  " Param: " <<fit_functions_all_[pt_i].at(Corr_i)->GetParameter(par_i) << " Error: " <<values_fit_params_in_bins_in_X_.back()<< endl;

	    }
	  values_fit_params_in_X_.push_back(values_fit_params_in_bins_in_X_);
	  }
	  else offset_for_X = pt_i;
	}
      values_fit_params_in_X_all_.push_back(values_fit_params_in_X_);

    }
  cout <<"sieht gut aus..."<< endl;


  //"Const","X0","B*(X-X0)","C*(X-X_0)^2"
//   param_fit_labels_.push_back("Const");
//   param_fit_labels_.push_back("Const");
//   param_fit_labels_.push_back("X0");
//   param_fit_labels_.push_back("X0");
//   param_fit_labels_.push_back("B_X-X0_");
//   param_fit_labels_.push_back("B_X-X0_");
//   param_fit_labels_.push_back("C_X-X0_2");
//   param_fit_labels_.push_back("C_X-X0_2");

//[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),2)+[4]*pow(log(x),2)
  for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
    {
      std::vector < TGraphErrors* > fit_params_in_X_;
      cout <<"sieht gut aus1..."<< endl;
      //      std::vector < TGraphErrors* > fit_params_in_X_GMP_;
      for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
	{
	  cout <<"sieht gut aus2..."<< Corr_labels_[Corr_i].first <<endl;
	  //TGraphErrors* base_corr::make_graph(std::vector < std::vector < Double_t > >  Double_gauss_, Int_t y_para, TString title, Int_t X_par)
	  
	  
	  fit_params_in_X_.push_back(make_graph( param_fit_labels_,values_fit_params_in_X_all_[Corr_i], par_i, X_labels_[X_choice] +"vs_" + param_fit_labels_[par_i] + "_of_" + Corr_labels_[Corr_i].first, X_choice,offset_for_X));
	 
// 	   if(titleprint==0){
// 	     fit_params_in_X_.back()->SetTitle("");
// 	   }

 
	  //	      fit_params_in_X_.push_back(make_graph(values_fit_params_in_X_all_[Corr_i], 2*par_i, X_labels_[X_choice] +"vs." + param_fit_labels_[par_i] + "_of_" + Corr_labels_[Corr_i].first, X_choice));
	  cout <<"sieht gut aus3..."<< endl;
	  


	}
      fit_params_in_X_all_.push_back(fit_params_in_X_);
    }


      for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
	{
//   param_fit_labels_.push_back("Const");
//   param_fit_labels_.push_back("X0");
//   param_fit_labels_.push_back("B_X-X0_");
//   param_fit_labels_.push_back("C_X-X0_2");
//   param_fit_labels_.push_back("Sigma");


	  fit_param_functions_X0_.push_back(  new TF1(eta_region_labels_[eta_choice] + X_labels_[bin_choice] + "Fit_" + param_fit_labels_[1] +  "_of_" + Corr_labels_[Corr_i].first ,
	  "[0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)",20,800)) ;
	   if(titleprint==0){
	     fit_param_functions_X0_.back()->SetTitle("");
	   }

	  fit_param_functions_X0_.back()->SetParameters(1,2,3,4);
	  fit_param_functions_X0_.back()->SetParNames("a","b","c","d"); 
	  fit_param_functions_X0_.back()->SetLineColor(kRed);
	  cout << fit_param_functions_X0_.back()->GetName() << endl;

	  fit_params_in_X_all_[1].at(Corr_i)->Fit(fit_param_functions_X0_.back(),"ES","same",20,800);


	}




  std::vector < std::vector < std::vector < Double_t > > > new_Corr_coeffs_all_;

  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < std::vector < Double_t > > new_Corr_coeffs_;
      //      std::vector < std::vector < Double_t > > values_fit_params_in_X_GMP_;

      //      cout << "test" << endl;
      for(Int_t pt_i=0;pt_i<no_small_pt_bins_;pt_i++)
	{
	  cout << "TEST "<< pt_i <<" of " << fit_functions_all_.size() << " and " << Corr_i  <<" of " << fit_functions_all_[pt_i].size() << " " << fit_functions_all_[pt_i].at(Corr_i)->GetName() <<endl;
	  cout << Corr_coeffs_all_[pt_i].at(Corr_i).at(0) << " and " <<Corr_coeffs_all_[pt_i].at(Corr_i).at(1) << endl; 
       	  new_Corr_coeffs_.push_back(Corr_coeffs_all_[pt_i].at(Corr_i));


	}
      new_Corr_coeffs_all_.push_back(new_Corr_coeffs_);
    }

  cout << "funzt noch" << endl;
  correlation_labels_.push_back("Correlation_BC_");

   for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
     {
       Correlation_B_C_.push_back(make_graph( correlation_labels_,new_Corr_coeffs_all_[Corr_i],0, X_labels_[X_choice] +"vs_" + correlation_labels_[0] + "_of_" + Corr_labels_[Corr_i].first, X_choice));

// 	   if(titleprint==0){
// 	     Correlation_B_C_.back()->SetTitle("");
// 	   }

     }

}

void flex_param::Write_TGraphErrors()
{
TString root_resol_name_binning= "histos_parametrisations_" + eta_region_labels_[eta_choice] + X_labels_[bin_choice] + X_labels_[X_choice] + ".root";

   TFile *outf = new TFile(root_resol_name_binning,"UPDATE");

    TString saveFolder = "Parametrization_Plots";
    if(chdir(saveFolder) != 0){
      mkdir(saveFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(saveFolder);
    }
// void draw_TGraphErrors_save_PS(std::vector < TGraphErrors* > histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="ALP", TString resolution_fit = "no", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS") {


    for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
      {
	cout << fit_params_in_X_all_[par_i].at(0)->GetName() << endl;
	draw_TGraphErrors_save_PS("PS_Fit_params_"+param_fit_labels_[par_i] ,img_choice, fit_params_in_X_all_[par_i],eta_region_labels_[eta_choice] + X_labels_[bin_choice] +fit_params_in_X_all_[par_i].at(0)->GetName(),"nice","ALP","no","x1");
	draw_TGraphErrors_save_PS("PS_Fit_params_"+param_fit_labels_[par_i] ,img_choice, fit_params_in_X_all_[par_i],eta_region_labels_[eta_choice] + X_labels_[bin_choice] +fit_params_in_X_all_[par_i].at(0)->GetName(),"nice");
      }
    draw_TGraphErrors_save_PS("Correlation_B_C_" ,img_choice, Correlation_B_C_,eta_region_labels_[eta_choice] + X_labels_[bin_choice] +Correlation_B_C_.at(0)->GetName(),"nice","ALP","no","x1");
    draw_TGraphErrors_save_PS("Correlation_B_C_" ,img_choice, Correlation_B_C_,eta_region_labels_[eta_choice] + X_labels_[bin_choice] +Correlation_B_C_.at(0)->GetName(),"nice");




    cout <<"kurz vor schluss" << endl;
    for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
      {
	
	for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
	  {
	    	    fit_params_in_X_all_[par_i].at(Corr_i)->Write();
	  }
      }

    for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
      {
		Correlation_B_C_.at(Corr_i)->Write();
      }
    for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
      {
		fit_param_functions_X0_.at(Corr_i)->Write();
      }



 std::vector < TLegend * > leg_selected_;
 //   leg->SetHeader("Legende");
 for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
   {
     leg_selected_.push_back(new TLegend(0.22,0.21,0.49,0.41));
     leg_selected_.back()->SetFillColor(0);
  leg_selected_.back()->SetBorderSize(0);
  leg_selected_.back()->SetFillColor(0);
  leg_selected_.back()->SetTextFont(42);
  leg_selected_.back()->SetTextSize(0.04);
     for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
       {
	 leg_selected_.back()->AddEntry(fit_params_in_X_all_[par_i].at(Corr_selected_labels_[Sel_Corr_i]),Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+param_fit_labels_[par_i],"lep");
       }
   }

 std::vector < TLegend * > leg_params_;
 //   leg->SetHeader("Legende");
 for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
   {
     leg_params_.push_back(new TLegend(0.22,0.21,0.49,0.41));
     leg_params_.back()->SetFillColor(0);
  leg_params_.back()->SetBorderSize(0);
  leg_params_.back()->SetFillColor(0);
  leg_params_.back()->SetTextFont(42);
  leg_params_.back()->SetTextSize(0.04);
     for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
       {
	 leg_params_.back()->AddEntry(fit_params_in_X_all_[par_i].at(Corr_selected_labels_[Sel_Corr_i]),Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+param_fit_labels_[par_i],"lep");
       }
   }


 TString selection;
   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
     {
       selection+=Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first +"-";
     }
   cout << selection << endl;


  if(Corr_selected_labels_.size()>0)
{

     std::vector < TGraphErrors* > sel_all_BC_correlations_;
   std::vector < std::vector < TGraphErrors* > > sel_one_para_all_selec_in_one_all_;
   std::vector < std::vector < TGraphErrors* > > sel_all_paras_in_one_all_;

 for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
   {
        sel_all_BC_correlations_.push_back(Correlation_B_C_.at(Corr_selected_labels_[Sel_Corr_i]));


   std::vector < TGraphErrors* > sel_all_paras_in_one_;
    for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
      {
	sel_all_paras_in_one_.push_back(fit_params_in_X_all_[par_i].at(Corr_selected_labels_[Sel_Corr_i]));
      }
    sel_all_paras_in_one_all_.push_back(sel_all_paras_in_one_);
   }

 for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
   {
   std::vector < TGraphErrors* > sel_one_para_all_selec_in_one_;
   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
      {
	sel_one_para_all_selec_in_one_.push_back(fit_params_in_X_all_[par_i].at(Corr_selected_labels_[Sel_Corr_i]));
      }
   sel_one_para_all_selec_in_one_all_.push_back(sel_one_para_all_selec_in_one_);
   }


  for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
    {
     draw_graphs(sel_all_paras_in_one_all_[Sel_Corr_i], -5., 5., leg_selected_[Sel_Corr_i],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] +"PARAM_sel_all_params_in_one_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first + "_");

    }
  for(Int_t par_i=0;par_i<dof_in_fit_;par_i++)
    {
      draw_graphs(sel_one_para_all_selec_in_one_all_[par_i], param_fit_y_edges_[par_i].first, param_fit_y_edges_[par_i].second, leg_params_[par_i],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] +"PARAM_sel_one_para_all_selec_in_one_" + param_fit_labels_[par_i] + "_");
    }

   draw_graphs(sel_all_BC_correlations_, -1., 1., leg_params_[0],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] +"_PARAM_selected_BC_Correlation_");



}





  chdir("..");


   outf->Close();



}


void flex_param::Loop(Bool_t setvalues, Int_t par_bin_choice, Int_t par_X_choice, Int_t par_eta_choice, TString par_Corr_choice, TString par_img_choice)
{
//   In a ROOT session, you can do:
//      Root > .L flex_param.C
//      Root > flex_param t
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

   if(setvalues){
 bin_choice	= par_bin_choice;	
 X_choice	= par_X_choice;
 eta_choice	= par_eta_choice;	
 Corr_choice    = par_Corr_choice;
 img_choice	= par_img_choice;
   }



   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   Declare_Labels();

   cout << "Bitte Binningvariable wählen (ACHTUNG: I.A. sollte die Binningvariable hier die gleiche sein wie beim Rausschreiben der Korrekturplots (bisher Genpt)!!!):" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>bin_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << bin_choice << ", also " << X_labels_[bin_choice]<< ", ausgewaehlt. Danke..." << endl;


   cout << "Bitte Parametrisierung wählen:" <<endl;
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   if(!setvalues){
   cin >>X_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << X_choice << ", also " << X_labels_[X_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Pseudorapiditäts-Region auswählen:" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_eta_region_;i++)  cout << i << ": " << eta_region_labels_[i] << endl;
   cin >>eta_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << eta_choice << ", also " << eta_region_labels_[eta_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Korrekturen zum Vergleich auswaehlen:" <<endl;
   //   cout << "Test: " << (char)66 << endl;

   //   cout << Corr_labels_[2].first << endl;

   if(!setvalues){
   for(Int_t i=0;i<no_all_Corr_labels_;i++)  
     {
       cout << (char)(i+48) << ": " << Corr_labels_[i].first << endl;
     }
   cin >> Corr_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << Corr_choice << " ausgewaehlt. Danke... Gewaehlt wurde also: " << endl;


   for(Int_t i=0;i<no_all_Corr_labels_;i++)  
     {
       if(Corr_choice.Contains((char)(i+48)))
	 {
	   cout << (char)(i+48) << ": " << Corr_labels_[i].first << endl;
	   Corr_selected_labels_.push_back(i);
	 }
     }
   if(Corr_selected_labels_.size()==0)cout << "Es wurde nichts ausgewaehlt!!!" << endl;


   cout << "Bitte Exportformat für Bilder wählen ('A' für keinen Export)" <<endl;
   if(!setvalues){
   cin >>img_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << img_choice << " ausgewaehlt. Danke..." << endl;


   Import_Histos();

   Fit_Histos();

    Create_TGraphErrors();

   std::cout << "Funzt...auch hier :)" << std::endl;

   gStyle->SetOptStat(1);
   gStyle->SetOptFit(1);
   gStyle->SetPalette(1);

   test->Range(-1.5,-0.625,13.5,5.625);
   test->SetBorderSize(2);
   test->SetFrameFillColor(0);
   test->SetSelected(test);







   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!
         Write_Histos();
        Write_TGraphErrors();


   //   Book_Histos();




   





}
