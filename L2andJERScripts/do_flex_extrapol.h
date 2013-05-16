#ifndef do_flex_extrapol_h
#define do_flex_extrapol_h

class do_flex_extrapol {

public :
  do_flex_extrapol(): cov(3){};
  void     import_plots();
  //  TH1D*    import_rel_response_plots(TString eta_abseta, );
  void     define_cosmetics_and_cuts();
  void     define_eta_bins_and_labels();
  void     Loop(TString jet_type_, TString generatorone_, TString generatortwo_, TString image_ext_, TString root_export_, TString use_imported_kFSRAbs_, TString fine_coarse_, TString use_easy_mean_, TString use_fitted_kFSR_, TString corr_generation_="Spring10", TString ratio_of_mean_or_GM_="Means", Bool_t export_all_plots_=true, TString kFSR_eq_one_ ="kFSR_no_eq_one", TString MPF_or_rel_response_ ="rel_response");


   TString jet_type, generatorone, generatortwo, image_ext, root_export, use_imported_kFSRAbs, fine_coarse, use_easy_mean, use_fitted_kFSR, corr_generation, ratio_of_mean_or_GM;
   Bool_t export_all_plots;
   TString kFSR_eq_one, MPF_or_rel_response;

  TGraphErrors*  import_kFSR_vs_Abseta_val1;
  TGraphErrors*  import_kFSR_vs_Abseta_val2;
  TGraphErrors*  import_kFSR_vs_Abseta_res1;
  TGraphErrors*  import_kFSR_vs_Abseta_res2;
  TH1D*          import_kFSR_vs_Abseta_histo_res1;
  TF1*           import_kFSR_fit;
  TFitResult*  import_kFSR_fit_result;
  TMatrixDSym cov;

  std::vector<Int_t> markers_;
  std::vector<Int_t> colours_;
  std::vector<Int_t> line_styles_;
  std::vector <TString> ptthreecuts;
  std::vector <Double_t> ptthreecuts_Double_;



//  Double_t kostas_eta_binning[]={-5.191,-3.489,-3.139,-2.964,-2.853,-2.5,-2.411,-2.322,-1.93,-1.479,-1.305,-1.131,-0.957,-0.783,-0.522,-0.261,0.0,0.261,0.522,0.783,0.957,1.131,1.305,1.479,1.93,2.322,2.411,2.5,2.853,2.964,3.139,3.489,5.191};
//  //
//  Int_t  kostas_no_eta =32;
//  Int_t  kostas_zero_eta =16;
//
//  Double_t k_HFfix_eta_binning[]={-5.191,-2.964,-2.853,-2.5,-2.411,-2.322,-1.93,-1.479,-1.305,-1.131,-0.957,-0.783,-0.522,-0.261,0.0,0.261,0.522,0.783,0.957,1.131,1.305,1.479,1.93,2.322,2.411,2.5,2.853,2.964,5.191};
//  //
//  Int_t  k_HFfix_no_eta =28;
//  Int_t  k_HFfix_zero_eta =14;
//
//  Double_t eta_binning[]={-6.0,-4.0,-3.5,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.1,1.2,1.3,1.4,1.5,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.5,4.0,6.0};
//  //39
//  Int_t  no_eta =38;
//  Int_t  zero_eta =19;
//
//  Double_t trad_eta_binning[]={-6.0,-4.0,-3.5,-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.3,-1.1,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.1,1.3,1.5,1.8,2.1,2.4,2.7,3.0,3.5,4.0,6.0};
//  //29
//  Int_t  trad_no_eta =28;
//  Int_t  trad_zero_eta =14;
// 
  std::vector <Double_t> kostas_eta_binning;
  Int_t  kostas_no_eta, kostas_zero_eta;

  std::vector <Double_t> k_HFfix_eta_binning;
  Int_t  k_HFfix_no_eta, k_HFfix_zero_eta;

  std::vector <Double_t> JEC_Mikko_eta_binning;
  Int_t  JEC_Mikko_no_eta, JEC_Mikko_zero_eta;
  
  std::vector <Double_t> Fine_eta_binning;
  Int_t  Fine_no_eta, Fine_zero_eta;  

  std::vector <Double_t> eta_binning;
  Int_t  no_eta, zero_eta;

  std::vector <Double_t> trad_eta_binning;
  Int_t  trad_no_eta, trad_zero_eta;


  std::vector< std::pair <TString,TString> > eta_bins_labels_;
  std::vector< std::pair <Double_t,Double_t> > eta_bins_;
  std::vector< std::pair <TString,TString> > Abseta_bins_labels_;
  std::vector< std::pair <Double_t,Double_t> > Abseta_bins_;
  std::vector< std::pair <TString,TString> >  trad_eta_bins_labels_;
  std::vector< std::pair <Double_t,Double_t> >  trad_eta_bins_;
  std::vector< std::pair <TString,TString> >  trad_Abseta_bins_labels_;
  std::vector< std::pair <Double_t,Double_t> >  trad_Abseta_bins_;
  Int_t eta_bins;
  Int_t no_eta_bins;
  Int_t no_Abseta_bins;
  Int_t trad_eta_bins;
  //  Int_t trad_no_eta_bins= trad_eta_bins_.size();
  Int_t trad_no_Abseta_bins;



};



void do_flex_extrapol::import_plots(){

    TFile *inf;
    if(fine_coarse.Contains("k_HFfix")){
       inf = new TFile("../"+corr_generation+"/k_HFfix_use_easy_mean_TuneZ2_TuneZ2_PF_kFSR_histos.root","OPEN");
    }
    else if(fine_coarse.Contains("kostas")){

      //use corresponding PF kFSR
            inf = new TFile("../"+corr_generation+"/kostas_use_easy_mean_TuneZ2_TuneZ2_PF_kFSR_histos.root","OPEN");}
    else if(fine_coarse.Contains("Fine")){

      //use corresponding PF kFSR
            inf = new TFile("../"+corr_generation+"/Fine_use_easy_mean_TuneZ2_TuneZ2_PF_kFSR_histos.root","OPEN");}	    
    else if(fine_coarse.Contains("JEC_Mikko")){

      //use corresponding PF kFSR
            inf = new TFile("../"+corr_generation+"/JEC_Mikko_use_easy_mean_TuneZ2_TuneZ2_PF_kFSR_histos.root","OPEN");}
	        

	    
      //use individual kRad for each algorithm
      //      inf = new TFile("../"+corr_generation+"/kostas_use_easy_mean_TuneZ2_TuneZ2_"+ jet_type +"_kFSR_histos.root","OPEN");



      // use Full2011 MPF krad for each algorithm
      //      inf = new TFile("../2011Full2011_CORRF11DB_He_AK5_MC_F11Z2wPUsm_Y_f_kostas_MPF_AK5/kostas_use_easy_mean_TuneZ2_TuneZ2_"+ jet_type +"_kFSR_histos.root","OPEN");

      //      inf = new TFile("../2011Full2011_CORRF11DB_He_AK5_MC_F11Z2wPUsm_Y_f_kostas_AK5/kostas_use_easy_mean_TuneZ2_TuneZ2_"+ jet_type +"_kFSR_histos.root","OPEN");



      //      inf = new TFile("../"+corr_generation+"/kostas_use_easy_mean_TuneZ2_TuneZ2_"+ jet_type +"_kFSR_histos.root","OPEN");


      //      inf = new TFile("kostas_use_easy_mean_TuneZ2_TuneZ2_"+ jet_type+"_kFSR_histos.root,"OPEN");

    
    
    else inf = new TFile("coarse_use_easy_mean_TuneZ2_TuneZ2_"+ jet_type+"_kFSR_histos.root","OPEN");
    if (inf->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;

    import_kFSR_vs_Abseta_val1  = (TGraphErrors*)inf->Get("kFSR_vs_Abseta_val1");
    import_kFSR_vs_Abseta_val2  = (TGraphErrors*)inf->Get("kFSR_vs_Abseta_val2");
    import_kFSR_vs_Abseta_res1  = (TGraphErrors*)inf->Get("kFSR_vs_Abseta_res1");
    import_kFSR_vs_Abseta_res2  = (TGraphErrors*)inf->Get("kFSR_vs_Abseta_res2");
    import_kFSR_vs_Abseta_histo_res1  = (TH1D*)inf->Get("kFSR_vs_Abseta_histo_res1");
    import_kFSR_fit = (TF1* ) import_kFSR_vs_Abseta_histo_res1->GetFunction("kFSR_fit");
    import_kFSR_fit_result = (TFitResult*) inf->Get("TFitResult-kFSR_vs_Abseta_histo_res1-kFSR_fit");

    //    import_kFSR_fit_result->Print("V"); 
    //funzt bnoch nicht...
     TMatrixDSym temp = import_kFSR_fit_result->GetCovarianceMatrix(); 
     temp.Print();
     cov = temp;
    cov.Print();

}

void     do_flex_extrapol::define_cosmetics_and_cuts(){
  markers_.push_back(20);
  markers_.push_back(24);
  markers_.push_back(21);
  markers_.push_back(25);
  markers_.push_back(22);
  markers_.push_back(26);
  markers_.push_back(23);
  markers_.push_back(29);
  markers_.push_back(28);

  colours_.push_back(1);
  colours_.push_back(2);
  colours_.push_back(4);
  colours_.push_back(3);
  colours_.push_back(46);
  colours_.push_back(9);
  colours_.push_back(12);

  line_styles_.push_back(1);
  line_styles_.push_back(1);
  line_styles_.push_back(1);
  line_styles_.push_back(1);
  line_styles_.push_back(1);
  line_styles_.push_back(1);
  line_styles_.push_back(1);

if(fine_coarse.Contains("JEC_Mikko")){
  //  ptthreecuts.push_back("05");
    ptthreecuts.push_back("10");
    ptthreecuts.push_back("15");
    ptthreecuts.push_back("20");
  //  ptthreecuts.push_back("25");
    ptthreecuts.push_back("30");
  //  ptthreecuts.push_back("35");
  //ptthreecuts.push_back("40");

  //  ptthreecuts_Double_.push_back(0.05);
    ptthreecuts_Double_.push_back(0.10);
    ptthreecuts_Double_.push_back(0.15);
    ptthreecuts_Double_.push_back(0.20);
  //  ptthreecuts_Double_.push_back(0.25);
    ptthreecuts_Double_.push_back(0.30);
  //  ptthreecuts_Double_.push_back(0.35);
  //ptthreecuts_Double_.push_back(0.40);
}
else{
  //  ptthreecuts.push_back("05");
    ptthreecuts.push_back("10");
  //  ptthreecuts.push_back("15");
    ptthreecuts.push_back("20");
  //  ptthreecuts.push_back("25");
    ptthreecuts.push_back("30");
  //  ptthreecuts.push_back("35");
  ptthreecuts.push_back("40");

  //  ptthreecuts_Double_.push_back(0.05);
    ptthreecuts_Double_.push_back(0.10);
   // ptthreecuts_Double_.push_back(0.15);
    ptthreecuts_Double_.push_back(0.20);
  //  ptthreecuts_Double_.push_back(0.25);
    ptthreecuts_Double_.push_back(0.30);
  //  ptthreecuts_Double_.push_back(0.35);
  ptthreecuts_Double_.push_back(0.40);
}

}


void     do_flex_extrapol::define_eta_bins_and_labels(){

//    int data[] = {111, 112, 123, 134};
//    int size = sizeof( data ) / sizeof( data[0] );
//    vector<int> vec( data, &data[ size ] );


  Double_t kostas_eta_binning_[]={-5.191,-3.489,-3.139,-2.964,-2.853,-2.5,-2.411,-2.322,-1.93,-1.479,-1.305,-1.131,-0.957,-0.783,-0.522,-0.261,0.0,0.261,0.522,0.783,0.957,1.131,1.305,1.479,1.93,2.322,2.411,2.5,2.853,2.964,3.139,3.489,5.191};
  //
  int size_k = sizeof(kostas_eta_binning_)/sizeof(kostas_eta_binning_[0]);
  vector<Double_t> vec_k(kostas_eta_binning_, &kostas_eta_binning_[ size_k] );
  kostas_eta_binning = vec_k;
  kostas_no_eta =32;
  kostas_zero_eta =16;



  Double_t k_HFfix_eta_binning_[]={-5.191,-2.964,-2.853,-2.5,-2.411,-2.322,-1.93,-1.479,-1.305,-1.131,-0.957,-0.783,-0.522,-0.261,0.0,0.261,0.522,0.783,0.957,1.131,1.305,1.479,1.93,2.322,2.411,2.5,2.853,2.964,5.191};
  //
  int size_kHF = sizeof(k_HFfix_eta_binning_)/sizeof(k_HFfix_eta_binning_[0]);
  vector<Double_t> vec_kHF(k_HFfix_eta_binning_, &k_HFfix_eta_binning_[ size_kHF] );
  k_HFfix_eta_binning = vec_kHF;
  k_HFfix_no_eta =28;
  k_HFfix_zero_eta =14;

  Double_t JEC_Mikko_eta_binning_[]={-5.191,-3.2,-2.964,-2.5,-1.93,-1.305,0.0,1.305,1.93,2.5,2.964,3.2,5.191};
  //
  int size_JM = sizeof(JEC_Mikko_eta_binning_)/sizeof(JEC_Mikko_eta_binning_[0]);
  vector<Double_t> vec_JM(JEC_Mikko_eta_binning_, &JEC_Mikko_eta_binning_[ size_JM] );
  JEC_Mikko_eta_binning = vec_JM;
  JEC_Mikko_no_eta =12;
  JEC_Mikko_zero_eta =6;

  Double_t Fine_eta_binning_[]={-5.191,-3.489,-3.139,-2.964,-2.853,-2.5,-2.322,-2.172,-2.043,-1.930,-1.830,-1.740,-1.653,-1.566,-1.479,-1.392,-1.305,-1.218,-1.131,-1.044,-0.957,-0.879,-0.783,-0.696,-0.609,-0.522,-0.435,-0.348,-0.261,-0.174,-0.087,0.000,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.879,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.853,2.964,3.139,3.489,5.191};
  //
  int size_F = sizeof(Fine_eta_binning_)/sizeof(Fine_eta_binning_[0]);
  vector<Double_t> vec_F(Fine_eta_binning_, &Fine_eta_binning_[ size_F] );
  Fine_eta_binning = vec_F;
  Fine_no_eta =62;
  Fine_zero_eta =31;

  Double_t eta_binning_[]={-6.0,-4.0,-3.5,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.5,-1.4,-1.3,-1.2,-1.1,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.1,1.2,1.3,1.4,1.5,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.5,4.0,6.0};
  //39
  int size = sizeof(eta_binning_)/sizeof(eta_binning_[0]);
  vector<Double_t> vec(eta_binning_, &eta_binning_[ size] );
  eta_binning = vec;
  no_eta =38;
  zero_eta =19;

  Double_t trad_eta_binning_[]={-6.0,-4.0,-3.5,-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.3,-1.1,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.1,1.3,1.5,1.8,2.1,2.4,2.7,3.0,3.5,4.0,6.0};
  int size_tra = sizeof(trad_eta_binning_)/sizeof(trad_eta_binning_[0]);
  vector<Double_t> vec_tra(trad_eta_binning_, &trad_eta_binning_[ size_tra] );
  trad_eta_binning = vec_tra;
  //29
  trad_no_eta =28;
  trad_zero_eta =14;




  if(fine_coarse.Contains("coarse")){
    cout <<"coarse detected... adapt eta-binning" << endl;
    no_eta=trad_no_eta;
    zero_eta=trad_zero_eta;
  for(Int_t eta_i=0;eta_i<=no_eta;eta_i++){
    eta_binning[eta_i]=trad_eta_binning[eta_i];
    cout << eta_binning[eta_i] << " ";
  }
  cout << endl;
  cout << no_eta << " zero_eta " << zero_eta << endl;
  cout << trad_no_eta << " zero_eta " << trad_zero_eta << endl;
  }
  else if(fine_coarse.Contains("kostas")){
    no_eta=kostas_no_eta;
    zero_eta=kostas_zero_eta;
  for(Int_t eta_i=0;eta_i<=no_eta;eta_i++){
    eta_binning[eta_i]=kostas_eta_binning[eta_i];
    cout << eta_binning[eta_i] << " ";
  }
  }
  else if(fine_coarse.Contains("k_HFfix")){
    no_eta=k_HFfix_no_eta;
    zero_eta=k_HFfix_zero_eta;
  for(Int_t eta_i=0;eta_i<=no_eta;eta_i++){
    eta_binning[eta_i]=k_HFfix_eta_binning[eta_i];
    cout << eta_binning[eta_i] << " ";
  }
  }
  else if(fine_coarse.Contains("JEC_Mikko")){
    no_eta=JEC_Mikko_no_eta;
    zero_eta=JEC_Mikko_zero_eta;
  for(Int_t eta_i=0;eta_i<=no_eta;eta_i++){
    eta_binning[eta_i]=JEC_Mikko_eta_binning[eta_i];
    cout << eta_binning[eta_i] << " ";
  }
  }
  else if(fine_coarse.Contains("Fine")){
    no_eta=Fine_no_eta;
    zero_eta=Fine_zero_eta;
  for(Int_t eta_i=0;eta_i<=no_eta;eta_i++){
    eta_binning[eta_i]=Fine_eta_binning[eta_i];
    cout << eta_binning[eta_i] << " ";
  }
  }  
  for(Int_t eta_i=0;eta_i<no_eta;eta_i++){
    char buffer [50];
    char buffer2 [50];
    sprintf (buffer, "%.1f", eta_binning[eta_i]);
    
        cout << buffer;
    
    sprintf (buffer2, "%.1f", eta_binning[eta_i+1]);
    
        cout << ", " << buffer2 << endl;
    
    eta_bins_labels_.push_back(make_pair(buffer,buffer2));
    eta_bins_.push_back(make_pair(eta_binning[eta_i],eta_binning[eta_i+1]));
    if(eta_i>=zero_eta){
    Abseta_bins_labels_.push_back(make_pair(buffer,buffer2));
    Abseta_bins_.push_back(make_pair(eta_binning[eta_i],eta_binning[eta_i+1]));
    //    cout << "ABS:" << buffer<< ", " << buffer2 << endl;
    }

  }
  for(Int_t eta_i=0;eta_i< trad_no_eta;eta_i++){
    char buffer [50];
    char buffer2 [50];
    sprintf (buffer, "%.1f",  trad_eta_binning[eta_i]);
    
    //        cout << buffer;
    
    sprintf (buffer2, "%.1f",  trad_eta_binning[eta_i+1]);
    
    //        cout << ", " << buffer2 << endl;
    
    trad_eta_bins_labels_.push_back(make_pair(buffer,buffer2));
    trad_eta_bins_.push_back(make_pair(trad_eta_binning[eta_i],trad_eta_binning[eta_i+1]));
    if(eta_i>=trad_zero_eta){
    trad_Abseta_bins_labels_.push_back(make_pair(buffer,buffer2));
    trad_Abseta_bins_.push_back(make_pair(trad_eta_binning[eta_i],trad_eta_binning[eta_i+1]));
    //        cout << "ABS:" << buffer<< ", " << buffer2 << endl;
    }

  }


  eta_bins=eta_bins_.size();
  no_eta_bins= eta_bins_.size();
  no_Abseta_bins= Abseta_bins_.size();
  trad_eta_bins=trad_eta_bins_.size();
  //  Int_t trad_no_eta_bins= trad_eta_bins_.size();
  trad_no_Abseta_bins= trad_Abseta_bins_.size();





}



#endif 


