#include <vector>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TAxis.h>

void setStyle()
{
  gStyle->SetPalette(1);

  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(800); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //Position on screen
  gStyle->SetCanvasDefY(0);

  // For the frame
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetFrameLineStyle(0);
  gStyle->SetFrameLineWidth(1);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the histo:
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);

  // For the statistics box:
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("neMR");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.92);              
  gStyle->SetStatY(0.86);              
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.22);

  // For the legend
  gStyle->SetLegendBorderSize(1);

  //  Margins
  // -------------------------------------------
  gStyle->SetPadTopMargin(0.16);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.19);
  gStyle->SetPadRightMargin(0.16);

  // For the Global title:
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42,"");
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleFontSize(0.12);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.515);
  gStyle->SetTitleH(0.06);
  gStyle->SetTitleXOffset(0);
  gStyle->SetTitleYOffset(0);
  gStyle->SetTitleBorderSize(0);

  // For the axis labels:
  //  For the axis labels and titles
  // -------------------------------------------
  gStyle->SetTitleColor(1,"XYZ");
  gStyle->SetLabelColor(1,"XYZ");
  // For the axis labels:
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007,"XYZ");
  gStyle->SetLabelSize(0.045,"XYZ");
  
  // For the axis titles:
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.5);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}



void parametrization(Bool_t setvalues=0, Int_t par_bin_choice=0, Int_t par_X_choice=0, Int_t par_eta_choice=0, Int_t par_GMP_choice=1){

  setStyle();
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.02);

  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0);


  std::vector < std::pair < Double_t,Double_t > > eta_region_;
  Int_t no_eta_region_;
  Int_t eta_choice;
  std::vector <  TString > eta_region_labels_;


  std::vector <  TString > GMP_choice_;
  Int_t GMP_choice;
  Int_t no_GMP_choice_;
  GMP_choice_.push_back("GMP_");//GM<P_JW_raw
  GMP_choice_.push_back("");//JW_raw
  no_GMP_choice_=GMP_choice_.size();
   TString img = ".eps";

   //Kalibri-parameters

   //b
//
//Double_t b_paras_kalibri[4][5]=
//  {
//    {-4.22844,4.15415,2.282,1.21118,-7.43724},
//    {-4.44824,4.83274,-9.69512,1.09251,-7.6847},
//    {0.167434,-4.6283,25.6509,1.40079,-16.5885},
//    {1.87382,-10.8619,58.8731,3056.24,-69.0435}
//  };
//

Double_t b_paras_kalibri[4][5]=
{
  {3.56585,-2.94696,-4.14987,-1.29651,-8.37254},
  {4.23050,-4.64314,10.8783,-1.20905,-8.62383},
  {0.278775,2.84047,2.78045,-3.75321,-23.5930},
  {-1.26146,8.91478,-41.4318,-98.7111,-46.1337}
};




   //c

//Double_t c_paras_kalibri[4][5]=
//  {
//    {0.294896,2.70661,-1.37199,-9.92943,-7.78206},
//    {2.61907,-3.13131,1.4052,-7.06145,-5.50906},
//    {-0.998786,5.90247,-2.71171,-74597.9,-106.392},
//    {9.12808,8.18951,2.62724,-11.3193,0.0046}
//  };
//

Double_t c_paras_kalibri[4][5]=
  {
    {2.37207,-7.57939,2.13660,141.489,-28.0435},
    {2.84147,-8.14719,2.33879,52815.8,-103.224},
    {-7.35154,-4.59322,5.01315,8.19059,-0.048455},
    {3.09711,-12.6187,17.9219,0.393632,-1.69317}
  };



   //




  //GMP_JW_raw

  std::vector <  TString > X_labels_;
  std::vector <  TString > X_labels_img_;
  Int_t bin_choice;
  Int_t X_choice;
  Int_t no_X_labels_;
  X_labels_.push_back("p_{T}^{gen}");
  X_labels_.push_back("E_{T}^{gen}");
  X_labels_.push_back("E^{Gen}");
  X_labels_.push_back("p_{T}^{calo}");
  X_labels_.push_back("E_{T}^{calo}");
  X_labels_.push_back("E^{Calo}");
  X_labels_.push_back("p_{T}^{L2L3}");
  no_X_labels_=X_labels_.size();

  X_labels_img_.push_back("PTgen");
  X_labels_img_.push_back("ETgen");
  X_labels_img_.push_back("EGen");
  X_labels_img_.push_back("PTcalo");
  X_labels_img_.push_back("ETcalo");
  X_labels_img_.push_back("ECalo");
  X_labels_img_.push_back("PTL2L3");


  eta_region_labels_.push_back("barrel_0_0_1_305");
  eta_region_labels_.push_back("tr_endcap_1_305_2_65");
  eta_region_labels_.push_back("HE_2_65_2_964_");
  eta_region_labels_.push_back("Cherenkov_2_964_5_191");

  eta_region_.push_back(make_pair(0.0,1.305)); 
  eta_region_.push_back(make_pair(1.305,2.65)); 
  eta_region_.push_back(make_pair(2.65,2.964)); 
  eta_region_.push_back(make_pair(2.964,5.191)); 
  no_eta_region_=eta_region_.size();

   if(setvalues){
 bin_choice	= par_bin_choice;	
 X_choice	= par_X_choice;
 eta_choice	= par_eta_choice;	
 // img_choice	= par_img_choice;
 GMP_choice     = par_GMP_choice;
   }




   cout << "Bitte Binningvariable wählen (ACHTUNG: I.A. sollte die Binningvariable hier die gleiche sein wie beim Rausschreiben der Korrekturplots (bisher Genpt)!!!):" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>bin_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << bin_choice << ", also " << X_labels_[bin_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Parametrisierung wählen:" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
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

   cout << "Bitte wählen, ob Mean oder Gaussian Mean genommen werden soll:" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_GMP_choice_;i++)  cout << i << ": " << GMP_choice_[i] << endl;
   cin >>GMP_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << GMP_choice << ", also " << GMP_choice_[GMP_choice]<< ", ausgewaehlt. Danke..." << endl;
   //   img = "_" + GMP_choice_[GMP_choice] + img;


TString root_resol_name_binning= "histos_parametrisations_" + eta_region_labels_[eta_choice] + X_labels_[bin_choice] + X_labels_[X_choice]+ ".root";

TFile *_file = TFile::Open(root_resol_name_binning);

 std::vector < Double_t > test;

  std::vector < TCanvas* > printcanvases_;

   TCanvas *X0 = new TCanvas("X0", "X0",0,0,500,500);
   X0->SetLogx();
   printcanvases_.push_back(X0);
   TCanvas *B = new TCanvas("B", "B",0,0,500,500);
   B->SetLogx();
   printcanvases_.push_back(B);
   TCanvas *C = new TCanvas("C", "C",0,0,500,500);
   C->SetLogx();
   printcanvases_.push_back(C);
   TCanvas *Sigma = new TCanvas("Sigma", "Sigma",0,0,500,500);
   Sigma->SetLogx();
   printcanvases_.push_back(Sigma);

   // TCanvas *Sigma = new TCanvas("Sigma", "Sigma",341,360,700,530);

   TGraphErrors*   X_vs_X0;
   TGraphErrors*   X_vs_X0_red;
   TGraphErrors*   X_vs_X0_all;
   TGraphErrors*   X_vs_B;
   TGraphErrors*   X_vs_B_red;
   TGraphErrors*   X_vs_B_all;
   TGraphErrors*   X_vs_C;
   TGraphErrors*   X_vs_C_red;
   TGraphErrors*   X_vs_C_all;
   TGraphErrors*   X_vs_Sigma; 
   TGraphErrors*   X_vs_Sigma_red;
   TGraphErrors*   X_vs_Sigma_all;

   X_vs_X0 = (TGraphErrors*) _file->Get(  X_labels_[X_choice] + "vs_X0_of_" + GMP_choice_[GMP_choice] + "JW_raw");
   X_vs_B = (TGraphErrors*) _file->Get(  X_labels_[X_choice] + "vs_B_X-X0__of_" + GMP_choice_[GMP_choice] + "JW_raw");
   X_vs_C = (TGraphErrors*) _file->Get(  X_labels_[X_choice] + "vs_C_X-X0_2_of_" + GMP_choice_[GMP_choice] + "JW_raw");
   X_vs_Sigma = (TGraphErrors*) _file->Get(  X_labels_[X_choice] + "vs_Sigma_of_" + GMP_choice_[GMP_choice] + "JW_raw");
   //GMP_JW_raw

   X_vs_X0->SetMarkerStyle(21);
   X_vs_B->SetMarkerStyle(21);
   X_vs_C->SetMarkerStyle(21);
   X_vs_Sigma->SetMarkerStyle(21);

   X_vs_X0->SetMarkerSize(2);
   X_vs_B->SetMarkerSize(2);
   X_vs_C->SetMarkerSize(2);
   X_vs_Sigma->SetMarkerSize(2);

   cout << X_vs_X0->GetName() << endl;
   cout << "import successful" << endl;

   Double_t fit_low_loose = 40.;
   Double_t fit_low_tight = 10.;
   Double_t fit_hig = 990.;//improve this later(after import)...



   fit_low_loose = X_vs_X0->GetX()[0] + 40.;
   fit_low_tight = X_vs_X0->GetX()[0] - 1.;
   fit_hig = X_vs_X0->GetX()[X_vs_X0->GetN()-2] + 5.;
   cout << "upper fit bopundary set to " << fit_hig << endl;


   X0->cd();
    TF1* X0_param = new TF1("X0_param","([0]/10.)+([1]/100.)*log(x)+([2]/10000.)*x",fit_low_tight,fit_hig);
    X0_param->SetParameters(0.2,-0.02,-0.0002);

    //   X_vs_X0->SetOptFit(0000);


   X_vs_X0->Draw("ALP");
   X_vs_X0->GetYaxis()->SetTitle("#mu");
   X_vs_X0->GetXaxis()->SetTitle(X_labels_[X_choice]+"[GeV]");
   X_vs_X0->Fit(X0_param,"","same",fit_low_loose,fit_hig);
   X_vs_X0->Print();
   X0_param->Print();

   TCanvas *X0_red = new TCanvas("X0_red", "X0_red",600,0,500,500);
   X0_red->SetLogx();

    TF2* X0_minus = new TF2("Fit_X0_minus","y - (([0]/10.)+([1]/100.)*log(x)+([2]/10000.)*x)",fit_low_loose,fit_hig);
   X0_minus->SetParameters(X0_param->GetParameters());
   X_vs_X0_red = (TGraphErrors*) X_vs_X0->Clone();
    X_vs_X0_red->Apply(X0_minus);
    X_vs_X0_red->Draw("ALP");

    TF1* X0_param_low = new TF1("Fit_X0_low"," (([0]/10.))*exp( ([1]/10.)*x)",fit_low_tight,fit_hig);
    X0_param_low->SetParameters(0.2,-0.22);
    X_vs_X0_red->Fit(X0_param_low,"M","same",fit_low_tight,fit_hig);

    X_vs_X0_all = (TGraphErrors*) X_vs_X0->Clone(); //(TGraphErrors*) _file->Get( X_labels_[X_choice] + "vs_X0_of_JW_raw");//(TGraphErrors*) X_vs_X0->Clone();

    TF1* X0_param_all = new TF1("Fit_X0_all"," (([0]/10.)+([1]/100.)*log(x)+([2]/10000)*x) + (([3]/10.))*exp( ([4]/10.)*x)",fit_low_tight,fit_hig);
    X0_param_all->SetLineColor(kRed);
    for(Int_t i=0;i<3;i++)
      {
	X0_param_all->SetParameter(i,X0_param->GetParameter(i));
	if(i<2)X0_param_all->SetParameter(i+3,X0_param_low->GetParameter(i));
      }
   X0_param->Print();

   //   X0_param_all->FixParameter(4,0);

    TCanvas *Xall = new TCanvas("Xall", "Xall",0,400,500,500);
   Xall->SetLogx();

   TPaveStats *st = (TPaveStats*)X_vs_X0_all->GetListOfFunctions()->FindObject("stats");
   st->Delete();

   X_vs_X0_all->Draw("ALP");
   X_vs_X0_all->Fit(X0_param_all,"EM","same",fit_low_tight,fit_hig);
   //   X_vs_X0_all->SetStats(0);  



   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS X0!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS X0!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS X0!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS X0!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;

   cout <<"{";
   for(Int_t i=0; i<5;i++)cout << X0_param_all->GetParameter(i) <<",";
   cout << "}" << endl;

   //////////////////////////BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
   B->cd();
    TF1* B_param = new TF1("B_param","[0]+([1]/10.)*log(x)+([2]/10000.)*x",fit_low_tight,fit_hig);
    B_param->SetParameters(0.2,-0.02,-0.0002);

   X_vs_B->Draw("ALP");
   X_vs_B->GetYaxis()->SetTitle("b");
   X_vs_B->GetXaxis()->SetTitle(X_labels_[X_choice]+"[GeV]");
   X_vs_B->Fit(B_param,"","same",fit_low_loose,fit_hig);

   TCanvas *B_red = new TCanvas("B_red", "B_red",600,0,500,500);
   B_red->SetLogx();

    TF2* B_minus = new TF2("Fit_B_minus","y - ([0]+([1]/10.)*log(x)+([2]/10000.)*x)",fit_low_loose,fit_hig);
   B_minus->SetParameters(B_param->GetParameters());
   X_vs_B_red = (TGraphErrors*) X_vs_B->Clone();
    X_vs_B_red->Apply(B_minus);
    X_vs_B_red->Draw("ALP");


    TF1* B_param_low = new TF1("Fit_B_low"," (([0]*10.))*exp( ([1]/100.)*x)",fit_low_tight,fit_hig);
    B_param_low->SetParameters(0.2,-0.22);
    X_vs_B_red->Fit(B_param_low,"M","same",fit_low_tight,fit_hig);

    X_vs_B_all = (TGraphErrors*) X_vs_B->Clone(); //(TGraphErrors*) _file->Get( X_labels_[X_choice] + "vs_B_X-X0__of_JW_raw");//(TGraphErrors*) X_vs_B->Clone();

    TF1* B_param_all = new TF1("Fit_B_all"," ([0]+([1]/10.)*log(x)+([2]/10000.)*x) + (([3]*10.))*exp( ([4]/100.)*x)",fit_low_tight,fit_hig);
    B_param_all->SetLineColor(kRed);
    for(Int_t i=0;i<3;i++)
      {
	B_param_all->SetParameter(i,B_param->GetParameter(i));
	if(i<2)B_param_all->SetParameter(i+3,B_param_low->GetParameter(i));
      }
   B_param->Print();
   //   B_param_all->FixParameter(4,0);


    TCanvas *Ball = new TCanvas("Ball", "Ball",0,400,500,500);
   Ball->SetLogx();

   X_vs_B_all->Draw("ALP");
   X_vs_B_all->Fit(B_param_all,"EM","same",fit_low_tight,fit_hig);

    TF1* B_param_kalibri = new TF1("B_param_kalibri"," ([0]+([1]/10.)*log(x)+([2]/10000.)*x) + (([3]*10.))*exp( ([4]/100.)*x)",fit_low_tight,fit_hig);
    B_param_kalibri->SetLineColor(kBlue);
    for(Int_t i=0;i<5;i++)
      {
	B_param_kalibri->SetParameter(i, b_paras_kalibri[eta_choice][i]);
      }
    TCanvas *Ball_kalibri = new TCanvas("Ball_kalibri", "Ball_kalibri",0,400,500,500);
   Ball_kalibri->SetLogx();

   X_vs_B_all->Draw("ALP");
   X_vs_B_all->Fit(B_param_all,"EM","same",fit_low_tight,fit_hig);
    if(bin_choice==0)B_param_kalibri->Draw("same");

   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS BBBBBBBB!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS BBBBBBBB!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS BBBBBBBB!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS BBBBBBBB!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;

   cout <<"";
   for(Int_t i=0; i<5;i++)cout << B_param_all->GetParameter(i) <<" ";
   cout << "" << endl;




   //////////////////////////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   C->cd();
    TF1* C_param = new TF1("C_param","([0]*10.)+[1]*log(x)+([2]/100.)*x",fit_low_tight,fit_hig);
    C_param->SetParameters(0.2,-0.02,-0.0002);

   X_vs_C->Draw("ALP");
   X_vs_C->GetYaxis()->SetTitle("c");
   X_vs_C->GetXaxis()->SetTitle(X_labels_[X_choice]+"[GeV]");
   X_vs_C->Fit(C_param,"","same",fit_low_loose,fit_hig);

   TCanvas *C_red = new TCanvas("C_red", "C_red",600,0,500,500);
   C_red->SetLogx();

    TF2* C_minus = new TF2("Fit_C_minus","y - (([0]*10.)+[1]*log(x)+([2]/100.)*x)",fit_low_loose,fit_hig);
   C_minus->SetParameters(C_param->GetParameters());
   X_vs_C_red = (TGraphErrors*) X_vs_C->Clone();
    X_vs_C_red->Apply(C_minus);
    X_vs_C_red->Draw("ALP");


    TF1* C_param_low = new TF1("Fit_C_low"," (([0]*10.))*exp( ([1]/100.)*x)",fit_low_tight,fit_hig);
    C_param_low->SetParameters(0.2,-0.22);
    X_vs_C_red->Fit(C_param_low,"M","same",fit_low_tight,fit_hig);

    X_vs_C_all = (TGraphErrors*) X_vs_C->Clone();// (TGraphErrors*) _file->Get( X_labels_[X_choice] + "vs_C_X-X0_2_of_JW_raw");//(TGraphErrors*) X_vs_C->Clone();

    TF1* C_param_all = new TF1("Fit_C_all"," (([0]*10.)+[1]*log(x)+([2]/100.)*x) + (([3]*10.))*exp( ([4]/100.)*x)",fit_low_tight,fit_hig);
    C_param_all->SetLineColor(kRed);
    for(Int_t i=0;i<3;i++)
      {
	C_param_all->SetParameter(i,C_param->GetParameter(i));
	if(i<2)C_param_all->SetParameter(i+3,C_param_low->GetParameter(i));
      }
   C_param->Print();
   //   C_param_all->FixParameter(4,0);


    TCanvas *Call = new TCanvas("Call", "Call",0,400,500,500);
   Call->SetLogx();

   X_vs_C_all->Draw("ALP");
   X_vs_C_all->Fit(C_param_all,"EM","same",fit_low_tight,fit_hig);

    TF1* C_param_kalibri = new TF1("C_param_kalibri"," ([0]+([1]/10.)*log(x)+([2]/10000.)*x) + (([3]*10.))*exp( ([4]/100.)*x)",fit_low_tight,fit_hig);
    C_param_kalibri->SetLineColor(kBlue);
    for(Int_t i=0;i<5;i++)
      {
	C_param_kalibri->SetParameter(i, c_paras_kalibri[eta_choice][i]);
      }
    //        if(bin_choice==0)C_param_kalibri->Draw("same");
    TCanvas *Call_kalibri = new TCanvas("Call_kalibri", "Call_kalibri",0,400,500,500);
   Call_kalibri->SetLogx();

   X_vs_C_all->Draw("ALP");
   X_vs_C_all->Fit(C_param_all,"EM","same",fit_low_tight,fit_hig);
    if(bin_choice==0)C_param_kalibri->Draw("same");

   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS CCCCC!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS CCCCC!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS CCCCC!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS CCCCC!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;

   cout <<"";
   for(Int_t i=0; i<5;i++)cout << C_param_all->GetParameter(i) <<" ";
   cout << "" << endl;



   //////////////////////////SigmaSigmaBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
   Sigma->cd();
    TF1* Sigma_param = new TF1("Sigma_param","([0]/100.)+([1]/1000.)*log(x)+([2]/1000000.)*x",fit_low_tight,fit_hig);
    Sigma_param->SetParameters(0.2,-0.02,-0.0002);

   X_vs_Sigma->Draw("ALP");
   X_vs_Sigma->GetYaxis()->SetTitle("#sigma");
   X_vs_Sigma->GetXaxis()->SetTitle(X_labels_[X_choice]+"[GeV]");
   X_vs_Sigma->Fit(Sigma_param,"","same",fit_low_loose,fit_hig);

   TCanvas *Sigma1 = new TCanvas("XSigma", "Sigma1",600,0,500,500);
   Sigma1->SetLogx();

    TF2* Sigma_minus = new TF2("Fit_Sigma_minus","y - (([0]/100.)+([1]/1000.)*log(x)+([2]/1000000.)*x)",fit_low_loose,fit_hig);
   Sigma_minus->SetParameters(Sigma_param->GetParameters());
   X_vs_Sigma_red = (TGraphErrors*) X_vs_Sigma->Clone();
    X_vs_Sigma_red->Apply(Sigma_minus);
    X_vs_Sigma_red->Draw("ALP");


    TF1* Sigma_param_low = new TF1("Fit_Sigma_low"," (([0]/100.))*exp( ([1]/10.)*x)",fit_low_tight,fit_hig);
    Sigma_param_low->SetParameters(0.2,-0.22);
    X_vs_Sigma_red->Fit(Sigma_param_low,"M","same",fit_low_tight,fit_hig);

    X_vs_Sigma_all = (TGraphErrors*) X_vs_Sigma->Clone();// (TGraphErrors*) _file->Get( X_labels_[X_choice] + "vs_Sigma_of_JW_raw");//(TGraphErrors*) X_vs_Sigma->Clone();

    cout <<"bsi hierhin ok..." << endl;

    TF1* Sigma_param_all = new TF1("Fit_Sigma_all"," (([0]/100.)+([1]/1000.)*log(x)+([2]/1000000.)*x) + (([3]/100.))*exp( ([4]/10.)*x)",fit_low_tight,fit_hig);
    Sigma_param_all->SetLineColor(kRed);
    for(Int_t i=0;i<3;i++)
      {
	Sigma_param_all->SetParameter(i,Sigma_param->GetParameter(i));
	if(i<2)Sigma_param_all->SetParameter(i+3,Sigma_param_low->GetParameter(i));
      }
   Sigma_param->Print();
   //   Sigma_param_all->FixParameter(4,0);


    TCanvas *Sigmaall = new TCanvas("Sigmaall", "Sigmaall",0,400,500,500);
   Sigmaall->SetLogx();

   X_vs_Sigma_all->Draw("ALP");
   X_vs_Sigma_all->Fit(Sigma_param_all,"EM","same",fit_low_tight,fit_hig);


   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS Sigma!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS Sigma!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS Sigma!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS Sigma!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;

   cout <<"{";
   for(Int_t i=0; i<5;i++)cout << Sigma_param_all->GetParameter(i) <<",";
   cout << "}" << endl;




   printcanvases_.push_back(X0_red);
   printcanvases_.push_back(Xall);
   printcanvases_.push_back(B_red);
   printcanvases_.push_back(Ball);
   printcanvases_.push_back(Ball_kalibri);
   printcanvases_.push_back(C_red);
   printcanvases_.push_back(Call);
   printcanvases_.push_back(Call_kalibri);
   printcanvases_.push_back(Sigma1);
   printcanvases_.push_back(Sigmaall);


   TString SaveTo = "Fitted_paras_";
     if(chdir(SaveTo) != 0){ 
       mkdir(SaveTo, S_IRWXU|S_IRWXG|S_IRWXO); 
       chdir(SaveTo); 
     } 


   for(unsigned int
	 i=0;i<printcanvases_.size();i++)printcanvases_[i]->Print(X_labels_img_[bin_choice]+"_"+X_labels_img_[X_choice]+"_"+
   eta_region_labels_[eta_choice]+"_"+"_" + GMP_choice_[GMP_choice] +
   printcanvases_[i]->GetName()+ img);

   chdir("..");
   TString fit_para_name= "Fitted_paras_" + eta_region_labels_[eta_choice]  + X_labels_[bin_choice] + ".root";
   TFile *outf = new TFile(fit_para_name,"RECREATE");

   Sigma_param_all->Write();
   X0_param_all->Write();
   B_param_all->Write();
   C_param_all->Write();

   outf->Close();

   cout << "The chosen eta-region was: " << eta_region_labels_[eta_choice] << endl;
   cout <<"" << endl;

   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS X0!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;
   cout <<"{";
   for(Int_t i=0; i<5;i++)cout << X0_param_all->GetParameter(i) <<",";
   cout << "}" << endl;


   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS Sigma!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;

   cout <<"{";
   for(Int_t i=0; i<5;i++)cout << Sigma_param_all->GetParameter(i) <<",";
   cout << "}" << endl;

   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS BBBBBBBB!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;

   cout <<"";
   for(Int_t i=0; i<5;i++)cout << B_param_all->GetParameter(i) <<" ";
   cout << "" << endl;

   cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOO   !!!!!PARAMETERS CCCCC!!!! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<< endl;

   cout <<"";
   for(Int_t i=0; i<5;i++)cout << C_param_all->GetParameter(i) <<" ";
   cout << "" << endl;



}
