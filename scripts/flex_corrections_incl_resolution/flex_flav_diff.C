#define flex_flav_diff_cxx
#include "THelpers.h"
#include "base_corr.C"
#include "flex_flav_diff.h"

void flex_flav_diff::Import_Histos()
{
TString  root_resol_name_binning;




 for(unsigned int eta_i=0;eta_i<sel_eta_.size();eta_i++)  
   {
     std::vector <TFile*> in_PDG_; 
     std::vector <std::vector <TGraphErrors*> > tlj_Response_Graphs_mean_selec_PDG_;
     for(unsigned int PDG_i=0;PDG_i<sel_PDG_.size();PDG_i++)  
       {
	 root_resol_name_binning= "histos_resolution_" + PDG_labels_[sel_PDG_[PDG_i]]  +"_" + eta_region_labels_[sel_eta_[eta_i]] +"_"+ X_labels_[bin_choice] + ".root";

	 in_PDG_.push_back(new TFile(root_resol_name_binning,"OPEN"));
	 std::vector <TGraphErrors*> tlj_Response_Graphs_mean_selec_PDG_eta_;
	 for(Int_t Corr_i =0; Corr_i<no_Corr_labels_;Corr_i++)
	   {
	     tlj_Response_Graphs_mean_selec_PDG_eta_.push_back((TGraphErrors*)in_PDG_.back()->Get(X_labels_[X_choice] +"vs.Mean_" + Corr_labels_[Corr_i].first));
	   }
	     tlj_Response_Graphs_mean_selec_PDG_.push_back(tlj_Response_Graphs_mean_selec_PDG_eta_);

       }
     infiles_.push_back(in_PDG_);
     tlj_Response_Graphs_mean_selec_.push_back(tlj_Response_Graphs_mean_selec_PDG_);
   }


 //tlj_Response_Graphs_mean_selec_[eta_i].at(PDG_i).at(Corr_i)

     
    TString saveFolder = "Flavour_Comparison";
    if(chdir(saveFolder) != 0){
      mkdir(saveFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(saveFolder);
    }

    cout <<"created folder"<< endl;


       //Personal selection Mean, default X
    //       tlj_Response_Graphs_mean_selec_.push_back(make_graph(Double_gauss_all_[Corr_selected_labels_[Sel_Corr_i]], 0, "Selected_" + X_labels_[X_choice] +"vs.Mean_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first, X_choice));

}


void flex_flav_diff::Book_Histos()
{


}

void flex_flav_diff::Fit_Histos()
{
}

void flex_flav_diff::Write_Histos()
{


  std::vector < TGraphErrors* > flav_diff_selected_;
  std::vector < TGraphErrors* > flav_mean_selected_flav_1_;
  std::vector < TGraphErrors* > flav_mean_selected_flav_2_;

  for(unsigned int Sel_Corr_i=0; Sel_Corr_i< Corr_selected_labels_.size();Sel_Corr_i++)
    {
      flav_diff_selected_.push_back((TGraphErrors*)tlj_Response_Graphs_mean_selec_[0].at(sel_PDG_.at(0)).at(Corr_selected_labels_[Sel_Corr_i])->Clone());
      flav_mean_selected_flav_2_.push_back((TGraphErrors*)tlj_Response_Graphs_mean_selec_[0].at(sel_PDG_.at(1)).at(Corr_selected_labels_[Sel_Corr_i])->Clone());
      flav_mean_selected_flav_1_.push_back((TGraphErrors*)tlj_Response_Graphs_mean_selec_[0].at(sel_PDG_.at(0)).at(Corr_selected_labels_[Sel_Corr_i])->Clone());
      cout << "no. of  points " << flav_diff_selected_.back()->GetN() << endl;

//       flav_diff_selected_.back()->SetTitle("");
//       flav_mean_selected_flav_1_.back()->SetTitle("");
//       flav_mean_selected_flav_2_.back()->SetTitle("");

      flav_diff_selected_.back()->SetLineWidth(3);
      //      flav_diff_selected_.back()->SetTitle("");
      flav_mean_selected_flav_1_.back()->SetLineWidth(3);
      flav_mean_selected_flav_2_.back()->SetLineWidth(3);

      for(Int_t p_i=0; p_i< flav_diff_selected_.back()->GetN();p_i++)
	{
	  flav_diff_selected_.back()->SetPoint(p_i,flav_diff_selected_.back()->GetX()[p_i], - flav_mean_selected_flav_2_.back()->GetY()[p_i] + flav_diff_selected_.back()->GetY()[p_i] );
	}
    }

      draw_TGraphErrors_save_PS("PS_Flav_diff_mean",img_choice,flav_diff_selected_,PDG_labels_[sel_PDG_.at(0)] + "-" + PDG_labels_[sel_PDG_.at(1)] + eta_region_labels_[sel_eta_[0]] + X_labels_[bin_choice] + "Response-difference","nice","ALP","no","x1",0,-1,-0.5,0.5);


      draw_TGraphErrors_save_PS("TEST_PS_Flav_diff_mean",img_choice,flav_mean_selected_flav_2_,PDG_labels_[sel_PDG_.at(0)] + "-" + PDG_labels_[sel_PDG_.at(1)] + eta_region_labels_[sel_eta_[0]] + X_labels_[bin_choice] + "Response_flavour2","nice","ALP","no","x1",0,-1,0.5,1.5);
      draw_TGraphErrors_save_PS("TEST_PS_Flav_diff_mean",img_choice,flav_mean_selected_flav_1_,PDG_labels_[sel_PDG_.at(0)] + "-" + PDG_labels_[sel_PDG_.at(1)] + eta_region_labels_[sel_eta_[0]] + X_labels_[bin_choice] + "Response_flavour1","nice","ALP","no","x1",0,-1,0.5,1.5);


 TLegend *leg_selected;
 leg_selected = new TLegend(0.7,0.7-0.03*Corr_selected_labels_.size(),0.95,0.85);
 leg_selected->SetFillColor(kWhite);
  leg_selected->SetBorderSize(0);
  leg_selected->SetFillColor(0);
  leg_selected->SetTextFont(42);
  leg_selected->SetTextSize(0.04);

 //   leg->SetHeader("Legende");
 for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
   {
   leg_selected->AddEntry(flav_diff_selected_[Sel_Corr_i],Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].second,"lep");
   }

 // draw_graphs( flav_diff_selected_, 0., 0.15, leg_selected,  PDG_labels_[PDG_choice] + eta_region_labels_[eta_choice] + X_labels_[bin_choice] +"RESPONSE_Flavour_difference",0,1000);
 Double_t ylow =-0.03;
 Double_t yhigh = 0.20;
 Double_t xlow = 0.;
 Double_t xhigh = 1000;


 TString PDF_PNG_name = PDG_labels_[PDG_choice] + eta_region_labels_[eta_choice] + X_labels_[bin_choice] +"RESPONSE_Flavour_difference";

 for(unsigned int i=0;i<Corr_selected_labels_.size();i++)
   {
     PDF_PNG_name +=Corr_labels_[Corr_selected_labels_[i]].first;
     PDF_PNG_name +="_";
   }



if(flav_diff_selected_.size()>=1)
   {
     flav_diff_selected_[0]->Draw("ALP");
     flav_diff_selected_[0]->SetTitle("");
     flav_diff_selected_[0]->GetYaxis()->SetTitle("R_{uds}-R_{gluon}");
     flav_diff_selected_[0]->GetXaxis()->SetTitle(flav_diff_selected_[0]->GetXaxis()->GetTitle() + (TString) "[GeV]");
     flav_diff_selected_[0]->GetYaxis()->SetRangeUser(ylow,yhigh);
     flav_diff_selected_[0]->GetXaxis()->SetRangeUser(xlow,xhigh);
     flav_diff_selected_[0]->Draw("ALP");
   }
 int markers[11]={20,24,21,25,22,26,29,30,23,28,34};

 for(unsigned int a_i=0; a_i<flav_diff_selected_.size();a_i++)
   {
     flav_diff_selected_[a_i]->SetLineColor(a_i+1);
     flav_diff_selected_[a_i]->SetMarkerColor(a_i+1);
     flav_diff_selected_[a_i]->SetMarkerSize(1.8);
     flav_diff_selected_[a_i]->SetLineStyle(a_i+1);
     flav_diff_selected_[a_i]->SetLineWidth(2);
     if(a_i<11)flav_diff_selected_[a_i]->SetMarkerStyle(markers[a_i]);
     flav_diff_selected_[a_i]->Draw("LP same"); 
   }
 leg_selected->Draw();

 TLine *line = new TLine(xlow,0.0,xhigh,0.0);
 line->SetLineStyle(2);
 line->SetLineColor(1);
 line->Draw();

 test->Print(PDF_PNG_name+".root");
 test->Print(PDF_PNG_name+".png");
 test->Print(PDF_PNG_name+".eps");

 test->SetLogx();

 test->Print(PDF_PNG_name+"_logX_" +".root");
 test->Print(PDF_PNG_name+"_logX_" +".png");
 test->Print(PDF_PNG_name+"_logX_" +".eps");

 test->SetLogx(0);


  //void SetPoint(Int_t i, Double_t  x, Double_t  y)



//   for(Int_t Corr_i =0; Corr_i<Corr_selected_labels_.size();Corr_i++)
//     {
      
//     }


//   char buffer[5];
//    gStyle->SetOptStat(0);
//   gStyle->SetOptFit(0000);

//   TString Sel_Response = "Corrected_Response_selected";
//   test->SetLogy(1);
//   for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
//     {
//       tlj_corrected_response_barrel_all_[Corr_selected_labels_[Corr_selected_labels_.size()-1]].at(pt_i)->Draw();  
//       tlj_corrected_response_barrel_all_[Corr_selected_labels_[Corr_selected_labels_.size()-1]].at(pt_i)->SetStats(0);  
//       for(unsigned int Sel_Corr_i=0; Sel_Corr_i< Corr_selected_labels_.size()-1;Sel_Corr_i++)
// 	{
// 	  tlj_corrected_response_barrel_all_[Corr_selected_labels_[Sel_Corr_i]].at(pt_i)->SetLineColor(Sel_Corr_i+2);  
// 	  tlj_corrected_response_barrel_all_[Corr_selected_labels_[Sel_Corr_i]].at(pt_i)->SetStats(0);  
// 	  tlj_corrected_response_barrel_all_[Corr_selected_labels_[Sel_Corr_i]].at(pt_i)->Draw("same");  
// 	}
//       sprintf (buffer, "%f", pt_bins_[pt_i].first);
//       test->Print((  Sel_Response +"_"+   buffer    ) + ".pdf");
//     }
//   test->SetLogy(0);
//    gStyle->SetOptStat(1);
//    gStyle->SetOptFit(1);



}




void flex_flav_diff::Create_TGraphErrors()
{

}

void flex_flav_diff::Write_TGraphErrors()
{

}


void flex_flav_diff::Loop(Bool_t setvalues, Int_t par_bin_choice, Int_t par_X_choice, TString par_eta_choice, TString par_Corr_choice, TString par_img_choice, TString par_PDG_choice)
{
//   In a ROOT session, you can do:
//      Root > .L flex_flav_diff.C
//      Root > flex_flav_diff t
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
 string_eta_choice	= par_eta_choice;	
 Corr_choice    = par_Corr_choice;
 img_choice	= par_img_choice;
 string_PDG_choice     = par_PDG_choice;
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
   if(!setvalues){
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>X_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << X_choice << ", also " << X_labels_[X_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Pseudorapiditäts-Regionen auswählen:" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_eta_region_;i++)  cout << i << ": " << eta_region_labels_[i] << endl;
   cin >>string_eta_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << string_eta_choice << " ausgewaehlt. Danke... Gewaehlt wurde also:" << endl;


   for(Int_t i=0;i<no_eta_region_;i++)  
     {
       if(string_eta_choice.Contains((char) (i+48)))
	 {
	   cout << i << ": " << eta_region_labels_[i] << endl;
	   sel_eta_.push_back(i);
	 }
     }
   if( sel_eta_.size()==0)cout << "Es wurde nichts ausgewaehlt!!!" << endl;


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


   cout << "Bitte PDG-Auswahl treffen wählen:" <<endl;
   if(!setvalues){
   for(Int_t i=0;i< no_PDG_labels_;i++)  cout << i << ": " <<  PDG_labels_[i] << endl;
   cin >>string_PDG_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << string_PDG_choice <<" ausgewaehlt. Danke... Gewaehlt wurde also:" << endl;


   for(Int_t i=0;i<no_PDG_labels_;i++)  
     {
       //       cout << (char) (i+48) << endl;
       if(string_PDG_choice.Contains((char) (i+48)))
	 {
	   cout << i << ": " << PDG_labels_[i] << endl;
	   sel_PDG_.push_back(i);
	 }
     }
   if( sel_PDG_.size()==0)cout << "Es wurde nichts ausgewaehlt!!!" << endl;



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
