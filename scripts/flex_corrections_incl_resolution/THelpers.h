#ifndef THelpers_h
#define THelpers_h

#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>



double delta_r(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2)
{
Double_t delta = TMath::Sqrt(TMath::Power((eta1-eta2),2)+TMath::Power((phi1-phi2),2));
 return delta;
}

template <class T> 
T deltaPhi (T phi1, T phi2) { 
  T result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

template <class T>
T deltaR2 (T eta1, T phi1, T eta2, T phi2) {
  T deta = eta1 - eta2;
  T dphi = deltaPhi (phi1, phi2);
  return deta*deta + dphi*dphi;
}

template <class T>
T deltaR (T eta1, T phi1, T eta2, T phi2) {
  return sqrt (deltaR2 (eta1, phi1, eta2, phi2));
}

template<typename T1, typename T2>
double deltaR2( const T1 & t1, const T2 & t2 ) {
  return deltaR2( t1.eta(), t1.phi(), t2.eta(), t2.phi() );
} 

template<typename T1, typename T2>
double deltaR( const T1 & t1, const T2 & t2 ) {
  return deltaR( t1.eta(), t1.phi(), t2.eta(), t2.phi() );
} 










void draw_TH1D_save_with_Gauss_Fit_PS(std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  TCanvas* temp_canvas;
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

if(format == "six") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
 else if(format == "square") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 501, 501 );
 else if(format == "nice") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 800, 600 );
 else {cout << "No valid format defined, use default 'six' instead" << endl; temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
    }


 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   TH1D *histo = histos_.at(i_histos);


    TString filename = PS_name;
    if(filename.Contains("/")) {
      filename.ReplaceAll("/","_");
    }

    filename=filename + "_"+draw_options;

 if(Logxyz.Contains("x1")){cout << "set logx" << endl; temp_canvas->SetLogx();filename=filename + "_logx_";}
 if(Logxyz.Contains("y1")){cout << "set logy" << endl; temp_canvas->SetLogy();filename=filename + "_logy_";}
 if(Logxyz.Contains("z1")){cout << "set logz" << endl; temp_canvas->SetLogz();filename=filename + "_logz_";}

 if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_zlow!=0&&ru_zhig!=-1){histo->GetZaxis()->SetRangeUser(ru_zlow,ru_zhig);testout << filename <<"_ruz_" << ru_zlow <<"_"<<ru_zhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
    
    TString pngFolder = "PS";
    if(chdir(pngFolder) != 0){
      mkdir(pngFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(pngFolder);
    }
 
    if(format.Contains("norm")){
      cout << filename <<endl;
      histo->Draw(draw_options);
      cout << "Doublegauss Anzahl der Eintraege: " << histo->GetEntries() << endl;
      if(histo->GetEntries()>100)histo->GetFunction("gaus")->Draw("same");
      if(i_histos==0)temp_canvas->Print(filename + "_norm.ps[");
      temp_canvas->Print(filename + "_norm.ps");
      if(i_histos==histos_.size()-1)temp_canvas->Print(filename + "_norm.ps]");
      
      
    }
    else{
      histo->Draw(draw_options);
      if(histo->GetEntries()>150) histo->GetFunction("gaus")->Draw("same");
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");
      
      
    }

    //    temp_canvas->Delete();

  chdir("..");
  
   }
    temp_canvas->Close();
}



void draw_TH1D_save_PS(std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT") {

  TCanvas* temp_canvas;
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

  //  if(histos_.size()<2){cout << "Fehler" <<endl;}
  //  else
  //    {
  //  for(Int_t a_i=0;a_i<no_pt_bins;a_i++)

// for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
//   cout << "bla" << endl;

if(format == "six") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
 else if(format == "square") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 501, 501 );
 else if(format == "nice") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 800, 600 );
 else {cout << "No valid format defined, use default 'six' instead" << endl; temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
    }



 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   TH1D *histo = histos_.at(i_histos);


    TString filename = PS_name;
    if(filename.Contains("/")) {
      filename.ReplaceAll("/","_");
    }

    filename=filename + "_"+draw_options;

 if(Logxyz.Contains("x1")){cout << "set logx" << endl; temp_canvas->SetLogx();filename=filename + "_logx_";}
 if(Logxyz.Contains("y1")){cout << "set logy" << endl; temp_canvas->SetLogy();filename=filename + "_logy_";}
 if(Logxyz.Contains("z1")){cout << "set logz" << endl; temp_canvas->SetLogz();filename=filename + "_logz_";}

 if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_zlow!=0&&ru_zhig!=-1){histo->GetZaxis()->SetRangeUser(ru_zlow,ru_zhig);testout << filename <<"_ruz_" << ru_zlow <<"_"<<ru_zhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
    
    TString pngFolder = "PS";
    if(chdir(pngFolder) != 0){
      mkdir(pngFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(pngFolder);
    }
 
    if(format.Contains("norm")){
      cout << filename <<endl;
      histo->Draw(draw_options);
      if(i_histos==0)temp_canvas->Print(filename + "_norm.ps[");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName() + ("_norm" + img_exp));
      temp_canvas->Print(filename + "_norm.ps");
      if(i_histos==histos_.size()-1)temp_canvas->Print(filename + "_norm.ps]");
      
      
    }
    else{
      histo->Draw(draw_options);
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName() + img_exp);
      if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");
      
      
    }

    //    temp_canvas->Delete();

  chdir("..");
  
   }
    temp_canvas->Close();
}



void draw_TH1D_save_PS( TString img_exp, std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1){
draw_TH1D_save_PS(histos_,PS_name, format, draw_options, Logxyz, ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow,ru_zhig, img_exp); 

}



void draw_TH2D_save_PS(std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT", Bool_t with_correlation_factors = false) {

  TCanvas* temp_canvas;
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

if(format == "six") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
 else if(format == "square") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 501, 501 );
 else if(format == "nice") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 800, 600 );
 else {cout << "No valid format defined, use default 'six' instead" << endl; temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
    }

TPaveText pt(.7,.55,.9,.65,"NDC ARC br");
pt.SetFillColor(kWhite);
pt.SetBorderSize(5);
pt.SetLineColor(kBlack);
pt.SetLabel("Correlation:");





 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   TH2D *histo = histos_.at(i_histos);


    TString filename = PS_name;
    if(filename.Contains("/")) {
      filename.ReplaceAll("/","_");
    }

    filename=filename + "_"+draw_options;

 if(Logxyz.Contains("x1")){cout << "set logx" << endl; temp_canvas->SetLogx();filename=filename + "_logx_";}
 if(Logxyz.Contains("y1")){cout << "set logy" << endl; temp_canvas->SetLogy();filename=filename + "_logy_";}
 if(Logxyz.Contains("z1")){cout << "set logz" << endl; temp_canvas->SetLogz();filename=filename + "_logz_";}

 if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_zlow!=0&&ru_zhig!=-1){histo->GetZaxis()->SetRangeUser(ru_zlow,ru_zhig);testout << filename <<"_ruz_" << ru_zlow <<"_"<<ru_zhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 TProfile* histo_profile;
 histo_profile = histo->ProfileX("histo_profile");
 histo_profile->SetLineWidth(4);    

    TString pngFolder = "PS";
    if(chdir(pngFolder) != 0){
      mkdir(pngFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(pngFolder);
    }




    if(format.Contains("norm")){
cout << filename <<endl;
      histo->Draw(draw_options);
      histo_profile->SetStats(0);
      histo_profile->Draw("same");

      if(with_correlation_factors)
	{
	  char buffer [10];
	  sprintf (buffer, "is %f", histo->GetCorrelationFactor());
	  pt.Clear();
	  pt.AddText( buffer);
	  pt.Draw();
	}
      if(i_histos==0)temp_canvas->Print(filename + "_norm.ps[");
      temp_canvas->Print(filename + "_norm.ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName() + ("_norm" + img_exp));
      if(i_histos==histos_.size()-1)temp_canvas->Print(filename + "_norm.ps]");
      
      
    }
    else{
      histo->Draw(draw_options);
      histo_profile->SetStats(0);
      histo_profile->Draw("same");
      if(with_correlation_factors)
	{
	  char buffer [10];
	  sprintf (buffer, "is %f", histo->GetCorrelationFactor());
	  pt.Clear();
	  pt.AddText( buffer);
	  pt.Draw();
	}
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName() + img_exp);
     if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");
      
      
    }

    //    temp_canvas->Delete();

  chdir("..");
  
   }
    temp_canvas->Close();

}

void draw_TH2D_save_PS(TString img_exp, std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_TH2D_save_PS(histos_, PS_name,format, draw_options,Logxyz,ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow, ru_zhig, img_exp);

}

void draw_TH2D_save_PS(Bool_t with_correlation_factors, TString img_exp, std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_TH2D_save_PS(histos_, PS_name,format, draw_options,Logxyz,ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow, ru_zhig, img_exp,with_correlation_factors);

}

void draw_Colz_Prof_GMP_save_PS(std::vector<TH2D*> histos_colz_,std::vector<TProfile*> histos_prof_, std::vector<TH1D*> histos_GMP_, TString PS_name ="DEFAULT_PS", TString format = "nice", TString draw_options ="colz", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT") {

  TCanvas* temp_canvas;
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

if(format == "six") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
 else if(format == "square") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 501, 501 );
 else if(format == "nice") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 800, 600 );
 else {cout << "No valid format defined, use default 'six' instead" << endl; temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
    }


 for (unsigned int i_histos=0;i_histos<histos_colz_.size();i_histos++)
   //    cout << " " << *it;
   {
   TH2D *histo_colz = histos_colz_.at(i_histos);
   TProfile *histo_prof = histos_prof_.at(i_histos);
   TH1D *histo_GMP = histos_GMP_.at(i_histos);
   histo_prof->SetLineWidth(4);    
   histo_GMP->SetLineWidth(4);    
   histo_GMP->SetLineColor(28);    




    TString filename = PS_name;
    if(filename.Contains("/")) {
      filename.ReplaceAll("/","_");
    }

    filename=filename + "_"+draw_options;

 if(Logxyz.Contains("x1")){cout << "set logx" << endl; temp_canvas->SetLogx();filename=filename + "_logx_";}
 if(Logxyz.Contains("y1")){cout << "set logy" << endl; temp_canvas->SetLogy();filename=filename + "_logy_";}
 if(Logxyz.Contains("z1")){cout << "set logz" << endl; temp_canvas->SetLogz();filename=filename + "_logz_";}

 if(ru_xlow!=0&&ru_xhig!=-1){histo_colz->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_ylow!=0&&ru_yhig!=-1){histo_colz->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_zlow!=0&&ru_zhig!=-1){histo_colz->GetZaxis()->SetRangeUser(ru_zlow,ru_zhig);testout << filename <<"_ruz_" << ru_zlow <<"_"<<ru_zhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}

    TString pngFolder = "PS";
    if(chdir(pngFolder) != 0){
      mkdir(pngFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(pngFolder);
    }

    if(format.Contains("norm")){
cout << filename <<endl;
      histo_colz->Draw(draw_options);
      histo_prof->SetStats(0);
      histo_prof->Draw("same");
      histo_GMP->SetStats(0);
      histo_GMP->Draw("same");
      if(i_histos==0)temp_canvas->Print(filename + "_norm.ps[");
      temp_canvas->Print(filename + "_norm.ps");
       if(img_exp.Contains("."))temp_canvas->Print(histo_colz->GetName() + ("_norm" + ("_2D_prof_GMP" + img_exp)));
     if(i_histos==histos_colz_.size()-1)temp_canvas->Print(filename + "_norm.ps]");
      
      
    }
    else{
      histo_colz->Draw(draw_options);
      histo_prof->SetStats(0);
      histo_prof->Draw("same");
      histo_GMP->SetStats(0);
      histo_GMP->Draw("same");
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo_colz->GetName() +("_2D_prof_GMP" +  img_exp));
      if(i_histos==histos_colz_.size()-1)      temp_canvas->Print(filename + ".ps]");
      
    }

    //    temp_canvas->Delete();

  chdir("..");
  
   }
    temp_canvas->Close();
}


void draw_Colz_Prof_GMP_save_PS(TString img_exp, std::vector<TH2D*> histos_colz_,std::vector<TProfile*> histos_prof_, std::vector<TH1D*> histos_GMP_, TString PS_name ="DEFAULT_PS", TString format = "nice", TString draw_options ="colz", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_Colz_Prof_GMP_save_PS(histos_colz_,histos_prof_,histos_GMP_,PS_name,format,draw_options,Logxyz,  ru_xlow, ru_xhig,ru_ylow, ru_yhig, ru_zlow, ru_zhig,img_exp);


}



void draw_TH2D_save_Png(TH2D* histo, TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {
  TCanvas* temp_canvas;
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);


if(format == "six") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
 else if(format == "square") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 501, 501 );
 else if(format == "nice") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 800, 600 );
 else {cout << "No valid format defined, use default 'six' instead" << endl; temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
    }

    TString filename = histo->GetName();
    if(filename.Contains("/")) {
      filename.ReplaceAll("/","_");
    }

    filename=filename + "_"+draw_options;

 if(Logxyz.Contains("x1")){cout << "set logx" << endl; temp_canvas->SetLogx();filename=filename + "_logx_";}
 if(Logxyz.Contains("y1")){cout << "set logy" << endl; temp_canvas->SetLogy();filename=filename + "_logy_";}
 if(Logxyz.Contains("z1")){cout << "set logz" << endl; temp_canvas->SetLogz();filename=filename + "_logz_";}

 if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_zlow!=0&&ru_zhig!=-1){histo->GetZaxis()->SetRangeUser(ru_zlow,ru_zhig);testout << filename <<"_ruz_" << ru_zlow <<"_"<<ru_zhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 TProfile* histo_profile;
 histo_profile = histo->ProfileX("histo_profile");
    
    TString pngFolder = "pngs";
    if(chdir(pngFolder) != 0){
      mkdir(pngFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(pngFolder);
    }

    if(format.Contains("norm")){
cout << filename <<endl;
      histo->Draw(draw_options);
      histo_profile->Draw("same");
      temp_canvas->Print(filename + "_norm.png");
    }
    else{
      histo->Draw(draw_options);
      histo_profile->Draw("same");
      temp_canvas->Print(filename + ".png");   
    }

    temp_canvas->Close();
    //    temp_canvas->Delete();

  chdir("..");
}


void draw_TGraphErrors_save_PS(std::vector < TGraphErrors* > histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="ALP", TString resolution_fit = "no", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, TString img_exp="ENTER_WITH_POINT") {


  TCanvas* temp_canvas;
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

  //  if(histos_.size()<2){cout << "Fehler" <<endl;}
  //  else
  //    {
  //  for(Int_t a_i=0;a_i<no_pt_bins;a_i++)

// for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
//   cout << "bla" << endl;


TF1 *resolution = new TF1("resolution","TMath::Sqrt(TMath::Power([0]/TMath::Sqrt(x),2)+TMath::Power([1]/x,2)+TMath::Power([2],2))",20,1000);
 resolution->SetParName(0,"S");
 resolution->SetParName(1,"N");
 resolution->SetParName(2,"C");
 Double_t fit_limit=1000.;


if(format == "six") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
 else if(format == "square") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 501, 501 );
 else if(format == "nice") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 800, 600 );
 else {cout << "No valid format defined, use default 'six' instead" << endl; temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
    }


 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   TGraphErrors *histo = histos_.at(i_histos);


    TString filename = PS_name;
    if(filename.Contains("/")) {
      filename.ReplaceAll("/","_");
    }

    filename=filename + "_"+draw_options;

 if(Logxyz.Contains("x1")){cout << "set logx" << endl; temp_canvas->SetLogx();filename=filename + "_logx_";}
 if(Logxyz.Contains("y1")){cout << "set logy" << endl; temp_canvas->SetLogy();filename=filename + "_logy_";}
 if(Logxyz.Contains("z1")){cout << "set logz" << endl; temp_canvas->SetLogz();filename=filename + "_logz_";}

 if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}
 if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();}

    TString pngFolder = "PS";
    if(chdir(pngFolder) != 0){
      mkdir(pngFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(pngFolder);
    }

    if(format.Contains("norm")){
cout << filename <<endl;
      histo->Draw(draw_options);
     if(resolution_fit.Contains("yes")) histo->Fit("resolution","EM","same",30,fit_limit);

      if(i_histos==0)temp_canvas->Print(filename + "_norm.ps[");
      temp_canvas->Print(filename + "_norm.ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName() +("norm_TGraphErr" +  img_exp));
      if(i_histos==histos_.size()-1)temp_canvas->Print(filename + "_norm.ps]");
      
      
    }
    else{
      histo->Draw(draw_options);
      if(resolution_fit.Contains("yes"))histo->Fit("resolution","EM","same",30,fit_limit);
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName() +("_TGraphErr" +  img_exp));
      if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");
      
      
    }

    //    temp_canvas->Delete();

  chdir("..");
  
   }
    temp_canvas->Close();
}

void draw_TGraphErrors_save_PS( TString img_exp, std::vector < TGraphErrors* > histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="ALP", TString resolution_fit = "no", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1) {

  draw_TGraphErrors_save_PS(histos_, PS_name ,format, draw_options,  resolution_fit,Logxyz, ru_xlow, ru_xhig, ru_ylow,ru_yhig,img_exp);
}




#endif
