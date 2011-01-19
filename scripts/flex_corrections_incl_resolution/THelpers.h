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
#include <TObject.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TH1.h>
#include <THStack.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>


//Helper functions for analysis
//
//A: Delta_X-like functions
//B: Define_Canvas
//C: manipulate_filenames
//D: Export multiple TH1D/TH2D/TGraphErrors to any image extension
//   or multipage PS-files
//   - options including PS_name, canvas-format, draw_options, Logxyz of Canvas, 
//     SetRangeUser (x,y(,z where appropriate)), image type, Subfolder
//   - TH2D are drawn with an overlayed profile by default
//
//
//
//
//
//


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






std::pair <Int_t, Int_t> define_Canvas(TString format)
{
  std::pair <Int_t, Int_t> dimensions;
if(format == "six") dimensions = make_pair(420,460);
 else if(format == "square") dimensions = make_pair(501,501);
 else if(format == "nice") dimensions = make_pair(1024,768);
 else if(format == "flat") dimensions = make_pair(1024,740);
 else {cout << "No valid format defined, use default 'six' instead" << endl; dimensions = make_pair(420,460);
    }
 return dimensions;
}



TString manipulate_filenames(TCanvas* draw_canvas, TH1* histo, TString filename, TString SaveTo, TString draw_options, TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1)
{
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

     if(filename.Contains("/")) { 
       filename.ReplaceAll("/","_"); 
     } 

     filename=filename + "_"+draw_options; 

  if(Logxyz.Contains("x1")){cout << "set logx" << endl; draw_canvas->SetLogx();filename=filename + "_logx_";} 
  if(Logxyz.Contains("y1")){cout << "set logy" << endl; draw_canvas->SetLogy();filename=filename + "_logy_";} 
  if(Logxyz.Contains("z1")){cout << "set logz" << endl; draw_canvas->SetLogz();filename=filename + "_logz_";} 

  if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();} 
  if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();} 
  if(ru_zlow!=0&&ru_zhig!=-1){histo->GetZaxis()->SetRangeUser(ru_zlow,ru_zhig);testout << filename <<"_ruz_" << ru_zlow <<"_"<<ru_zhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();} 
  
     if(chdir(SaveTo) != 0){ 
       mkdir(SaveTo, S_IRWXU|S_IRWXG|S_IRWXO); 
       chdir(SaveTo); 
     } 

     return filename;
}

TString manipulate_filenames(TCanvas* draw_canvas, THStack* histo, TString filename, TString SaveTo, TString draw_options, TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1)
{
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);
     if(filename.Contains("/")) { 
       filename.ReplaceAll("/","_"); 
     } 

     filename=filename + "_"+draw_options; 

  if(Logxyz.Contains("x1")){cout << "set logx" << endl; draw_canvas->SetLogx();filename=filename + "_logx_";} 
  if(Logxyz.Contains("y1")){cout << "set logy" << endl; draw_canvas->SetLogy();filename=filename + "_logy_";} 
  //  if(Logxyz.Contains("z1")){cout << "set logz" << endl; draw_canvas->SetLogz();filename=filename + "_logz_";} 

  histo->Draw();//needs to be drawed, otherwise axes do not exist properly


  if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();} 
  if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();} 
  
     if(chdir(SaveTo) != 0){ 
       mkdir(SaveTo, S_IRWXU|S_IRWXG|S_IRWXO); 
       chdir(SaveTo); 
     } 

     return filename;
}


TString manipulate_filenames(TCanvas* draw_canvas, TGraphErrors* histo, TString filename, TString SaveTo, TString draw_options, TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1)
{
  std::string s;
  std::stringstream testout;
  testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

     if(filename.Contains("/")) { 
       filename.ReplaceAll("/","_"); 
     } 

     filename=filename + "_"+draw_options; 

  if(Logxyz.Contains("x1")){cout << "set logx" << endl; draw_canvas->SetLogx();filename=filename + "_logx_";} 
  if(Logxyz.Contains("y1")){cout << "set logy" << endl; draw_canvas->SetLogy();filename=filename + "_logy_";} 
  if(Logxyz.Contains("z1")){cout << "set logz" << endl; draw_canvas->SetLogz();filename=filename + "_logz_";} 

  if(ru_xlow!=0&&ru_xhig!=-1){histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);testout << filename <<"_rux_" << ru_xlow <<"_"<<ru_xhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();} 
  if(ru_ylow!=0&&ru_yhig!=-1){histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);testout << filename <<"_ruy_" << ru_ylow <<"_"<<ru_yhig <<"_";s=testout.str();filename = s.c_str();cout << filename <<endl;testout.str("");testout.clear();} 
  
     if(chdir(SaveTo) != 0){ 
       mkdir(SaveTo, S_IRWXU|S_IRWXG|S_IRWXO); 
       chdir(SaveTo); 
     } 

     return filename;
}


void draw_TH1D_save_with_Gauss_Fit_PS(std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS") {

  std::pair <Int_t, Int_t> dimensions = define_Canvas(format);
  TCanvas* temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, dimensions.first, dimensions.second);
  temp_canvas->SetRightMargin(0.12);
  temp_canvas->SetTopMargin(0.13);
  //temp_canvas->SetGridy(); 
  //temp_canvas->SetGridx();
 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   TH1D *histo = histos_.at(i_histos);


    TString filename = PS_name;
    filename = manipulate_filenames(temp_canvas,histo,filename,SaveTo,draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig,ru_zlow,ru_zhig);

      histo->Draw(draw_options);
      if(histo->GetEntries()>150) histo->GetFunction("gaus")->Draw("same");
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName()+ Logxyz + img_exp);
      if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");

  chdir("..");
  
   }
    temp_canvas->Close();
}
void draw_TH1D_save_with_Gauss_Fit_PS(TString SaveTo, TString img_exp,std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {
  draw_TH1D_save_with_Gauss_Fit_PS( histos_, PS_name, format, draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig,ru_zlow,ru_zhig, img_exp, SaveTo);
}


void draw_TH1D_save_PS(std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS") {

  std::pair <Int_t, Int_t> dimensions = define_Canvas(format);
  TCanvas* temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, dimensions.first, dimensions.second);
  temp_canvas->SetRightMargin(0.12);
  temp_canvas->SetTopMargin(0.13);
  //temp_canvas->SetGridy(); 
  //temp_canvas->SetGridx();

 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   TH1D *histo = histos_.at(i_histos);


    TString filename = PS_name;
    filename = manipulate_filenames(temp_canvas,histo,filename,SaveTo,draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig,ru_zlow,ru_zhig);
 
      histo->Draw(draw_options);
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName()+ Logxyz + img_exp);
      if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");

  chdir("..");
  
   }
    temp_canvas->Close();
}



void draw_TH1D_save_PS( TString img_exp, std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1){
draw_TH1D_save_PS(histos_,PS_name, format, draw_options, Logxyz, ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow,ru_zhig, img_exp); 

}

void draw_TH1D_save_PS(TString SaveTo, TString img_exp, std::vector<TH1D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1){
draw_TH1D_save_PS(histos_,PS_name, format, draw_options, Logxyz, ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow,ru_zhig, img_exp, SaveTo); 

}

void draw_THStack_save_PS(TLegend *leg, std::vector<THStack*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="hist leg", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS") {

  std::pair <Int_t, Int_t> dimensions = define_Canvas(format);
  TCanvas* temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, dimensions.first, dimensions.second);
  temp_canvas->SetRightMargin(0.12);
  temp_canvas->SetTopMargin(0.13);
  //temp_canvas->SetGridy(); 
  //temp_canvas->SetGridx();

 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   THStack *histo = histos_.at(i_histos);


    TString filename = PS_name;
    filename = manipulate_filenames(temp_canvas,histo,filename,SaveTo,draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig);

    gPad->SetTheta(3.77);
    gPad->SetPhi(2.9);

    if(draw_options.Contains("hist"))
      {
      histo->Draw("hist");
      if(ru_xlow!=0&&ru_xhig!=-1)histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);
      if(ru_ylow!=0&&ru_yhig!=-1)histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);
      if(draw_options.Contains("legend"))leg->DrawClone();
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName()+ Logxyz + img_exp);
      if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");
      }


       if(draw_options.Contains("nostack"))
      {
      histo->Draw("nostack,e1p");
      if(ru_xlow!=0&&ru_xhig!=-1)histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);
      if(ru_ylow!=0&&ru_yhig!=-1)histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);
      histo->Draw("nostack,e1p");
      if(draw_options.Contains("legend"))leg->DrawClone();
      if(i_histos==0)      temp_canvas->Print((filename +"no_stack")+ ".ps[");
      temp_canvas->Print((filename +"no_stack") + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName()+ Logxyz+("no_stack" + img_exp));
      if(i_histos==histos_.size()-1)      temp_canvas->Print((filename +"no_stack") + ".ps]");
      }

       if(draw_options.Contains("norm"))
      {

	TIter next(histo->GetHists());
	TObject *obj;// =  next();
	while ((obj =  next()))
	  {
	    Double_t value_integral;
	    TH1D* temp=(TH1D*)obj;
	    value_integral = temp->Integral();
	    //	    value_integral = temp->Integral();
	    temp->Scale(1/value_integral);
	  }
//
      histo->Draw("nostack,e1p");
      if(ru_xlow!=0&&ru_xhig!=-1)histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);
      if(ru_ylow!=0&&ru_yhig!=-1)histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);
      histo->GetYaxis()->SetTitle("arbitrary units");
      histo->Draw("nostack,e1p");
      if(draw_options.Contains("legend"))leg->DrawClone();
      if(i_histos==0)      temp_canvas->Print((filename +"no_stack")+ ".ps[");
      temp_canvas->Print((filename +"no_stack") + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName()+ Logxyz+("no_stack_norm" + img_exp));
      if(i_histos==histos_.size()-1)      temp_canvas->Print((filename +"no_stack_norm") + ".ps]");
      }



       if(draw_options.Contains("lego"))
      {
      //      histo->Draw("nostack,e1p");
      histo->Draw("lego1");
      if(ru_xlow!=0&&ru_xhig!=-1)histo->GetXaxis()->SetRangeUser(ru_xlow,ru_xhig);
      if(ru_ylow!=0&&ru_yhig!=-1)histo->GetYaxis()->SetRangeUser(ru_ylow,ru_yhig);
      histo->Draw("lego1");
      if(draw_options.Contains("legend"))leg->DrawClone();
      if(i_histos==0)      temp_canvas->Print((filename +"lego")+ ".ps[");
      temp_canvas->Print((filename +"lego") + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName()+ Logxyz+("lego" + img_exp));
      if(i_histos==histos_.size()-1)      temp_canvas->Print((filename +"lego") + ".ps]");
      }

  chdir("..");

   }
    temp_canvas->Close();
}

void draw_THStack_save_PS(std::vector<THStack*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="hist", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS") {

  TLegend *temp_leg = new TLegend(0.5,0.85,0.51,0.86);
  draw_THStack_save_PS(temp_leg, histos_, PS_name, format,draw_options, Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig,img_exp,SaveTo);
}


void draw_TH2D_save_PS(std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT", Bool_t with_correlation_factors = false, TString SaveTo="PS_test", TString w_profile ="true") {

  std::pair <Int_t, Int_t> dimensions = define_Canvas(format);
  TCanvas* temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, dimensions.first, dimensions.second);
  temp_canvas->SetRightMargin(0.12);
  temp_canvas->SetTopMargin(0.13);
  //temp_canvas->SetGridy(); 
  //temp_canvas->SetGridx();

  /*  TCanvas* temp_canvas;
if(format == "six") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
 else if(format == "square") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 501, 501 );
 else if(format == "nice") temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 1024, 768 );
 else {cout << "No valid format defined, use default 'six' instead" << endl; temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, 420, 460 );
    }
  */

TPaveText pt(.7,.55,.9,.65,"NDC ARC br");
pt.SetFillColor(kWhite);
pt.SetBorderSize(5);
pt.SetLineColor(kBlack);
pt.SetLabel("Correlation:");



 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++) {
   TH2D *histo = histos_.at(i_histos);
   
   
   TString filename = PS_name;
   filename = manipulate_filenames(temp_canvas,histo,filename,SaveTo,draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig,ru_zlow,ru_zhig);
   
   TProfile* histo_profile;
   histo_profile = histo->ProfileX("histo_profile");
   histo_profile->SetLineWidth(4);    
   
   histo->Draw(draw_options);
   histo_profile->SetStats(0);
   if(w_profile.Contains("true"))histo_profile->Draw("same");//"true"
   if(w_profile.Contains("marker")){//"marker"
     histo_profile->Draw("same");
     histo_profile->Rebin(3);
     //     cout << "blatest" << endl;
     histo_profile->SetMarkerStyle(29);
     histo_profile->SetMarkerColor(1);
     histo_profile->SetMarkerSize(4);
     histo_profile->Draw("same");
     //     histo_profile->Draw("same");
   }
   if(with_correlation_factors) {
     char buffer [10];
     sprintf (buffer, "is %f", histo->GetCorrelationFactor());
     pt.Clear();
     pt.AddText( buffer);
     pt.Draw();
   }
   if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
   temp_canvas->Print(filename + ".ps");
   if(img_exp.Contains("."))temp_canvas->Print(histo->GetName()+ Logxyz + img_exp);
   if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");
   
   chdir("..");
   delete histo_profile;
 }
 temp_canvas->Close();
 
}

void draw_TH2D_save_PS(TString img_exp, std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_TH2D_save_PS(histos_, PS_name,format, draw_options,Logxyz,ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow, ru_zhig, img_exp);

}


void draw_TH2D_save_PS(Bool_t with_correlation_factors, TString img_exp, std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_TH2D_save_PS(histos_, PS_name,format, draw_options,Logxyz,ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow, ru_zhig, img_exp,with_correlation_factors);

}
void draw_TH2D_save_PS(TString SaveTo,Bool_t with_correlation_factors, TString img_exp, std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_TH2D_save_PS(histos_, PS_name,format, draw_options,Logxyz,ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow, ru_zhig, img_exp,with_correlation_factors,SaveTo);

}

void draw_TH2D_save_PS(TString w_profile, TString SaveTo,Bool_t with_correlation_factors, TString img_exp, std::vector<TH2D*> histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_TH2D_save_PS(histos_, PS_name,format, draw_options,Logxyz,ru_xlow, ru_xhig, ru_ylow, ru_yhig, ru_zlow, ru_zhig, img_exp,with_correlation_factors,SaveTo, w_profile);

}

void draw_Colz_Prof_GMP_save_PS(std::vector<TH2D*> histos_colz_,std::vector<TProfile*> histos_prof_, std::vector<TH1D*> histos_GMP_, TString PS_name ="DEFAULT_PS", TString format = "nice", TString draw_options ="colz", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS") {

  std::pair <Int_t, Int_t> dimensions = define_Canvas(format);
  TCanvas* temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, dimensions.first, dimensions.second);
  temp_canvas->SetRightMargin(0.12);
  temp_canvas->SetTopMargin(0.13);
  //temp_canvas->SetGridy(); 
  //temp_canvas->SetGridx();

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
    filename = manipulate_filenames(temp_canvas,histo_colz,filename,SaveTo,draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig,ru_zlow,ru_zhig);

      histo_colz->Draw(draw_options);
      histo_prof->SetStats(0);
      histo_prof->Draw("same");
      histo_GMP->SetStats(0);
      histo_GMP->Draw("same");
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo_colz->GetName() +("_2D_prof_GMP"  + Logxyz+ img_exp));
      if(i_histos==histos_colz_.size()-1)      temp_canvas->Print(filename + ".ps]");

  chdir("..");
  
   }
    temp_canvas->Close();
}


void draw_Colz_Prof_GMP_save_PS(TString img_exp, std::vector<TH2D*> histos_colz_,std::vector<TProfile*> histos_prof_, std::vector<TH1D*> histos_GMP_, TString PS_name ="DEFAULT_PS", TString format = "nice", TString draw_options ="colz", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_Colz_Prof_GMP_save_PS(histos_colz_,histos_prof_,histos_GMP_,PS_name,format,draw_options,Logxyz,  ru_xlow, ru_xhig,ru_ylow, ru_yhig, ru_zlow, ru_zhig,img_exp);


}

void draw_Colz_Prof_GMP_save_PS(TString SaveTo, TString img_exp, std::vector<TH2D*> histos_colz_,std::vector<TProfile*> histos_prof_, std::vector<TH1D*> histos_GMP_, TString PS_name ="DEFAULT_PS", TString format = "nice", TString draw_options ="colz", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {

  draw_Colz_Prof_GMP_save_PS(histos_colz_,histos_prof_,histos_GMP_,PS_name,format,draw_options,Logxyz,  ru_xlow, ru_xhig,ru_ylow, ru_yhig, ru_zlow, ru_zhig,img_exp,SaveTo);


}



void draw_TH2D_save_Png(TH2D* histo, TString format = "six", TString draw_options ="", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, Double_t ru_zlow=0, Double_t ru_zhig=-1) {
  std::pair <Int_t, Int_t> dimensions = define_Canvas(format);
  TCanvas* temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, dimensions.first, dimensions.second);
  temp_canvas->SetRightMargin(0.12);
  temp_canvas->SetTopMargin(0.13);
  //temp_canvas->SetGridy(); 
  //temp_canvas->SetGridx();

    TString filename = histo->GetName();
    filename = manipulate_filenames(temp_canvas,histo,filename,"png",draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig,ru_zlow,ru_zhig);
 TProfile* histo_profile;
 histo_profile = histo->ProfileX("histo_profile");
    
    TString pngFolder = "pngs";
    if(chdir(pngFolder) != 0){
      mkdir(pngFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(pngFolder);
    }

      histo->Draw(draw_options);
      histo_profile->Draw("same");
      temp_canvas->Print(filename + ".png");   

    temp_canvas->Close();
    //    temp_canvas->Delete();

  chdir("..");
}


void draw_TGraphErrors_save_PS(std::vector < TGraphErrors* > histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="ALP", TString resolution_fit = "no", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1, TString img_exp="ENTER_WITH_POINT", TString SaveTo="PS") {

  std::pair <Int_t, Int_t> dimensions = define_Canvas(format);
  TCanvas* temp_canvas = new TCanvas("temp_canvas", "temp_canvas", 1, 1, dimensions.first, dimensions.second);
  temp_canvas->SetRightMargin(0.12);
  temp_canvas->SetTopMargin(0.13);
  //temp_canvas->SetGridy(); 
  //temp_canvas->SetGridx();

  TF1 *resolution = new TF1("resolution","TMath::Sqrt(TMath::Power([0]/TMath::Sqrt(x),2)+TMath::Power([1]/x,2)+TMath::Power([2],2))",30,1000); //was used before...
 resolution->SetParName(0,"S");
 resolution->SetParName(1,"N");
 resolution->SetParName(2,"C");
//TF1 *resolution = new TF1("resolution","TMath::Sqrt(TMath::Power([0]/TMath::Sqrt(x),2)+TMath::Power([1],2))",30,1000);
// resolution->SetParName(0,"S");
// resolution->SetParName(1,"C");
 resolution->SetLineColor(kRed);
 Double_t fit_limit=1000.;

 for (unsigned int i_histos=0;i_histos<histos_.size();i_histos++)
   //    cout << " " << *it;
   {
   TGraphErrors *histo = histos_.at(i_histos);


    TString filename = PS_name;
    filename = manipulate_filenames(temp_canvas,histo,filename,SaveTo,draw_options,Logxyz,ru_xlow,ru_xhig,ru_ylow,ru_yhig);
      histo->Draw(draw_options);
      if(resolution_fit.Contains("yes"))histo->Fit("resolution","EM","same",30,fit_limit);
      if(i_histos==0)      temp_canvas->Print(filename + ".ps[");
      temp_canvas->Print(filename + ".ps");
      if(img_exp.Contains("."))temp_canvas->Print(histo->GetName() +("_TGraphErr" + Logxyz+  img_exp));
      if(i_histos==histos_.size()-1)      temp_canvas->Print(filename + ".ps]");
      

  chdir("..");
  
   }
    temp_canvas->Close();
}

void draw_TGraphErrors_save_PS( TString img_exp, std::vector < TGraphErrors* > histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="ALP", TString resolution_fit = "no", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1) {

  draw_TGraphErrors_save_PS(histos_, PS_name ,format, draw_options,  resolution_fit,Logxyz, ru_xlow, ru_xhig, ru_ylow,ru_yhig,img_exp);
}

void draw_TGraphErrors_save_PS(TString SaveTo, TString img_exp, std::vector < TGraphErrors* > histos_, TString PS_name ="DEFAULT_PS", TString format = "six", TString draw_options ="ALP", TString resolution_fit = "no", TString Logxyz="x0_y0_z0", Double_t ru_xlow=0, Double_t ru_xhig=-1, Double_t ru_ylow=0, Double_t ru_yhig=-1) {

  draw_TGraphErrors_save_PS(histos_, PS_name ,format, draw_options,  resolution_fit,Logxyz, ru_xlow, ru_xhig, ru_ylow,ru_yhig,img_exp,SaveTo);
}




#endif
