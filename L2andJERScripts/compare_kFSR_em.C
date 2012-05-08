#include "tdrstyle_mod.C"
#include "MakeDateDir.h"

void compare_kFSR_em(TString dir1, TString dir2, TString algo="PF",
  TString dir_prefix="KOSTAS_L1_on_pt_plain_", TString label_dir1="", TString label_dir2="",TString binning_select ="kostas"){
  setTDRStyle();
  TCanvas* c = new TCanvas("c","",600,600);
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  Double_t kFSR_extrapol_ymin=0.95;
  Double_t kFSR_extrapol_ymax=1.1;


    TFile *inf_Dir1;
          inf_Dir1 = new TFile(dir_prefix+dir1+"/"+binning_select+"_use_easy_mean_TuneZ2_TuneZ2_" + algo + "_kFSR_histos.root","OPEN");
    //    inf_Dir1 = new TFile(dir_prefix+dir1+"/"+binning_select+"_TuneZ2_TuneZ2_" + algo + "_kFSR_histos.root","OPEN");

    TFile *inf_Dir2;
        inf_Dir2 = new TFile(dir_prefix+dir2+"/"+binning_select+"_use_easy_mean_TuneZ2_TuneZ2_" + algo + "_kFSR_histos.root","OPEN");
	//    inf_Dir2 = new TFile(dir_prefix+dir2+"/"+binning_select+"_TuneZ2_TuneZ2_" + algo + "_kFSR_histos.root","OPEN");
    if(label_dir1!="")dir1=label_dir1;
    if(label_dir2!="")dir2=label_dir2;
 
    Double_t chisquared1, ndf1, chisquared2, ndf2;
    char buffer [50];

    import_kFSR_vs_Abseta_histo_res1_Dir1  = (TH1D*)inf_Dir1->Get("kFSR_vs_Abseta_histo_res1");
    import_kFSR_vs_Abseta_histo_res1_Dir2  = (TH1D*)inf_Dir2->Get("kFSR_vs_Abseta_histo_res1");
    import_kFSR_vs_Abseta_histo_res1_Dir1->SetStats(0);
    import_kFSR_vs_Abseta_histo_res1_Dir2->SetStats(0);
    import_kFSR_vs_Abseta_histo_res1_Dir1->GetYaxis()->SetTitle("k_{rad} correction");
    import_kFSR_vs_Abseta_histo_res1_Dir2->GetYaxis()->SetTitle("k_{rad} correction");



    import_kFSR_vs_Abseta_histo_res1_Dir1->Draw("");
    import_kFSR_vs_Abseta_histo_res1_Dir1->GetYaxis()->SetRangeUser(kFSR_extrapol_ymin,kFSR_extrapol_ymax);
    
    TF1 *kFSR_fit = new TF1("kFSR_fit","[0]+[1]*cosh(x)/(1+cosh(x)*[2])",import_kFSR_vs_Abseta_histo_res1_Dir1->GetXaxis()->GetXmin(),import_kFSR_vs_Abseta_histo_res1_Dir1->GetXaxis()->GetXmax()); //was used before...
    import_kFSR_vs_Abseta_histo_res1_Dir1->GetFunction("kFSR_fit")->GetParameters();
    kFSR_fit->SetParameters(import_kFSR_vs_Abseta_histo_res1_Dir1->GetFunction("kFSR_fit")->GetParameters());
    import_kFSR_vs_Abseta_histo_res1_Dir1->GetFunction("kFSR_fit")->SetLineColor(1);
    kFSR_fit->SetLineColor(2);
    kFSR_fit->SetLineWidth(2);
    TFitResultPtr r = import_kFSR_vs_Abseta_histo_res1_Dir1->Fit(kFSR_fit,"S E");
    TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
    chisquared1= r->Chi2();
    ndf1= r->Ndf();

    TGraphErrors* band1= get_error_band(cov, kFSR_fit);
    band1->SetFillColor(15);
    //    band1->SetFillStyle(3001);
    
    band1->Draw("same 3");
    import_kFSR_vs_Abseta_histo_res1_Dir1->Draw("e same");
    cmsPrel(intLumi=36, false);

    TLine *line = new TLine(import_kFSR_vs_Abseta_histo_res1_Dir1->GetXaxis()->GetXmin(),1.,import_kFSR_vs_Abseta_histo_res1_Dir1->GetXaxis()->GetXmax(),1.);
    line->Draw();
    line->SetLineStyle(2);
    line->SetLineColor(1);
    sprintf (buffer, "Fit uncertainty (#chi^{2}/ndf = %.1f / %d )", chisquared1, ndf1);
    cout << buffer << endl;
    TLegend *leg_kFSR1;
    leg_kFSR1 = new TLegend(0.2,0.70,0.65,0.85);
    leg_kFSR1->SetFillColor(kWhite);
    leg_kFSR1->SetTextFont(42);
      //   leg->SetHeader("Legende");
    leg_kFSR1->AddEntry(kFSR_fit,"p_{0}+ #frac{p_{1}cosh(|#eta|)}{1+p_{2}cosh(|#eta|)}","l");
    leg_kFSR1->AddEntry(band1,buffer,"f");
    //    leg_kFSR->AddEntry(import_kFSR_vs_Abseta_histo_res1_Dir1,"Data","lep");


    leg_kFSR1->Draw();
    //"kFSR_comp_"+
    cmsPrel(intLumi=36, false);
    c->SaveAs(GetDateDir()+"/kFSR_comp_"+dir1+"_kFSR_"+ algo +"_em.eps");

    import_kFSR_vs_Abseta_histo_res1_Dir2->Draw("e");
    import_kFSR_vs_Abseta_histo_res1_Dir2->GetYaxis()->SetRangeUser(kFSR_extrapol_ymin,kFSR_extrapol_ymax);
    r = import_kFSR_vs_Abseta_histo_res1_Dir2->Fit(kFSR_fit,"S");
    cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
    TGraphErrors* band2= get_error_band(cov, kFSR_fit);
    band2->SetFillColor(15);
    //    band2->SetFillStyle(3001);
    band2->Draw("same 3");
    chisquared2= r->Chi2();
    ndf2= r->Ndf();
    cout << "CHI2: " << chisquared2 << " ndf: " << ndf2 << endl;
    import_kFSR_vs_Abseta_histo_res1_Dir2->Draw("e same");
    sprintf (buffer, "Fit uncertainty (#chi^{2}/ndf = %.1f / %d )", chisquared2, ndf2);
    cout << buffer << endl;
    TLegend *leg_kFSR2;
    leg_kFSR2 = new TLegend(0.2,0.70,0.65,0.85);
    leg_kFSR2->SetFillColor(kWhite);
    leg_kFSR2->SetTextFont(42);
      //   leg->SetHeader("Legende");
    leg_kFSR2->AddEntry(kFSR_fit,"p_{0}+ #frac{p_{1}cosh(|#eta|)}{1+p_{2}cosh(|#eta|)}","l");
    leg_kFSR2->AddEntry(band2,buffer,"f");
    //    leg_kFSR->AddEntry(import_kFSR_vs_Abseta_histo_res1_Dir1,"Data","lep");

    leg_kFSR2->Draw();
    line->Draw();
    cmsPrel(intLumi=36, false);
    c->SaveAs(GetDateDir()+"/kFSR_comp_"+dir2+"_kFSR_"+ algo +"_em.eps");





    TLegend *leg_kFSR_comb;
    leg_kFSR_comb = new TLegend(0.2,0.60,0.65,0.85);
    leg_kFSR_comb->SetFillColor(kWhite);
    leg_kFSR_comb->SetTextFont(42);
      //   leg->SetHeader("Legende");
//    leg_kFSR_comb->AddEntry(import_kFSR_vs_Abseta_histo_res1_Dir1,"p_{0}+ #frac{p_{1}cosh(|#eta|)}{1+p_{2}cosh(|#eta|)} ("+dir1+")","lp");
//    leg_kFSR_comb->AddEntry(band1,"Fit uncertainty ("+dir1+")","f");
//    leg_kFSR_comb->AddEntry(import_kFSR_vs_Abseta_histo_res1_Dir2,"p_{0}+ #frac{p_{1}cosh(|#eta|)}{1+p_{2}cosh(|#eta|)} ("+dir2+")","lp");
//    leg_kFSR_comb->AddEntry(band2,"Fit uncertainty ("+dir2+")","f");
    leg_kFSR_comb->AddEntry(import_kFSR_vs_Abseta_histo_res1_Dir1,dir1,"lp");
    //    leg_kFSR_comb->AddEntry(band1,"Fit uncertainty","f");
    leg_kFSR_comb->AddEntry(import_kFSR_vs_Abseta_histo_res1_Dir2,dir2,"lp");
    //    leg_kFSR_comb->AddEntry(band2,"Fit uncertainty","f");



    import_kFSR_vs_Abseta_histo_res1_Dir2->SetLineColor(2);
    import_kFSR_vs_Abseta_histo_res1_Dir2->SetMarkerColor(2);
    import_kFSR_vs_Abseta_histo_res1_Dir2->SetMarkerStyle(21);
    import_kFSR_vs_Abseta_histo_res1_Dir1->GetFunction("kFSR_fit")->SetLineColor(1);
    import_kFSR_vs_Abseta_histo_res1_Dir2->GetFunction("kFSR_fit")->SetLineColor(2);
    //    band2->SetLineColor();


//
//    import_kFSR_vs_Abseta_histo_res1_Dir1->Draw("e");
//    band1->Draw("same 3");
//    import_kFSR_vs_Abseta_histo_res1_Dir2->Draw("same e");
//    band2->SetFillColor(42);
//    band2->SetFillStyle(3001);
//    band2->Draw("same 3");
//
//    import_kFSR_vs_Abseta_histo_res1_Dir1->Draw("e same");
//    import_kFSR_vs_Abseta_histo_res1_Dir2->Draw("e same");
//
//    leg_kFSR_comb->Draw();
//    cmsPrel(intLumi=36, false);
//

    import_kFSR_vs_Abseta_histo_res1_Dir1->Draw("hist e");
    //    band1->Draw("same 3");
    import_kFSR_vs_Abseta_histo_res1_Dir2->Draw("hist same e");
    band2->SetFillColor(42);
    band2->SetFillStyle(3001);
    //    band2->Draw("same 3");

    import_kFSR_vs_Abseta_histo_res1_Dir1->Draw("hist e same");
    import_kFSR_vs_Abseta_histo_res1_Dir2->Draw("hist e same");

    leg_kFSR_comb->Draw();
    cmsPrel(intLumi=36, false);
    latex.DrawLatex(.25,.9,algo+" Jets");


    c->SaveAs(GetDateDir()+"/kFSR_comp_"+dir1+"_"+dir2+"_kFSR_"+ algo +"_em.eps");

}


TGraphErrors* get_error_band(TMatrixDSym cov, TF1 *kFSR_fit){

  //define array for x-values and fill it
  Double_t xmin= kFSR_fit->GetXaxis()->GetXmin();
  Double_t xmax= kFSR_fit->GetXaxis()->GetXmax();
  
  Int_t granularity =100;
  Double_t step = (xmax-xmin)/granularity;
  Double_t x_values[100];
  Double_t ex_values[100];
  Double_t y_values[100];
  Double_t ey_values[100];

  for(Int_t i=0;i<granularity;i++){
  //calc x value
    x_values[i]=xmin+i*step;
    ex_values[i]=0;
  //calc y value
    y_values[i]=kFSR_fit->Eval(x_values[i]);
    ey_values[i]=0.001;
  
    Double_t par_a=kFSR_fit->GetParameter(0);
    Double_t par_b=kFSR_fit->GetParameter(1);
    Double_t par_c=kFSR_fit->GetParameter(2);


  //calc partial derivatives and fill into vectors
    Double_t dqdz[3];
    dqdz[0]=1;
    dqdz[1]=cosh(x_values[i])/(1+par_c*cosh(x_values[i]));
    dqdz[2]=(-par_b*cosh(x_values[i])*cosh(x_values[i]))/TMath::Power((1+par_c*cosh(x_values[i])),2);


    TMatrixD result(1,1);
    TMatrixD part_left(1,3);
    TMatrixD part_right(3,1);
    for(Int_t j=0;j<3;j++){
      part_left[0][j]=dqdz[j];
      part_right[j][0]=dqdz[j];
    }
    
    
    //    part_left.Print();
    //    part_right.Print();
    //    cov.Print();
    
    part_right=cov*part_right;
    result= part_left*part_right;
    Double_t ey = TMath::Sqrt(result[0][0]);
    //    cout << "Error in y: " << ey << endl;
    //multiply v*m*v
    
    ey_values[i]=ey;
    
  }
  

  //build tgrapherror and return
  TGraphErrors* done= new	TGraphErrors(granularity, x_values, y_values, ex_values ,  ey_values);
  done->SetFillColor(15); 
  return done;


}
