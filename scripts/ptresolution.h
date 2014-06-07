// pt_resolution file from Mikko

#ifndef __ptresolution_h__
#define __ptresolution_h__


#include "TMath.h"

namespace ptResolutionForSmearing {

// Hauke's resolutions
int _nres = 8;

float Matthias_eta_bins [] =  {0, 0.5, 1.1, 1.7, 2.3, 5.2};
float JME_eta_bins [] =  {0, 1.1, 1.7, 2.3, 5.2}; //jinst
float Kristin_eta_bins [] = {0,0.5,1.1,1.7,2.3,2.8,3.2,5.2};



float *pointerToEtaBins;
float *pointerToScaleFactors;

int find_ieta(float eta){
  int ieta=-1;
  for(int i=0;i<_nres-1;i++){
    if(std::abs(eta)>=pointerToEtaBins[i]&&std::abs(eta)<pointerToEtaBins[i+1])ieta=i;
  }
  return ieta;
}

//// k-scale for PFJets from Matthias, July 14, Hamburg Calib
float _pscale_ak5pf_M[] = {1.052, 1.057, 1.096, 1.134, 1.288};

// downward fluctuation values
float _pscale_ak5pf_M_d[] =
{
1.052-TMath::Sqrt(pow(0.012,2)+pow(0.061,2)),
1.057-TMath::Sqrt(pow(0.012,2)+pow(0.055,2)), 
1.096-TMath::Sqrt(pow(0.017,2)+pow(0.062,2)), 
1.134-TMath::Sqrt(pow(0.035,2)+pow(0.085,2)), 
1.288-TMath::Sqrt(pow(0.127,2)+pow(0.153,2))};

// upward fluctuation values
float _pscale_ak5pf_M_u[] =
{
1.052+TMath::Sqrt(pow(0.012,2)+pow(0.061,2)),
1.057+TMath::Sqrt(pow(0.012,2)+pow(0.055,2)), 
1.096+TMath::Sqrt(pow(0.017,2)+pow(0.062,2)), 
1.134+TMath::Sqrt(pow(0.035,2)+pow(0.085,2)), 
1.288+TMath::Sqrt(pow(0.127,2)+pow(0.153,2))};

//k-scale for PFjets from Kristin, April 14 2014
float _pscale_ak5pf_K[] = {1.077, 1.100, 1.119, 1.205, 1.191};

//downward variation values
float _pscale_ak5pf_K_d[] = {1.077-0.026, 1.100-0.028, 1.119-0.029, 1.205-0.045, 1.191-0.079};

//upward variation values
float _pscale_ak5pf_K_u[] = {1.077+0.026, 1.100+0.028, 1.119+0.029, 1.205+0.045, 1.191+0.079};

//k-scale for PFjets from Kristin with forward extension, June 03 2014
float _pscale_ak5pf_Kfe[] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};

//downward variation values
float _pscale_ak5pf_Kfe_d[] = {1.079-0.026, 1.099-0.028, 1.121-0.029, 1.208-0.046, 1.254-0.062, 1.395-0.063, 1.056-0.191};

//upward variation values
float _pscale_ak5pf_Kfe_u[] = {1.079+0.026, 1.099+0.028, 1.121+0.029, 1.208+0.046, 1.254+0.062, 1.395+0.063, 1.056+0.191};


  //according to JINST paper/JER twiki
float _pscale_ak5calo_M_JME[] ={ 1.088, 1.139, 1.082, 1.065};//Calo
float _pscale_ak5jpt_M_JME[] ={ 1.087, 1.213, 1.018, 1.068};//JPT


void configureSmearfactor(TString chooseScaleFactors) {
  pointerToScaleFactors=0;
  pointerToEtaBins=0;

  //choose eta-binning
  if(chooseScaleFactors.Contains("Matthias")){
    pointerToEtaBins=Matthias_eta_bins;
    _nres=6;
  }
  else if(chooseScaleFactors.Contains("KristinFE")){
    pointerToEtaBins=Kristin_eta_bins;
    _nres=8;
  }
  else if(chooseScaleFactors.Contains("Kristin")){
    pointerToEtaBins=Matthias_eta_bins;
    _nres=6;
  }
  else if(chooseScaleFactors.Contains("JME")){
    pointerToEtaBins=JME_eta_bins;
    _nres=5; 
  }
  std::cout<<_nres<<std::endl;
  //choose scale-factors
  if(chooseScaleFactors=="PF_Matthias")pointerToScaleFactors=_pscale_ak5pf_M;
  else if(chooseScaleFactors=="PF_Matthias_u")pointerToScaleFactors=_pscale_ak5pf_M_u;
  else if(chooseScaleFactors=="PF_Matthias_d")pointerToScaleFactors=_pscale_ak5pf_M_d;

  else if(chooseScaleFactors=="PF_Kristin")pointerToScaleFactors=_pscale_ak5pf_K;
  else if(chooseScaleFactors=="PF_Kristin_u")pointerToScaleFactors=_pscale_ak5pf_K_u;
  else if(chooseScaleFactors=="PF_Kristin_d")pointerToScaleFactors=_pscale_ak5pf_K_d;
  
  else if(chooseScaleFactors=="PF_KristinFE")pointerToScaleFactors=_pscale_ak5pf_Kfe;
  else if(chooseScaleFactors=="PF_KristinFE_u")pointerToScaleFactors=_pscale_ak5pf_Kfe_u;
  else if(chooseScaleFactors=="PF_KristinFE_d")pointerToScaleFactors=_pscale_ak5pf_Kfe_d;

  else if(chooseScaleFactors=="Calo_JME")pointerToScaleFactors=_pscale_ak5calo_M_JME;
  else if(chooseScaleFactors=="JPT_JME")pointerToScaleFactors=_pscale_ak5calo_M_JME;
  else std::cout << "scale factors not properly defined... this will definitely fail" << std::endl;

  if(pointerToEtaBins==0 || pointerToScaleFactors ==0) std::cout << "Pointers should not be null-pointers by now. This will fail!" << std::endl;

  std::cout << "!!!!!Using " << _nres-1 << " bins configured via " << chooseScaleFactors << " for smearing."<< std::endl;
  std::cout << "Scalefactors: "<< std::endl;
  for(int i =0;i<_nres-1;i++)std::cout << pointerToEtaBins[i] <<"<"<<pointerToEtaBins[i+1] <<": \t" << pointerToScaleFactors[i] <<std::endl;
}

double getScaleFactor(float eta){
  int ieta = find_ieta(eta);
  return pointerToScaleFactors[ieta];
  
}


}

#endif // __ptresolution_h__
