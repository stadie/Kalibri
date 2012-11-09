// pt_resolution file from Mikko

#ifndef __ptresolution_h__
#define __ptresolution_h__


#include "TMath.h"


// Hauke's resolutions
int _nres = 6;

float Matthias_eta_bins [] =  {0, 0.5, 1.1, 1.7, 2.3, 5.2};
float JME_eta_bins [] =  {0, 1.1, 1.7, 2.3, 5.2}; //jinst



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


  //according to JINST paper/JER twiki
float _pscale_ak5calo_M_JME[] ={ 1.088, 1.139, 1.082, 1.065};//Calo
float _pscale_ak5jpt_M_JME[] ={ 1.087, 1.213, 1.018, 1.068};//JPT


void configureSmearfactor(TString chooseScaleFactors) {
  

  if(chooseScaleFactors.Contains("Matthias")){
    pointerToEtaBins=Matthias_eta_bins;
    _nres=6;
  }
  else if(chooseScaleFactors.Contains("JME")){
    pointerToEtaBins=JME_eta_bins;
    _nres=5;
  }
  if(chooseScaleFactors=="PF_Matthias")pointerToScaleFactors=_pscale_ak5pf_M;
  if(chooseScaleFactors=="PF_Matthias_u")pointerToScaleFactors=_pscale_ak5pf_M_u;
  if(chooseScaleFactors=="PF_Matthias_d")pointerToScaleFactors=_pscale_ak5pf_M_d;

  if(chooseScaleFactors=="Calo_JME")pointerToScaleFactors=_pscale_ak5calo_M_JME;
  if(chooseScaleFactors=="JPT_JME")pointerToScaleFactors=_pscale_ak5calo_M_JME;

  std::cout << "Using " << _nres-1 << " bins configured via " << chooseScaleFactors << " for smearing."<< std::endl;
  std::cout << "Scalefactors: "<< std::endl;
  for(int i =0;i<_nres-1;i++)std::cout << pointerToEtaBins[i] <<"<"<<pointerToEtaBins[i+1] <<": \t" << pointerToScaleFactors[i] <<std::endl;
}

double getScaleFactor(float eta){
  int ieta = find_ieta(eta);
  return pointerToScaleFactors[ieta];
  
}



#endif // __ptresolution_h__
