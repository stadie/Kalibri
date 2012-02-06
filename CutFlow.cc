#include "CutFlow.h"
#include "ConfigFile.h"
#include "CalibData.h"
#include "TwoJetsPtBalanceEvent.h"
#include "TString.h"
#include "TVector2.h"

CutFlow::CutFlow(const std::string& configfile){
//  : config_(configFile){
  //CutFlow::CutFlow(){
  nEvents_=0;
  nMaxEtaCut_=0;
  config_ = new ConfigFile(configfile.c_str());


}


void CutFlow::setEvent(Event* event){
  event_=event;
  nEvents_++;
}

void CutFlow::setDiJetEvent(Event* event){
  if((*event).type() != PtBalance) std::cout<< "WARNING: event is no PtBalance-event. Will break..." << std::endl;
  diJetEvent_=dynamic_cast<TwoJetsPtBalanceEvent*>(event);
  nEvents_++;
}


void CutFlow::setAllSuppDiJetCuts(){
  std::string cutPrefix = "SuppDiJet ";
  minDeltaPhi_       = config_->read<double>(cutPrefix+"Min Delta Phi",2.5);
  minRel3rdJetEt_    = config_->read<double>(cutPrefix+"Min cut on relative n+1 Jet Et",0.);
  maxRel3rdJetEt_    = config_->read<double>(cutPrefix+"Max cut on relative n+1 Jet Et",1.);
  nMinDeltaPhi_=0;
  nMinCutOn3rdJet_=0;
  nMaxCutOn3rdJet_=0;
  useMinDeltaPhi_=true;
  useMinRel3rdJetEt_=true;
  useMaxRel3rdJetEt_=true;
}


bool CutFlow::doAllSuppDiJetCuts(){
  bool eventSurvives=true;
  Jet * j1 = diJetEvent_->getJet1();
  Jet * j2 = diJetEvent_->getJet2();
  Jet * j3 = diJetEvent_->getJet3();

  if( std::abs(TVector2::Phi_mpi_pi(j1->phi() - j2->phi())) < minDeltaPhi_ ) {
    nMinDeltaPhi_++;
    eventSurvives=false;
  }
  if(diJetEvent_->hasJet3()  && std::abs(j3->corFactors().getL2L3() * j3->pt()) < minRel3rdJetEt_*diJetEvent_->ptDijetCorrL2L3() ) {
    nMinCutOn3rdJet_++;
    eventSurvives=false;
  }
  if(diJetEvent_->hasJet3()  && std::abs(j3->corFactors().getL2L3() * j3->pt()) > maxRel3rdJetEt_*diJetEvent_->ptDijetCorrL2L3() ) {
    nMaxCutOn3rdJet_++;
    eventSurvives=false;
  }
  return eventSurvives;
}


bool CutFlow::doMaxEtaCut(){
    Jet * j1 = diJetEvent_->getJet1();
    if(std::abs(j1->eta())<1.3) return true;
    else   {
      nMaxEtaCut_++;
      return false;
    }
}


void CutFlow::printCutFlow(){
  std::cout << nEvents_ << " events before additional cuts " << std::endl;
  //  std::cout << "  " << nEvents_-nMaxEtaCut_  << " events after eta-cut: "<< std::endl;
  if(useMinDeltaPhi_){
    std::cout << "  " << (nEvents_-=nMinDeltaPhi_)  << std::flush;
    std::cout << " dijet events with DeltaPhi > " << minDeltaPhi_ << std::endl;
  }
  if(useMinRel3rdJetEt_){
    std::cout << "  " << (nEvents_-=nMinCutOn3rdJet_) << std::flush;
    std::cout << " dijet events with " << minRel3rdJetEt_ << " < pt(jet3) / ptAve \n";
  }
  if(useMaxRel3rdJetEt_){
  std::cout << "  " << (nEvents_-=nMaxCutOn3rdJet_) << std::flush;
  std::cout << " dijet events with pt(jet3) / ptAve < " << maxRel3rdJetEt_ << "\n";
  }
}

CutFlow::~CutFlow()
{
}
