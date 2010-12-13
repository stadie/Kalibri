//
//    Class for jet bins
//
//    first version: Hartmut Stadie 2010/05/10
//    $Id: JetBin.cc,v 1.5 2010/11/01 15:47:43 stadie Exp $
//   


#include "JetBin.h"
#include "Jet.h"
#include "CorFactors.h"

void JetBin::addJet(float Et, float EmEt, float HadEt ,float OutEt, float E,
		    float eta,float phi, float phiphi, float etaeta, 
		    float genPt, float dR, const CorFactors& corFactors)
{
  //std::cout << "jet added for par eta bin:" << f_.parIndex() << std::endl;
  sumMess_.pt += Et;
  sumMess_.EMF += EmEt;
  sumMess_.HadF += HadEt;
  sumMess_.OutF += OutEt;
  sumMess_.E += E;
  sumMess_.eta += eta;
  //sumMess_.phi += phi;
  sumMess_.phiphi += phiphi;
  sumMess_.etaeta += etaeta;
  sumPt2_ += Et * Et;

  sumGenPt_ += genPt;
  sumGenPt2_ += genPt * genPt;
  sumdR_ += dR;
  
  sumL1_ += corFactors.getL1() * Et;
  sumL2_ += corFactors.getL2() * Et;
  sumL3_ += corFactors.getL3() * corFactors.getL2() * Et;
  sumLres_ += corFactors.getLRes() * corFactors.getL3() * corFactors.getL2() * Et;
  sumL4_ += corFactors.getL4() * corFactors.getLRes() * corFactors.getL3() * corFactors.getL2() * Et;
  sumL5_ += corFactors.getL5() * corFactors.getL4() * corFactors.getLRes() * corFactors.getL3() * corFactors.getL2() * Et;
  sumJPT_ += corFactors.getJPT() * Et;
  sumJPTL2L3_ += corFactors.getJPTL2L3() * Et;

  ++njets_;
}

Jet* JetBin::createJet() const {
  if(! njets_) return 0;
  float w = 1.0/njets_;
  //std::cout << "Jet:" << sumMess_.pt*w << " , " << sumGenPt_ * w 
  //	    << ", " << sumL3_ * w << std::endl;
  Jet* j = new Jet(sumMess_.pt*w, sumMess_.EMF*w, sumMess_.HadF*w, 
		   sumMess_.OutF*w, sumMess_.E*w,
		   sumMess_.eta*w, sumMess_.phi*w, sumMess_.phiphi*w, 
		   sumMess_.etaeta*w, Jet::unknown,
		   sumGenPt_*w, sumdR_*w, 
		   new CorFactors(sumL1_*w,sumL2_/sumMess_.pt,sumL3_/sumL2_,sumLres_/sumL3_,
				  sumL4_/sumLres_,sumL5_/sumL4_,sumJPT_/sumMess_.pt,
				  sumJPTL2L3_/sumMess_.pt),
		   *f_,errf_,*gf_);
  float err = sqrt(w*sumPt2_ - j->pt() * j->pt() + w*sumGenPt2_ - j->genPt() * j->genPt());
  j->setError(err);
  //std::cout << j->pt() << ":" << j->error() << " ; " << err << '\n';
  return j;
}
