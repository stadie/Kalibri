#include "TwoJetsInvMassEvent.h"

#include "TLorentzVector.h"
double TwoJetsInvMassEvent::chi2() const
{
  double et1 = jet1->correctedEt(jet1->Et());
  double c1 = et1/ jet1->Et();
   
  double et2 = jet2->correctedEt(jet2->Et());
  double c2 = et2/ jet2->Et();
  

  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1->eta(),jet1->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2->eta(),jet2->phi(),0);
  
  double err2inv = c1 * jet1->Error();
  err2inv *= err2inv;
  double err2 = c2 * jet2->Error();
  err2inv += err2 * err2;
  err2inv = 1/err2inv;
  
  double chi2 = truth - (p1+p2).M();
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  return chi2;
}
 
double  TwoJetsInvMassEvent::correctedMass() const {
  double et1 = jet1->correctedEt(jet1->Et());
  double et2 = jet2->correctedEt(jet2->Et());
  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1->eta(),jet1->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2->eta(),jet2->phi(),0);
  return (p1+p2).M();
}


double TwoJetsInvMassEvent::chi2_fast_simple(double * temp_derivative1, 
					     double * temp_derivative2, double const epsilon) const
{
  double et1 = jet1->correctedEt(jet1->Et());
  double c1 = et1/ jet1->Et();
   
  double et2 = jet2->correctedEt(jet2->Et());
  double c2 = et2/ jet2->Et();
  

  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1->eta(),jet1->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2->eta(),jet2->phi(),0);
  
  double err2inv = c1 * jet1->expectedError(et1);
  err2inv *= err2inv;
  double err2 = c2 * jet2->expectedError(et2);
  err2inv += err2 * err2;
  err2inv = 1/err2inv;
  
  double chi2 = truth - (p1+p2).M();
  chi2 *= chi2 * err2inv; 
  if(chi2 != chi2) {//check for NAN
    std::cout <<et1 << ", " << et2 << ", " <<  jet1->Et() << ", " << jet2->Et() << ", " << chi2 << '\n';
  }
  chi2 = weight * TData::ScaleResidual(chi2);
  double temp1,temp2;
  const Jet::VariationColl& varcoll1 = jet1->varyParsDirectly(epsilon);
  const Jet::VariationColl& varcoll2 = jet2->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i1 = varcoll1.begin() ; i1 != varcoll1.end() ; ++i1) {
    p1.SetPtEtaPhiM(i1->lowerEt,jet1->eta(),jet1->phi(),0);
    c1 = i1->lowerEt/jet1->Et();
    Jet::VariationCollIter i2 = find(varcoll2.begin(),varcoll2.end(),i1->parid);
    if(i2 != varcoll2.end()) {
      p2.SetPtEtaPhiM(i2->lowerEt,jet2->eta(),jet2->phi(),0);
      c2 = i2->lowerEt/jet2->Et();
      err2 = c2 * i2->lowerError;
    } else {
      p2.SetPtEtaPhiM(et2,jet2->eta(),jet2->phi(),0);
      c2 = et2/ jet2->Et();
      err2 = c2 * jet2->expectedError(et2);
    }
    temp1 = truth - (p1+p2).M();
    err2inv = c1 * i1->lowerError;
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    p1.SetPtEtaPhiM(i1->upperEt,jet1->eta(),jet1->phi(),0);
    c1 = i1->upperEt/jet1->Et();
    if(i2 != varcoll2.end()) {
      p2.SetPtEtaPhiM(i2->upperEt,jet2->eta(),jet2->phi(),0);
      c2 = i2->upperEt/ jet2->Et();
      err2 = c2 * i2->upperError;
    } 
    temp2 = truth - (p1+p2).M();
    err2inv = c1 * i1->upperError;
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  for(Jet::VariationCollIter i2 = varcoll2.begin() ; i2 != varcoll2.end() ; ++i2) {
    p2.SetPtEtaPhiM(i2->lowerEt,jet2->eta(),jet2->phi(),0);
    c2 = i2->lowerEt/jet2->Et();
    Jet::VariationCollIter i1 = find(varcoll1.begin(),varcoll1.end(),i2->parid);
    if(i1 != varcoll1.end()) {
      continue;
    } else {
      p1.SetPtEtaPhiM(et1,jet1->eta(),jet1->phi(),0);
      c1 = et1/ jet1->Et();
      err2 = c1 * jet1->expectedError(et1);
    }
    temp1 = truth - (p1+p2).M();
    err2inv = c2 * i2->lowerError;
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    p2.SetPtEtaPhiM(i2->upperEt,jet2->eta(),jet2->phi(),0);
    c2 = i2->upperEt/jet2->Et();

    temp2 = truth - (p1+p2).M();
    err2inv = c2 * i2->upperError;
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;

}
