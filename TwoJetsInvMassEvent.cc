//
//    Class for all events with two jets constraint to one invariant mass
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: TwoJetsInvMassEvent.cc,v 1.8 2009/11/24 16:52:59 stadie Exp $
//     
#include "TwoJetsInvMassEvent.h"

#include "TLorentzVector.h"
#include <algorithm>

double TwoJetsInvMassEvent::chi2() const
{
  double chi2 = chi2_fast(0, 0, 0);
  return chi2;
}
 
double  TwoJetsInvMassEvent::correctedMass() const {
  double et1 = jet1_->correctedEt(jet1_->Et());
  double et2 = jet2_->correctedEt(jet2_->Et());
  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
  return (p1+p2).M();
}


double TwoJetsInvMassEvent::chi2_fast_simple(double * temp_derivative1, 
					     double * temp_derivative2, const double* epsilon) const
{
  const double et1 = jet1_->correctedEt(jet1_->Et());
  double c1 = et1/ jet1_->Et();
   
  const double et2 = jet2_->correctedEt(jet2_->Et());
  double c2 = et2/ jet2_->Et();
  

  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
  
  // alpha = 1 - cos(angle(p1,p2))
  const double alpha = 1 - cos(p1.Angle(p2.Vect()));
  double m = sqrt(2 * p1.P()*p2.P() * alpha);
  double dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
  double dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
  
  double err2inv = dmdp1 * c1 * jet1_->Error();
  err2inv *= err2inv;
  double err2 = dmdp2 * c2 * jet2_->Error();
  err2inv += err2 * err2;
  err2inv = 1/err2inv;
 
  double chi2 = truth_ - m;
  chi2 *= chi2 * err2inv; 
  if(chi2 != chi2) {//check for NAN
    std::cout <<et1 << ", " << et2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = weight_ * Event::scaleResidual(-log(err2inv) + chi2);

  if(!temp_derivative1) return chi2;

  double temp1,temp2;
  const Jet::VariationColl& varcoll1 = jet1_->varyParsDirectly(epsilon);
  const Jet::VariationColl& varcoll2 = jet2_->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i1 = varcoll1.begin() ; i1 != varcoll1.end() ; ++i1) {
    p1.SetPtEtaPhiM(i1->lowerEt,jet1_->eta(),jet1_->phi(),0);
    c1 = i1->lowerEt/jet1_->Et();
    Jet::VariationCollIter i2 = find(varcoll2.begin(),varcoll2.end(),i1->parid);
    if(i2 != varcoll2.end()) {
      assert(i1->parid == i2->parid);
      p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
      c2 = i2->lowerEt/jet2_->Et();
      err2 = c2 * jet2_->Error();
    } else {
      p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
      c2 = et2/ jet2_->Et();
      err2 = c2 * jet2_->Error();
    } 
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    err2 *= dmdp2;
    temp1 = truth_ - m;
    err2inv = dmdp1 * c1 * jet1_->Error();
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(-log(err2inv) + temp1);
    p1.SetPtEtaPhiM(i1->upperEt,jet1_->eta(),jet1_->phi(),0);
    c1 = i1->upperEt/jet1_->Et();
    if(i2 != varcoll2.end()) {
      p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
      c2 = i2->upperEt/ jet2_->Et();
      err2 = c2 * jet2_->Error();
    } else {
      err2 = c2 * jet2_->Error();
    }
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    temp2 = truth_ - m;
    err2 *= dmdp2;
    err2inv = dmdp1 * c1 * jet1_->Error();
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(-log(err2inv) + temp2);
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    //if(i2 != varcoll2.end()) {
    //  temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    //  temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    //}
  }
  for(Jet::VariationCollIter i2 = varcoll2.begin() ; i2 != varcoll2.end() ; ++i2) {
    p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
    c2 = i2->lowerEt/jet2_->Et();
    Jet::VariationCollIter i1 = find(varcoll1.begin(),varcoll1.end(),i2->parid);
    if(i1 != varcoll1.end()) {
      continue;
    } else {
      p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
      c1 = et1/ jet1_->Et();
      err2 = c1 * jet1_->Error();
    }
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    temp1 = truth_ - m;
    err2inv = dmdp2 * c2 * jet2_->Error();
    err2inv *= err2inv;
    err2 *= dmdp1;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(-log(err2inv) + temp1);
    p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
    c2 = i2->upperEt/jet2_->Et();
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    temp2 = truth_ - m;
    err2 = dmdp1 * c1 * jet1_->Error();
    err2inv = dmdp2 * c2 * jet2_->Error();
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(-log(err2inv) + temp2);
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;

}

double TwoJetsInvMassEvent::chi2_fast_const_error(double * temp_derivative1, 
						  double * temp_derivative2, const double* epsilon) const
{
  const double et1 = jet1_->correctedEt(jet1_->Et());
  const double et2 = jet2_->correctedEt(jet2_->Et());

  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
  
  // alpha = 1 - cos(angle(p1,p2))
  const double alpha = 1 - cos(p1.Angle(p2.Vect()));
  double m = sqrt(2 * p1.P()*p2.P() * alpha);
  double chi2 = truth_ - m;
  double err2inv = 0.01;
  chi2 *= chi2 * err2inv; 
  if(chi2 != chi2) {//check for NAN
    std::cout <<et1 << ", " << et2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = weight_ * Event::scaleResidual(chi2);

  if(!temp_derivative1) return chi2;

  double temp1,temp2;
  const Jet::VariationColl& varcoll1 = jet1_->varyParsDirectly(epsilon);
  const Jet::VariationColl& varcoll2 = jet2_->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i1 = varcoll1.begin() ; i1 != varcoll1.end() ; ++i1) {
    p1.SetPtEtaPhiM(i1->lowerEt,jet1_->eta(),jet1_->phi(),0);
    Jet::VariationCollIter i2 = find(varcoll2.begin(),varcoll2.end(),i1->parid);
    if(i2 != varcoll2.end()) {
      assert(i1->parid == i2->parid);
      p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
    } else {
      p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
    } 
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    temp1 = truth_ - m;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(temp1);
    p1.SetPtEtaPhiM(i1->upperEt,jet1_->eta(),jet1_->phi(),0);
    if(i2 != varcoll2.end()) {
      p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
    } 
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    temp2 = truth_ - m;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(temp2);
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    //if(i2 != varcoll2.end()) {
    //  temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    //  temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    //}
  }
  for(Jet::VariationCollIter i2 = varcoll2.begin() ; i2 != varcoll2.end() ; ++i2) {
    p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
    Jet::VariationCollIter i1 = find(varcoll1.begin(),varcoll1.end(),i2->parid);
    if(i1 != varcoll1.end()) {
      continue;
    } else {
      p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
    }
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    temp1 = truth_ - m;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(temp1);
    p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    temp2 = truth_ - m;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(temp2);
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;

}

double TwoJetsInvMassEvent::chi2_fast_scaled(double * temp_derivative1, 
					     double * temp_derivative2, const double* epsilon) const
{
  double deltaE = 0.001;
  double et1 = jet1_->correctedEt(jet1_->Et());
  double et1prime = (jet1_->correctedEt(jet1_->Et() + deltaE) - 
		     jet1_->correctedEt(jet1_->Et() - deltaE))/2/deltaE;
  et1prime /= jet1_->Et();
  double et2 = jet2_->correctedEt(jet2_->Et());
  double et2prime =  (jet2_->correctedEt(jet2_->Et() + deltaE) - 
		      jet2_->correctedEt(jet2_->Et() - deltaE))/2/deltaE;
  et2prime /= jet2_->Et();
  double c1 = et1/jet1_->Et() + et1prime;
  double c2 = et2/jet2_->Et() + et2prime;
 
  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
  
  // alpha = 1 - cos(angle(p1,p2))
  double alpha = 1 - cos(p1.Angle(p2.Vect()));
  double m = sqrt(2 * p1.P()*p2.P() * alpha);
  double dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
  double dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
  
  double err2inv = dmdp1 * c1 * jet1_->Error();
  err2inv *= err2inv;
  double err2 = dmdp2 * c2 * jet2_->Error();
  //std::cout << "m:" << m << " = " << (p1+p2).M() << " sigma:" << sqrt(err2inv) << ", " << err2 << '\n';
  err2inv += err2 * err2;
  err2inv = 1/err2inv;
 
  double chi2 = truth_ - m;
  chi2 *= chi2 * err2inv; 
  if(chi2 != chi2) {//check for NAN
    std::cout <<et1 << ", " << et2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = weight_ * Event::scaleResidual(-log(err2inv) + chi2);

  if(!temp_derivative1) return chi2;

  double temp1,temp2;
  const Jet::VariationColl& varcoll1 = jet1_->varyParsDirectly(epsilon);
  const Jet::VariationColl& varcoll2 = jet2_->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i1 = varcoll1.begin() ; i1 != varcoll1.end() ; ++i1) {
    p1.SetPtEtaPhiM(i1->lowerEt,jet1_->eta(),jet1_->phi(),0);
    c1 = (i1->lowerEt + i1->lowerEtDeriv)/jet1_->Et();
    Jet::VariationCollIter i2 = find(varcoll2.begin(),varcoll2.end(),i1->parid);
    if(i2 != varcoll2.end()) {
      p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
      c2 = (i2->lowerEt + i2->lowerEtDeriv)/jet2_->Et();
      err2 = c2 * jet2_->Error();
    } else {
      p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
      c2 = et2/jet2_->Et() + et2prime;
      err2 = c2 * jet2_->Error();
    } 
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    err2 *= dmdp2;
    temp1 = truth_ - m;
    err2inv = dmdp1 * c1 * jet1_->Error();
    //std::cout << "m:" << m << " sigma:" << err2inv << ", " << err2 << '\n';
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(-log(err2inv) + temp1);
    if(std::abs(temp1 - chi2) > 1.0) {
      std::cout << "Et1:" << jet1_->Et() << " c1:" << et1prime << ", " << et1/jet1_->Et() << '\n';
      std::cout << "Et2:" << jet2_->Et() << " c1:" << et2prime << ", " << et2/jet2_->Et() << '\n';
      std::cout << "chi2:" << chi2 << "  temp1:" << temp1 << ", " << c1 << ", " << c2 << '\n';
      assert(std::abs(temp1 - chi2) < 1.0);
    }
    p1.SetPtEtaPhiM(i1->upperEt,jet1_->eta(),jet1_->phi(),0);
    c1 = (i1->upperEt + i1->upperEtDeriv)/jet1_->Et();
    if(i2 != varcoll2.end()) {
      p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
      c2 = (i2->upperEt + i2->upperEtDeriv)/jet2_->Et();
      err2 = c2 * jet2_->Error();
    }  else {
      err2 = c2 * jet2_->Error();
    }
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    temp2 = truth_ - m;
    err2 *= dmdp2;
    err2inv = dmdp1 * c1 * jet1_->Error();
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(-log(err2inv) + temp2); 
    //std::cout << "chi2:" << chi2 << "  temp2:" << temp2 << '\n';
    assert(std::abs(temp2 - chi2) < 1.0);
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  for(Jet::VariationCollIter i2 = varcoll2.begin() ; i2 != varcoll2.end() ; ++i2) {
    p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
    c2 = (i2->lowerEt + i2->lowerEtDeriv) / jet2_->Et();
    Jet::VariationCollIter i1 = find(varcoll1.begin(),varcoll1.end(),i2->parid);
    if(i1 != varcoll1.end()) {
      continue;
    } else {
      p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
      c1 = et1/jet1_->Et() + et1prime;
      err2 = c1 * jet1_->Error();
    }
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    temp1 = truth_ - m;
    err2inv = dmdp2 * c2 * jet2_->Error();
    err2inv *= err2inv;
    err2 *= dmdp1;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(-log(err2inv) + temp1); 
    //std::cout << "chi2:" << chi2 << "  temp1:" << temp1 << '\n';
    assert(std::abs(temp1 - chi2) <1.0);
    p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
    c2 = (i2->upperEt + i2->upperEtDeriv) / jet2_->Et();
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    dmdp1 = p2.P() * alpha / m * jet1_->E() / jet1_->Et();
    dmdp2 = p1.P() * alpha / m * jet2_->E() / jet2_->Et();
    temp2 = truth_ - m;
    err2 = dmdp1 * c1 * jet1_->Error();
    err2inv = dmdp2 * c2 * jet2_->Error();
    err2inv *= err2inv;
    err2inv += err2 * err2;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(-log(err2inv) + temp2);  
    //std::cout << "chi2:" << chi2 << "  temp2:" << temp2 << '\n';
    assert(std::abs(temp2 - chi2) < 1.0);
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}

double TwoJetsInvMassEvent::chi2_fast_inv(double * temp_derivative1, 
					  double * temp_derivative2, const double* epsilon) const
{
  if(flagged_bad_) return chi2plots_;
  const double et1 = jet1_->correctedEt(jet1_->Et());
  const double et2 = jet2_->correctedEt(jet2_->Et());

  TLorentzVector p1,p2;
  p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
  p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
  
  // alpha = 1 - cos(angle(p1,p2))
  double alpha = 1 - cos(p1.Angle(p2.Vect()));
  double m = sqrt(2 * p1.P() * p2.P() * alpha);

  double K = truth_ / m;
  double truth1 = et1 * K; 
  double truth2 = et2 * K; 
  
  double expErr1, expErr2;
  double expEt1 = jet1_->expectedEt(truth1,truth1,expErr1);  
  double expEt2 = jet2_->expectedEt(truth2,truth2,expErr2);
  if((expEt1 < 0) || (expEt2 < 0)) {
    flagged_bad_ = true;
    return chi2plots_;;
  }

  //std::cout << "m:" << m << " ..." <<  sqrt(2 * truth1 *jet1_->E() / jet1_->Et() * truth2 * jet2_->E() / jet2_->Et() * alpha)
  //	    << " == " << truth << '\n';;
  
  double chi21 = expEt1 - jet1_->Et();
  double chi22 = expEt2 - jet2_->Et();
  double err21 = expErr1;
  double err22 = expErr2;
  chi21 *= chi21;
  chi22 *= chi22;
  err21 *= err21;
  err22 *= err22;

  double chi2 = chi21/err21 + chi22/err22 + log(err21) + log(err22);
  if(chi2 != chi2) {//check for NAN
    std::cout <<et1 << ", " << et2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = weight_ * Event::scaleResidual(chi2);

  if(!temp_derivative1) return chi2;

  double temp1,temp2;
  const Jet::VariationColl& varcoll1 = jet1_->varyParsDirectly(epsilon);;
  const Jet::VariationColl& varcoll2 = jet2_->varyParsDirectly(epsilon);;
  for(Jet::VariationCollIter i1 = varcoll1.begin() ; i1 != varcoll1.end() ; ++i1) { 
    p1.SetPtEtaPhiM(i1->lowerEt,jet1_->eta(),jet1_->phi(),0);
    Jet::VariationCollIter i2 = find(varcoll2.begin(),varcoll2.end(),i1->parid);
    if(i2 != varcoll2.end()) {  
      p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
    } else {
      p2.SetPtEtaPhiM(et2,jet2_->eta(),jet2_->phi(),0);
    }  
    m = sqrt(2 * p1.P() * p2.P() * alpha);
    K = truth_ / m;
    truth1 = p1.Pt() * K; 
    truth2 = p2.Pt() * K;
    double oldpar = par_[i1->parid];
    par_[i1->parid] -= epsilon[i1->parid];
    expEt1 = jet1_->expectedEt(truth1,truth1,err21);  
    expEt2 = jet2_->expectedEt(truth2,truth2,err22); 
    if((expEt1 < 0) || (expEt2 < 0)) {
      flagged_bad_ = true;
      par_[i1->parid] = oldpar;
      return chi2plots_;
    }
    chi21 = expEt1 - jet1_->Et();
    chi22 = expEt2 - jet2_->Et();
    chi21 *= chi21;
    chi22 *= chi22;
    err21 *= err21;
    err22 *= err22;
    temp1 = chi21/err21 + chi22/err22 + log(err21) + log(err22);
    temp1 = weight_ * Event::scaleResidual(temp1);
    p1.SetPtEtaPhiM(i1->upperEt,jet1_->eta(),jet1_->phi(),0);
    if(i2 != varcoll2.end()) {
      p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
    } 
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    K = truth_ / m;
    truth1 = p1.Pt() * K; 
    truth2 = p2.Pt() * K;
    par_[i1->parid] = oldpar + epsilon[i1->parid]; 
    expEt1 = jet1_->expectedEt(truth1,truth1,err21);  
    expEt2 = jet2_->expectedEt(truth2,truth2,err22); 
    if((expEt1 < 0) || (expEt2 < 0)) {
      flagged_bad_ = true;
      par_[i1->parid] = oldpar;
      return chi2plots_;
    }
    chi21 = expEt1 - jet1_->Et();
    chi22 = expEt2 - jet2_->Et();
    par_[i1->parid] = oldpar;
    chi21 *= chi21;
    chi22 *= chi22;
    err21 *= err21;
    err22 *= err22;
    temp2 = chi21/err21 + chi22/err22 + log(err21) + log(err22);
    temp2 = weight_ * Event::scaleResidual(temp2);
    //std::cout << temp2 << " - " << temp1 << " : " << chi2 << '\n';
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  for(Jet::VariationCollIter i2 = varcoll2.begin() ; i2 != varcoll2.end() ; ++i2) {
    p2.SetPtEtaPhiM(i2->lowerEt,jet2_->eta(),jet2_->phi(),0);
    Jet::VariationCollIter i1 = find(varcoll1.begin(),varcoll1.end(),i2->parid);
    if(i1 != varcoll1.end()) {
      continue;
    } 
    p1.SetPtEtaPhiM(et1,jet1_->eta(),jet1_->phi(),0);
    m = sqrt(2 * p1.P() * p2.P() * alpha);
    K = truth_ / m;
    truth1 = p1.Pt() * K; 
    truth2 = p2.Pt() * K;
    double oldpar = par_[i2->parid];
    par_[i2->parid] -= epsilon[i2->parid];   
    expEt1 = jet1_->expectedEt(truth1,truth1,err21);  
    expEt2 = jet2_->expectedEt(truth2,truth2,err22); 
    if((expEt1 < 0) || (expEt2 < 0)) {
      flagged_bad_ = true;
      par_[i1->parid] = oldpar;
      return chi2plots_;
    }
    chi21 = expEt1 - jet1_->Et();
    chi22 = expEt2 - jet2_->Et();
    chi21 *= chi21;
    chi22 *= chi22;
    err21 *= err21;
    err22 *= err22;
    temp1 = chi21/err21 + chi22/err22 + log(err21) + log(err22);
    temp1 = weight_ * Event::scaleResidual(temp1);
    p2.SetPtEtaPhiM(i2->upperEt,jet2_->eta(),jet2_->phi(),0);
    m = sqrt(2 * p1.P()*p2.P() * alpha);
    K = truth_ / m;
    truth1 = p1.Pt() * K; 
    truth2 = p2.Pt() * K;
    par_[i2->parid] = oldpar + epsilon[i2->parid]; 
    expEt1 = jet1_->expectedEt(truth1,truth1,err21);  
    expEt2 = jet2_->expectedEt(truth2,truth2,err22); 
    if((expEt1 < 0) || (expEt2 < 0)) {
      flagged_bad_ = true;
      par_[i1->parid] = oldpar;
      return chi2plots_;
    }
    chi21 = expEt1 - jet1_->Et();
    chi22 = expEt2 - jet2_->Et();
    par_[i2->parid] = oldpar;
    chi21 *= chi21;
    chi22 *= chi22;
    err21 *= err21;
    err22 *= err22;
    temp2 = chi21/err21 + chi22/err22 + log(err21) + log(err22);
    temp2 = weight_ * Event::scaleResidual(temp2);
    //std::cout << i2->lowerEt << " - " << i2->upperEt << ", " << et2 << '\n';
    //std::cout << temp2 << " - " << temp1 << " : " << chi2 << '\n';
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;

}
