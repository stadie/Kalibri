// $Id: TwoJetsPtBalanceEvent.cc,v 1.12 2012/09/10 15:44:05 kirschen Exp $

#include "TwoJetsPtBalanceEvent.h"
#include "TVector2.h"

#include <iomanip>
#include <algorithm>

//!  \brief Calculates \f$ \chi^{2} \f$ from pt difference
// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_simple(double * temp_derivative1, 
					       double * temp_derivative2, const double* epsilon) const
{
  // Corrected jet pt
  double pt1 = parametrizedMess();
  double pt2 = parametrizedMess2();

  // Residual
  double res = chi2_fast_simple_res(pt1,pt2);

  // Squared error on residual
  double dRes2 = chi2_fast_simple_dRes2(pt1, pt2);

  // Likelihood
  double chi2 = res * res / dRes2;
  if(chi2 != chi2) {//check for NAN
    std::cout <<pt1 << ", " << pt2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = weight() * Event::scaleResidual(chi2);

  if(!temp_derivative1) return chi2;

  // Derivative calculation
  const Parameters::VariationColl& varColl1 = jet1_->varyParsDirectly(epsilon);
  const Parameters::VariationColl& varColl2 = jet2_->varyParsDirectly(epsilon);

  // Variation of parameters of first jet
  for(Parameters::VariationCollIter i1 = varColl1.begin() ; i1 != varColl1.end() ; ++i1) {
    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = i1->lowerEt;
    double pt2tmp = 0.;
    Parameters::VariationCollIter i2 = find(varColl2.begin(),varColl2.end(),i1->parid);
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->lowerEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = weight() * Event::scaleResidual(temp1);


    // Corrected pt and derivative in case of 
    // upper parameter variation
    pt1tmp = i1->upperEt;
    pt2tmp = 0.;
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->upperEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = weight() * Event::scaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of first jet

  // Variation of parameters of second jet
  for(Parameters::VariationCollIter i2 = varColl2.begin() ; i2 != varColl2.end() ; ++i2) {
    Parameters::VariationCollIter i1 = find(varColl1.begin(),varColl1.end(),i2->parid);
    if(i1 != varColl1.end()) continue; // Parameter already covered by jet 1 par variation

    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = pt1;
    double pt2tmp = i2->lowerEt;

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = weight() * Event::scaleResidual(temp1);


    // Corrected pt in case of 
    // upper parameter variation
    pt2tmp = i2->upperEt;

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_simple_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_simple_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = weight() * Event::scaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of second jet

  return chi2;
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_simple_res(double pt1, double pt2) const {
  return pt1 - pt2;
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_simple_dRes2(double pt1, double pt2) const {
  // Derivativ d(corr pt) / d(pt); for simplification
  // assuming constant correction function
  double dPt1 = pt1 / jet1_->Et();
  double dPt2 = pt2 / jet2_->Et();

  double dRes2 = dPt1 * dPt1 * error1_ * error1_;
  dRes2       += dPt2 * dPt2 * error2_ * error2_;
  
  return dRes2;
}



//!  \brief Calculates \f$ \chi^{2} \f$ from pt balance
// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_balance(double * temp_derivative1, 
						double * temp_derivative2, const double* epsilon) const
{
  // Corrected jet pt
  double pt1 = parametrizedMess();
  double pt2 = parametrizedMess2();

  // Residual
  double res = chi2_fast_balance_res(pt1,pt2);

  // Squared error on residual
  double dRes2 = chi2_fast_balance_dRes2(pt1, pt2);

  // Likelihood
  double chi2 = res * res / dRes2;
  if(chi2 != chi2) {//check for NAN
    std::cout <<pt1 << ", " << pt2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  chi2 = weight() * Event::scaleResidual(chi2);

  if(!temp_derivative1) return chi2;

  // Derivative calculation
  const Parameters::VariationColl& varColl1 = jet1_->varyParsDirectly(epsilon);
  const Parameters::VariationColl& varColl2 = jet2_->varyParsDirectly(epsilon);

  // Variation of parameters of first jet
  for(Parameters::VariationCollIter i1 = varColl1.begin() ; i1 != varColl1.end() ; ++i1) {
    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = i1->lowerEt;
    double pt2tmp = 0.;
    Parameters::VariationCollIter i2 = find(varColl2.begin(),varColl2.end(),i1->parid);
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->lowerEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = weight() * Event::scaleResidual(temp1);


    // Corrected pt and derivative in case of 
    // upper parameter variation
    pt1tmp = i1->upperEt;
    pt2tmp = 0.;
    if(i2 != varColl2.end()) { // Parameter also covered by jet 2
      assert(i1->parid == i2->parid);
      pt2tmp = i2->upperEt;
    } else {
      pt2tmp = pt2;
    } 

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = weight() * Event::scaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i1->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i1->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of first jet

  // Variation of parameters of second jet
  for(Parameters::VariationCollIter i2 = varColl2.begin() ; i2 != varColl2.end() ; ++i2) {
    Parameters::VariationCollIter i1 = find(varColl1.begin(),varColl1.end(),i2->parid);
    if(i1 != varColl1.end()) continue; // Parameter already covered by jet 1 par variation

    // Corrected pt in case of 
    // lower parameter variation
    double pt1tmp = pt1;
    double pt2tmp = i2->lowerEt;

    // Likelihood for lower parameter variation
    double temp1 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp1 = temp1 * temp1 / dRes2;
    temp1 = weight() * Event::scaleResidual(temp1);


    // Corrected pt in case of 
    // upper parameter variation
    pt2tmp = i2->upperEt;

    // Likelihood for upper parameter variation
    double temp2 = chi2_fast_balance_res(pt1tmp,pt2tmp);
    dRes2 = chi2_fast_balance_dRes2(pt1tmp, pt2tmp);
    temp2 = temp2 * temp2 / dRes2;
    temp2 = weight() * Event::scaleResidual(temp2);

    // Contribution to global derivative
    temp_derivative1[i2->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i2->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }  // End variation of parameters of second jet

  return chi2;
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_balance_res(double pt1, double pt2) const {
  return 2. * (pt1 - pt2) / (pt1 + pt2);
}



// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_fast_balance_dRes2(double pt1, double pt2) const {
  double ptDijet = (pt1 + pt2) / 2.;
  double ptBal = (pt1 - pt2) / (pt1 + pt2);

  // Derivativ d(corr pt) / d(pt); for simplification
  // assuming constant correction function
  double dPt1 = pt1 / jet1_->Et();
  double dPt2 = pt2 / jet2_->Et();

  double dRes1 = (1. - ptBal) * dPt1 * error1_;
  double dRes2 = (1. + ptBal) * dPt2 * error2_;
  
  return ( dRes1*dRes1 + dRes2*dRes2 ) / ptDijet / ptDijet;
}



//!  \brief Return absolute uncorrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbs(double pt1, double pt2) const {
  double x = pt1 * cos(jet1_->phi()) + pt2 * cos(jet2_->phi());
  double y = pt1 * sin(jet1_->phi()) + pt2 * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute uncorrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbs() const {
  double x = jet1_->Et() * cos(jet1_->phi()) + jet2_->Et() * cos(jet2_->phi());
  double y = jet1_->Et() * sin(jet1_->phi()) + jet2_->Et() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute generated pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbsGen() const {
  double x = jet1_->genPt() * cos(jet1_->phi()) + jet2_->genPt() * cos(jet2_->phi());
  double y = jet1_->genPt() * sin(jet1_->phi()) + jet2_->genPt() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute corrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbsCorr() const {
  double x = jet1_->correctedEt() * cos(jet1_->phi()) + jet2_->correctedEt() * cos(jet2_->phi());
  double y = jet1_->correctedEt() * sin(jet1_->phi()) + jet2_->correctedEt() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}



//!  \brief Return absolute L2L3 corrected pt sum
//!         \f$ |\vec{p}^{1}_{T} + \vec{p}^{2}_{T}| \f$
// --------------------------------------------------
double TwoJetsPtBalanceEvent::ptSumAbsCorrL2L3() const {
  double x = jet1_->corFactors().getL2L3() * jet1_->Et() * cos(jet1_->phi())
    + jet2_->corFactors().getL2L3() * jet2_->Et() * cos(jet2_->phi());
  double y = jet1_->corFactors().getL2L3() * jet1_->Et() * sin(jet1_->phi())
    + jet2_->corFactors().getL2L3() * jet2_->Et() * sin(jet2_->phi());

  return sqrt( x*x + y*y );
}


//!  \brief Calculates \f$ \chi^{2} \f$ from pt difference
// --------------------------------------------------
double TwoJetsPtBalanceEvent::chi2_relative(double * temp_derivative1, 
					    double * temp_derivative2, 
					    const double* epsilon) const
{
  // first jet: probe, 2nd jet: tag
  // Corrected jet pt
  double pt1 = jet1_->Et();
  double pt2 = jet2_->Et();
  double ptave = 0.5*(pt1 + pt2);
  
  double var1 = jet1_->expectedError(ptave);
  var1 *= var1;
  double var2 = jet2_->expectedError(ptave);
  var2 *= var2;
  /*
  ptave = (pt1/var1 + pt2/var2) /(1/var1+1/var2);
  var1 = jet1_->expectedError(ptave);
  var1 *= var1;
  var2 = jet2_->expectedError(ptave);
  var2 *= var2;
  ptave = (pt1/var1 + pt2/var2) /(1/var1+1/var2);
  */
  double L = jet1_->correctedEt(ptave)/ptave; 
  // Residual
  double res = L * pt1 - ptave - residual_;
  
  // Squared error on residual 
  const double deltaE = 1e-7 * jet1_->Et();
  double etprime  = (jet1_->correctedEt(ptave + deltaE) - 
  		     jet1_->correctedEt(ptave - deltaE))/2/deltaE;
  //var1 = etprime * jet1_->expectedError(ptave);
  double dr1 = 0.5 * ( (etprime -L)/ptave * pt1  + 2 * L  - 1);
  double dr2 = 0.5 * ( (etprime -L)/ptave * pt1 - 1);

  double var = dr1 * dr1 * var1 + dr2 * dr2 * var2 + varresidual_;
  // Likelihood
  double chi2 = res * res / var;

  if(chi2 != chi2) {//check for NAN
    std::cout << "TwoJetsPtBalanceEvent::chi2_relative: " << pt1 << ", " << pt2 << ", " <<  jet1_->Et() << ", " << jet2_->Et() << ", " << chi2 << '\n';
  }
  double scale = Event::scaleResidual(chi2)/chi2;
  chi2 = weight() * scale * (log(var) + chi2);

  if(!temp_derivative1) return chi2;

  // Derivative calculation
  const Parameters::VariationColl& varColl1 = jet1_->varyParsDirectly(epsilon,true,ptave);

  // Variation of parameters of first jet
  for(Parameters::VariationCollIter i = varColl1.begin() ; i != varColl1.end() ; ++i) {
    // Likelihood for lower parameter variation
    L = i->lowerEt/ptave;
    etprime = i->lowerEtDeriv;
    // Residual
    res = L * pt1 - ptave - residual_;
    // Squared error on residual 
    dr1 = 0.5 * ( (etprime -L)/ptave * pt1  + 2 * L  - 1);
    dr2 = 0.5 * ( (etprime -L)/ptave * pt1 - 1);

    var = dr1 * dr1 * var1 + dr2 * dr2 * var2 + varresidual_;
    double temp1 = res * res / var; 
    scale = Event::scaleResidual(temp1)/temp1;
    temp1 = weight() * scale * (log(var) + temp1);

    // Likelihood for upper parameter variation   
    L = i->upperEt/ptave;
    etprime = i->upperEtDeriv;
    // Residual
    res = L * pt1 - ptave - residual_;
    // Squared error on residual 
    dr1 = 0.5 * ( (etprime -L)/ptave * pt1  + 2 * L  - 1);
    dr2 = 0.5 * ( (etprime -L)/ptave * pt1 - 1);

    var = dr1 * dr1 * var1 + dr2 * dr2 * var2 + varresidual_;
    double temp2 = res * res / var;
    scale = Event::scaleResidual(temp2)/temp2;
    temp2 = weight() * scale * (log(var) + temp2);

    //std::cout << temp1 << ", " << temp2 << '\n';
    // Contribution to global derivative
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  } 
  return chi2;
}


//!  \brief Define low thresholds for Jet1 Pt used to calculate relPtJet3*
//!
// --------------------------------------------------
const float TwoJetsPtBalanceEvent::lowThresholdJet1_=10;

//!  \brief Define low thresholds for Jet3 Pt used to calculate relPtJet3*
//!
// --------------------------------------------------
const float TwoJetsPtBalanceEvent::lowThresholdJet3_=5;


//!  \brief Returns fraction of momentum of third jet and average dijet momentum  
//!
// --------------------------------------------------
double TwoJetsPtBalanceEvent::relPtJet3() const {
  if (!  hasJet3()) return 0;
  if( getJet1()->pt() < lowThresholdJet1_) return 0;
  double pJ3 =  getJet3()->pt();
  if(pJ3 < lowThresholdJet3_) pJ3=lowThresholdJet3_;
  return pJ3/ ptDijet();
}

//!  \brief Returns fraction of momentum of third jet and average dijet momentum (L2L3Corr)
//!
// --------------------------------------------------
double TwoJetsPtBalanceEvent::relPtJet3CorrL2L3() const {
  if (!  hasJet3()) return 0;
  if( getJet1()->pt() < lowThresholdJet1_) return 0;
  double pJ3 =  getJet3()->corFactors().getL2L3() *  getJet3()->pt();
  if(pJ3 < lowThresholdJet3_) pJ3=lowThresholdJet3_;
  return pJ3/ ptDijetCorrL2L3();
}



//!  \brief Returns fraction of momentum of third jet projected to dijet axis 
//!
// --------------------------------------------------
double TwoJetsPtBalanceEvent::relPtJet3Projection() const {
  if (! hasJet3()) return 0;
  if(getJet1()->pt() < lowThresholdJet1_) return 0;
  // Phi of dijet axis
  double pPhi = TVector2::Phi_mpi_pi(0.5*(getJet1()->phi() + getJet2()->phi()) + M_PI/2.);
  double pJ3 =  getJet3()->pt() * cos(TVector2::Phi_mpi_pi(pPhi-getJet3()->phi()));
  if(pJ3 < lowThresholdJet3_) pJ3=lowThresholdJet3_;
  return pJ3/ptDijet();
}



//!  \brief Returns fraction of momentum of third jet projected to dijet axis (L2L3Corr)
//!
// --------------------------------------------------
double TwoJetsPtBalanceEvent::relPtJet3ProjectionCorrL2L3() const {
  if (! hasJet3()) return 0;
  if(getJet1()->pt() < lowThresholdJet1_) return 0;
  // Phi of dijet axis
  double pPhi = TVector2::Phi_mpi_pi(0.5*(getJet1()->phi() + getJet2()->phi()) + M_PI/2.);
  double pJ3 = getJet3()->corFactors().getL2L3() * getJet3()->pt() * cos(TVector2::Phi_mpi_pi(pPhi-getJet3()->phi()));
  if(pJ3 < lowThresholdJet3_) pJ3=lowThresholdJet3_;
  return pJ3/ptDijetCorrL2L3();
}

//!  \brief Returns L2L3 corrected trigger pt-variable
//!
// --------------------------------------------------
double TwoJetsPtBalanceEvent::triggerPtVariableL2L3(bool useSingleJetTriggers){
 
  if(!useSingleJetTriggers)return ptDijetCorrL2L3();
  else return jet2_->corFactors().getL2L3() * jet2_->pt();


}
