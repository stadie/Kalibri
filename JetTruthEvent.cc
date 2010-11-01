//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.cc,v 1.23 2010/10/20 11:28:09 stadie Exp $
//   

#include "JetTruthEvent.h"

int JetTruthEvent::nflagged_ = 0;

JetTruthEvent::~JetTruthEvent() 
{ 
  delete jet_;
}

double JetTruthEvent::chi2() const
{
  double diff = (jet_->correctedEt(jet_->Et()) - truth_)/jet_->error();
  return weight_ * Event::scaleResidual(diff*diff);
}

double JetTruthEvent::chi2_fast_blobel(double * temp_derivative1, 
				       double * temp_derivative2, 
				       const double* epsilon) const
{
  double et = jet_->correctedEt(jet_->Et());
  double err2inv = jet_->error();
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth_ - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight_ * Event::scaleResidual(chi2);
  double temp1,temp2;
  const Parameters::VariationColl& varcoll = jet_->varyParsDirectly(epsilon);
  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth_ - i->lowerEt;
    //err2inv = jet_->expectedError(i->lowerEt);;
    //err2inv *= err2inv;
    //err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(temp1);
    temp2 = truth_ - i->upperEt; 
    //err2inv = jet_->expectedError(i->upperEt);
    //err2inv *= err2inv;
    //err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}

double JetTruthEvent::chi2_fast_simple_scaled(double * temp_derivative1, 
					      double * temp_derivative2, 
					      const double* epsilon) const
{
  double et = jet_->correctedEt(jet_->Et());
  double c = et / jet_->Et();
  //if(c <= 0) c = 1.0;
  double err2inv = c * jet_->expectedError(truth_/c);
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth_ - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight_ * Event::scaleResidual(chi2);
  if(chi2 != chi2) {//check for NAN
    std::cout << truth_ << ", " << et << ", " <<  jet_->Et() << ", " << c << ", " << chi2 << '\n';
  }
  double temp1,temp2;
  const Parameters::VariationColl& varcoll = jet_->varyParsDirectly(epsilon);
  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth_ - i->lowerEt;
    c = i->lowerEt / jet_->Et();  
    //if(c <= 0) c = 1.0;
    err2inv = c*jet_->expectedError(truth_/c );
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(temp1);
    assert(temp1 == temp1);
    temp2 = truth_ - i->upperEt;
    c = i->upperEt / jet_->Et();  
    //if(c <= 0) c = 1.0;
    err2inv = c * jet_->expectedError(truth_/c); 
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(temp2);
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative   
  }
  assert(chi2 == chi2);
  return chi2;
}


double JetTruthEvent::chi2_fast_scaled(double * temp_derivative1, 
				       double * temp_derivative2, 
				       const double* epsilon) const
{
  if(flagged_bad_) return 0;
  float et = jet_->correctedEt(jet_->Et());
  float c = et / jet_->Et();
  const float deltaE = 1e-4 * jet_->Et();
  float etprime  = (jet_->correctedEt(jet_->Et() + deltaE) - 
		     jet_->correctedEt(jet_->Et() - deltaE))/2/deltaE;
  if(etprime < 0.1) {
    //std::cout <<"low deriv:" << etprime << " for " << jet_->Et() << " " 
    //	      << jet_->correctedEt(jet_->Et() + deltaE) 
    //	      << " - " << jet_->correctedEt(jet_->Et() - deltaE) << " step:" << deltaE << "\n";
    flagged_bad_ = true;
    ++nflagged_;
    return 0;
  }
  if((etprime == 0) || (c <= 0)) {
    //std::cout << "warning: deriv zero!\n";
    flagged_bad_ = true;
    ++nflagged_;
    return 0;
  }
  double err2 = etprime * jet_->expectedError(jet_->Et());
  err2 *= err2;
  double chi2 = truth_ - et;
  chi2 *= chi2 / err2;
  //chi2 = weight_ * Event::scaleResidual(log(err2) + chi2);
  chi2 = weight_ *  Event::scaleResidual(chi2);
  if(chi2 != chi2) {//check for NAN
    std::cout << truth_ << ", " << et << ", " <<  jet_->Et() << ", " << c << ", " << log(err2) << ", " << etprime << ", " << chi2 << '\n';
    assert(chi2 == chi2);
  }
  assert(! std::isinf(chi2));
  if( chi2 / weight_ > 1e04) {
    std::cout << "too large:"<< truth_ << ", " << et << ", " <<  jet_->Et() << ", " << c << ", " << log(err2) << ", " << etprime << ", " << chi2 << '\n';
  }
  double temp1,temp2;
  const Parameters::VariationColl& varcoll = jet_->varyParsDirectly(epsilon);
  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth_ - i->lowerEt;
    c = i->lowerEt / jet_->Et();
    //if(c <= 0) c = 1.0;
    err2 = i->lowerEtDeriv * jet_->expectedError(jet_->Et());
    if((i->lowerEtDeriv < 0.1)|| (c <= 0)) {
      //std::cout << "warning: deriv too low!\n"; 
      flagged_bad_ = true;
      ++nflagged_;
      return 0;
      //err2 = c * jet_->expectedError(jet_->Et());
    }
    err2 *= err2;
    temp1 *= temp1 / err2; 
    //temp1 = weight * (log(err2) + temp1);
    //temp1 = weight_ * Event::scaleResidual(log(err2) + temp1);
    temp1 = weight_ * Event::scaleResidual(temp1);
    assert(temp1 == temp1);
    assert(! std::isinf(temp1));
    temp2 = truth_ - i->upperEt;
    c = i->upperEt / jet_->Et();  
    //if(c <= 0) c = 1.0;
    err2 = i->upperEtDeriv * jet_->expectedError(jet_->Et());
    if((i->upperEtDeriv < 0.1) || (c <= 0)) {
      //std::cout << "warning: deriv too low!\n";
      flagged_bad_ = true; 
      ++nflagged_;
      return 0;
      //err2 = c * jet_->expectedError(jet_->Et());
    }
    err2 *= err2;
    temp2 *= temp2 / err2;
    //temp2 = weight_ * Event::scaleResidual(log(err2) + temp2);
    temp2 = weight_ * Event::scaleResidual(temp2);
    assert(temp2 == temp2);
    assert(! std::isinf(temp2));
    if( temp1 / weight_ > 1e04) {
    std::cout << "temp1 too large:" << truth_ << ", " << i->lowerEt << ", " <<  jet_->Et() << ", " << c << ", " << log(err2) << ", " << i->lowerEtDeriv << ", " << temp1 << '\n';
    }
    if( temp2 / weight_ > 1e04) {
    std::cout << "temp2 too large:" << truth_ << ", " << i->upperEt << ", " <<  jet_->Et() << ", " << c << ", " << log(err2) << ", " <<  i->upperEtDeriv << ", " << temp2 << '\n';
    }
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 - 2 * chi2 + temp1); // for 2nd derivative   
  }
  assert(chi2 == chi2);
  return chi2;
}

double JetTruthEvent::chi2_fast_simple(double * temp_derivative1, 
				       double * temp_derivative2, 
				       const double* epsilon) const
{
  double et = jet_->correctedEt(jet_->Et());
  double c = et/jet_->Et(); 
  double err2inv = jet_->expectedError(truth_/c);
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth_ - et;
  chi2 *= chi2 * err2inv;  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth_ << ", " << et << ", " <<  jet_->Et() << ", " << c << ", " << chi2 << '\n';
  }
  chi2 = weight_ * Event::scaleResidual(chi2);
  double temp1,temp2;
  const Parameters::VariationColl& varcoll = jet_->varyParsDirectly(epsilon);
  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth_ - i->lowerEt;
    c = i->lowerEt/jet_->Et(); 
    err2inv = jet_->expectedError(truth_/c);
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(temp1);
    temp2 = truth_ - i->upperEt;
    c = i->upperEt/jet_->Et(); 
    err2inv = jet_->expectedError(truth_/c); 
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}

double JetTruthEvent::chi2_fast_invert(double * temp_derivative1, 
				double * temp_derivative2, 
				const double* epsilon) const
{
  if(flagged_bad_) return 0;
  //find expected measurement of jet Et 
  float err2inv;
  double expectedEt = jet_->expectedEt(truth_,truth_,err2inv);
  if(expectedEt < 0) {
    flagged_bad_ = true;
    ++nflagged_;
    //std::cout << "Inversion failed: flag event bad.\n";
    return 0;
    //return chi2_fast_simple_scaled(temp_derivative1,temp_derivative2,epsilon);
  }
  //calculate chi2
  //double err2inv = jet_->expectedError(expectedEt);
  //std::cout << "Jet Et:" << expectedEt << "," << jet_->Et() << " error:" << err2inv 
  //  	    << "  true Et:" << truth_ << '\n';
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = jet_->Et() - expectedEt;
  chi2 *= chi2 * err2inv;
  chi2 = weight_ * Event::scaleResidual(chi2);
  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth_ << ", " << expectedEt << ", " << chi2 << '\n';
    assert(false);
  }
  
  //calculate chi2 for derivatives
  double temp1,temp2;
  const Parameters::VariationColl& varcoll = jet_->varyPars(epsilon,truth_,truth_);
  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if(( i->lowerEt < 0 ) || ( i->upperEt < 0 )) {
      flagged_bad_ = true;
      ++nflagged_;
      return 0;
    }
  }

  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    //std::cout << i->parid << ":" << i->lowerEt << ";" << i->upperEt << " diff:" 
    //	      << i->upperEt - i->lowerEt << ", " << i->upperError <<
    //  ", " << i->lowerError << '\n';
    if((std::abs((i->lowerEt - expectedEt)/expectedEt) > 0.01) || 
       (std::abs((i->upperEt - expectedEt)/expectedEt) > 0.01)) {
      std::cout << "strange extrapolation result modifying par:" << i->parid << ":" 
		<< expectedEt << "  - " << i->lowerEt << "  + " << i->upperEt 
		<< "  uncor  jet Et:" << jet_->Et() << " truth:" << truth_ << std::endl;
      //continue;
    }   
    temp1 = i->lowerEt - jet_->Et();
    err2inv = i->lowerError;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight_ * Event::scaleResidual(temp1);
    assert(temp1 == temp1);
    temp2 = i->upperEt - jet_->Et(); 
    err2inv = i->upperError;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight_ * Event::scaleResidual(temp2);
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    assert( temp_derivative1[i->parid] ==  temp_derivative1[i->parid] );
    assert( temp_derivative2[i->parid] ==  temp_derivative2[i->parid] );
  }
  assert(chi2 == chi2);
  if(chi2/weight_ > 1000) {
    std::cout << "from invert: Jet Et:" << jet_->Et() << "  expected Et:" << expectedEt << " error:" << sqrt(1/err2inv)  
	      << "  true Et:" << truth_ << "  chi2:" << chi2 << '\n';
  }
  return chi2;
}



//!  \brief Use correct likelihood
//!
//!  Calculates the summand of the negative log-likelihood
//!  for non-constant error:
//!  \f[
//!   -2\ln L = c +
//!   \sum_{i}\left(\ln\sigma^{2}_{i} + \left(\frac{x_{i}-\mu_{i}}{\sigma_{i}}\right)^{2}\right)
//!  \f]
//!
//!  \param temp_derivative1 Summand of this event to first derivative
//!  \param temp_derivative2 Summand of this event to second derivative
//!  \param epsilon Step width in numerical derivative calculation
//!  \return Summand of this event to negative log-likelihood
// ------------------------------------------------------------
double JetTruthEvent::chi2_log_fast_invert(double * temp_derivative1, 
					   double * temp_derivative2, 
					   const double* epsilon) const
{
  if(flagged_bad_) return 0;
  //find expected measurement of jet Et 
  float err2;
  float expectedEt = jet_->expectedEt(truth_,truth_,err2);
  if(expectedEt < 0) {
    flagged_bad_ = true;
    ++nflagged_;
    return 0;
  }
  err2 *= err2;
  double chi2 = jet_->Et() - expectedEt;
  chi2 *= chi2 / err2;
  chi2 = weight_ * Event::scaleResidual(log(err2) + chi2);
  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth_ << ", " << expectedEt << ", " << chi2 << '\n';
  }

  //calculate chi2 for derivatives
  double temp1,temp2;
  const Parameters::VariationColl& varcoll = jet_->varyPars(epsilon,truth_,expectedEt);
  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if(( i->lowerEt < 0 ) || ( i->upperEt < 0 )) {
      flagged_bad_ = true;
      ++nflagged_;
      return 0;
    }
  }

  for(Parameters::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if((std::abs((i->lowerEt - expectedEt)/expectedEt) > 0.05) || 
       (std::abs((i->upperEt - expectedEt)/expectedEt) > 0.05)) {
      std::cout << "strange extrapolation result modifying par:" << i->parid << ":" 
		<< expectedEt << "  - " << i->lowerEt << "  + " << i->upperEt 
		<< "  uncor  jet Et:" << jet_->Et() << " truth:" << truth_ << std::endl;
      //flagged_bad_ = true;
      //++nflagged_;
      //return 0;
      continue;
    }   
    temp1 = i->lowerEt - jet_->Et();
    err2 = i->lowerError;
    err2 *= err2;
    temp1 *= temp1 / err2;
    temp1 = weight_ * Event::scaleResidual(log(err2) + temp1);
    assert(temp1 == temp1);
    temp2 = i->upperEt - jet_->Et(); 
    err2 = i->upperError;
    err2 *= err2;
    temp2 *= temp2 / err2;
    temp2 = weight_ * Event::scaleResidual(log(err2) + temp2);
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    assert( temp_derivative1[i->parid] ==  temp_derivative1[i->parid] );
    assert( temp_derivative2[i->parid] ==  temp_derivative2[i->parid] );
  }
  assert(chi2 == chi2);
  if(chi2/weight_ > 1000) {
    std::cout << "from invert: Jet Et:" << jet_->Et() << "  expected Et:" << expectedEt << " error:" << sqrt(err2)  
	      << "  true Et:" << truth_ << "  chi2:" << chi2 << '\n';
  }
  return chi2;
}
 

void JetTruthEvent::printStats()
{
  std::cout << "JetTruthEvent: " << nflagged_ << " events flagged bad.\n";
}
