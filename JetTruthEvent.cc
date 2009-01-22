//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.cc,v 1.6 2009/01/16 08:46:40 stadie Exp $
//   

#include "JetTruthEvent.h"


double JetTruthEvent::chi2() const
{
  double diff = (jet->correctedEt(jet->Et()) - truth)/jet->Error();
  return weight * diff*diff;
}

double JetTruthEvent::chi2_fast_blobel(double * temp_derivative1, 
				       double * temp_derivative2, 
				       double const epsilon) const
{
  double et = jet->correctedEt(jet->Et());
  double err2inv = jet->expectedError(et);
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    err2inv = jet->expectedError(i->lowerEt);;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = truth - i->upperEt; 
    err2inv = jet->expectedError(i->upperEt);
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}

double JetTruthEvent::chi2_fast_simple_scaled(double * temp_derivative1, 
					      double * temp_derivative2, 
					      double const epsilon) const
{
  double et = jet->correctedEt(jet->Et());
  double s = jet->Et()/et;
  double err2inv = s * jet->expectedError(truth * s);
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << jet->Et() << ", " << s << ", " << jet->expectedError(truth/s) << ", " << chi2 << '\n';
  }
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    s = jet->Et() / i->lowerEt;
    err2inv = s * jet->expectedError(truth * s);
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = truth - i->upperEt; 
    s = jet->Et() / i->upperEt;
    err2inv = s * jet->expectedError(truth * s);
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}


double JetTruthEvent::chi2_fast_simple(double * temp_derivative1, 
				       double * temp_derivative2, 
				       double const epsilon) const
{
  double et = jet->correctedEt(jet->Et());
  double err2inv = jet->expectedError(truth * jet->Et() / et);
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    err2inv = jet->expectedError(truth * jet->Et() / i->lowerEt);
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = truth - i->upperEt; 
    err2inv = jet->expectedError(truth * jet->Et() / i->upperEt);
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}

double JetTruthEvent::chi2_fast_invert(double * temp_derivative1, 
				double * temp_derivative2, 
				double const epsilon) const
{
  //find expected measurement of jet Et 
  double scale = truth;
  double expectedEt = jet->expectedEt(truth,scale);
  if(expectedEt < 0) {//return 1.5;
    //return chi2 from modified measurement instead!
    double et = jet->correctedEt(jet->Et());
    double s= et/jet->Et();
    double err2inv = s * jet->Error();
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    double chi2 = truth - et;
    chi2 *= chi2 * err2inv;
    chi2 = weight * TData::ScaleResidual(chi2);
    return chi2;
  }
  //calculate chi2
  double err2inv = jet->expectedError(expectedEt);
  //std::cout << "Jet Et:" << expectedEt << "," << jet->Et() << " error:" << err2inv << '\n';
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = jet->Et() - expectedEt;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << scale << ", " << expectedEt << ", " << chi2 << '\n';
  }
  
  //calculate chi2 for derivatives
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyPars(epsilon,truth,scale);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if((std::abs((i->lowerEt - expectedEt)/scale) > 0.01) || 
       (std::abs((i->upperEt - expectedEt)/scale) > 0.01)) {
      std::cout << "strange extrapolation result modifying par:" << i->parid << ":" 
		<< expectedEt << "  - " << i->lowerEt << "  + " << i->upperEt 
		<< "  uncor  jet Et:" << jet->Et() << " truth:" << truth << std::endl;
      continue;
    }
    temp1 = i->lowerEt - jet->Et();
    err2inv = i->lowerError;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = i->upperEt - jet->Et(); 
    err2inv = i->upperError;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  assert(chi2 == chi2);
  return chi2;
}
