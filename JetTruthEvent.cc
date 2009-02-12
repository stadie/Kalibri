//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.cc,v 1.9 2009/02/10 08:47:26 stadie Exp $
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
  double err2inv = jet->Error();
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    //err2inv = jet->expectedError(i->lowerEt);;
    //err2inv *= err2inv;
    //err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = truth - i->upperEt; 
    //err2inv = jet->expectedError(i->upperEt);
    //err2inv *= err2inv;
    //err2inv = 1/err2inv;
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
  double c = (et - jet->EmEt() - jet->OutEt()) / jet->HadEt();
  if(c == 0) c = 1.0;
  double err2inv = c * jet->expectedError((truth - jet->EmEt() - jet->OutEt())/c + jet->EmEt() + jet->OutEt() );
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << et << ", " <<  jet->Et() << ", " << c << ", " << chi2 << '\n';
  }
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    c = (i->lowerEt - jet->EmEt() - jet->OutEt()) / jet->HadEt();  
    if(c == 0) c = 1.0;
    err2inv = c*jet->expectedError((truth - jet->EmEt() - jet->OutEt())/c + jet->EmEt() + jet->OutEt() );
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    assert(temp1 == temp1);
    temp2 = truth - i->upperEt;
    c = (i->upperEt - jet->EmEt() - jet->OutEt()) / jet->HadEt();  
    if(c == 0) c = 1.0;
    err2inv = c * jet->expectedError((truth - jet->EmEt() - jet->OutEt())/c + jet->EmEt() + jet->OutEt() ); 
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative   
  }
  assert(chi2 == chi2);
  return chi2;
}

double JetTruthEvent::chi2_fast_simple(double * temp_derivative1, 
				       double * temp_derivative2, 
				       double const epsilon) const
{
  double et = jet->correctedEt(jet->Et());
  double c = (et - jet->EmEt() - jet->OutEt()) / jet->HadEt(); 
  double err2inv = jet->expectedError((truth - jet->EmEt() - jet->OutEt())/c + jet->EmEt() + jet->OutEt() );
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = truth - et;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyParsDirectly(epsilon);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    temp1 = truth - i->lowerEt;
    c = (i->lowerEt - jet->EmEt() - jet->OutEt()) / jet->HadEt(); 
    err2inv = jet->expectedError((truth - jet->EmEt() - jet->OutEt())/c + jet->EmEt() + jet->OutEt() );
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = truth - i->upperEt;
    c = (i->upperEt - jet->EmEt() - jet->OutEt()) / jet->HadEt(); 
    err2inv = jet->expectedError((truth - jet->EmEt() - jet->OutEt())/c + jet->EmEt() + jet->OutEt() ); 
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
  double expectedEt = jet->expectedEt(truth,truth);
  if(expectedEt < 0) {
    return chi2_fast_simple_scaled(temp_derivative1,temp_derivative2,epsilon);
  }
  //calculate chi2
  double err2inv = jet->expectedError(expectedEt);
  //std::cout << "Jet Et:" << expectedEt << "," << jet->Et() << " error:" << err2inv 
  //  	    << "  true Et:" << truth << '\n';
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = jet->Et() - expectedEt;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  
  if(chi2 != chi2) {//check for NAN
    std::cout << truth << ", " << expectedEt << ", " << chi2 << '\n';
  }
  
  //calculate chi2 for derivatives
  double temp1,temp2;
  const Jet::VariationColl& varcoll = jet->varyPars(epsilon,truth,expectedEt);
  for(Jet::VariationCollIter i = varcoll.begin() ; i != varcoll.end() ; ++i) {
    if((std::abs((i->lowerEt - expectedEt)/expectedEt) > 0.01) || 
       (std::abs((i->upperEt - expectedEt)/expectedEt) > 0.01)) {
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
    assert(temp1 == temp1);
    temp2 = i->upperEt - jet->Et(); 
    err2inv = i->upperError;
    err2inv *= err2inv;
    err2inv = 1/err2inv;
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2); 
    assert(temp2 == temp2);
    temp_derivative1[i->parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
    assert( temp_derivative1[i->parid] ==  temp_derivative1[i->parid] );
    assert( temp_derivative2[i->parid] ==  temp_derivative2[i->parid] );
  }
  assert(chi2 == chi2);
  if(chi2 > 1000) {
    std::cout << "from invert: Jet Et:" << jet->Et() << "  expected Et:" << expectedEt << " error:" << sqrt(1/err2inv)  
	      << "  true Et:" << truth << "  chi2:" << chi2 << '\n';
  }
  return chi2;
}
