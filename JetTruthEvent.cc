//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: JetTruthEvent.cc,v 1.2 2008/12/17 09:38:36 stadie Exp $
//   

#include "JetTruthEvent.h"


double JetTruthEvent::chi2() const
{
  double diff = (jet->correctedEt(jet->Et()) - truth)/jet->Error();
  return weight * diff*diff;
}


double JetTruthEvent::chi2_fast(double * temp_derivative1, double * temp_derivative2, 
				double const epsilon) const
{
  //find expected measurement of jet Et 
  double scale = truth;
  double expectedEt = jet->expectedEt(truth,scale);
  //calculate chi2
  double err2inv = jet->Error();
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double chi2 = jet->Et() - expectedEt;
  chi2 *= chi2 * err2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  //calculate chi2 for derivatives
  int npar = jet->nPar();
  double et1,et2,temp1,temp2;
  for(int i = 0 ; i < npar ; ++i) {
    int parid = jet->varyPar(i,epsilon,truth,scale,et2,et1);
    //std::cout << truth << ", " << expectedEt << "," << et1 << ", " << et2 << std::endl; 
    if((std::abs((et1 - expectedEt)/scale) > 0.05) || 
       (std::abs((et2 - expectedEt)/scale) > 0.05)) {
      std::cout << "strange extrapolation result:" << expectedEt << "; " 
		<< et1 << "; " << et2 << " jet:" << jet->Et() 
		<< " truth:" << truth << std::endl;
    }
    temp1 = et1 - jet->Et(); 
    temp1 *= temp1 * err2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = et2 - jet->Et();
    temp2 *= temp2 * err2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}
