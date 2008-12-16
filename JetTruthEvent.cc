//
//    Class for all events with one jet and truth informatio
//
//    first version: Hartmut Stadie 2008/12/14
//    $Id: EventReader.h,v 1.1 2008/12/12 13:43:15 stadie Exp $
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
  double expectedEt = 0;
  {
    //find root of truth - jet->correctedEt(expectedEt)
    int i = 0;
    double x1 = truth * 0.5,x2 = 2 * truth;
    double y1 = truth - jet->correctedEt(x1), y2 = truth - jet->correctedEt(x2);
    //std::cout << y1 << "," << y2 << '\n';
    assert(y2 < 0);
    assert(y1 > 0);
    double xm,ym;
    while(std::abs(x1 - x2) > 0.1) { 
      xm  = 0.5 * (x1 + x2);
      ym = jet->correctedEt(xm);
      if(ym < 0) {
	y2 = ym;
	x2 = xm;
      } else {
	x1 = xm;
	y1 = ym;
      }
      ++i;
      if(i < 20) break;
    }
    //std::cout << "I:" << i << ", " <<  y1 << "," << y2 << '\n';
    expectedEt = 0.5 * (x1 + x2);
  }
  //expectedEt = jet->Et();
  //calculate chi2
  double err2inv = jet->Error();
  err2inv *= err2inv;
  err2inv = 1/err2inv;
  double s = jet->correctedEt(expectedEt)/expectedEt;
  double et = s * jet->Et();
  double s2inv =1/s;
  s2inv *= s2inv;
  double chi2 = et - truth;
  chi2 *= chi2 * err2inv * s2inv;
  chi2 = weight * TData::ScaleResidual(chi2);
  //calculate chi2 for derivatives
  int npar = jet->nPar();
  double et1,et2,temp1,temp2;
  for(int i = 0 ; i < npar ; ++i) {
    int parid = jet->varyPar(i,epsilon,expectedEt,et2,et1);
    et1 = et1/expectedEt * jet->Et();
    et2 = et2/expectedEt * jet->Et();
    temp1 = et1 - truth; 
    s2inv = jet->Et()/et1;
    s2inv *= s2inv;
    temp1 *= temp1 * err2inv * s2inv;
    temp1 = weight * TData::ScaleResidual(temp1);
    temp2 = et2 - truth;
    s2inv = jet->Et()/et2;
    s2inv *= s2inv;
    temp2 *= temp2 * err2inv * s2inv;
    temp2 = weight * TData::ScaleResidual(temp2);
    temp_derivative1[parid] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[parid] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative
  }
  return chi2;
}
