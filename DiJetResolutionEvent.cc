// $Id: $

#include "DiJetResolutionEvent.h"

//#include <iomanip>


// --------------------------------------------------
DiJetResolutionEvent::DiJetResolutionEvent(Jet* jet1, Jet* jet2, double deltaPhi12, double pPhi,
					   double ptJet3, double ptJet4, double pJ3, double pSJ, double ptRef,
					   double ptHat, double weight, const ResolutionFunction& pdf,
					   double min, double max, double eps, int niter)
  : Event(weight,ptHat,0),
    pdf_(&pdf), jet1_(jet1), jet2_(jet2),
    deltaPhi12_(deltaPhi12), pPhi_(pPhi),
    ptRef_(ptRef), ptJet3_(ptJet3), ptJet4_(ptJet4), pJ3_(pJ3), pSJ_(pSJ),
    kMaxNIter_(niter),kEps_(eps),kMin_(min),kMax_(max) {};



// --------------------------------------------------
DiJetResolutionEvent::~DiJetResolutionEvent() { 
  delete jet1_;
  delete jet2_;
}


// --------------------------------------------------
void DiJetResolutionEvent::setParameters(Parameters* param) { 
  jet1_->setParameters(param);
  jet2_->setParameters(param);
  pdf_ = &(param->function(*pdf_));
}


//!  \brief Get the negative log-likelihood and the
//!         contributions to the 1. and 2. derivatives
//!
//!  Calculates the negative log-likelihood \f$ -2\ln L \f$
//!  of this event. Moreover the contribution to the 
//!  first and second derivative ('temp_derivative1',
//!  'temp_derivative2') of the global \f$ -\ln L \f$
//!  function is calculated numerically and returned
//!  by reference, where
//!  \f[
//!   -2\frac{\partial (\ln L) }{\partial p}
//!   = 2\sum \frac{\textrm{temp\_derivative1}}{2\epsilon}
//!  \f]
//!  and
//!  \f[
//!   -2\frac{\partial^{2} (\ln L) }{\partial p^{2}}
//!   = 2\sum \frac{\textrm{temp\_derivative2}}{\epsilon^{2}}
//!  \f]
//!
//!  \param temp_derivative1 Pointer to first derivative contribution
//!  \param temp_derivative2 Pointer to second derivative contribution
//!  \param epsilon Step sizes \f$ \epsilon \f$ for derivative calculation
//!  \return The negative log-likelihood of this event
// --------------------------------------------------
double DiJetResolutionEvent::chi2_fast(double * temp_derivative1,
			    double * temp_derivative2,
			    const double* epsilon) const {
  double f = chi2();

  double   oldpar;
  double   temp1;
  double   temp2;

  // Vary parameters of response pdf
  for(int i = 0; i < pdf_->nPars(); i++) {
    oldpar = pdf_->firstPar()[i];

    pdf_->firstPar()[i] += epsilon[pdf_->parIndex()+i];
    temp1 = chi2();

//     std::cout << std::endl;
//     std::cout << std::setprecision(10) << i << ": " << oldpar << " -- > " << f << std::endl;
//     std::cout << "   " << pdf_->firstPar()[i] << " -- > " << temp1 << std::endl;


    pdf_->firstPar()[i] -= 2.*epsilon[pdf_->parIndex()+i];
    temp2 = chi2();
    
    //    std::cout << "   " << pdf_->firstPar()[i] << " -- > " << temp1 << std::endl;


    pdf_->firstPar()[i] = oldpar;

    //    std::cout << i << ": " << pdf_->parIndex()+i << " += " << temp1 - temp2 << std::endl;


    temp_derivative1[pdf_->parIndex()+i] += temp1 - temp2;
    temp_derivative2[pdf_->parIndex()+i] += temp1 + temp2 - 2*f;

    if( !(f == f) ) {
      std::cerr << "ERROR in DiJetResolutionEvent::chi2_fast(): nan error\n";
      std::cerr << "  " << i << ": " << pdf_->firstPar()[i] << std::endl;
    }
  }

  return f;
}



//! \brief Using pdf without spectrum (delta function)
// --------------------------------------------------
double DiJetResolutionEvent::chi2Simple() const {
   double pdf = pdfPtMeas(jet1()->pt(),jet2()->pt(),0.);
   if( pdf != pdf ) {
     pdf = 0.;
     std::cerr << "WARINGING: pdf = " << pdf << std::endl;
   } else if( pdf < 0. ) {
     pdf = 0.;
     std::cerr << "WARINGING: pdf = " << pdf << std::endl;
   } else if( pdf > 0. ) {    
     pdf = -2.*log(pdf);
     if( pdf <= 0. ) std::cerr << "WARNING: pdf = " << pdf << std::endl;
     if( pdf != pdf) std::cout << ">> pdf = " << pdf << std::endl;
   } else {
     pdf = 0.;
   }
 return weight()*pdf;
}



//!  \brief Get the negative log-likelihood of this event
//!
//!  Calculates the probability of this event configuration
//!  \f[
//!   P(m_{1},m_{2}) = \int\;dt\;f(t)r(m_{1}/t)r(m_{2}/t),
//!  \f]
//!  where \f$ f \f$ is the truth pdf and \f$ r \f$ is the
//!  response pdf. \f$ r \f$ is normalized to unity,
//!  \f$ f \f$ is normalized such that \f$ P \f$ is 
//!  normalized to unity (see truthPDF(t)). The method
//!  returns \f$ -\ln(P) \f$.
//!
//!  The integral is calculated numerically using the
//!  Simpson's 3/8 rule.
//!
//!  \return The negative log-likelihood of this event
// --------------------------------------------------
double DiJetResolutionEvent::chi2Spectrum() const {

  double h = kMax_ - kMin_;     // Integration interval
  double pint = 0.;              // Current value of integral over response pdfs
  double pint_old  = 1.;              // Value of integral over response pdfs from previous iteration
  double eps = 1.;
  int nIter = 0;               // Current iteration in interval splitting
  std::vector<double> pp;         // Product of current function values of response pdfs
  std::vector<double> pp_old;     // Product of function values of response pdfs from previous iteration

  // Iterate until precision or max. number iterations reached
  while( eps > kEps_ && nIter < kMaxNIter_ ) {
    pint_old = pint;
    pint     = 0;
    pp_old   = pp;
    pp.clear();
    h       /= 3.;    // In each iteration, split h into 3 new intervals
    
    // Loop over nodes xi i.e. interval borders
    for(int i = 0; i <= pow(3.0,nIter+1); ++i){
      double t = kMin_ + i * h;  // Pt at node
      
      // Calculate probability only at new nodes
      if(nIter == 0 || i % 3 != 0) {
	pp.push_back(pdfPtMeas(jet1()->pt(),jet2()->pt(),t)*pdfPtTrue(t));
      } else {
	pp.push_back(pp_old.at(i/3));       // Store product of pdfs previously calcluated
      }
    }

    // Sum up weighted function values
    for(unsigned int i = 0; i < pp.size(); i++)	{
      double w = 1.;                       // Weight w from Simpson's rule
      if( i > 0 && i < (pp.size() - 1) ) { // w = 1 for x0 and last node
	if( i % 3 == 0 ) {                 // w = 2 for x3, x6, ...
	  w = 2.;
	} else {
	  w = 3.;
	}
      }
      pint += w * (pp.at(i));              // Sum up weighted function values
    }
    pint *= (3. * h / 8.);                 // Apply overall normalization
    nIter++;

    if( pint_old ) eps = std::abs((pint - pint_old) / pint_old);
  }
  if( !(pint == pint) ) {
    std::cerr << "ERROR in DiJetResolutionEvent::chi2Spectrum(): pint = nan" << std::endl;
    std::cerr << "  " << jet1()->genPt() << ",  " << jet1()->pt() << std::endl;
    std::cerr << "  " << jet2()->genPt() << ",  " << jet2()->pt() << std::endl;
  }
  if( pint <= 0. ) return 0.;
  
  return  weight()*(-2.*log(pint));
}



//!  \brief Print event parameters
// --------------------------------------------------
void DiJetResolutionEvent::printInitStats() const {
  std::cout << "Event type: " << type() << "\n";
  std::cout << "Integration parameters\n";
  std::cout << "  niter: " << kMaxNIter_ << "\n";
  std::cout << "  eps:   " << kEps_ << "\n";
  std::cout << "  range: " << kMin_ << " < truth < " << kMax_ << " (GeV)\n";
}


