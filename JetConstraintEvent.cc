//
//    Class for constraints on the jet correction
//
//    first version: Hartmut Stadie 2009/07/23
//    $Id: JetConstraintEvent.cc,v 1.10 2010/11/01 15:47:43 stadie Exp $
//   


#include "JetConstraintEvent.h"


//interface to Data
JetConstraintEvent::~JetConstraintEvent()
{
  for(unsigned int i = 0, njets = jets_.size() ; i < njets ; ++i) {
    delete jets_[i];
  }
  jets_.clear();
}

void JetConstraintEvent::addJet(double truePt, const Jet* j, 
				Function* globalFunc) 
{
  if((truePt > maxpt_) || (truePt < minpt_) || 
     (std::abs(j->eta()) > maxeta_) || (std::abs(j->eta()) < mineta_)) return; 
  
  Jet* jet = j->clone();
  if(globalFunc) jet->setGlobalFunction(*globalFunc);
  error_ = error_ * jets_.size() + jet->error();
  jets_.push_back(jet); 
  trusum2_ += jet->Et() * jet->Et();
  trusum_ += jet->Et();
  error_ = sqrt(trusum2_ -  trusum_ *  trusum_/jets_.size())/jets_.size();
  ptHat_ = truth();
}


void JetConstraintEvent::setParameters(Parameters* param) { 
  for(unsigned int i = 0, njets = jets_.size() ; i < njets ; ++i) {
    jets_[i]->setParameters(param);
  }
}
 
double JetConstraintEvent::chi2() const
{
  double sumet = 0;
  unsigned int njets = jets_.size();
  if(! njets) return 0;
  for(unsigned int i = 0 ; i < njets ; ++i) {
    sumet += jets_[i]->correctedEt(jets_[i]->Et());
  }
  double chi2 = (sumet - trusum_) / error_;
  return weight_ * chi2 * chi2 / jets_.size();
}

double JetConstraintEvent::chi2_fast(double * temp_derivative1, double * temp_derivative2,  double * temp_derivative3, double * temp_derivative4, const double* epsilon) const
{
  double sumet = 0;
  unsigned int njets = jets_.size();
  if(! njets) {
    chi2plots_ = 0;
    return 0;
  }
  for(unsigned int i = 0 ; i < njets ; ++i) {
    sumet += jets_[i]->correctedEt(jets_[i]->Et());
  }
  double chi2 = (sumet - trusum_) / error_;
  chi2 *= chi2;
  chi2 *= weight_;
  std::cout << "JetConstraintEvent::chi2_fast:  Et:" << trusum_/jets_.size() 
     	    << " chi2:" << chi2 << " sum Et:" << sumet << " error:" << error_ << std::endl;
  //derivative....
  for(VarMap::iterator i = varmap_.begin() ; i != varmap_.end() ; ++i) {
    i->second.uppersum_ = sumet;
    i->second.lowersum_ = sumet;
  }
  for(unsigned int i = 0 ; i < njets ; ++i) {
    const Parameters::VariationColl& varcoll = jets_[i]->varyParsDirectly(epsilon,false);
    double et = jets_[i]->correctedEt(jets_[i]->Et());
    for(Parameters::VariationCollIter j = varcoll.begin() ; j != varcoll.end() ; ++j) {
      VarMap::iterator k = varmap_.find(j->parid);
      if(k == varmap_.end()) {
	k = varmap_.insert(k,std::make_pair(j->parid,Variation()));
	k->second.uppersum_ = sumet;
	k->second.lowersum_ = sumet;
      }
      k->second.uppersum_ += j->upperEt - et;
      k->second.lowersum_ += j->lowerEt - et;
    }
  }
  for(VarMap::iterator i = varmap_.begin() ; i != varmap_.end() ; ++i) {
    double temp1 = (i->second.lowersum_ - trusum_) / error_;
    temp1 *= temp1;
    temp1 *= weight_;
    double temp2 = (i->second.uppersum_ - trusum_) / error_;
    temp2 *= temp2;
    temp2 *= weight_;
//     std::cout << "JetConstraintEvent::chi2_fast:  Et:" << jets[0]->Et() 
// 	      << " temp1:" << temp1 << " sum Et:" <<  i->second.lowersum 
// 	      << std::endl;
//     std::cout << "JetConstraintEvent::chi2_fast:  Et:" << jets[0]->Et() 
// 	      << " temp2:" << temp2 << " sum Et:" <<  i->second.uppersum 
// 	      << std::endl;
    temp_derivative1[i->first] += (temp2 - temp1); // for 1st derivative
    temp_derivative2[i->first] += (temp2 + temp1 - 2 * chi2); // for 2nd derivative   
  }
  chi2plots_ = chi2;
  return chi2;
}
