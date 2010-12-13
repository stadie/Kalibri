//!  \brief   Container class for jet correction factors
//
//    $Id: CorFactors.h,v 1.6 2010/10/20 11:28:16 stadie Exp $
//   
#ifndef CORFACTORS_H
#define CORFACTORS_H

class CorFactors
{
 public :
  CorFactors(float L1, float L2, float L3, float Lres, float L4, 
	     float L5, float JPT, float JPTL2L3) 
    : l1_(L1),l2_(L2),l3_(L3),lres_(Lres),l4_(L4),l5_(L5),jpt_(JPT),jptL2L3_(JPTL2L3) {};
  float getL1()  const { return l1_; }    //!< Return L1 correction factor (zero-suppression)
  float getL2()  const { return l2_; }    //!< Return L2 correction factor (relative, in eta)
  float getL3()  const { return l3_; }    //!< Return L3 correction factor (absolute, in pt)
  float getLRes()  const { return lres_; }//!< Return residual correction factor (in eta and pt)  
  float getL4()  const { return l4_; }    //!< Return L4 correction factor (electromagnetic fraction)
  float getL5()  const { return l5_; }    //!< Return L5 correction factor (flavor)
  float getJPT() const { return jpt_; }   //!< Return Jet+Track correction factor
  float getL2L3() const { return l2_*l3_; }   //!< Return product of L2 and L3 correction factors
  float getL2L3Res() const { return l2_*l3_*lres_; }   //!< Return product of L2 and L3 correction factors
  float getL2L3L4() const { return l2_*l3_*lres_*l4_; }   //!< Return product of L2,L3 and L4 correction factors
  float getJPTL2L3() const { return jptL2L3_; }   //!< Return product of L2 and L3 correction factors for Jet+Track
  float getToL2() const { return l1_*l2_; }         //!< Return factor needed to get L2 corrected from raw jets: L1*L2
  float getToL3() const { return getToL2()*l3_; }   //!< Return factor needed to get L3 corrected from raw jets: L1*L2*L3  
  float getToLRes() const { return getToL3()*lres_; }   //!< Return factor needed to get L3 corrected from raw jets: L1*L2*L3
  float getToL4() const { return getToLRes()*l4_; }   //!< Return factor needed to get L4 corrected from raw jets: L1*L2*L3*L4
  float getToL5() const { return getToL4()*l5_; }   //!< Return factor needed to get L5 corrected from raw jets: L1*L2*L3*L4*L5
  float getToJPTL3() const { return jpt_*l1_*jptL2L3_; }   //!< Return factor needed to get L3 corrected from raw jets for JPT: JPT*L1*JPTL2L3
 private :
  float l1_;      //!< Level 1 correction factor (zero-suppression)
  float l2_;      //!< Level 2 correction factor (relative, in eta)
  float l3_;      //!< Level 3 correction factor (absolute, in pt)
  float lres_;    //!< residual correction in data
  float l4_;      //!< Level 4 correction factor (electromagnetic fraction, JW)
  float l5_;      //!< Level 5 correction factor (flavor)
  float jpt_;     //!< Jet+Track correction factor
  float jptL2L3_; //!< Product of level 2 and level 3 correction factors for Jet+Track
}; 

#endif
