//
// $Id: CalibData.h,v 1.82 2010/10/20 11:28:16 stadie Exp $
//
#ifndef CalibData_h
#define CalibData_h

#include <iostream>
#include <vector> 
#include <cmath>
#include <cassert>

class Parameters;

//!  \brief Type of data
//!
//!  \sa TAbstractData 
enum DataType {Default, TrackTower, GammaJet, TrackCluster, MessMess, PtBalance,
               InvMass, typeTowerConstraint, ParLimit, TypeSmearPhotonJet, TypeSmearDiJet, JetConstraint };

//!  \brief Base class of a measurement
//!
//!  A measurement can represent a tower, a track, or a jet.
//!
//!  \note: The parametrized (pt-)measurement, wich will be compared to the
//!      'truth' will remain a single 'double' value (i.e. the 
//!       same type as 'truth')!
//!
//!  \sa Jet, Tower, Track, JetWithTowers, JetWithTracks
//!
//!  \author Christian Autermann, Hartmut Stadie
//!  $Id: CalibData.h,v 1.82 2010/10/20 11:28:16 stadie Exp $
class Measurement
{
public:
 Measurement() :
  pt(0.),EMF(0.),HadF(0.),OutF(0.),E(0.),eta(0.),phi(0.),phiphi(0),
   etaeta(0)
    {
    }
 Measurement(float Et,float EmEt,float HadEt,float OutEt,float E,
	     float eta,float phi, float nphiphi = 0,float netaeta = 0)
   : pt(Et),EMF(EmEt),HadF(HadEt),OutF(OutEt),E(E),eta(eta),phi(phi),
    phiphi(nphiphi),etaeta(netaeta) 
  {
    assert(phiphi == phiphi);
  }
  virtual ~Measurement() {};
  //all common variables
  float pt;     //!< Total transverse momentum (pt = EMF + HadF + OutF)
  float EMF;    //!< Pt from the ECAL part of the tower(s)		
  float HadF;   //!< Pt from the HCAL part of the towers(s)		
  float OutF;   //!< Pt fromt the HO part of the tower(s)		
  float E;      //!< Total energy					
  float eta;    //!< Pseudorapidity eta				
  float phi;    //!< Polar angle phi  
  float phiphi; //!< Phi-Phi moment (width in phi) 
  float etaeta; //!< Eta-Eta moment (width in eta) 
};




//!  \brief A track measurement
//!
//!  \sa Measurement, TJet, TTower, Jet, JetWithTowers
//!
//!  \todo Document members
//!
//!  \author Jan Thomsen
//!  $Id: CalibData.h,v 1.82 2010/10/20 11:28:16 stadie Exp $
class TTrack : public Measurement
{
public:
  TTrack():Measurement(){};
  TTrack(float Et, float EmEt, float HadEt ,float OutEt, float E,float eta,
	 float phi,int TrackId, int TowerId, float DR, float DRout, 
	 float etaOut, float phiOut, float EM1, float EM5, float Had1, 
	 float Had5, float TrackChi2, int NValidHits, bool TrackQualityT, 
	 float MuDR, float MuDE, float Efficiency) 
    : Measurement(Et,EmEt,HadEt,OutEt,E,eta,phi),TrackId(TrackId),TowerId(TowerId),
    DR(DR),DRout(DRout),etaOut(etaOut),phiOut(phiOut),EM1(EM1),EM5(EM5),Had1(Had1),
    Had5(Had5),TrackChi2(TrackChi2),NValidHits(NValidHits),TrackQualityT(TrackQualityT),
    MuDR(MuDR),MuDE(MuDE),Efficiency(Efficiency) {}
  virtual ~TTrack() {}
//variables specific only to Tracks
  int TrackId;
  int TowerId;
  float DR;
  float DRout;
  float etaOut;
  float phiOut;
  float EM1;
  float EM5;
  float Had1;
  float Had5;
  float TrackChi2;
  int NValidHits;
  bool TrackQualityT;
  float MuDR;
  float MuDE;
  float Efficiency;
};


//!  \brief Interface to the data 
//!
//!  A Event object represents one event. It holds the measured
//!  quantities of that event (see Measurement) and allows
//!  access to the corrected measurement. Moreover, the normalized,
//!  weighted, squared, and squared residual \f$ z^{2} \f$ of this
//!  event, which enters the global \f$ \chi^{2} = \sum z^{2} \f$
//!  function, is calculated.
//!
//!  Event is a virtual base class. The derived interfaces are
//!  specific for a certain type of data.
//!
//!  There are currently two different calibration schemes, resulting
//!  in two different sets of data classes derived from Event:
//!  -# Original calibration scheme ("Correction of the measurement")
//!     There is a second base class for this scheme, TAbstractData,
//!     deriving from Event. All interfaces for specific data types
//!     derive from TAbstractData in this scheme.
//!  -# New calibration scheme ("Variation of the truth") 
//!     The available data types are:
//!  \author Christian Autermann
//!  \date Wed Jul 18 13:54:50 CEST 2007
//! $Id: CalibData.h,v 1.82 2010/10/20 11:28:16 stadie Exp $
class Event
{
public:
 Event(double w = 0, double pthat = 0) : weight_(w),ptHat_(pthat) {}
  virtual ~Event() {}
  virtual Measurement *mess() const = 0;                           //!< Get Measurement object
  virtual double truth() const = 0;                                 //!< Get truth of measurement
  virtual double parametrizedMess() const = 0;                      //!< Get corrected measurement
  virtual void setParameters(Parameters* param) = 0;                  //!< Set Parameters
  virtual DataType type() const = 0;                                //!< Get DataType
  double weight() const { return weight_;}                          //!< Get weight
  void   setWeight(double w)  {weight_ = w;}                           //!< Set weight
  double ptHat() const { return ptHat_; }                              //!< Get event scale


  //!  \brief Get the normalized, squared residual \f$ z^{2} \f$ of this event
  //!
  //!  The normalized, squared residual \f$ z^{2} \f$ of
  //!  this event is calculated. It is weighted with
  //!  GetWeight(), and scaled with ScaleResidual.
  //!  It enters the global \f$ \chi^{2} = \sum z^{2} \f$
  //!  function.
  //!
  //!  \return The normalized, squared residual\f$ z^{2} \f$ of this event
  virtual double chi2() const = 0;


  //!  \brief Chi2 value from last iteration
  //!  \return Chi2 value from last iteration
  // ------------------------------------------
  virtual double chi2_plots() const = 0;


  //!  \brief Get the normalized, squared residual\f$ z^{2} \f$ of this event
  //!         and calculate the first and second derivatives
  //!
  //!  The normalized, squared residual \f$ z^{2} \f$ of
  //!  this event is calculated. It is weighted with
  //!  GetWeight(), and scaled with ScaleResidual.
  //!  It enters the global \f$ \chi^{2} = \sum z^{2} \f$
  //!  function.
  //!
  //!  Moreover, the contribution of this event to the 
  //!  first and second derivative ('temp_derivative1',
  //!  'temp_derivative2') of the global \f$ \chi^{2} \f$
  //!  function is calculated numerically and returned
  //!  by reference, where
  //!  \f[
  //!   \frac{\partial \chi^{2} }{\partial p}
  //!   = \sum \frac{\textrm{temp\_derivative1}}{2\epsilon}
  //!  \f]
  //!  and
  //!  \f[
  //!   \frac{\partial^{2} \chi^{2} }{\partial p^{2}}
  //!   = \sum \frac{\textrm{temp\_derivative2}}{\epsilon^{2}}
  //!  \f]
  //!
  //!  \param temp_derivative1 Pointer to first derivative contribution
  //!  \param temp_derivative2 Pointer to second derivative contribution
  //!  \param epsilon Step sizes \f$ \epsilon \f$  for derivative calculation
  //!  \return The normalized, squared residual\f$ z^{2} \f$ of this event
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, const double *epsilon) const = 0;


  virtual void updateError() = 0;  //!< Update error terms using current corrected energies


  //!  \brief Scale residual for outlier treatment
  //!
  //!  Points to one of the following functions to
  //!  scale the squared, normalized, and weighted
  //!  residual
  //!  \f$ z^{2} = \chi^{2}/\textrm{weight} \f$:
  //!   - scaleNone(double z2)
  //!   - scaleCauchy(double z2)
  //!   - scaleHuber(double z2)
  //!   - scaleTukey(double z2)
  //!
  //!  \param z2 Normalized and squared residual
  //!  \return Scaled residual
  static double (*scaleResidual)(double z2);


  //!  \brief No scaling of residuals
  //!
  //!  \note This is the default
  //!
  //!  \param z2 Normalized and squared residual
  //!  \return Scaled residual
  static double scaleNone(double z2){ return z2; }


  static double scaleCauchy(double z2);  //!< Scaling of residual with Cauchy function
  static double scaleHuber(double z2);   //!< Scaling of residual with Huber function  
  
  //!  \brief Cut on residuals
  //!
  //!  discards events with \f$ |residual| > 1.5 \sigma \f$
  //!
  //!  \param z2 Normalized and squared residual
  //!  \return Scaled residual
  static double scaleTukey(double z2);  //!< Scaling of residual a la  Tukey


 protected:
  double weight_;
  double ptHat_;  
};




#endif
