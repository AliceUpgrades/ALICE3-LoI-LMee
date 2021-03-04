/// \author R+Preghenella
/// \email  preghenella@bo.infn.it
/// \date   March 2020

#ifndef __GENERICDETECTOR_H__
#define __GENERICDETECTOR_H__

#include "PID.h"

class genericDetector : public PID
{
public:
  genericDetector() = default;
  virtual ~genericDetector() = default;

  enum EDetector_t {kBarrel, kForward};
  
  /** setters **/
  void setType(EDetector_t val) { mType = val; };
  void setName(string val)      { mName = val; };
  void setLength(double val)    { mLength = val; };
  void setRadius(double val)    { mRadius = val; };
  void setPositionZ(double val) { mPositionZ = val; };
  void setRadiusIn(double val)  { mRadiusIn = val; };
  void setRadiusOut(double val) { mRadiusOut = val; };
  void setMagneticField(double val) { mMagneticField = val; };
  
  /** methods to override **/
  virtual bool valid (double eta, double p) override { return isHit(eta, p); };
  virtual double numSigma (double eta, double p, PID::type PID) override = 0;
  virtual double maxP (double eta, double numSigma, PID::type PID) override = 0;
  virtual double minP (double eta, double numSigma, PID::type PID) override = 0;
  string name () override { return mName; };
  void description() override { return mDescription; };

  double maxPt (double eta, double numSigma, PID::type PID) { return maxP(eta, numSigma, PID) / cosh(eta); };
  double minPt (double eta, double numSigma, PID::type PID) { return ptMin(); };
  
  double etaMin();
  double etaMax();
  double ptMin();
  double pMin(double eta) { return ptMin() * cosh(eta); };
  bool   isHit(double eta, double p) { return eta > etaMin() && eta < etaMax() && p > pMin(eta); };
  double trackLength(double eta);
  
protected:

  string mName          = "genericDetector"; 
  string mDescription   = "Detector description"; 
  EDetector_t mType     = kBarrel;
  double mLength        = 200.; // [cm]
  double mRadius        = 200.; // [cm]
  double mPositionZ     = 200.; // [cm]
  double mRadiusIn      = 20.;  // [cm]
  double mRadiusOut     = 200.; // [cm]
  double mMagneticField = 2.;   // [T]

  const double mLightSpeed   = 29.9792458;    // speed of light [cm/ns]
  const double mMassElectron = 0.00051099891; // electron mass [GeV]
  const double mMassPion     = 0.13957018;    // pion mass [GeV]
  const double mMassKaon     = 0.493677;      // kaon mass [GeV]
  const double mMassProton   = 0.93827208816; // proton mass [GeV]

};

double
genericDetector::etaMin()
{
  switch (mType) {
  case kBarrel:
    return -log( tan( atan2(mRadius, -mLength) * 0.5 ) );
  case kForward:
    return -log( tan( atan2(mRadiusOut, mPositionZ) * 0.5 ) );
  }
  return 0.;
}

double
genericDetector::etaMax()
{
  switch (mType) {
  case kBarrel:
    return -log( tan( atan2(mRadius, mLength) * 0.5 ) );
  case kForward:
    return -log( tan( atan2(mRadiusIn, mPositionZ) * 0.5 ) );
  }
  return 0.;
}

double
genericDetector::ptMin()
{
  switch (mType) {
  case kBarrel:
    return 0.003 * mMagneticField * 0.5 * mRadius;
  case kForward:
    return 0.003 * mMagneticField * 0.5 * mRadiusIn;
  }
  return 0.;
}

double
genericDetector::trackLength(double eta)
{
  auto theta = 2. * atan( exp(-eta) );

  switch (mType) {
  case kBarrel:
    return mRadius / sin(theta);
  case kForward:
    return mPositionZ / cos(theta);
  }
  return 0.;
}

#endif /* __GENERICDETECTOR_H__ */
