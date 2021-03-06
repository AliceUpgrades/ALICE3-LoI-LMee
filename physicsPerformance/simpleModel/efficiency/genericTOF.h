/// \author R+Preghenella
/// \email  preghenella@bo.infn.it
/// \date   March 2020

#ifndef __GENERICTOF_H__
#define __GENERICTOF_H__
	
#include "genericDetector.h"
	
class genericTOF : public genericDetector
{
public:
  genericTOF() = default;
  virtual ~genericTOF() = default;

  /** setters **/
  void setSigma(double val) { mSigma = val; };

  /** getters **/
  double getSigma() const { return mSigma; };
  
  /** methods to override **/
  double numSigma (double eta, double p, PID::type PID) override;
  double maxP (double eta, double numSigma, PID::type PID) override;
  double minP (double eta, double numSigma, PID::type PID) override { return pMin(eta); };
  
  double timeOfFlight(double l, double m, double p) { return l / (p / sqrt(m * m + p * p) * mLightSpeed); };
  
protected:

  // TOF parameters
  double mSigma; // time-of-flight resolution [ns]

};
	
double genericTOF::numSigma(double eta, double p, PID::type PID)
{
  double mass1, mass2;
  switch (PID) {
  case e_pi:
    mass1 = mMassElectron;
    mass2 = mMassPion;
    break;
  case pi_k:
    mass1 = mMassPion;
    mass2 = mMassKaon;
    break;
  case k_p:
    mass1 = mMassKaon;
    mass2 = mMassProton;
    break;
  }
  double length = trackLength(eta);
  return (timeOfFlight(length, mass2, p) - timeOfFlight(length, mass1, p)) / mSigma;
}

double genericTOF::maxP(double eta, double numSigma, PID::type PID)
{
  double mass1, mass2;
  switch (PID) {
  case e_pi:
    mass1 = mMassElectron;
    mass2 = mMassPion;
    break;
  case pi_k:
    mass1 = mMassPion;
    mass2 = mMassKaon;
    break;
  case k_p:
    mass1 = mMassKaon;
    mass2 = mMassProton;
    break;
  }
  double length = trackLength(eta);

  /** let's do it numerically, starting from the threshold **/
  double p = pMin(eta);
  while (timeOfFlight(length, mass2, p) - timeOfFlight(length, mass1, p) > numSigma * mSigma) {
    p += 0.001;
  }
  return p;
}

#endif /* __GENERICTOF_H__ */
