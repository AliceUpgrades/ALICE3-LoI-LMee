#ifndef __PID_H__
#define __PID_H__
	
#include <string>
	
class PID
{
public:
  PID() {}
  virtual ~PID() {}
	
  enum type {e_pi, pi_k , k_p};

  virtual bool   valid    (double eta, double p                      )=0;
  virtual double numSigma (double eta, double p,        PID::type PID)=0;
  virtual double maxP     (double eta, double numSigma, PID::type PID)=0;
  virtual double minP     (double eta, double numSigma, PID::type PID)=0;
  virtual std::string name()=0;
  virtual void description()=0;
		
protected:
	
};
	
#endif /* __PID_H__ */
