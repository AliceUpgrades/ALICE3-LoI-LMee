//R__LOAD_LIBRARY(libDelphes)
//R__LOAD_LIBRARY(libDelphesO2)

#include "anaEEstudy.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"
// #include "pid.h"
#include <algorithm>
using namespace std;


// class GenParticle;
// class Track;
// class TClonesArray;
// #include <GenParticle.h>
// #include <Track.h>
// #include <TClonesArray.h>
// #include <DelphesClasses.h>




bool bSmear        = true;
bool bUseiTOF      = false;
bool bSemiCentral  = false;
bool bUseWeights  = true;

double Bz = 5;            // in kG,  becomes overwritten by the generateEfficiencies.sh skript
double eMass = 0.000511;



// --- binning ---
// Double_t pt_binning[]  = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0};
Double_t pt_binning[]  = {0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0};
// Double_t eta_binning[]  = {-5.0,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.0,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0};
Int_t nEtaBins = 100;    Double_t maxEta = 5.;
Int_t nPhiBins = 180;    Double_t maxPhi = 2*TMath::Pi();

// const Int_t n_ptee_bin_c = 520;
// Double_t ptee_bin_c[n_ptee_bin_c+1] = {};
// Double_t mee_bin_c[2];
// Int_t n_mee_bin_c = 1;
// Double_t dca_bin_c[2];
// int n_dca_bin_c = 1;

// ptee
// const Int_t n_ptee_bin_c = 202;
const Int_t n_ptee_bin_c = 105;
Double_t ptee_bin_c[n_ptee_bin_c+1] = {};
// mee
Double_t mee_bin_c[121];
Int_t n_mee_bin_c = 120;
// DCA
Double_t dca_bin_c[101];
int n_dca_bin_c = 100;


// outer TOF (100 cm)
double tof_radius = 100.; // [cm]
double tof_length = 280.; // [cm]
double tof_sigmat = 0.02; // [ns]
double tof_sigma0 = 0.20; // [ns]

// inner TOF (19 cm)
double inner_tof_radius = 20.;  // [cm]
double inner_tof_length = 56.;  // [cm]
double inner_tof_sigmat = 0.05; // [ns]
double inner_tof_sigma0 = 0.20; // [ns]


// TOF ele pt acceptance
double tof_EleAccep_p_cut = 0.6;  // [GeV/c]
double itof_EleAccep_p_cut = 0.3;  // [GeV/c]
// TOF additional Pion/Muon rejection
double tof_PionRej_p_cut = 0.4; // [GeV/c] above this value the nSigmaRICHEle_forTOFPID cut is apllied. Thus the TOF is only accepting particles within nSigmaRICHEle and nSigmaTOFele
double nSigmaRICHEle_forTOFPID = 3.; //


// RICH
double rich_radius = 100.;
double rich_length = 280.;

// RICH ele pt acceptance
// double rich_Ele_pt_cut = 1.8; // [GeV/c]
double rich_PionRejection_p_cut = 1.0; // [GeV/c]

// Prefilter PID properties
bool   pf_PID_useTOF      = kTRUE;
bool   pf_PID_useRICH     = kTRUE;
bool   pf_PID_usePS       = kTRUE;
double pf_PID_NSigTOFele  = 3.;
double pf_PID_NSigTOFpi   = 3.0;
double pf_PID_NSigRICHele = 3.;
double pf_PID_NSigRICHpi  = 4.0;
double pf_PID_PSptacc     = 1.5;



// PID Scenarios            =   /*    woPF,   TOFRICH,TOFPS,TOFRICHPS PID             */ /*  30-50% large eta      */ /*  0-10% additional Prefilter studies    */       compare eta acceptance with MC PID                               different PF pt cuts                        /*     TOF, RICH,w wo PID */        different pt cuts           different PID scenarios                 using iTOF layer         only TOF, RICH, both           use TOF   use only PS   use TOF+PS
bool bUseTrackEff           =   /*                   kTRUE                            */ /*            kTRUE       */ /*            kTRUE                       */ /*                       kFALSE                         */  /*                      kTRUE                      */  /*         kFALSE         */                                                                                                                                 /*   kTRUE */    kTRUE     /*   kTRUE */  ;
bool useMCPID[]             = { /*kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE */ /* kFALSE, kFALSE, kFALSE */ /* kFALSE, kFALSE, kFALSE, kFALSE, kFALSE */ /*  kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE  */  /* kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE  */  /*  kTRUE,  kTRUE,  kTRUE */                                                                                                                                     kFALSE,     kFALSE,        kFALSE    };
bool usePreFilter[]         = { /*kFALSE,  kFALSE,  kFALSE,   kTRUE,   kTRUE,   kTRUE */ /* kFALSE,  kTRUE,  kTRUE */ /* kFALSE,  kTRUE,  kTRUE,  kTRUE,  kTRUE */ /*  kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE  */  /* kFALSE,  kTRUE,  kTRUE,  kTRUE,  kTRUE,  kTRUE  */  /* kFALSE, kFALSE, kFALSE */                                                                                                                                     kFALSE,     kFALSE,        kFALSE    };
bool useTOFPID[]            = { /* kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE */ /*  kTRUE,  kTRUE,  kTRUE */ /*  kTRUE,  kTRUE,  kTRUE,  kTRUE,  kTRUE */ /* kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE  */  /*  kTRUE,  kTRUE,  kTRUE,  kTRUE,  kTRUE,  kTRUE  */  /*  kTRUE, kFALSE, kFALSE */      /*kTRUE,   kTRUE,     */  /*kTRUE,  kTRUE ,  kTRUE,   kTRUE */    /* kTRUE  ,   kTRUE */   /*kTRUE , kFALSE,*/ /*kTRUE  */      kTRUE,     kFALSE,         kTRUE    };
bool useRICHPID[]           = { /* kTRUE,  kFALSE,   kTRUE,   kTRUE,  kFALSE,   kTRUE */ /*  kTRUE,  kTRUE,  kTRUE */ /*  kTRUE,  kTRUE,  kTRUE,  kTRUE,  kTRUE */ /* kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE  */  /*  kTRUE,  kTRUE,  kTRUE,  kTRUE,  kTRUE,  kTRUE  */  /* kFALSE,  kTRUE, kFALSE */      /*kFALSE,  kFALSE,    */  /*kFALSE, kTRUE ,  kTRUE,   kTRUE */    /* kFALSE ,  kFALSE */   /*kFALSE, kTRUE, */ /*kTRUE  */     kFALSE,     kFALSE,        kFALSE    };
bool usePreShPID[]          = { /*kFALSE,   kTRUE,   kTRUE,  kFALSE,   kTRUE,   kTRUE */ /* kFALSE, kFALSE, kFALSE */ /* kFALSE, kFALSE, kFALSE, kFALSE, kFALSE */ /* kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE  */  /* kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE  */  /* kFALSE, kFALSE, kFALSE */      /*kFALSE,  kFALSE,    */  /*kFALSE, kFALSE, kFALSE,  kFALSE */    /* kFALSE ,  kFALSE */   /*kFALSE, kFALSE,*/ /*kFALSE */     kFALSE,      kTRUE,         kTRUE    };
double nSigmaTOFEle[]       = { /*   3.0,     3.0,     3.0,     3.0,     3.0,     3.0 */ /*    3.0,    3.0,    3.0 */ /*    3.0,    3.0,    3.0,    3.0,    3.0 */ /*    3.0,     3.0,     3.0,     3.0,     3.0,     3.0  */  /*    3.0,    3.0,    3.0,    3.0,    3.0,    3.0  */  /*    3.0,    3.0,    3.0 */      /*3.0,     3.0,       */  /*  3.0,    3.0 ,    3.0,     3.0 */    /*  3.0   ,    3.0  */   /*  3.0 ,   3.0, */ /*  3.0  */        3.0,        3.0,           3.0    };
double nSigmaTOFPi[]        = { /*   3.0,     3.0,     3.0,     3.0,     3.0,     3.0 */ /*    3.0,    3.0,    3.0 */ /*    3.0,    3.0,    3.0,    3.0,    3.0 */ /*    3.0,     3.0,     3.0,     3.0,     3.0,     3.0  */  /*    3.0,    3.0,    3.0,    3.0,    3.0,    3.0  */  /*    3.0,    3.0,    3.0 */      /*3.0,     3.0,       */  /*  3.0,    3.0 ,    3.0,     3.0 */    /*  3.0   ,    3.0  */   /*  3.0 ,   3.0, */ /*  3.0  */        3.0,        3.0,           3.0    };
double nSigmaRICHEle[]      = { /*   3.0,     3.0,     3.0,     3.0,     3.0,     3.0 */ /*    3.0,    3.0,    3.0 */ /*    3.0,    3.0,    3.0,    3.0,    3.0 */ /*    3.0,     3.0,     3.0,     3.0,     3.0,     3.0  */  /*    3.0,    3.0,    3.0,    3.0,    3.0,    3.0  */  /*    3.0,    3.0,    3.0 */      /*3.0,     3.0,       */  /*  3.0,    3.0 ,    3.0,     3.0 */    /*  0.0   ,    0.0  */   /*  3.0 ,   3.0, */ /*  3.0  */        3.0,        3.0,           3.0    };
double nSigmaRICHPi[]       = { /*   4.0,     4.0,     4.0,     4.0,     4.0,     4.0 */ /*    4.0,    4.0,    4.0 */ /*    4.0,    4.0,    4.0,    4.0,    4.0 */ /*    4.0,     4.0,     4.0,     4.0,     4.0,     4.0  */  /*    4.0,    4.0,    4.0,    4.0,    4.0,    4.0  */  /*    4.0,    4.0,    4.0 */      /*4.0,     4.0,       */  /*  3.0,    3.0 ,    3.5,     4.0 */    /*  99.   ,    99.  */   /*  4.0 ,   4.0, */ /*  4.0  */        4.0,        4.0,           4.0    };
double PtCut02[]            = { /*  0.08,    0.08,    0.08,    0.08,    0.08,    0.08 */ /*   0.08,   0.08,   0.08 */ /*   0.08,   0.08,   0.08,   0.08,   0.08 */ /*   0.08,    0.08,    0.08,    0.08,    0.08,    0.08  */  /*   0.08,   0.08,   0.08,   0.08,   0.08,   0.08  */  /*   0.08,   0.08,   0.08 */      /*0.04,    0.08,      */  /*  0.03,   0.03,   0.03,     0.3 */    /*  0.04  ,   0.0   */   /*  0.08,   0.08,*/ /*  0.08 */       0.08,       0.08,          0.08    };
double PtCut05[]            = { /*   0.2,     0.2,     0.2,     0.2,     0.2,     0.2 */ /*    0.2,    0.2,    0.2 */ /*    0.2,    0.2,    0.2,    0.2,    0.2 */ /*    0.2,     0.2,     0.2,     0.2,     0.2,     0.2  */  /*    0.2,    0.2,    0.2,    0.2,    0.2,    0.2  */  /*    0.2,    0.2,    0.2 */      /*0.2,     0.08,      */  /*  0.03,   0.03,   0.03,     0.3 */    /*  0.08  ,   0.0   */   /*  0.2 ,   0.2, */ /*  0.2  */        0.2,        0.2,           0.2    };
double PtCut20[]            = { /*   0.3,     0.3,     0.3,     0.3,     0.3,     0.3 */ /*    0.3,    0.3,    0.3 */ /*    0.3,    0.3,    0.3,    0.3,    0.3 */ /*    0.3,     0.3,     0.3,     0.3,     0.3,     0.3  */  /*    0.3,    0.3,    0.3,    0.3,    0.3,    0.3  */  /*    0.3,    0.3,    0.3 */      /*0.2,     0.08,      */  /*  0.03,   0.03,   0.03,     0.3 */    /*  0.08  ,   0.0   */   /*  0.2 ,   0.2, */ /*  0.3  */        0.3,        0.3,           0.3    };
double EtaCut[]             = { /*   0.8,     0.8,     0.8,     0.8,     0.8,     0.8 */ /*   1.75,   1.75,   1.75 */ /*    0.8,    0.8,    0.8,    0.8,    0.8 */ /*    0.8,    1.25,    1.75,     2.5,     3.0,     4.0  */  /*    0.8,    0.8,    0.8,    0.8,    0.8,    0.8  */  /*    0.8,    0.8,    0.8 */                                                                                                                                        0.8,        0.8,           0.8    };
bool busePFpid[]            = { /*kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE,  kFALSE */                                                                                                                                                                                                                                                                                                                                                               kFALSE,     kFALSE,        kFALSE    };
bool usePFtrackVec[]        = { /* kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE,   kTRUE */ /* kFALSE, kFALSE,  kTRUE */ /* kFALSE, kFALSE,  kTRUE,  kTRUE,  kTRUE */                                                             /* kFALSE, kFALSE,  kTRUE,  kTRUE,  kTRUE,  kTRUE  */  /* kFALSE, kFALSE, kFALSE */                                                                                                                                     kFALSE,     kFALSE,        kFALSE    };
double pf_pt_cut[]          = { /*  0.02,    0.02,    0.02,    0.02,    0.02,    0.02 */ /*   0.08,   0.08,   0.02 */ /*   0.08,   0.08,   0.02,   0.02,   0.02 */ /*   0.08,    0.08,    0.08,    0.08,    0.08,    0.08  */  /*   0.00,   0.08,   0.02,   0.02,   0.02,   0.02  */  /*   0.00,   0.00,   0.00 */              /* use small pt cut first in rising order */                                                                             0.02,       0.02,          0.02    };
double pf_opAngle_cut[]     = { /*   0.1,     0.1,     0.1,     0.1,     0.1,     0.1 */ /*    0.1,    0.1,    0.1 */ /*    0.1,    0.1,   0.08,    0.1,   0.12 */                                                             /*   0.00,   0.08,   0.08,    0.1,   0.12,   0.12  */  /*   0.00,   0.00,   0.00 */                                                                                                                                        0.1,        0.1,           0.1    };
double pf_mass_cut[]        = { /*  0.05,    0.05,    0.05,    0.05,    0.05,    0.05 */ /*   0.05,   0.05,   0.05 */ /*   0.05,   0.05,   0.05,   0.05,   0.05 */                                                             /*   0.00,   0.05,   0.05,   0.05,   0.05,   0.06  */  /*   0.00,   0.00,   0.00 */                                                                                                                                       0.05,       0.05,          0.05    };
double ps_EleAccep_pt_cut[] = { /*   0.0,     0.4,     1.5,     0.0,     0.4,     1.5 */                                                                                                                                                                                                                                                                                                                                                                  0.0,        0.0,           0.0    };

// TOF pte > 0.04 B = 0.2 T (highest priority) or TOF pte > 0.08 B = 0.5 T
// TOF RICH pte > 0.2 B = 0.5 T (highest priority) or TOF RICH pte > 0.08  B = 0.5 T

// Eta cut (pseudorapidity)
// double EtaCut = 1.1;





bool doPID(Track * track, bool useTOF, bool useRICH, double p_tofMaxAcc, double p_tofPionRej, double p_richPionRej, double nSigmaTOFele, double nSigmaTOFpi, double nSigmaRICHele, double nSigmaRICHpi, std::array<float, 5> PIDnsigmaTOF, std::array<float, 5> PIDnsigmaRICH){
  double p = track->P;

  bool TOFpid = kFALSE;
  bool RICHpid = kFALSE;
  //TOF PID
  if(useTOF && (toflayer.hasTOF(*track)) && (p < p_tofMaxAcc)) {
    if (p > p_tofPionRej && richdetector.hasRICH(*track)) {
      if((fabs(PIDnsigmaTOF[0]) < nSigmaTOFele) && (fabs(PIDnsigmaRICH[0]) < nSigmaRICHEle_forTOFPID))
        TOFpid = true; // is within 3 sigma of the electron band (TOF)
    } else if(p <= p_tofPionRej){
      if(fabs(PIDnsigmaTOF[0]) < nSigmaTOFele)
        TOFpid = true; // is within 3 sigma of the electron band (TOF)
    }
    else TOFpid = false; // This is rejecting all the heavier particles which do not have a RICH signal in the pt area of 0.4-0.6 GeV/c


    if(fabs(PIDnsigmaTOF[2]) < nSigmaTOFpi)
      TOFpid = false; // is within 3 sigma of the pion band (TOF)
  }
  else TOFpid = false;

  //RICH PID
  if(useRICH && richdetector.hasRICH(*track)) {
    if(fabs(PIDnsigmaRICH[0]) < nSigmaRICHele) RICHpid = true; // is within 3 sigma of the electron band (RICH)
    if((fabs(PIDnsigmaRICH[2]) < nSigmaRICHpi) && (p > p_richPionRej) ) RICHpid = false; // is within 3 sigma of the pion band (RICH)
  }
  else RICHpid = false;


  return RICHpid || TOFpid;

  // if (!(RICHpid || TOFpid) /*&& !PreShpid*/) return false; // check if TOF or RICH signal is true.
  // else return true;
  // ################## end of PID selection ##################
}

bool doTOFPID(Track * track, bool useTOF, double p_tofMaxAcc, double p_tofPionRej, double nSigmaTOFele, double nSigmaTOFpi, std::array<float, 5> PIDnsigmaTOF){
  double p = track->P;
  bool TOFpid = false;
  //TOF PID
  if(useTOF && (toflayer.hasTOF(*track)) /*&& (p < p_tofMaxAcc)*/) {
    if(fabs(PIDnsigmaTOF[0]) < nSigmaTOFele)
      TOFpid = true; // is within 3 sigma of the electron band (TOF)
    else TOFpid = false;


    if(fabs(PIDnsigmaTOF[2]) < nSigmaTOFpi)
      TOFpid = false; // is within 3 sigma of the pion band (TOF)
  }
  else TOFpid = false;

  return TOFpid;
  // ################## end of PID selection ##################
}

bool hasStrangeAncestor(GenParticle *particle, TClonesArray *particles)
{
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle *)particles->At(imother);
  auto pid = mother->PID;

  switch (abs(pid)) {
  case 310:  // K0s
  case 3122: // Lambda
  case 3112: // Sigma-
  case 3222: // Sigma+
  case 3312: // Xi-
  case 3322: // Xi0  was not included by Roberto in his example
  case 3334: // Omega-
    return true;
  }

  return hasStrangeAncestor(mother, particles);
}


// eta cut for tracks
bool etaCut(Track *tr, Int_t iSce){
  // evaluate as true if criterion is passed
  bool eta = abs(tr->Eta) < EtaCut[iSce];
  return (eta);
}
// eta cut for particles
bool etaCut(GenParticle *pa, Int_t iSce){
  // evaluate as true if criterion is passed
  bool eta = abs(pa->Eta) < EtaCut[iSce];
  return (eta);
}
// eta cut for LorentzVector
bool etaCut(TLorentzVector LV, Int_t iSce){
  // evaluate as true if criterion is passed
  bool eta = abs(LV.Eta()) < EtaCut[iSce];
  return (eta);
}

// calculate rapidity
double etatorap(double pt, double eta, double mass){
  double nominator = TMath::Sqrt(mass*mass + pt*pt * TMath::CosH(eta)*TMath::CosH(eta)) + pt*TMath::SinH(eta);
  double denominator = TMath::Sqrt(mass*mass + pt*pt);
  double y = TMath::Log(nominator/denominator);
  return y;
}

// kinematic cuts for tracks
bool kineCuts(Track *tr, Int_t iSce){
  // check pt and eta for track
  // evaluate as true if criterion is passed
  bool pt;
  if(Bz == 2) pt = tr->PT > PtCut02[iSce];
  else if(Bz == 5) pt = tr->PT > PtCut05[iSce];
  else if(Bz == 20) pt = tr->PT > PtCut20[iSce];
  bool eta = abs(tr->Eta) < EtaCut[iSce];
  // all have to be true
  return (pt && eta);
}

// kinematic cuts for generated particles
bool kineCuts(GenParticle *pa, Int_t iSce){
  // check pt and eta for particle
  // evaluate as true if criterion is passed
  bool pt;
  if(Bz == 2) pt = pa->PT > PtCut02[iSce];
  else if(Bz == 5) pt = pa->PT > PtCut05[iSce];
  else if(Bz == 20) pt = pa->PT > PtCut20[iSce];
  bool eta = abs(pa->Eta) < EtaCut[iSce];
  // all have to be true
  return (pt && eta);
}

// kinematic cuts for generated smeared particles
bool kineCuts(TLorentzVector LV, Int_t iSce){
  // check pt and eta for particle
  // evaluate as true if criterion is passed
  bool pt;
  if(Bz == 2) pt = LV.Pt() > PtCut02[iSce];
  else if(Bz == 5) pt = LV.Pt() > PtCut05[iSce];
  else if(Bz == 20) pt = LV.Pt() > PtCut20[iSce];
  bool eta = abs(LV.Eta()) < EtaCut[iSce];
  // all have to be true
  return (pt && eta);
}

bool hasHeavyAncestor(GenParticle *particle, TClonesArray *particles){
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle *)particles->At(imother);
  auto pid = 0;
  if(mother) pid = mother->PID;

  if (((abs(pid)>=400) && (abs(pid)<=439)) || ((abs(pid)>=4000) && (abs(pid)<=4399))
   || ((abs(pid)>=500) && (abs(pid)<=549)) || ((abs(pid)>=5000) && (abs(pid)<=5499))) {
        return true;
  }
  // switch (abs(pid)) {
  // case 411:  // D+
  // case 421:  // D0
  // case 431:  // Ds+
  // case 4122: // Lambdac+
  // case 4132: // Xic0
  // case 4232: // Xic+
  // case 511:  // B0
  // case 521:  // B+
  // case 531:  // Bs0
  // case 541:  // Bc+
  // case 5122: // Lambdab0
  // case 5132: // Xib-
  // case 5232: // Xib0
  // case 5332: // Omegab-
  //   return true;
  // }
  return hasHeavyAncestor(mother, particles);
}

bool isBeauty(int pid){
  if (((abs(pid)>=500) && (abs(pid)<=549)) || ((abs(pid)>=5000) && (abs(pid)<=5499))) {
        return true;
  }
  // switch (abs(pid)) {
  // case 511:  // B0
  // case 521:  // B+
  // case 531:  // Bs0
  // case 541:  // Bc+
  // case 5122: // Lambdab0
  // case 5132: // Xib-
  // case 5232: // Xib0
  // case 5332: // Omegab-
  //   return true;
  // }
  return false;
}

bool hasBeautyAncestor(GenParticle *particle, TClonesArray *particles){
  //check until beauty ancestor is found
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle *)particles->At(imother);
  auto pid = mother->PID;
  if (isBeauty(pid)) return true;
  return hasBeautyAncestor(mother, particles);
}

bool isCharm(int pdg) {
  if (((abs(pdg)>=400) && (abs(pdg)<=439)) || ((abs(pdg)>=4000) && (abs(pdg)<=4399))) {
        return true;
  }
  // switch (abs(pdg)) {
  // case 411:  // D+
  // case 421:  // D0
  // case 431:  // Ds+
  // case 4122: // Lambdac+
  // case 4132: // Xic0
  // case 4232: // Xic+
  //   return true;
  // }
  return false;
}

bool hasCharmAncestor(GenParticle *particle, TClonesArray *particles)
{
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle *)particles->At(imother);
  auto pid = mother->PID;
  auto bCharm = false;
  bCharm = isCharm(pid);

  if (hasBeautyAncestor(mother, particles)) {bCharm = false;}
  return bCharm;
}


Bool_t IsStable(Int_t pdg, Bool_t &ischarged)
{
  //
  // Decide whether particle (pdg) is stable
  //


  // All ions/nucleons are considered as stable
  // Nuclear code is 10LZZZAAAI
  if(pdg>1000000000)return kTRUE;

  const Int_t kNstable = 19;
  Int_t i;

  Int_t pdgStable[kNstable] = {
     22,             // Photon
     11,          // Electron
     -13,          // Muon
     211,            // Pion
     321,             // Kaon
     310,           // K0s
     130,            // K0l
     2212,            // Proton
     2112,           // Neutron
     3122,           // Lambda_0
     3212,         // Sigma0
     3112,        // Sigma Minus
     3222,         // Sigma Plus
     3312,               // Xsi Minus
     3322,               // Xsi
     3334,               // Omega
     12,               // Electron Neutrino
     14,              // Muon Neutrino
     16              // Tau Neutrino
  };

  const Int_t kNotCharge = 11;

  Int_t pdgNotCharge[kNotCharge] = {
    22,             // Photon
    310,           // K0s
    130,            // K0l
    2112,           // Neutron
    3122,           // Lambda_0
    3212,         // Sigma0
    3322,               // Xsi
    3334,               // Omega
    12,               // Electron Neutrino
    14,              // Muon Neutrino
    16              // Tau Neutrino
 };

  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstable; i++) {
   if (TMath::Abs(pdg) == TMath::Abs(pdgStable[i])) {
     isStable = kTRUE;
     break;
   }
  }

  ischarged = kTRUE;
 for (i = 0; i < kNotCharge; i++) {
  if (TMath::Abs(pdg) == TMath::Abs(pdgNotCharge[i])) {
    ischarged = kFALSE;
    break;
  }
 }

  return isStable;
}

// function to get the corresponding weight for a track
Double_t GetWeight(TH1* h, Double_t pt){
  Int_t iBin = h->GetXaxis()->FindBin(pt);
  Double_t weight = 1.;
  weight = h->GetBinContent(iBin);
  return weight;
}




//_____________________________________________________________________________________________
TLorentzVector  ApplySmearing(TObjArray* fArrResoPt, TObjArray* fArrResoEta, TObjArray* fArrResoPhi_Pos, TObjArray* fArrResoPhi_Neg, const TLorentzVector& vec, Double_t ch)
{
  //
  // Smearing in pt, eta, phi: Run 2 method
  //

  if(!fArrResoPt || !fArrResoEta || !fArrResoPhi_Pos || !fArrResoPhi_Neg){
    printf("No smearing of generated tracks for efficiencies\n");
    return vec;
  }

  //Double_t theta, phi, pt, p, mass, eta;
  Double_t phi, pt, mass, eta;
  TLorentzVector resvec;

  mass = eMass;
  pt = vec.Pt();
  phi = vec.Phi();
  eta = vec.Eta();

  // smear pt
  Int_t ptbin     = ((TH2D *)(fArrResoPt->At(0)))->GetXaxis()->FindBin(pt);
  Int_t ptbin_max = ((TH2D *)(fArrResoPt->At(0)))->GetXaxis()->GetNbins();
  // make sure that no underflow or overflow bins are used
  if (ptbin < 1)
    ptbin = 1;
  else if (ptbin > ptbin_max)
    ptbin = ptbin_max;

  Double_t sPt = pt;
  if (((TH1D *)(fArrResoPt->At(ptbin)))->Integral() > 0.000001){
    Double_t smearing = ((TH1D *)(fArrResoPt->At(ptbin)))->GetRandom() * pt;
    sPt = pt - smearing;
  }

  // smear eta
  ptbin     = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->FindBin(pt);
  ptbin_max = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->GetNbins();
  if (ptbin < 1)
    ptbin = 1;
  else if (ptbin > ptbin_max)
    ptbin = ptbin_max;
  Double_t sEta = eta;
  if (((TH1D *)(fArrResoEta->At(ptbin)))->Integral() > 0.000001){
    Double_t smearing = ((TH1D *)(fArrResoEta->At(ptbin)))->GetRandom();
    sEta = eta - smearing;
  }

  // smear phi
  ptbin     = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->FindBin(pt);
  ptbin_max = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->GetNbins();
  if (ptbin < 1)
    ptbin = 1;
  if (ptbin > ptbin_max)
    ptbin = ptbin_max;

  Double_t sPhi = phi;
  if (ch > 0) {
    if (((TH1D *)(fArrResoPhi_Pos->At(ptbin)))->Integral() > 0.000001){
      Double_t smearing = ((TH1D *)(fArrResoPhi_Pos->At(ptbin)))->GetRandom();
      sPhi = phi - smearing;
    }
  } else if (ch < 0) {
    if (((TH1D *)(fArrResoPhi_Neg->At(ptbin)))->Integral() > 0.000001){
      Double_t smearing = ((TH1D *)(fArrResoPhi_Neg->At(ptbin)))->GetRandom();
      sPhi = phi - smearing;
    }
  }

  // printf(" Original Pt = %f Phi %f Eta %f -> final pt = %f Phi %f Eta %f
  // \n",pt,phi,eta,sPt,sPhi,sEta);

  const Double_t sPx = sPt * cos(sPhi);
  const Double_t sPy = sPt * sin(sPhi);
  const Double_t sPz = sPt * sinh(sEta);
  const Double_t sP = sPt * cosh(sEta);
  const Double_t sE = sqrt(sP * sP + mass * mass);

  resvec.SetPxPyPzE(sPx, sPy, sPz, sE);
  return resvec;
}



// Pair histograms ULS
std::string title2d = ";#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c})";
std::string title3d = ";#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});DCA_{ee} (mm)";




//#####################################
//#####################################
//         begin of analysis
//#####################################
//#####################################

void anaEEstudy(
    const char *inputFile = "delphes100k.root", // one fo the LF and charm part. Has 1M charm+beauty events
    const char *outputFile = "anaEEstudyLFcc.root" // merge output files after analysis was run to keep file size moderate
  )
{
  cout << "    Bz = " << Bz  << " kG"<< endl;
  // Stopwatch
  TStopwatch* watch = new TStopwatch;
  watch->Start();

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  auto treeReader = new ExRootTreeReader(&chain);
  auto numberOfEntries = treeReader->GetEntries();
  auto numberOfEntriesTree = numberOfEntries;


// binning
  Int_t nPtBins  = sizeof(pt_binning)/sizeof(*pt_binning) -1;
  Double_t eta_binning[nEtaBins+1];
  Double_t phi_binning[nPhiBins+1];
  for (int iBin = 0; iBin < nEtaBins+1; iBin++) eta_binning[iBin] = -maxEta + iBin * 2.* maxEta/nEtaBins;
  for (int iBin = 0; iBin < nPhiBins+1; iBin++) phi_binning[iBin] = 0. + iBin * maxPhi/nPhiBins;


  // // ptee
  // for(Int_t i=0  ;i<40   ;i++) {
  //   ptee_bin_c[i] = 0.005 * i;        //from 0 to 0.2 GeV/c, every 0.005 GeV/c
  //   //printf("bin %f for %d\n",ptee_bin_c[i],i);
  // }
  // for(Int_t i=40  ;i<=520   ;i++) {
  //   ptee_bin_c[i] = 0.01 * (i-  40) +  0.2;        //from 0.2 to 5 GeV/c, every 0.01 GeV/c
  //   //printf("bin %f for %d\n",ptee_bin_c[i],i);
  // }
  // for(Int_t i=0  ;i<2   ;i++) {
  //   mee_bin_c[i] = 4. *i;//from 0 to 2. GeV/c, every 0.02 GeV/c
  // }
  //
  // // DCA
  // for(int k = 0 ; k < 2; k++){
  //   dca_bin_c[k] = k*10.;
  // }

  // ptee
  for(Int_t i=0  ;i<5   ;i++) {
    ptee_bin_c[i] = 0.02 * i;        //from 0 to 0.1 GeV/c, every 0.02 GeV/c
    //printf("bin %f for %d\n",ptee_bin_c[i],i);
  }
  for(Int_t i=5  ;i<=105   ;i++) {
    ptee_bin_c[i] = 0.1 * (i-  5) +  0.1;        //from 0.1 to 10 GeV/c, every 0.1 GeV/c
    //printf("bin %f for %d\n",ptee_bin_c[i],i);
  }
  // mee
  for(Int_t i=0  ;i<100   ;i++) {
    mee_bin_c[i] = 0.02 * (i-  0) +  0.0;//from 0 to 2. GeV/c, every 0.02 GeV/c
  }
  for(Int_t i=100  ;i<=120   ;i++) {
    mee_bin_c[i] = 0.1 * (i-  100) +  2.0;//from 2 to 4. GeV/c, every 0.1 GeV/c
  }
  // DCA
  for(int k = 0 ; k < 101; k++){
    dca_bin_c[k] = k*0.1;
  }


  // ptee
  // for(Int_t i=0  ;i<4   ;i++) {
  //   ptee_bin_c[i] = 0.025 * (i-  0) +  0.0;//from 0 to 0.075 GeV/c, every 0.025 GeV/c
  //   //printf("bin %f for %d\n",ptee_bin_c[i],i);
  // }
  // for(Int_t i=4 ;i<=202  ;i++) {
  //   ptee_bin_c[i] = 0.05  * (i- 4) +  0.1;//from 0.1 to 10 GeV/c, evety 0.05 GeV/c
  //   //printf("bin %f for %d\n",ptee_bin_c[i],i);
  // }

  // for(Int_t i=0  ;i<10   ;i++) {
  //   ptee_bin_c[i] = 0.01 * (i-  0) +  0.0;//from 0 to 0.1 GeV/c, every 0.01 GeV/c
  //   //printf("bin %f for %d\n",ptee_bin_c[i],i);
  // }
  // for(Int_t i=10 ;i<20  ;i++) {
  //   ptee_bin_c[i] = 0.02  * (i- 10) +  0.1;//from 0.1 to 0.3 GeV/c, evety 0.02 GeV/c
  //   //printf("bin %f for %d\n",ptee_bin_c[i],i);
  // }
  // for(Int_t i=20 ;i<215  ;i++) {
  //   ptee_bin_c[i] = 0.05  * (i- 20) +  0.3;//from 0.3 to 10 GeV/c, evety 0.05 GeV/c
  //   //printf("bin %f for %d\n",ptee_bin_c[i],i);
  // }
  //
  // // mee
  // for(Int_t k = 0 ; k < 401; k++){
  //   mee_bin_c[k] = k*0.01; // 4./400.
  // }
  // // DCA
  // for(int k = 0 ; k < 201; k++){
  //   dca_bin_c[k] = k*0.1;
  // }


  // Get pointers to branches used in this analysis
  auto events = treeReader->UseBranch("Event");
  auto tracks = treeReader->UseBranch("Track");
  // auto innertracks = treeReader->UseBranch("InnerTrack");
  particles = treeReader->UseBranch("Particle");


  // Smearing for generated tracks
  smearingfile = "resolution.root";
  TFile *fRes = TFile::Open(smearingfile.Data());
  // ReadResoFile(fRes);

  // Set resolutions: default Run 2
  TObjArray *fArrResoPt = 0x0;
  TObjArray *fArrResoEta = 0x0;
  TObjArray *fArrResoPhi_Pos = 0x0;
  TObjArray *fArrResoPhi_Neg = 0x0;
  if(fRes && fRes->IsOpen()){
    // Oton
    fArrResoPt = (TObjArray *)fRes->Get("RelPtResArrCocktail");
    fArrResoEta = (TObjArray *)fRes->Get("EtaResArrVsPt");
    fArrResoPhi_Pos = (TObjArray *)fRes->Get("PhiPosResArrVsPt");
    fArrResoPhi_Neg = (TObjArray *)fRes->Get("PhiEleResArrVsPt");
  }
  else {
    printf("no resolution file given!");
  }


  // read file with weighthing corrections for leight and heavy flavor TrackSmearer
  TFile* fOutputWeightsHF = TFile::Open("hfe_weights.root", "READ");
  hweightHFtoe_lowerPt = (TH1F*) fOutputWeightsHF->Get("hftoe_weights_cocktail_data_ratio");   // weights for tracks with pt below 0.6 GeV/c
  hweightHFtoe_higherPt = (TH1F*) fOutputWeightsHF->Get("hftoe_weights_paper_data_ratio");  // weights for tracks with pt equal higher than 0.6 GeV/c

  TFile* fOutputWeightsLF = TFile::Open("lfe_weights.root", "READ");
  hweightLFtoe = (TH1F*) fOutputWeightsLF->Get("lftoe_weights_paper_data_ratio");


  // // inner TOF layer
  // innertoflayer.setup(inner_tof_radius, inner_tof_length, inner_tof_sigmat,inner_tof_sigma0);

  // outer TOF layer
  toflayer.setup(tof_radius, tof_length, tof_sigmat,tof_sigma0);

  // RICH detector
  richdetector.setup(rich_radius, rich_length);
  richdetector.setIndex(1.03);
  richdetector.setRadiatorLength(2.);
  richdetector.setEfficiency(0.4);
  richdetector.setSigma(7.e-3);

  // PreShower detector
  preshower.setup();

  bool bUseTOFPID = kFALSE;
  bool bUseRICHPID = kFALSE;
  if (std::all_of(std::begin(useTOFPID), std::end(useTOFPID), [](bool i){return i;}))  bUseTOFPID = kTRUE;
  if (std::all_of(std::begin(useRICHPID), std::end(useRICHPID), [](bool i){return i;}))  bUseRICHPID = kTRUE;
  cout << "    bUseTOFPID = " << bUseTOFPID << endl;
  cout << "    bUseRICHPID = " << bUseRICHPID << endl;

  // smearer
  // look up tables are selected in generateEfficiencies.sh script
  o2::delphes::TrackSmearer smearer;
  smearer.loadTable(11, "./lutCovm.el.dat");
  smearer.loadTable(13, "./lutCovm.mu.dat");
  smearer.loadTable(211, "./lutCovm.pi.dat");
  smearer.loadTable(321, "./lutCovm.ka.dat");
  smearer.loadTable(2212, "./lutCovm.pr.dat");

  smearer.useEfficiency(bUseTrackEff);



  int nPIDscenarios = sizeof(nSigmaTOFEle) / sizeof(nSigmaTOFEle[0]);
  TList* listGenTracks = new TList();
  TList* listRecTracks = new TList();
  listGenTracks->SetOwner();
  listRecTracks->SetOwner();
  listGenTracks->SetName("GenTracks");
  listRecTracks->SetName("RecTracks");


  // Event histograms
  auto nTracks = new TH1F("nTracks",";Tracks",10000,0,10000);
  auto nParticles = new TH1F("nParticles",";Particles",1000,0,1000);
  auto nTracksGen = new TH1F("nTracksGen",";Tracks",8000,0,16000);
  auto nTracksCent = new TH1F("nTracksCent",";Tracks",10000,0,10000);
  auto nParticlesFIT = new TH1F("nParticlesFIT",";Particles",15000,0,15000);
  auto nParticlesFITCent = new TH1F("nParticlesFITCent",";Particles",15000,0,15000);
  auto nParticlesMidRapidity = new TH1F("nParticlesMidRapidity",";Particles",15000,0,15000);
  auto nParticlesMidRapidityCent = new TH1F("nParticlesMidRapidityCent",";Particles",15000,0,15000);

  // Get number of Tracks selected by centrality within mid rapidity range
  auto hdNdeta_midrap_gen = new TH1F("hdNdeta_midrap_gen",";#eta",1,-0.5,0.5);
  auto hdNdeta_midrap_rec = new TH1F("hdNdeta_midrap_rec",";#eta",1,-0.5,0.5);

  // Check pt of tracks with and without hasTOF(), hasRICH() signals
  auto hTrackP = new TH1F("hTrackP",";#it{p} (GeV/#it{c})",200,0.,1.);
  auto hTrackP_hasTOForRICH = new TH1F("hTrackP_hasTOForRICH",";#it{p} (GeV/#it{c})",200,0.,0.2);
  auto hTrackPt = new TH1F("hTrackPt",";#it{p}_{T} (GeV/#it{c})",200,0.,1.);
  auto hTrackPt_hasTOForRICH = new TH1F("hTrackPt_hasTOForRICH",";#it{p}_{T} (GeV/#it{c})",200,0.,0.2);




  // Check smearing of generated tracks
  auto hSmearing_For_Eff_pt = new TH2F("hSmearing_For_Eff_pt","(#it{p}_{T,gen}-#it{p}_{T,rec})/#it{p}_{T,gen};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.5,0.5);
  auto hSmearing_For_Eff_phi_pos = new TH2F("hSmearing_For_Eff_phi_pos","#phi_{gen}-#phi_{rec};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.4,0.1);
  auto hSmearing_For_Eff_phi_neg = new TH2F("hSmearing_For_Eff_phi_neg","#phi_{gen}-#phi_{rec};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.1,0.4);
  auto hSmearing_For_Eff_eta = new TH2F("hSmearing_For_Eff_eta","#eta_{gen}-#eta_{rec};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.1,0.1);

  // Track histograms
  auto hBeforeSmearing_Pt_Eta_Phi_rec = new TH3F("hBeforeSmearing_Pt_Eta_Phi_rec", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  auto hAfterSmearing_Pt_Eta_Phi_rec = new TH3F("hAfterSmearing_Pt_Eta_Phi_rec", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hAfterKineCuts_Pt_Eta_Phi_rec[nPIDscenarios];
  TH3F* hTrack_ElePos_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hNegTrack_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hPosTrack_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Muon_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Pion_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Kaon_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Proton_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hAllTracks_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_Muon_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_Pion_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_Kaon_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_Proton_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_ElePos_Rec_Pt_Eta_Phi_afterPSPID[nPIDscenarios];
  TH3F* hTrack_Pion_Rec_Pt_Eta_Phi_afterPSPID[nPIDscenarios];
  TH3F* hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_beforePID[nPIDscenarios];
  TH3F* hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_afterPID[nPIDscenarios];
  TH3F* hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_afterPID[nPIDscenarios];
  TH3F* hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_afterPID[nPIDscenarios];
  TH2F* hMeeOpAngle_RecULS_bPF_pi0Dalitz[nPIDscenarios];
  TH2F* hMeeOpAngle_RecULS_aPF_pi0Dalitz[nPIDscenarios];
  TH2F* hMeeOpAngle_RecULS_bPF_bkg[nPIDscenarios];
  TH2F* hMeeOpAngle_RecULS_aPF_bkg[nPIDscenarios];
  TH2F* hMeeOpAngle_RecULS_bPF_HFPi0[nPIDscenarios];
  TH2F* hMeeOpAngle_RecULS_aPF_HFPi0[nPIDscenarios];
  TH1F* hTracks_ElePos_Gen_inPFvector_Pt[nPIDscenarios];
  TH1F* hTracks_ElePos_Rec_inPFvector_Pt[nPIDscenarios];

  for (int j = 0; j < nPIDscenarios; ++j){
    hAfterKineCuts_Pt_Eta_Phi_rec[j] = new TH3F(Form("hAfterKineCuts_Pt_Eta_Phi_rec_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hNegTrack_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hNegTrack_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hPosTrack_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hPosTrack_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hAllTracks_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hAllTracks_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Muon_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Muon_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Pion_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pion_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Kaon_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Kaon_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Proton_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Proton_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_ElePos_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Muon_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_Muon_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Pion_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_Pion_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Kaon_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_Kaon_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Proton_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_Proton_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_Rec_Pt_Eta_Phi_afterPSPID[j]  = new TH3F(Form("hTrack_ElePos_Rec_Pt_Eta_Phi_afterPSPID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Pion_Rec_Pt_Eta_Phi_afterPSPID[j]  = new TH3F(Form("hTrack_Pion_Rec_Pt_Eta_Phi_afterPSPID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_beforePID[j]  = new TH3F(Form("hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_beforePID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_afterPID[j] = new TH3F(Form("hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_afterPID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_afterPID[j] = new TH3F(Form("hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_afterPID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_afterPID[j] = new TH3F(Form("hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_afterPID_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTracks_ElePos_Rec_inPFvector_Pt[j] = new TH1F(Form("hTracks_ElePos_Rec_inPFvector_Pt_sce%i",j+1),";#it{p}_{T} (GeV/#it{c})",n_ptee_bin_c,ptee_bin_c);
  }



  // Pair histograms
  std::vector<TH2F*> vecPFhistos[nPIDscenarios];
  for (int j = 0; j < nPIDscenarios; ++j){
  // Pair histograms ULS
  hMPtDCA_ULS_gen[j]         = new TH3F(Form("hMPtDCA_ULS_gen_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec[j]         = new TH3F(Form("hMPtDCA_ULS_rec_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle[j]         = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle_charmTOe[j]  = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_charmTOe_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle_beautyTOe[j]  = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_beautyTOe_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle_hfTOe[j]   = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_hfTOe_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle_lfTOee[j] = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_lfTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle_ccTOee[j] = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_ccTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle_bbTOee[j] = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_bbTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_ULS_rec_MCpidEle_hfTOee[j] = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_hfTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  // Pair histograms LS
  hMPtDCA_LS_gen[j]         = new TH3F(Form("hMPtDCA_LS_gen_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec[j]         = new TH3F(Form("hMPtDCA_LS_rec_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_MCpidEle[j]         = new TH3F(Form("hMPtDCA_LS_rec_MCpidEle_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  // hMPtDCA_LS_rec_misIDoneLeg[j]      = new TH3F(Form("hMPtDCA_LS_rec_misIDoneLeg_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  // hMPtDCA_LS_rec_misIDtwoLeg[j]      = new TH3F(Form("hMPtDCA_LS_rec_misIDtwoLeg_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  // hMPtDCA_LS_rec_misIDPion[j]        = new TH3F(Form("hMPtDCA_LS_rec_misIDPion_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  // hMPtDCA_LS_rec_misIDhf[j]          = new TH3F(Form("hMPtDCA_LS_rec_misIDhf_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);

  hMPtDCA_LS_rec_charmTOe[j]        = new TH3F(Form("hMPtDCA_LS_rec_charmTOe_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_beautyTOe[j]       = new TH3F(Form("hMPtDCA_LS_rec_beautyTOe_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_hfTOe[j]           = new TH3F(Form("hMPtDCA_LS_rec_hfTOe_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_lfTOee[j]           = new TH3F(Form("hMPtDCA_LS_rec_lfTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_lfTOee_selectPDG[j] = new TH3F(Form("hMPtDCA_LS_rec_lfTOee_selectPDG_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_ccTOee[j]           = new TH3F(Form("hMPtDCA_LS_rec_ccTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_bbTOee[j]           = new TH3F(Form("hMPtDCA_LS_rec_bbTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  hMPtDCA_LS_rec_hfTOee[j]           = new TH3F(Form("hMPtDCA_LS_rec_hfTOee_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // before/after prefilter plots
    hMeeOpAngle_RecULS_bPF_pi0Dalitz[j] = new TH2F(Form("hMeeOpAngle_RecULS_bPF_pi0Dalitz_sce%i",j+1), ";#it{m}_{ee} GeV/#it{c}^{2};#phi rad", 40, 0., 0.4, 28, 0. , 0.7);
    hMeeOpAngle_RecULS_aPF_pi0Dalitz[j] = new TH2F(Form("hMeeOpAngle_RecULS_aPF_pi0Dalitz_sce%i",j+1), ";#it{m}_{ee} GeV/#it{c}^{2};#phi rad", 40, 0., 0.4, 28, 0. , 0.7);
    hMeeOpAngle_RecULS_bPF_bkg[j] = new TH2F(Form("hMeeOpAngle_RecULS_bPF_bkg_sce%i",j+1), ";#it{m}_{ee} GeV/#it{c}^{2};#phi rad", 40, 0., 0.4, 28, 0. , 0.7);
    hMeeOpAngle_RecULS_aPF_bkg[j] = new TH2F(Form("hMeeOpAngle_RecULS_aPF_bkg_sce%i",j+1), ";#it{m}_{ee} GeV/#it{c}^{2};#phi rad", 40, 0., 0.4, 28, 0. , 0.7);
    hMeeOpAngle_RecULS_bPF_HFPi0[j] = new TH2F(Form("hMeeOpAngle_RecULS_bPF_HFPi0_sce%i",j+1), ";#it{m}_{ee} GeV/#it{c}^{2};#phi rad", 40, 0., 0.4, 28, 0. , 0.7);
    hMeeOpAngle_RecULS_aPF_HFPi0[j] = new TH2F(Form("hMeeOpAngle_RecULS_aPF_HFPi0_sce%i",j+1), ";#it{m}_{ee} GeV/#it{c}^{2};#phi rad", 40, 0., 0.4, 28, 0. , 0.7);
    vecPFhistos[j].push_back(hMeeOpAngle_RecULS_bPF_pi0Dalitz[j]);
    vecPFhistos[j].push_back(hMeeOpAngle_RecULS_aPF_pi0Dalitz[j]);
    vecPFhistos[j].push_back(hMeeOpAngle_RecULS_bPF_bkg[j]);
    vecPFhistos[j].push_back(hMeeOpAngle_RecULS_aPF_bkg[j]);
    vecPFhistos[j].push_back(hMeeOpAngle_RecULS_bPF_HFPi0[j]);
    vecPFhistos[j].push_back(hMeeOpAngle_RecULS_aPF_HFPi0[j]);
  }




  TList* lRecPIDscenario[nPIDscenarios];
  for (int iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
    lRecPIDscenario[iPID_scenario] = new TList();
    lRecPIDscenario[iPID_scenario]->SetName(Form("PIDscenario_%i", iPID_scenario+1));
    lRecPIDscenario[iPID_scenario]->SetOwner();
    lRecPIDscenario[iPID_scenario]->Add(hAfterKineCuts_Pt_Eta_Phi_rec[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Muon_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Pion_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Kaon_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Proton_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_ElePos_Rec_Pt_Eta_Phi_afterPSPID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Pion_Rec_Pt_Eta_Phi_afterPSPID[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hNegTrack_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hPosTrack_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Muon_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Pion_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Kaon_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTrack_Proton_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hAllTracks_Rec_Pt_Eta_Phi[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hTracks_ElePos_Rec_inPFvector_Pt[iPID_scenario]);

    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle_charmTOe[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle_beautyTOe[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle_hfTOe[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle_lfTOee[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle_ccTOee[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle_bbTOee[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_rec_MCpidEle_hfTOee[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_MCpidEle[iPID_scenario]);
    // lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_misIDoneLeg[iPID_scenario]);
    // lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_misIDtwoLeg[iPID_scenario]);
    // lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_misIDPion[iPID_scenario]);
    // lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_misIDhf[iPID_scenario]);

    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_charmTOe[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_beautyTOe[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_hfTOe[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_lfTOee[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_lfTOee_selectPDG[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_ccTOee[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_bbTOee[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_rec_hfTOee[iPID_scenario]);
    // add prefilter plots
    lRecPIDscenario[iPID_scenario]->Add(hMeeOpAngle_RecULS_bPF_pi0Dalitz[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMeeOpAngle_RecULS_aPF_pi0Dalitz[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMeeOpAngle_RecULS_bPF_bkg[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMeeOpAngle_RecULS_aPF_bkg[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMeeOpAngle_RecULS_bPF_HFPi0[iPID_scenario]);
    lRecPIDscenario[iPID_scenario]->Add(hMeeOpAngle_RecULS_aPF_HFPi0[iPID_scenario]);
    // listRecTracks->Add(lRecPIDscenario);
  }


  TH3F* hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts  = new TH3F("hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts  = new TH3F("hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_Pion_Gen_Pt_Rap_Phi = new TH3F("hTrack_Pion_Gen_Pt_Rap_Phi",";#it{p}_{T} (GeV/#it{c});y;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);


  TH3F* hTrack_All_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_LF_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_HF_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_MotherPi0_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_GenSmeared_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Muon_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Pion_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Kaon_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Proton_Gen_Pt_Eta_Phi[nPIDscenarios];


  for (int j = 0; j < nPIDscenarios; ++j) hTrack_All_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_All_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_LF_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_LF_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_HF_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_HF_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_MotherPi0_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_MotherPi0_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_GenSmeared_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_GenSmeared_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Muon_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Muon_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Pion_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pion_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Kaon_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Kaon_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Proton_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Proton_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", n_ptee_bin_c, ptee_bin_c, nEtaBins, eta_binning, nPhiBins, phi_binning);

  for (int j = 0; j < nPIDscenarios; ++j) hTracks_ElePos_Gen_inPFvector_Pt[j] = new TH1F(Form("hTracks_ElePos_Gen_inPFvector_Pt_sce%i",j+1),";#it{p}_{T} (GeV/#it{c})",n_ptee_bin_c,ptee_bin_c);


  TList* lGenPIDscenario[nPIDscenarios];
  for (int iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
    lGenPIDscenario[iPID_scenario] = new TList();
    lGenPIDscenario[iPID_scenario]->SetName(Form("PIDscenario_%i", iPID_scenario+1));
    lGenPIDscenario[iPID_scenario]->SetOwner();
    lGenPIDscenario[iPID_scenario]->Add(hTrack_All_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_ElePos_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_ElePos_LF_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_ElePos_HF_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_ElePos_MotherPi0_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_ElePos_GenSmeared_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_Muon_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_Pion_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_Kaon_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hTrack_Proton_Gen_Pt_Eta_Phi[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hMPtDCA_ULS_gen[iPID_scenario]);
    lGenPIDscenario[iPID_scenario]->Add(hMPtDCA_LS_gen[iPID_scenario]);

    lGenPIDscenario[iPID_scenario]->Add(hTracks_ElePos_Gen_inPFvector_Pt[iPID_scenario]);
    // listGenTracks->Add(lGenPIDscenario);
  }

  std::vector<TList*> vecTListGenPIDscenarios;
  for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
    vecTListGenPIDscenarios.push_back(lGenPIDscenario[iPID_scenario]);
  }



  auto hM_Pt_sig_gen = new TH2F("hM_Pt_sig_gen",title2d.c_str(),300.,0,3.,200,0.,20.);

  TH1F* hTime0 = new TH1F("hTime0_sce%i",";t_{0} (ns)", 1000, -1., 1.);



  std::vector<TList*> vecTListRecPIDscenarios;
  for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
    vecTListRecPIDscenarios.push_back(lRecPIDscenario[iPID_scenario]);
  }

  // Setting Track and Particle vectors
  //-----------------------------------
  std::vector<Track *> vecElectron[nPIDscenarios],vecPositron[nPIDscenarios];
  std::vector<Track *> vecNegTracks[nPIDscenarios],vecPosTracks[nPIDscenarios];
  std::vector<Track *> vecNegTracks_pfPID[nPIDscenarios],vecPosTracks_pfPID[nPIDscenarios];
  std::vector<GenParticle *> vecNegPFTracks[nPIDscenarios],vecPosPFTracks[nPIDscenarios];
  // std::vector<Track *> vecNegPFTracks[nPIDscenarios],vecPosPFTracks[nPIDscenarios];
  std::vector<GenParticle *> vecElectronGen[nPIDscenarios],vecPositronGen[nPIDscenarios], vecElectronGenSmeared[nPIDscenarios],vecPositronGenSmeared[nPIDscenarios];
  std::vector<GenParticle *> vecGen;
  std::vector<Track *> vecPIDtracks; // , vecTOFtracks

  int nTotalTracks = 0;
  int nTotalParticles = 0;


  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {
    // if(ientry % (numberOfEntries/100) == 0) printf("events processed: %i of %i \n", ientry,numberOfEntries);

    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);
    nTracksGen->Fill(tracks->GetEntries());
    // nTracksGen->Fill(innertracks->GetEntries());
    nParticles->Fill(particles->GetEntries());
    // nTotalParticles = nTotalParticles + particles->GetEntries();
    // nTotalTracks = nTotalTracks + tracks->GetEntries();




    //##################################################
    //############## generated Track Loop ##############
    //##################################################
    TLorentzVector LV; // for smearing of generated particles
    TLorentzVector LV_smeared; // for smearing of generated particles
    Int_t numParticlesFIT=0;            // (2.2<eta&&eta< 5.0) || (-3.4<eta&&eta<-2.3))
    Int_t numParticlesMidRapidity=0;    // |eta|<0.5


    // loop over particles and fill vectors for generated electrons and positrons
    for (Int_t iparticle = 0; iparticle < particles->GetEntries(); ++iparticle) {

      // get particle
      auto particle = (GenParticle *)particles->At(iparticle);

      if(particle->Charge != 0) vecGen.push_back(particle);

      // Check isprimary and ischarged
      Bool_t isphysicalprimary = kTRUE;
      Bool_t ischarged = kTRUE;
      Bool_t ischargedmother = kTRUE;

      auto pdg = 0;
      pdg = particle->PID;
      if(!IsStable(pdg,ischarged)) isphysicalprimary = kFALSE;
      // mother
      auto imother  = particle->M1;
      auto mother   = imother != -1 ? (GenParticle *)particles->At(imother) : (GenParticle *)nullptr;
      auto mPid     = 0;
      if(mother) mPid = mother->PID;
      if(IsStable(pdg,ischarged) && IsStable(mPid,ischargedmother)) {
        //printf("pdg and mpdg are %d and %d, not counted\n",pdg,mPid);
        isphysicalprimary = kFALSE;
      }

      if(isphysicalprimary && ischarged) {
        // increase numer of particle if its in acceptance of  FIT detecotr
        double eta = particle->Eta;
        if( ((2.2<eta&&eta< 5.0) || (-3.4<eta&&eta<-2.3)) && particle->Charge != 0) numParticlesFIT++;
        if( (fabs(eta) < 0.5) && particle->Charge != 0) numParticlesMidRapidity++;
      }
    }






    // centrality selection
    //---------------------
    // accept event if event has more then X particles within the FIT detector acceptance (2.2<eta&&eta< 5.0) || (-3.4<eta&&eta<-2.3)
    nParticlesFIT->Fill(numParticlesFIT);
    nParticlesMidRapidity->Fill(numParticlesMidRapidity);
    if( (!bSemiCentral && numParticlesFIT < 5179 )                                 // PbPb: 5179,   pp: 58    // centrality of 0-10%
    // if( (numParticlesFIT < 584) || (numParticlesFIT > 1692) ) {       // PbPb: 584 - 1692,  // centrality of 40-60%
       || (bSemiCentral && ((numParticlesFIT < 1045) || (numParticlesFIT > 2549)) ) ){   // PbPb: 1045 - 2549,    // centrality of 30-50%
      vecGen.clear();
    continue;
    }
    nParticlesFITCent->Fill(numParticlesFIT);
    nParticlesMidRapidityCent->Fill(numParticlesMidRapidity);
    hdNdeta_midrap_gen->Fill(0.,numParticlesMidRapidity);






    // loop over particles within centrality
    //--------------------------------------
    for (auto particle : vecGen) {
      // if( (particle->Eta > -0.5 && particle->Eta<0.5) && particle->Charge != 0 ) hdNdeta_midrap_gen->Fill(particle->Eta);

      // get particle pdg
      auto pPID = particle->PID;
      // cout << " particle PID = " << pPID << endl;

      //Get mother
      auto imother  = particle->M1;
      auto mother   = imother != -1 ? (GenParticle *)particles->At(imother) : (GenParticle *)nullptr;
      auto mPid     = 0;
      if(mother) mPid = mother->PID;
      //Get grandmother
      auto igmother = 0;
      if(mother) igmother = mother->M1;
      auto gmother  = igmother != -1 ? (GenParticle *)particles->At(igmother) : (GenParticle *)nullptr;
      auto gmPid    = 0;
      if(gmother) gmPid = gmother->PID;

      // cout << " mother PID  = " << mPid << endl;
      // cout << " gmother PID = " << gmPid << endl << endl;

      // get weight to weight the contribution of LF and HF tracks
      Double_t weight = 1.;
      if(bUseWeights){
        if ( abs(pPID) == 11 && hasHeavyAncestor(particle, particles) ) {
          if (particle->PT < 0.1)                             weight = 1.;
          else if (particle->PT >= 0.1 && particle->PT < 0.6) weight = GetWeight(hweightHFtoe_lowerPt, particle->PT);
          else                                                weight = GetWeight(hweightHFtoe_higherPt, particle->PT);
        }
        else if ( abs(pPID) == 11 && !hasHeavyAncestor(particle, particles) ) {
          if(mother->PT < 0.1) weight = 0.7;
          else                 weight = GetWeight(hweightLFtoe, mother->PT);
          // printf("weight: %f, pt: %f \n", weight, mother->PT);
        }
        else {
          weight = GetWeight(hweightLFtoe, particle->PT);
          if(weight == 0.) weight = 1.;
        }
      }


      double phiGen = TVector2::Phi_0_2pi(particle->Phi);

      // Fill positve and negative vecotrs with tracks that don't have tracking efficiency, no smearing, true MC PID
      for (Int_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
        if( usePreFilter[iPID_scenario] && (abs(particle->PID) == 11) && (particle->PT >= pf_pt_cut[iPID_scenario]) && (fabs(particle->Eta) < 1.75) ){
          if     (particle->Charge < 0) vecNegPFTracks[iPID_scenario].push_back(particle);
          else if(particle->Charge > 0) vecPosPFTracks[iPID_scenario].push_back(particle);
          hTracks_ElePos_Gen_inPFvector_Pt[iPID_scenario]->Fill(particle->PT);
        }
      }

      // Smearing generated particles
      //-----------------------------
      LV.SetPtEtaPhiM(particle->PT,particle->Eta,phiGen,eMass);
      Double_t charge = -particle->PID;
      LV_smeared = ApplySmearing(fArrResoPt,fArrResoEta,fArrResoPhi_Pos,fArrResoPhi_Neg,LV,charge);

      //ing suport histograms to see effect of smearing
      if(charge > 0.)  hSmearing_For_Eff_phi_pos->Fill(particle->PT, LV.Phi() - LV_smeared.Phi(),weight);
      if(charge < 0.)  hSmearing_For_Eff_phi_neg->Fill(particle->PT, LV.Phi() - LV_smeared.Phi(),weight);
                       hSmearing_For_Eff_eta->Fill(particle->PT, LV.Eta() - LV_smeared.Eta(),weight);
      if(LV.Pt() > 0.) hSmearing_For_Eff_pt->Fill(particle->PT, (LV.Pt() - LV_smeared.Pt())/LV.Pt(),weight);

      double phiGenSmear = TVector2::Phi_0_2pi(LV_smeared.Phi());

      // Filling TH3 before kinematic cuts are applied
      if (abs(particle->PID) == 11 )  hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts       ->Fill(particle->PT,particle->Eta,phiGen,weight); // Pt Eta Phi of generated electrons + positrons
      if (abs(particle->PID) == 11 )  hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts->Fill(LV_smeared.Pt(),LV_smeared.Eta(),phiGenSmear,weight);     // Pt Eta Phi of generated electrons + positrons


      if (abs(particle->PID) == 211 )  {
        double rap = etatorap(particle->PT,particle->Eta, 0.13957018); // use charged pion mass
        if(abs(rap) <= 0.5) hTrack_Pion_Gen_Pt_Rap_Phi->Fill(particle->PT, rap, phiGen, weight);
      }

      // loop over iPID_scenario
      for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {

        // apply eta acceptance cut
        if(etaCut(particle, iPID_scenario)){

        // // kinematic cuts on particles
        // if (kineCuts(particle,iPID_scenario)){

          // Filling TH3 histograms for generated  Tracks
          //---------------------------------------------
                                            hTrack_All_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight); // Pt Eta Phi of all generated tracks
          if (abs(particle->PID) == 11 )    hTrack_ElePos_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 13 )    hTrack_Muon_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 211 )   hTrack_Pion_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 321 )   hTrack_Kaon_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 2212 )  hTrack_Proton_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight); // Pt Eta Phi of generated electrons + positrons
          if((abs(particle->PID) == 11) && (abs(mPid) == 111 || abs(mPid) == 221 || abs(mPid) == 331 || abs(mPid) == 223 || abs(mPid) == 333 || abs(mPid) == 113) ) hTrack_ElePos_LF_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight);
          if((abs(particle->PID) == 11) && (hasHeavyAncestor(particle, particles))) hTrack_ElePos_HF_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight);
          if((abs(particle->PID) == 11) && (abs(mPid) == 111))  hTrack_ElePos_MotherPi0_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen,weight);  // Pt Eta Phi of generated electrons + positrons from Pi0
          // cout << "PID: " << abs(particle->PID) << ", hasHFancestor = " << hasHeavyAncestor(particle, particles) << ", pt = " << particle->PT << ", Mother pt = " << mother->PT << ", weight = " << weight << endl;

        }

        // fill generated vectors
        if (kineCuts(particle,iPID_scenario)){
          if (particle->PID == 11 ) vecElectronGen[iPID_scenario].push_back(particle);       // vector filled with generated electrons
          else if (particle->PID == -11 ) vecPositronGen[iPID_scenario].push_back(particle); // vector filled with generated positrons
        }


        // // applying kinematic cuts on smeared LV vector
        // if (!kineCuts(LV_smeared,iPID_scenario)) continue;
        if (etaCut(LV_smeared, iPID_scenario)) {
                          // cout << "gen pt eta phi: " << particle->PT << " " << particle->Eta << " " << phiGen << endl;
                          // cout << "gen smeared pt eta phi: " << LV_smeared.Pt() << " " << LV_smeared.Eta() << " " << phiGenSmear << endl;


          // Filling TH3 histograms for generated smeared Tracks
          //----------------------------------------------------
          // look only at electrons
          if (abs(particle->PID) == 11 )  hTrack_ElePos_GenSmeared_Pt_Eta_Phi[iPID_scenario]->Fill(LV_smeared.Pt(),LV_smeared.Eta(),phiGenSmear,weight); // Pt Eta Phi of generated smeared electrons + positrons
        }

        // fill generated smeared vectors
        if (kineCuts(LV_smeared,iPID_scenario)){
          if (particle->PID == 11 )       vecElectronGenSmeared[iPID_scenario].push_back(particle);        // vector filled with generated smeared electrons
          else if (particle->PID == -11 ) vecPositronGenSmeared[iPID_scenario].push_back(particle);        // vector filled with generated smeared positrons
        }

      }
    }
    vecGen.clear();






    //##################################################
    //########## reconstructed Track Loop ##############
    //##################################################




    // loop over tracks and fill vectors for electrons and positrons
    for (Int_t itrack = 0; itrack < tracks->GetEntries(); ++itrack)
    {
      // get track and corresponding particle
      auto track = (Track *)tracks->At(itrack);
      auto particle = (GenParticle *)track->Particle.GetObject();

      //Get mother
      auto imother  = particle->M1;
      auto mother   = imother != -1 ? (GenParticle *)particles->At(imother) : (GenParticle *)nullptr;
      auto mPid     = mother->PID;
      //Get grandmother
      auto igmother = mother->M1;
      auto gmother  = igmother != -1 ? (GenParticle *)particles->At(igmother) : (GenParticle *)nullptr;
      auto gmPid    = 0;
      if(gmother) gmPid = gmother->PID;

      double phiRec = TVector2::Phi_0_2pi(track->Phi);


      // get weight to weight the contribution of LF and HF tracks
      Double_t weight = 1.;
      if(bUseWeights){
        if ( abs(particle->PID) == 11 && hasHeavyAncestor(particle, particles) ) {
          if (track->PT < 0.1)                          weight = 1.;
          else if (track->PT >= 0.1 && track->PT < 0.6) weight = GetWeight(hweightHFtoe_lowerPt, track->PT);
          else                                          weight = GetWeight(hweightHFtoe_higherPt, track->PT);
        }
        else if ( abs(particle->PID) == 11 && !hasHeavyAncestor(particle, particles) ) {
          if(mother->PT < 0.1) weight = 0.7;
          else                 weight = GetWeight(hweightLFtoe, mother->PT);
        }
        else {
          weight = GetWeight(hweightLFtoe, track->PT);
          if(weight == 0.) weight = 1.;
        }
      }

      // Fill positve and negative vecotrs with tracks that don't have tracking efficiency, no smearing, true MC PID
      for (Int_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
        if( usePreFilter[iPID_scenario] && (abs(particle->PID) == 11) && (track->PT >= pf_pt_cut[iPID_scenario]) && (fabs(track->Eta) < 1.75) ){
          // if     (particle->Charge < 0) vecNegPFTracks[iPID_scenario].push_back(track);
          // else if(particle->Charge > 0) vecPosPFTracks[iPID_scenario].push_back(track);
          hTracks_ElePos_Rec_inPFvector_Pt[iPID_scenario]->Fill(track->PT);
        }
      }


      // smear track if requested
      if(abs(particle->PID) == 11 ) hBeforeSmearing_Pt_Eta_Phi_rec->Fill(track->PT,track->Eta,phiRec,weight);
      // cout << "track Pt,Eta,Phi: " << track->PT << ", " <<  track->Eta << ", " << track->Phi << endl;
      if (bSmear) if (!smearer.smearTrack(*track)) continue; // strange syntax, but works
      if(TMath::IsNaN(track->Phi)){
          // cout << "track: " << track << endl;
          // cout << "track Pt,Eta,Phi: " << track->PT << ", " <<  track->Eta << ", " << track->Phi << endl;
          continue;
      }

      // cut away tracks that are way off.
      if (fabs(track->D0) > 0.4) continue; // adopt to just stay in the beampipe?
      if (fabs(track->DZ) > 3.) continue;

      // check if has TOF
      // if (toflayer.hasTOF(*track)) vecTOFtracks.push_back(track);


      // push all tracks into a vector.
      vecPIDtracks.push_back(track);
    }






    nTracks->Fill(vecPIDtracks.size());
    // if (vecPIDtracks.size() < 3750) {vecPIDtracks.clear(); continue;} // dirty dirty centrality
    nTracksCent->Fill(vecPIDtracks.size());

    // std::array<float, 2> tzero;
    // toflayer.eventTime(vecTOFtracks, tzero);
    // hTime0->Fill(tzero[0]);
    // vecTOFtracks.clear();





    // loop over tracks and fill vectors for reconstructed electrons and positrons
    // for (Int_t itrack = 0; itrack < tracks->GetEntries(); ++itrack) {
    // for(auto track : vecPIDtracks){
    for(Int_t itrack = 0; itrack < vecPIDtracks.size(); itrack++){

      // get track and corresponding particle
      // auto track = (Track *)tracks->At(itrack);
      auto track = (Track *)vecPIDtracks.at(itrack);
      auto particle = (GenParticle *)track->Particle.GetObject();

      //Get mother
      auto imother  = particle->M1;
      auto mother   = imother != -1 ? (GenParticle *)particles->At(imother) : (GenParticle *)nullptr;
      auto mPid     = mother->PID;
      // //Get grandmother
      // auto igmother = mother->M1;
      // auto gmother  = igmother != -1 ? (GenParticle *)particles->At(igmother) : (GenParticle *)nullptr;
      // auto gmPid    = 0;
      // if(gmother) gmPid = gmother->PID;

      // get weight to weight the contribution of LF and HF tracks
      Double_t weight = 1.;
      if(bUseWeights){
        if ( abs(particle->PID) == 11 && hasHeavyAncestor(particle, particles) ) {
          if (particle->PT < 0.1)                             weight = 1.;
          else if (particle->PT >= 0.1 && particle->PT < 0.6) weight = GetWeight(hweightHFtoe_lowerPt, track->PT);
          else                                                weight = GetWeight(hweightHFtoe_higherPt, track->PT);
        }
        else if ( abs(particle->PID) == 11 && !hasHeavyAncestor(particle, particles) ) {
          if(mother->PT < 0.1) weight = 0.7;
          else                 weight = GetWeight(hweightLFtoe, mother->PT);
        }
        else {
          weight = GetWeight(hweightLFtoe, track->PT);
          if(weight == 0.) weight = 1.;
        }
      }


      if( (particle->Eta > -0.5 && particle->Eta<0.5) && particle->Charge != 0 ) hdNdeta_midrap_rec->Fill(particle->Eta,weight);

      // if(track->PT < 0.1) continue; // debug code
      double phiRec = TVector2::Phi_0_2pi(track->Phi);
      if(TMath::IsNaN(phiRec)) continue;


      if(abs(particle->PID) == 11 ) hAfterSmearing_Pt_Eta_Phi_rec->Fill(track->PT,track->Eta,phiRec,weight);


      auto p = track->P;

      //make inner TOF PID
      std::array<float, 5> PIDdeltatInnerTOF, PIDnsigmaInnerTOF;
      // if(bUseiTOF) innertoflayer.makePID(*track, PIDdeltatInnerTOF, PIDnsigmaInnerTOF);

      //make outer TOF PID
      std::array<float, 5> PIDdeltatTOF, PIDnsigmaTOF;
      toflayer.makePID(*track, PIDdeltatTOF, PIDnsigmaTOF);
      //make RICH PID
      std::array<float, 5> PIDdeltaangleRICH, PIDnsigmaRICH;
      richdetector.makePID(*track, PIDdeltaangleRICH, PIDnsigmaRICH);

      // fill nsigma TOF   (before PID selection is applied)
      hTrackPt->Fill(track->PT,weight);
      hTrackP->Fill(track->P,weight);
      if (toflayer.hasTOF(*track) || richdetector.hasRICH(*track)){
        hTrackPt_hasTOForRICH->Fill(track->PT,weight);
        hTrackP_hasTOForRICH->Fill(track->P,weight);
      }


      for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {

        // apply eta acceptance cut
        if(!etaCut(track, iPID_scenario)) continue;


        if(abs(particle->PID) == 11 ) {
                                        hTrack_ElePos_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed electrons + positrons from Pi0 before PID
          if(abs(mPid) == 111)          hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed electrons + positrons from Pi0 before PID
          if(abs(mPid) == 111 || abs(mPid) == 221 || abs(mPid) == 331 || abs(mPid) == 223 || abs(mPid) == 333 || abs(mPid) == 113)
                                        hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);   // Pt Eta Phi of reconstructed electrons + positrons from LF before PID (pi0,eta,etaprime,omega,phi,rho)
          if(hasHeavyAncestor(particle, particles))
                                        hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);   // Pt Eta Phi of reconstructed electrons + positrons from HF before PID
        }
        else if (abs(particle->PID) == 13)    hTrack_Muon_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 211)   hTrack_Pion_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 321)   hTrack_Kaon_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 2212)  hTrack_Proton_Rec_Pt_Eta_Phi_beforePID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);

        if(abs(particle->PID) == 11 ) hAfterKineCuts_Pt_Eta_Phi_rec[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);

        //#########################
        //     prefilter PID
        //#########################
        // PID equivilant to standard PID, just with different Nsigma cuts and pt acceptance

        bool   i_useTOFPID;
        bool   i_useRICHPID;
        bool   i_usePSPID;
        double i_SigmaTOFEle;
        double i_SigmaTOFPi;
        double i_SigmaRICHEle;
        double i_SigmaRICHPi;
        double ps_pt_acc;
        if(busePFpid[iPID_scenario] == kTRUE){
          i_useTOFPID    = pf_PID_useTOF;
          i_useRICHPID   = pf_PID_useRICH;
          i_usePSPID     = pf_PID_usePS;
          i_SigmaTOFEle  = pf_PID_NSigTOFele;
          i_SigmaTOFPi   = pf_PID_NSigTOFpi;
          i_SigmaRICHEle = pf_PID_NSigRICHele;
          i_SigmaRICHPi  = pf_PID_NSigRICHpi;
          ps_pt_acc      = pf_PID_PSptacc;
        }
        else {
          i_useTOFPID    = useTOFPID[iPID_scenario];
          i_useRICHPID   = useRICHPID[iPID_scenario];
          i_usePSPID     = usePreShPID[iPID_scenario];
          i_SigmaTOFEle  = nSigmaTOFEle[iPID_scenario];
          i_SigmaTOFPi   = nSigmaTOFPi[iPID_scenario];
          i_SigmaRICHEle = nSigmaRICHEle[iPID_scenario];
          i_SigmaRICHPi  = nSigmaRICHPi[iPID_scenario];
          ps_pt_acc      = ps_EleAccep_pt_cut[iPID_scenario];
        }


        // Apply Preshower electron selection
        bool isElectronPreShower = false;
        if( i_usePSPID && preshower.hasPreShower(*track) ) isElectronPreShower = preshower.isElectron(*track,1800);
        if ( (abs(particle->PID) == 11)  && (isElectronPreShower == true) )   hTrack_ElePos_Rec_Pt_Eta_Phi_afterPSPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        if ( (abs(particle->PID) == 211) && (isElectronPreShower == true) )   hTrack_Pion_Rec_Pt_Eta_Phi_afterPSPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);

        bool bool_doPID = false;
        if(                       i_useRICHPID && !bUseiTOF && doPID(track, i_useTOFPID, i_useRICHPID, tof_EleAccep_p_cut, tof_PionRej_p_cut, rich_PionRejection_p_cut, i_SigmaTOFEle, i_SigmaTOFPi, i_SigmaRICHEle, i_SigmaRICHPi, PIDnsigmaTOF, PIDnsigmaRICH) ) bool_doPID = true;
        else if(  i_useTOFPID && !i_useRICHPID && !bUseiTOF && doTOFPID(track, i_useTOFPID, tof_EleAccep_p_cut, tof_PionRej_p_cut, i_SigmaTOFEle, i_SigmaTOFPi, PIDnsigmaTOF) ) bool_doPID = true;
        else if(  i_useTOFPID &&                   bUseiTOF && doTOFPID(track, i_useTOFPID, itof_EleAccep_p_cut, tof_PionRej_p_cut, i_SigmaTOFEle, i_SigmaTOFPi, PIDnsigmaTOF) ) bool_doPID = true;

        // ----------- PID strategy -----------
        // if:
        // check if any of the detectors is used
        // if so -> check wheather (track pt < PS pt acceptance and PS identifies track as electron)
        //          AND (track not identified as electron by TOF and/or RICH)
        // then  -> reject track
        //
        // else if: if no detector is used AND no monte carlo PID is used
        // then  -> reject track

        if ( (i_usePSPID || i_useTOFPID || i_useRICHPID) && ( ((isElectronPreShower == false) || (track->PT < ps_pt_acc) ) && (bool_doPID == false)) ) continue;
        else if ( !(i_usePSPID || i_useTOFPID || i_useRICHPID) && !useMCPID[iPID_scenario] ) continue;
        // if ( (isElectronPreShower == false) && (bool_doPID == false) &&  !useMCPID[iPID_scenario] ) continue;
        if (useMCPID[iPID_scenario] && (abs(track->PID) != 11)) continue;

        // fill track vector
        if (track->Charge < 0) vecNegTracks_pfPID[iPID_scenario].push_back(track);       //  vector of all reconstructed negative tracks
        else if (track->Charge > 0) vecPosTracks_pfPID[iPID_scenario].push_back(track);  //  vector of all reconstructed positve tracks

      }
    }
    vecPIDtracks.clear();






    //#########################
    //     Apply PreFilter
    //#########################
    for (int iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
      if(usePreFilter[iPID_scenario] && !usePFtrackVec[iPID_scenario]){
        PreFilter(&vecNegTracks[iPID_scenario], &vecPosTracks[iPID_scenario], iPID_scenario, vecPFhistos[iPID_scenario]);
      }
      else if(usePreFilter[iPID_scenario] && usePFtrackVec[iPID_scenario]){
        std::vector<Int_t> vecMotherID;
        std::vector<std::tuple<GenParticle*, GenParticle*>> vecGenTuples;
        PreFilter(&vecNegPFTracks[iPID_scenario], &vecPosTracks_pfPID[iPID_scenario], vecMotherID, vecGenTuples, iPID_scenario, vecPFhistos[iPID_scenario]);
        PreFilter(&vecPosPFTracks[iPID_scenario], &vecNegTracks_pfPID[iPID_scenario], vecMotherID, vecGenTuples, iPID_scenario, vecPFhistos[iPID_scenario]);
        // PreFilter(&vecElectron[iPID_scenario], &vecPositron[iPID_scenario], particles, iPID_scenario, hMeeOpAngle_RecULS_bPF[iPID_scenario], hMeeOpAngle_RecULS_aPF[iPID_scenario]);
        vecMotherID.clear();
        vecGenTuples.clear();
      }
      vecNegPFTracks[iPID_scenario].clear();
      vecPosPFTracks[iPID_scenario].clear();
    }




    //#########################
    //  Apply standard  PID       (is used only when separate PF PID is applied)
    //#########################
    for (int iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
      if(busePFpid[iPID_scenario]){
        // printf("iSce = %i\n", iPID_scenario);
        ApplyPID(&vecNegTracks_pfPID[iPID_scenario], &vecNegTracks[iPID_scenario], iPID_scenario);
        ApplyPID(&vecPosTracks_pfPID[iPID_scenario], &vecPosTracks[iPID_scenario], iPID_scenario);
      }
      else{
        vecNegTracks[iPID_scenario] = vecNegTracks_pfPID[iPID_scenario];
        vecPosTracks[iPID_scenario] = vecPosTracks_pfPID[iPID_scenario];
      }
      vecNegTracks_pfPID[iPID_scenario].clear();
      vecPosTracks_pfPID[iPID_scenario].clear();
    }




    // fill single track histograms
    for (int iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
      // loop over negative track vector
      for (int itrack = vecNegTracks[iPID_scenario].size()-1; itrack >= 0 ; itrack--) {

        Track* track = (Track *) vecNegTracks[iPID_scenario].at(itrack);
        auto particle = (GenParticle *)track->Particle.GetObject();

        //Get mother
        auto imother  = particle->M1;
        auto mother   = imother != -1 ? (GenParticle *)particles->At(imother) : (GenParticle *)nullptr;
        auto mPid     = mother->PID;

        // get weight to weight the contribution of LF and HF tracks
        Double_t weight = 1.;
        if(bUseWeights){
          if ( abs(particle->PID) == 11 && hasHeavyAncestor(particle, particles) ) {
            if (particle->PT < 0.1)                             weight = 1.;
            else if (particle->PT >= 0.1 && particle->PT < 0.6) weight = GetWeight(hweightHFtoe_lowerPt, track->PT);
            else                                                weight = GetWeight(hweightHFtoe_higherPt, track->PT);
          }
          else if ( abs(particle->PID) == 11 && !hasHeavyAncestor(particle, particles) ) {
            if(mother->PT < 0.1) weight = 0.7;
            else                 weight = GetWeight(hweightLFtoe, mother->PT);
          }
          else {
            weight = GetWeight(hweightLFtoe, track->PT);
            if(weight == 0.) weight = 1.;
          }
        }

        // if(track->PT < 0.1) continue; // debug code
        double phiRec = TVector2::Phi_0_2pi(track->Phi);


        // fill histograms
                                    hAllTracks_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight); // Pt Eta Phi of all type of reconstructed particles
        if     (track->Charge < 0)  hNegTrack_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed negative tracks
        else if(track->Charge > 0)  hPosTrack_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed positive tracks



        // Pt Eta Phi plots for different type of particles
        if(abs(particle->PID) == 11 ) {   // considering only electrons & positrons
                                               hTrack_ElePos_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed electrons + positrons
          if(abs(mPid) == 111)                 hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed electrons + positrons from Pi0
          if(abs(mPid) == 111 || abs(mPid) == 221 || abs(mPid) == 331 || abs(mPid) == 223 || abs(mPid) == 333 || abs(mPid) == 113)
                                               hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);   // Pt Eta Phi of reconstructed electrons + positrons from LF (pi0,eta,etaprime,omega,phi,rho)
          if(hasHeavyAncestor(particle, particles))
                                               hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);   // Pt Eta Phi of reconstructed electrons + positrons from HF
        }
        else if (abs(particle->PID) == 13)    hTrack_Muon_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 211)   hTrack_Pion_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 321)   hTrack_Kaon_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 2212)  hTrack_Proton_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else std::cout << "Particle not identified!    Particle PID: " << abs(particle->PID) << std::endl;

        //apply pt cuts before filling vectors for pairing
        if (!kineCuts(track, iPID_scenario)) {
            vecNegTracks[iPID_scenario].erase(vecNegTracks[iPID_scenario].begin()+itrack);
            continue;  // end itrack loop

        }
        if (particle->PID == 11 ) vecElectron[iPID_scenario].push_back(track);       // vector of reconstructed electrons
      }

      // loop over positve track vector
      for (int itrack = vecPosTracks[iPID_scenario].size()-1; itrack >= 0 ; itrack--) {
        Track* track = (Track *) vecPosTracks[iPID_scenario].at(itrack);
        auto particle = (GenParticle *)track->Particle.GetObject();

        //Get mother
        auto imother  = particle->M1;
        auto mother   = imother != -1 ? (GenParticle *)particles->At(imother) : (GenParticle *)nullptr;
        auto mPid     = mother->PID;

        // get weight to weight the contribution of LF and HF tracks
        Double_t weight = 1.;
        if(bUseWeights){
          if ( abs(particle->PID) == 11 && hasHeavyAncestor(particle, particles) ) {
            if (particle->PT < 0.1)                             weight = 1.;
            else if (particle->PT >= 0.1 && particle->PT < 0.6) weight = GetWeight(hweightHFtoe_lowerPt, track->PT);
            else                                                weight = GetWeight(hweightHFtoe_higherPt, track->PT);
          }
          else if ( abs(particle->PID) == 11 && !hasHeavyAncestor(particle, particles) ) {
            if(mother->PT < 0.1) weight = 0.7;
            else                 weight = GetWeight(hweightLFtoe, mother->PT);
          }
          else {
            weight = GetWeight(hweightLFtoe, track->PT);
            if(weight == 0.) weight = 1.;
          }
        }

        // if(track->PT < 0.1) continue; // debug code
        double phiRec = TVector2::Phi_0_2pi(track->Phi);


        // fill histograms
                                    hAllTracks_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight); // Pt Eta Phi of all type of reconstructed particles
        if     (track->Charge < 0)  hNegTrack_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed negative tracks
        else if(track->Charge > 0)  hPosTrack_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed positive tracks

        // Pt Eta Phi plots for different type of particles
        if(abs(particle->PID) == 11 ) {   // considering only electrons & positrons
                                          hTrack_ElePos_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed electrons + positrons
          if(abs(mPid) == 111)            hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);  // Pt Eta Phi of reconstructed electrons + positrons from Pi0
          if(abs(mPid) == 111 || abs(mPid) == 221 || abs(mPid) == 331 || abs(mPid) == 223 || abs(mPid) == 333 || abs(mPid) == 113)
                                          hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);   // Pt Eta Phi of reconstructed electrons + positrons from LF (pi0,eta,etaprime,omega,phi,rho)
          if(hasHeavyAncestor(particle, particles))
                                          hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_afterPID[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);   // Pt Eta Phi of reconstructed electrons + positrons from HF
        }
        else if (abs(particle->PID) == 13)    hTrack_Muon_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 211)   hTrack_Pion_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 321)   hTrack_Kaon_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else if (abs(particle->PID) == 2212)  hTrack_Proton_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec,weight);
        else std::cout << "Particle not identified!    Particle PID: " << abs(particle->PID) << std::endl;

        //apply pt cuts before filling vectors for pairing
        if (!kineCuts(track, iPID_scenario)) {
            vecPosTracks[iPID_scenario].erase(vecPosTracks[iPID_scenario].begin()+itrack);
            continue;  // end itrack loop
        }
        if (particle->PID == -11 ) vecPositron[iPID_scenario].push_back(track);       // vector of reconstructed positrons
      }
    }




    // cout << "Number of tracks in vecNegTracks[0] = " <<  vecNegTracks[0].size() << endl;

    //##################################################
    //############   ULS and LS pairing   ##############
    //############ for Gen and Rec tracks ##############
    //##################################################

    bool pairULS = kTRUE;
    bool MCpidEle = kTRUE;

    for (Int_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
      dileptonPairingRec(vecNegTracks[iPID_scenario], vecPosTracks[iPID_scenario],  pairULS, !MCpidEle, iPID_scenario);
      dileptonPairingRec(vecNegTracks[iPID_scenario], vecNegTracks[iPID_scenario], !pairULS, !MCpidEle, iPID_scenario);
      dileptonPairingRec(vecPosTracks[iPID_scenario], vecPosTracks[iPID_scenario], !pairULS, !MCpidEle, iPID_scenario);
      dileptonPairingRec(vecElectron[iPID_scenario],  vecPositron[iPID_scenario],   pairULS,  MCpidEle, iPID_scenario);
      dileptonPairingRec(vecElectron[iPID_scenario],  vecElectron[iPID_scenario],  !pairULS,  MCpidEle, iPID_scenario);
      dileptonPairingRec(vecPositron[iPID_scenario],  vecPositron[iPID_scenario],  !pairULS,  MCpidEle, iPID_scenario);

      vecNegTracks[iPID_scenario].clear();
      vecPosTracks[iPID_scenario].clear();
      vecElectron[iPID_scenario].clear();
      vecPositron[iPID_scenario].clear();


      dileptonPairingGen(vecElectronGen[iPID_scenario], vecPositronGen[iPID_scenario], pairULS, iPID_scenario);
      dileptonPairingGen(vecElectronGen[iPID_scenario], vecElectronGen[iPID_scenario], !pairULS, iPID_scenario);
      dileptonPairingGen(vecPositronGen[iPID_scenario], vecPositronGen[iPID_scenario], !pairULS, iPID_scenario);

      vecElectronGen[iPID_scenario].clear();
      vecPositronGen[iPID_scenario].clear();
      vecElectronGenSmeared[iPID_scenario].clear();
      vecPositronGenSmeared[iPID_scenario].clear();
    }

  }


  auto fout = TFile::Open(outputFile, "RECREATE");
  // TList* listGenerated;
  // listULS = new TList(); listULS->SetName("ULS"); listULS->SetOwner();
  // listULS->Add(nTracks);
  // listULS->Write();

  // fout->cd();
  // listGenerated->Write("generated",1);

  fout->cd();
  fout->mkdir("generated/");
  fout->cd("generated/");

  for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
    fout->mkdir(Form("generated/%s",vecTListGenPIDscenarios.at(iPID_scenario)->GetName()));
    fout->cd(Form("generated/%s",vecTListGenPIDscenarios.at(iPID_scenario)->GetName()));
    vecTListGenPIDscenarios.at(iPID_scenario)->Write();
    vecTListGenPIDscenarios.at(iPID_scenario)->Clear();
  }
  vecTListGenPIDscenarios.clear();

  fout->cd();
  fout->mkdir("reconstructed");
  fout->cd("reconstructed");

  for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
    fout->mkdir(Form("reconstructed/%s",vecTListRecPIDscenarios.at(iPID_scenario)->GetName()));
    fout->cd(Form("reconstructed/%s",vecTListRecPIDscenarios.at(iPID_scenario)->GetName()));
    vecTListRecPIDscenarios.at(iPID_scenario)->Write();
    vecTListRecPIDscenarios.at(iPID_scenario)->Clear();
  }
  vecTListRecPIDscenarios.clear();


  fout->cd();
  hSmearing_For_Eff_phi_pos->Write();
  hSmearing_For_Eff_phi_neg->Write();
  hSmearing_For_Eff_eta->Write();
  hSmearing_For_Eff_pt->Write();

  nParticles->Write();
  nTracksGen->Write();
  nTracks->Write();
  nTracksCent->Write();
  nParticlesFIT->Write();
  nParticlesFITCent->Write();
  nParticlesMidRapidity->Write();
  nParticlesMidRapidityCent->Write();
  hdNdeta_midrap_gen->Write();
  hdNdeta_midrap_rec->Write();
  hTrackP->Write();
  hTrackP_hasTOForRICH->Write();
  hTrackPt->Write();
  hTrackPt_hasTOForRICH->Write();

  hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts->Write();
  hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts->Write();
  hTrack_Pion_Gen_Pt_Rap_Phi->Write();
  hBeforeSmearing_Pt_Eta_Phi_rec->Write();
  hAfterSmearing_Pt_Eta_Phi_rec->Write();

  hM_Pt_sig_gen->Write();

  fout->Close();

  // stop watch
  watch->Stop();
  watch->Print();

}



// function to ULS or LS pair reconstructed particles
void dileptonPairingRec(std::vector<Track *> vec_track_neg,std::vector<Track *> vec_track_pos, bool pairULS, bool MCpidEle, Int_t iSce){
  //#########################################
  //  Pairing two "reconstructed" particles
  //#########################################
  // pairing reconstructed ULS or LS pairs
  double vec_track_neg_size = vec_track_neg.size();
  double vec_track_pos_size = vec_track_pos.size();
  if(!pairULS && (vec_track_neg_size == 1)) return;
  int    iPos_start = 0;
  int    iNeg_end = vec_track_neg_size;
  if(!pairULS) iNeg_end = vec_track_neg_size-1;
  TLorentzVector LV1,LV2,LV;
  double dca= 0. ,dca1 = 0.,dca2 = 0.;

  for (int iEle = 0; iEle < iNeg_end; iEle++){
    auto track1 = (Track *) vec_track_neg.at(iEle);
    auto particle1 = (GenParticle *)track1->Particle.GetObject();
    auto imother1 = particle1->M1;
    auto mother1 = imother1 != -1 ? (GenParticle *)particles->At(imother1) : (GenParticle *)nullptr;
    auto m1Pid = mother1->PID;
    LV1.SetPtEtaPhiM(track1->PT,track1->Eta,track1->Phi,eMass);

    if(!pairULS) iPos_start = iEle+1;
    for (auto iPos = iPos_start; iPos < vec_track_pos_size; iPos++){
      auto track2 = (Track *) vec_track_pos.at(iPos);
      auto particle2 = (GenParticle *)track2->Particle.GetObject();
      auto imother2 = particle2->M1;
      auto mother2 = imother2 != -1 ? (GenParticle *)particles->At(imother2) : (GenParticle *)nullptr;
      auto m2Pid = mother2->PID;
      LV2.SetPtEtaPhiM(track2->PT,track2->Eta,track2->Phi,eMass);
      LV = LV1 + LV2;

      dca1 = (track1->D0/track1->ErrorD0);
      dca2 = (track2->D0/track2->ErrorD0);
      dca = sqrt((dca1*dca1 + dca2*dca2) / 2);


      // get weights for both particles to weight the contribution of LF and HF tracks
      Double_t weight_neg = 1.;
      Double_t weight_pos = 1.;
      if(bUseWeights){
      // weight for negative particle
        if (abs(particle1->PID) == 11 && hasHeavyAncestor(particle1, particles)) {
          if (track1->PT < 0.1)                           weight_neg = 1.;
          else if (track1->PT >= 0.1 && track1->PT < 0.6) weight_neg = GetWeight(hweightHFtoe_lowerPt, track1->PT);
          else                                            weight_neg = GetWeight(hweightHFtoe_higherPt, track1->PT);
        }
        else if (abs(particle1->PID) == 11 && !hasHeavyAncestor(particle1, particles)) {
          if(mother1->PT < 0.1) weight_neg = 0.7;
          else                  weight_neg = GetWeight(hweightLFtoe, mother1->PT);
        }
        else {
          weight_neg = GetWeight(hweightLFtoe, track1->PT);
          if(weight_neg == 0.) weight_neg = 1.;
        }

        // weight for positive particle
        if ( abs(particle2->PID) == 11 && hasHeavyAncestor(particle2, particles) ) {
          if (track2->PT < 0.1)                           weight_pos = 1.;
          else if (track2->PT >= 0.1 && track2->PT < 0.6) weight_pos = GetWeight(hweightHFtoe_lowerPt, track2->PT);
          else                                            weight_pos = GetWeight(hweightHFtoe_higherPt, track2->PT);
        }
        else if (abs(particle2->PID) == 11 && !hasHeavyAncestor(particle2, particles) ) {
          if(mother2->PT < 0.1) weight_pos = 0.7;
          else                  weight_pos = GetWeight(hweightLFtoe, mother2->PT);
        }
        else {
          weight_pos = GetWeight(hweightLFtoe, track2->PT);
          if(weight_pos == 0.) weight_pos = 1.;
        }
      }

      std::vector<Int_t> vecLFpdgs = {111,221,331,223,333,113};
      if(pairULS && !MCpidEle)  hMPtDCA_ULS_rec[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if(!pairULS && !MCpidEle) hMPtDCA_LS_rec[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if(pairULS && MCpidEle)  hMPtDCA_ULS_rec_MCpidEle[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((pairULS && MCpidEle) && ( (hasCharmAncestor(particle1, particles)  ) || (hasCharmAncestor(particle2, particles)  )) )  hMPtDCA_ULS_rec_MCpidEle_charmTOe[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((pairULS && MCpidEle) && ( (hasBeautyAncestor(particle1, particles) ) || (hasBeautyAncestor(particle2, particles) )) )  hMPtDCA_ULS_rec_MCpidEle_beautyTOe[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((pairULS && MCpidEle) && ( (hasCharmAncestor(particle1, particles) || hasBeautyAncestor(particle1, particles) ) || (hasCharmAncestor(particle2, particles) || hasBeautyAncestor(particle2, particles) )) )  hMPtDCA_ULS_rec_MCpidEle_hfTOe[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((pairULS && MCpidEle) && (std::find(vecLFpdgs.begin(), vecLFpdgs.end(), m1Pid) != vecLFpdgs.end()) && (std::find(vecLFpdgs.begin(), vecLFpdgs.end(), m2Pid) != vecLFpdgs.end()) )  hMPtDCA_ULS_rec_MCpidEle_lfTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((pairULS && MCpidEle) && ( hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles) ) )  hMPtDCA_ULS_rec_MCpidEle_ccTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((pairULS && MCpidEle) && ( hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles) ) )  hMPtDCA_ULS_rec_MCpidEle_bbTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((pairULS && MCpidEle) && ( (hasCharmAncestor(particle1, particles) || hasBeautyAncestor(particle1, particles)) && (hasCharmAncestor(particle2, particles) || hasBeautyAncestor(particle2, particles)) ) )  hMPtDCA_ULS_rec_MCpidEle_hfTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);

      if(!pairULS && MCpidEle) hMPtDCA_LS_rec_MCpidEle[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      // if((!pairULS && !MCpidEle) && (  (fabs(particle1->PID) != 11)||(fabs(particle2->PID)!=11) ))  hMPtDCA_LS_rec_misIDoneLeg[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      // if((!pairULS && !MCpidEle) && (  (fabs(particle1->PID) != 11)&&(fabs(particle2->PID)!=11) ))  hMPtDCA_LS_rec_misIDtwoLeg[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      // if((!pairULS && !MCpidEle) && (  (fabs(particle1->PID) == 211)&&(fabs(particle2->PID)==211) ))  hMPtDCA_LS_rec_misIDPion[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      // if((!pairULS && !MCpidEle) && ( ((fabs(particle1->PID) != 11) && (hasCharmAncestor(particle1, particles) || hasBeautyAncestor(particle1, particles)) ) || ((fabs(particle2->PID)!=11) && (hasCharmAncestor(particle2, particles) || hasBeautyAncestor(particle2, particles)) )) )  hMPtDCA_LS_rec_misIDhf[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);

      if((!pairULS && MCpidEle) && ( (hasCharmAncestor(particle1, particles)  ) || (hasCharmAncestor(particle2, particles)  )) )  hMPtDCA_LS_rec_charmTOe[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((!pairULS && MCpidEle) && ( (hasBeautyAncestor(particle1, particles) ) || (hasBeautyAncestor(particle2, particles) )) )  hMPtDCA_LS_rec_beautyTOe[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((!pairULS && MCpidEle) && ( (hasCharmAncestor(particle1, particles) || hasBeautyAncestor(particle1, particles) ) || (hasCharmAncestor(particle2, particles) || hasBeautyAncestor(particle2, particles) )) )  hMPtDCA_LS_rec_hfTOe[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);

      if((!pairULS && MCpidEle) && ( hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles) ) )  hMPtDCA_LS_rec_ccTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((!pairULS && MCpidEle) && ( hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles) ) )  hMPtDCA_LS_rec_bbTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      // if((!pairULS && MCpidEle) && ( ((hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles)) || (hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles)) )) )  hMPtDCA_LS_rec_hfTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((!pairULS && MCpidEle) && ( (hasCharmAncestor(particle1, particles) || hasBeautyAncestor(particle1, particles)) && (hasCharmAncestor(particle2, particles) || hasBeautyAncestor(particle2, particles)) ) )  hMPtDCA_LS_rec_hfTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((!pairULS && MCpidEle) && ( ((!hasCharmAncestor(particle1, particles)&&!hasCharmAncestor(particle2, particles)) && (!hasBeautyAncestor(particle1, particles)&&!hasBeautyAncestor(particle2, particles)) )) )  hMPtDCA_LS_rec_lfTOee[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if((!pairULS && MCpidEle) && (std::find(vecLFpdgs.begin(), vecLFpdgs.end(), m1Pid) != vecLFpdgs.end()) && (std::find(vecLFpdgs.begin(), vecLFpdgs.end(), m2Pid) != vecLFpdgs.end()) ) hMPtDCA_LS_rec_lfTOee_selectPDG[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      // if((!pairULS && MCpidEle) && ( ((abs(m1Pid) == 111 && abs(m2Pid) == 111) || (abs(m1Pid) == 221 && abs(m2Pid) == 221) || (abs(m1Pid) == 331 && abs(m2Pid) == 331) || (abs(m1Pid) == 223 && abs(m2Pid) == 223) || (abs(m1Pid) == 333 && abs(m2Pid) == 333) || (abs(m1Pid) == 113 && abs(m2Pid) == 113))) ) hMPtDCA_LS_rec_lfTOee_selectPDG[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);

    }
  }
}

// function to ULS or LS pair geneated particles
void dileptonPairingGen(std::vector<GenParticle *> vec_track_neg,std::vector<GenParticle *> vec_track_pos, bool pairULS, Int_t iSce){
  //#########################################
  //    Pairing two "generated" particles
  //#########################################
  // pairing reconstructed ULS or LS pairs
  double vec_track_neg_size = vec_track_neg.size();
  double vec_track_pos_size = vec_track_pos.size();
  if(!pairULS && (vec_track_neg_size == 1)) return;
  int    iPos_start = 0;
  int    iNeg_end = vec_track_neg_size;
  if(!pairULS) iNeg_end = vec_track_neg_size-1;
  TLorentzVector LV1,LV2,LV;
  double dca= 0. ,dca1 = 0.,dca2 = 0.;

  for (int iEle = 0; iEle < iNeg_end; iEle++){
    auto track1 = (Track *) vec_track_neg.at(iEle);
    auto particle1 = (GenParticle *) vec_track_neg.at(iEle);
    auto imother1 = particle1->M1;
    auto mother1 = imother1 != -1 ? (GenParticle *)particles->At(imother1) : (GenParticle *)nullptr;
    auto m1Pid = mother1->PID;
    LV1.SetPtEtaPhiM(particle1->PT,particle1->Eta,particle1->Phi,eMass);

    if(!pairULS) iPos_start = iEle+1;
    for (int iPos = iPos_start; iPos < vec_track_pos_size; iPos++){
      auto track2 = (Track *) vec_track_pos.at(iPos);
      auto particle2 = (GenParticle *) vec_track_pos.at(iPos);
      auto imother2 = particle2->M1;
      auto mother2 = imother2 != -1 ? (GenParticle *)particles->At(imother2) : (GenParticle *)nullptr;
      auto m2Pid = mother2->PID;
      LV2.SetPtEtaPhiM(particle2->PT,particle2->Eta,particle2->Phi,eMass);
      LV = LV1 + LV2;

      dca1 = (track1->D0/track1->ErrorD0);
      dca2 = (track2->D0/track2->ErrorD0);
      dca = sqrt((dca1*dca1 + dca2*dca2) / 2);


      // get weights for both particles to weight the contribution of LF and HF tracks
      Double_t weight_neg = 1.;
      Double_t weight_pos = 1.;
      if(bUseWeights){
      // weight for negative particle
        if ( abs(particle1->PID) == 11 && hasHeavyAncestor(particle1, particles) ) {
          if (particle1->PT < 0.1)                              weight_neg = 1.;
          else if (particle1->PT >= 0.1 && particle1->PT < 0.6) weight_neg = GetWeight(hweightHFtoe_lowerPt, particle1->PT);
          else                                                  weight_neg = GetWeight(hweightHFtoe_higherPt, particle1->PT);
        }
        else if ( abs(particle1->PID) == 11 && !hasHeavyAncestor(particle1, particles) ) {
          if(mother1->PT < 0.1) weight_neg = 0.7;
          else                  weight_neg = GetWeight(hweightLFtoe, mother1->PT);
        }
        else {
          weight_neg = GetWeight(hweightLFtoe, track1->PT);
          if(weight_neg == 0.) weight_neg = 1.;
        }

        // weight for positive particle
        if (abs(particle2->PID) == 11 && hasHeavyAncestor(particle2, particles) ) {
          if (particle2->PT < 0.1)                              weight_pos = 1.;
          else if (particle2->PT >= 0.1 && particle2->PT < 0.6) weight_pos = GetWeight(hweightHFtoe_lowerPt, particle2->PT);
          else                                                  weight_pos = GetWeight(hweightHFtoe_higherPt, particle2->PT);
        }
        else if ( abs(particle2->PID) == 11 && !hasHeavyAncestor(particle2, particles) ) {
          if(mother2->PT < 0.1) weight_pos = 0.7;
          else                  weight_pos = GetWeight(hweightLFtoe, mother2->PT);
        }
        else {
          weight_pos = GetWeight(hweightLFtoe, track2->PT);
          if(weight_pos == 0.) weight_pos = 1.;
        }
      }


      if(pairULS)  hMPtDCA_ULS_gen[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);
      if(!pairULS) hMPtDCA_LS_gen[iSce]->Fill(LV.Mag(),LV.Pt(),dca,weight_neg*weight_pos);

    }
  }
}



// Apply a prefilter on invariant mass and opening angle
void PreFilter(std::vector<Track *>* vec_track_neg, std::vector<Track *>* vec_track_pos, Int_t iPID_scenario, std::vector<TH2F*> vecPFhistos){
  //#########################################
  //    PreFilter for pairing two tracks
  //#########################################
  // pairing reconstructed ULS or LS pairs
  std::vector<Int_t> vecDelPos, vecDelEle;
  TLorentzVector LV1,LV2,LV;
  // double dca= 0. ,dca1 = 0.,dca2 = 0.;

  for (int iEle = 0; iEle < vec_track_neg->size() ; iEle++){
    auto track1 = (Track *) vec_track_neg->at(iEle);
    auto particle1 = (GenParticle *)track1->Particle.GetObject();
    auto imother1 = particle1->M1;
    auto mother1 = imother1 != -1 ? (GenParticle *)particles->At(imother1) : (GenParticle *)nullptr;
    auto m1Pid = mother1->PID;
    LV1.SetPtEtaPhiM(track1->PT,track1->Eta,track1->Phi,eMass);

    if( LV1.Pt() < pf_pt_cut[iPID_scenario] ) {
      vec_track_neg->erase(vec_track_neg->begin()+iEle);
      iEle--;
      continue;
    }

    for (auto iPos = 0; iPos < vec_track_pos->size(); iPos++){
      auto track2 = (Track *) vec_track_pos->at(iPos);
      auto particle2 = (GenParticle *)track2->Particle.GetObject();
      auto imother2 = particle2->M1;
      auto mother2 = imother2 != -1 ? (GenParticle *)particles->At(imother2) : (GenParticle *)nullptr;
      auto m2Pid = mother2->PID;
      LV2.SetPtEtaPhiM(track2->PT,track2->Eta,track2->Phi,eMass);
      LV = LV1 + LV2;
      double opAngle = fabs(TVector2::Phi_mpi_pi(LV1.Angle(LV2.Vect())));
      if(TMath::IsNaN(opAngle)){
          continue;
      }

      if( (imother1 == imother2) && (abs(m1Pid)==111) ) vecPFhistos.at(0)->Fill(LV.M(), opAngle);

      // apply prefilter pt cut from the corresponding PIDscenario
      if( LV2.Pt() < pf_pt_cut[iPID_scenario] ) {
        vec_track_pos->erase(vec_track_pos->begin()+iPos);
        iPos--;
        continue;
      }


      if( !((LV.M() < pf_mass_cut[iPID_scenario]) && (opAngle < pf_opAngle_cut[iPID_scenario])) && ((imother1 == imother2) && (abs(m1Pid)==111)) ) vecPFhistos.at(1) ->Fill(LV.M(), opAngle);

      if( (LV.M() < pf_mass_cut[iPID_scenario]) && (opAngle < pf_opAngle_cut[iPID_scenario]) ) {
        vecDelEle.push_back(iEle);
        vecDelPos.push_back(iPos);
      }
    }
  }

  std::sort(vecDelEle.begin(), vecDelEle.end());
  std::sort(vecDelPos.begin(), vecDelPos.end());

  vecDelEle.erase(std::unique(vecDelEle.begin(), vecDelEle.end()), vecDelEle.end());
  vecDelPos.erase(std::unique(vecDelPos.begin(), vecDelPos.end()), vecDelPos.end());

  for (int i = vecDelEle.size()-1; i >= 0; i--){
    // erase element at position stored in index vector vecDelEle
    vec_track_neg->erase(vec_track_neg->begin() + vecDelEle.at(i));
  }
  for (int j = vecDelPos.size()-1; j >= 0; j--){
    // erase element at position stored in index vector vecDelPos
    vec_track_pos->erase(vec_track_pos->begin() + vecDelPos.at(j));
  }

  vecDelEle.clear();
  vecDelPos.clear();
}


void PreFilter(std::vector<GenParticle *>* vec_particle,std::vector<Track *>* vec_track, std::vector<Int_t> vecMotherID, std::vector<std::tuple<GenParticle*, GenParticle*>> vecTuples, Int_t iPID_scenario, std::vector<TH2F*> vecPFhistos){
  //#################################################
  //  PreFilter for pairing GenParticles and Tracks
  //#################################################
  // pairing reconstructed ULS or LS pairs
  std::vector<Int_t> vecDelPos, vecDelEle;
  TLorentzVector LV1,LV2,LV;
  // double dca= 0. ,dca1 = 0.,dca2 = 0.;

  for (int iParticle = 0; iParticle < vec_particle->size() ; iParticle++){
    auto particle1 = (GenParticle *) vec_particle->at(iParticle);
    auto imother1 = particle1->M1;
    auto mother1 = imother1 != -1 ? (GenParticle *)particles->At(imother1) : (GenParticle *)nullptr;
    auto m1Pid = mother1->PID;
    LV1.SetPtEtaPhiM(particle1->PT,particle1->Eta,particle1->Phi,eMass);

    if( LV1.Pt() < pf_pt_cut[iPID_scenario] ) {
      vec_particle->erase(vec_particle->begin()+iParticle);
      iParticle--;
      continue;
    }

    for (auto iTrack = 0; iTrack < vec_track->size(); iTrack++){
      auto track2 = (Track *) vec_track->at(iTrack);
      auto particle2 = (GenParticle *)track2->Particle.GetObject();
      auto imother2 = particle2->M1;
      auto mother2 = imother2 != -1 ? (GenParticle *)particles->At(imother2) : (GenParticle *)nullptr;
      auto m2Pid = mother2->PID;
      LV2.SetPtEtaPhiM(track2->PT,track2->Eta,track2->Phi,eMass);
      LV = LV1 + LV2;
      double opAngle = fabs(TVector2::Phi_mpi_pi(LV1.Angle(LV2.Vect())));
      if(TMath::IsNaN(opAngle)){
          continue;
      }

      bool firstMother = kFALSE;
      bool twoMothers = kFALSE;
      std::tuple<GenParticle*,GenParticle*> tup_genpart12 = {particle1,particle2};
      std::tuple<GenParticle*,GenParticle*> tup_genpart21 = {particle2,particle1};
      // fill histogram of Pi0 dalitz decays, before prefilter
      if( !(std::find(vecMotherID.begin(), vecMotherID.end(), imother1) != vecMotherID.end()) && (imother1 == imother2) && (abs(m1Pid)==111) ) {
        vecPFhistos.at(0)->Fill(LV.M(), opAngle);
        vecMotherID.push_back(imother1);
        firstMother = kTRUE;
      }
      // fill histogram of background, before prefilter
      else if( !(std::find(vecTuples.begin(), vecTuples.end(), tup_genpart21) != vecTuples.end()) && (imother1 != imother2) ){
        vecPFhistos.at(2)->Fill(LV.M(), opAngle);

        // fill histogram of pairs with HF and Pi0Dalitz track, before prefilter
        if( (hasHeavyAncestor(particle1,particles) && (abs(m2Pid)==111)) || ((abs(m1Pid)==111) && hasHeavyAncestor(particle2,particles)) )
        vecPFhistos.at(4)->Fill(LV.M(), opAngle);

        vecTuples.push_back(tup_genpart12);
        twoMothers = kTRUE;
      }


      // apply prefilter pt cut from the corresponding PIDscenario
      if( LV2.Pt() < pf_pt_cut[iPID_scenario] ) {
        vec_track->erase(vec_track->begin()+iTrack);
        iTrack--;
        continue;
      }


      if( !((LV.M() < pf_mass_cut[iPID_scenario]) && (opAngle < pf_opAngle_cut[iPID_scenario])) ){
        // fill histogram of Pi0 dalitz decays, after prefilter
        if( firstMother )
           vecPFhistos.at(1) ->Fill(LV.M(), opAngle);

        // fill histogram of background, after prefilter
        else if( twoMothers ){
          vecPFhistos.at(3) ->Fill(LV.M(), opAngle);

          // fill histogram of pairs with HF and Pi0Dalitz track, after prefilter
          if( (hasHeavyAncestor(particle1,particles) && (abs(m2Pid)==111)) || ((abs(m1Pid)==111) && hasHeavyAncestor(particle2,particles)) )
          vecPFhistos.at(5)->Fill(LV.M(), opAngle);
        }
      }

      if( (LV.M() < pf_mass_cut[iPID_scenario]) && (opAngle < pf_opAngle_cut[iPID_scenario]) ) {
        vecDelEle.push_back(iParticle);
        vecDelPos.push_back(iTrack);
      }
    }
  }

  std::sort(vecDelEle.begin(), vecDelEle.end());
  std::sort(vecDelPos.begin(), vecDelPos.end());

  vecDelEle.erase(std::unique(vecDelEle.begin(), vecDelEle.end()), vecDelEle.end());
  vecDelPos.erase(std::unique(vecDelPos.begin(), vecDelPos.end()), vecDelPos.end());

  for (int i = vecDelEle.size()-1; i >= 0; i--){
    // erase element at position stored in index vector vecDelEle
    vec_particle->erase(vec_particle->begin() + vecDelEle.at(i));
  }
  for (int j = vecDelPos.size()-1; j >= 0; j--){
    // erase element at position stored in index vector vecDelPos
    vec_track->erase(vec_track->begin() + vecDelPos.at(j));
  }

  vecDelEle.clear();
  vecDelPos.clear();
}



void ApplyPID(std::vector<Track *>* vecTracks_pfPID, std::vector<Track *>* vecTracks, Int_t iPID_scenario) {
  //#########################
  //    standard  PID
  //#########################
  // prefilter should already been applied
  // take provided vector and loop over all particles
  // loop over all the PID scenarios to apply differnet PID cuts


  for (int itrack = 0; itrack < vecTracks_pfPID->size(); itrack++) {
    // get Track
    Track* track = (Track *) vecTracks_pfPID->at(itrack);
    auto particle = (GenParticle *)track->Particle.GetObject();

    //make outer TOF PID
    std::array<float, 5> PIDdeltatTOF, PIDnsigmaTOF;
    toflayer.makePID(*track, PIDdeltatTOF, PIDnsigmaTOF);
    //make RICH PID
    std::array<float, 5> PIDdeltaangleRICH, PIDnsigmaRICH;
    richdetector.makePID(*track, PIDdeltaangleRICH, PIDnsigmaRICH);


    //setting variables for i-th PID scenario
    bool   i_useTOFPID    = useTOFPID[iPID_scenario];
    bool   i_useRICHPID   = useRICHPID[iPID_scenario];

    double i_SigmaTOFEle  = nSigmaTOFEle[iPID_scenario];
    double i_SigmaTOFPi   = nSigmaTOFPi[iPID_scenario];
    double i_SigmaRICHEle = nSigmaRICHEle[iPID_scenario];
    double i_SigmaRICHPi  = nSigmaRICHPi[iPID_scenario];
    double ps_pt_acc      = ps_EleAccep_pt_cut[iPID_scenario];


    // Apply Preshower electron selection
    bool isElectronPreShower = false;
    if(usePreShPID[iPID_scenario]){
      if (preshower.hasPreShower(*track) ) isElectronPreShower = preshower.isElectron(*track,1800);
    }

    bool bool_doPID = false;
    if(                       i_useRICHPID && !bUseiTOF && doPID(track, i_useTOFPID, i_useRICHPID, tof_EleAccep_p_cut, tof_PionRej_p_cut, rich_PionRejection_p_cut, i_SigmaTOFEle, i_SigmaTOFPi, i_SigmaRICHEle, i_SigmaRICHPi, PIDnsigmaTOF, PIDnsigmaRICH) ) bool_doPID = true;
    else if(  i_useTOFPID && !i_useRICHPID && !bUseiTOF && doTOFPID(track, i_useTOFPID, tof_EleAccep_p_cut, tof_PionRej_p_cut, i_SigmaTOFEle, i_SigmaTOFPi, PIDnsigmaTOF) ) bool_doPID = true;
    else if(  i_useTOFPID &&                   bUseiTOF && doTOFPID(track, i_useTOFPID, itof_EleAccep_p_cut, tof_PionRej_p_cut, i_SigmaTOFEle, i_SigmaTOFPi, PIDnsigmaTOF) ) bool_doPID = true;

    // ----------- PID strategy -----------
    // if:
    // check if any of the detectors is used
    // if so -> check wheather (track pt < PS pt acceptance and PS identifies track as electron)
    //          AND (track not identified as electron by TOF and/or RICH)
    // then  -> reject track
    //
    // else if: if no detector is used AND no monte carlo PID is used
    // then  -> reject track

    if ( (usePreShPID[iPID_scenario] || i_useTOFPID || i_useRICHPID) && ( ((isElectronPreShower == false) || (track->PT < ps_pt_acc) ) && (bool_doPID == false)) ) { /*printf("Continue in standard PID: Use min one detector, not accepted by any detector \n");*/ continue;}
    else if ( !(usePreShPID[iPID_scenario] || i_useTOFPID || i_useRICHPID) && !useMCPID[iPID_scenario] ) {/*printf( "Continue in standard PID: No detector used for PID, no MC PID used \n");*/ continue;}
    // if ( (isElectronPreShower == false) && (bool_doPID == false) &&  !useMCPID[iPID_scenario] ) continue;
    if (useMCPID[iPID_scenario] && (abs(track->PID) != 11)) { /*printf("Continue in standard PID: Use MC PID, not an electron \n");*/ continue;}

    // fill track vector
    vecTracks->push_back(track);       //  vector of all reconstructed negative tracks
  }
}
