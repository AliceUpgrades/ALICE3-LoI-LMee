R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libDelphesO2)

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




bool bSmear    = true;
bool bUsePreSh = false;
bool bUseTOF   = true;
bool bUseRICH  = true;
double Bz = 0.2;            // becomes overwritten by the generateEfficiencies.sh skript
double eMass = 0.000511;


// Cinematic cuts on tracks
// double PtCut = 0.0;     // open cuts
// double PtCut = 0.04;    // right now pt cut in cut variations below // if B=0.2T becomes overwritten by the generateEfficiencies.sh skript
// double PtCut = 0.08;    // right now pt cut in cut variations below // if B=0.5T becomes overwritten by the generateEfficiencies.sh skript
double EtaCut = 1.1;
// double EtaCut = 10.;    // open cuts

// --- binning ---
Double_t pt_binning[]  = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0};
// Double_t eta_binning[]  = {-5.0,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.0,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0};
Int_t nEtaBins = 100;    Double_t maxEta = 5.;
Int_t nPhiBins = 180;    Double_t maxPhi = 2*TMath::Pi();
// ptee
const Int_t n_ptee_bin_c = 202;
Double_t ptee_bin_c[n_ptee_bin_c+1] = {};
// mee
Double_t mee_bin_c[401];
Int_t n_mee_bin_c = 400;
// DCA
Double_t dca_bin_c[201];
int n_dca_bin_c = 200;


// TOF
double tof_radius = 100.; // [cm]
double tof_length = 200.; // [cm]
double tof_sigmat = 0.02; // [ns]
double tof_sigma0 = 0.20; // [ns]

// TOF ele pt acceptance
double tof_EleAccep_p_cut = 0.6;  // [GeV/c]
// TOF additional Pion/Muon rejection
double tof_PionRej_p_cut = 0.4; // [GeV/c] above this value the nSigmaRICHEle_forTOFPID cut is apllied. Thus the TOF is only accepting particles within nSigmaRICHEle and nSigmaTOFele
double nSigmaRICHEle_forTOFPID = 3.; //


// RICH
double rich_radius = 100.;
double rich_length = 200.;

// RICH ele pt acceptance
// double rich_Ele_pt_cut = 1.8; // [GeV/c]
double rich_PionRejection_p_cut = 1.0; // [GeV/c]



// PID Scenarios
bool useTOFPID[]        = {/*kTRUE,   kTRUE,     */  kTRUE,  kTRUE,  kTRUE,   kTRUE};
bool useRICHPID_B2[]    = {/*kFALSE,  kFALSE,    */  kFALSE, kTRUE,  kTRUE,   kTRUE};
bool useRICHPID_B5[]    = {/*kTRUE,   kTRUE,     */  kFALSE, kTRUE,  kTRUE,   kTRUE};
bool usePreShPID[]      = {/*kFALSE,  kFALSE,    */  kFALSE, kFALSE, kFALSE,  kFALSE};
// TOF cuts on tracks
double nSigmaTOFEle[]   = {/*3.0,     3.0,       */    3.0,    3.0,    3.0,     3.0};
double nSigmaTOFPi[]    = {/*3.0,     3.0,       */    3.0,    3.0,    3.0,     3.0};
// RICH cuts on tracks
double nSigmaRICHEle[]  = {/*3.0,     3.0,       */    3.0,    3.0,    3.0,     3.0};
double nSigmaRICHPi[]   = {/*4.0,     4.0,       */    3.0,    3.0,    3.5,     4.0};
// TOF pte > 0.04 B = 0.2 T (highest priority) or TOF pte > 0.08 B = 0.5 T
// TOF RICH pte > 0.2 B = 0.5 T (highest priority) or TOF RICH pte > 0.08  B = 0.5 T
double PtCut02[]        = {/*0.04,    0.08,      */   0.03,   0.03,   0.03,    0.03};
double PtCut05[]        = {/*0.2,     0.08,       */  0.03,   0.03,   0.03,    0.03};




// bool doPID(Track * tr, bool useTOF, bool useRICH, bool usePreSh, double p_tofMaxAcc, double p_tofPionRej, double p_richPionRej, double nSigmaTOFele, double nSigmaTOFpi, double nSigmaRICHele, double nSigmaRICHpi, o2::delphes::TOFLayer toflayer, o2::delphes::RICHdetector richdetector, std::array<float, 5> PIDnsigmaTOF, std::array<float, 5> PIDnsigmaRICH){
//   double p = tr->P;
//
//   bool TOFpid = kFALSE;
//   bool RICHpid = kFALSE;
//   //apply TOF PID cuts
//   // require hasTOF and p < 0.6
//   if(useTOF && (toflayer.hasTOF(*tr)) && (p < p_tofMaxAcc)) {
//     // electron acceptance
//     // If  p > 0.4 GeV/c
//     if (p > p_tofPionRej) {
//       if( (fabs(PIDnsigmaTOF[0]) < nSigmaTOFele) && (fabs(PIDnsigmaRICH[0]) < nSigmaRICHEle_forTOFPID) ) TOFpid = true; // is within 3 sigma of the electron band (TOF+RICH)
//     }
//     // If p < 0.4 GeV/c
//     else if(p <= p_tofPionRej){  // p_tofPionRej = 0.4
//       if(fabs(PIDnsigmaTOF[0]) < nSigmaTOFele) TOFpid = true; // is within 3 sigma of the electron band (TOF)
//     }
//     else cout << "!!! something is going wrong !!! " << endl;
//
//     // pion rejection
//     if(fabs(PIDnsigmaTOF[2]) < nSigmaTOFpi) TOFpid = false; // is within 3 sigma of the pion band (TOF)
//   }
//
//   //apply RICH PID cuts, require hasRICH
//   if(useRICH && richdetector.hasRICH(*tr)) {
//     // electron acceptance
//     if(fabs(PIDnsigmaRICH[0]) < nSigmaRICHele) RICHpid = true; // is within 3 sigma of the electron band (RICH)
//     // pion rejection (p > 1.0 GeV/c)
//     if( (fabs(PIDnsigmaRICH[2]) < nSigmaRICHpi) && (p > p_richPionRej) ) RICHpid = false; // is within 3 sigma of the pion band (RICH)
//   }
//
//   // if ((TOFpid || RICHpid) && (p > 0.6 && fabs(PIDnsigmaRICH[0] > 3)) ) {
//   //   cout << __LINE__ << " bool of TOF = " << TOFpid  << ", bool of RICH = " << RICHpid << endl;
//   //   cout << " track p = " << p << endl;
//   //   cout << " track pt = " << tr->PT << endl;
//   //   cout << " PIDnsigmaTOF[0] = " << PIDnsigmaTOF[0] << endl;
//   //   cout << " PIDnsigmaTOF[2] = " << PIDnsigmaTOF[2] << endl;
//   //   cout << " PIDnsigmaRICH[0] = " <<  PIDnsigmaRICH[0] << endl;
//   //   cout << " PIDnsigmaRICH[2] = " << PIDnsigmaRICH[2] << endl << endl;
//   // }
//
//   if (!(RICHpid || TOFpid) /*&& !PreShpid*/) return false; // check if TOF or RICH signal is true.
//   else return true;
//   // ################## end of PID selection ##################
// }



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


// kinematic cuts for tracks
bool etaCut(Track *tr){
  // evaluate as true if criterion is passed
  bool eta = abs(tr->Eta) < EtaCut;
  // all have to be true
  return (eta);
}

// kinematic cuts for tracks
bool kineCuts(Track *tr, Int_t iSce){
  // check pt and eta for track
  // evaluate as true if criterion is passed
  bool pt;
  if(Bz == 0.2) pt = tr->PT > PtCut02[iSce];
  else if(Bz == 0.5) pt = tr->PT > PtCut05[iSce];
  bool eta = abs(tr->Eta) < EtaCut;
  // all have to be true
  return (pt && eta);
}

// kinematic cuts for generated particles
bool kineCuts(GenParticle *pa, Int_t iSce){
  // check pt and eta for particle
  // evaluate as true if criterion is passed
  bool pt;
  if(Bz == 0.2) pt = pa->PT > PtCut02[iSce];
  else if(Bz == 0.5) pt = pa->PT > PtCut05[iSce];
  bool eta = abs(pa->Eta) < EtaCut;
  // all have to be true
  return (pt && eta);
}

// kinematic cuts for generated smeared particles
bool kineCuts(TLorentzVector LV, Int_t iSce){
  // check pt and eta for particle
  // evaluate as true if criterion is passed
  bool pt;
  if(Bz == 0.2) pt = LV.Pt() > PtCut02[iSce];
  else if(Bz == 0.5) pt = LV.Pt() > PtCut05[iSce];
  bool eta = abs(LV.Eta()) < EtaCut;
  // all have to be true
  return (pt && eta);
}

bool hasHeavyAncestor(GenParticle *particle, TClonesArray *particles){
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle *)particles->At(imother);
  auto pid = mother->PID;

  switch (abs(pid)) {
  case 411:  // D+
  case 421:  // D0
  case 431:  // Ds+
  case 4122: // Lambdac+
  case 4132: // Xic0
  case 4232: // Xic+
  case 511:  // B0
  case 521:  // B+
  case 531:  // Bs0
  case 541:  // Bc+
  case 5122: // Lambdab0
  case 5132: // Xib-
  case 5232: // Xib0
  case 5332: // Omegab-
    return true;
  }
  return hasHeavyAncestor(mother, particles);
}

bool isBeauty(int pid){
  switch (abs(pid)) {
  case 511:  // B0
  case 521:  // B+
  case 531:  // Bs0
  case 541:  // Bc+
  case 5122: // Lambdab0
  case 5132: // Xib-
  case 5232: // Xib0
  case 5332: // Omegab-
    return true;
  }
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
  switch (abs(pdg)) {
  case 411:  // D+
  case 421:  // D0
  case 431:  // Ds+
  case 4122: // Lambdac+
  case 4132: // Xic0
  case 4232: // Xic+
    return true;
  }
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


Bool_t IsStable(Int_t pdg)
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

  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstable; i++) {
   if (TMath::Abs(pdg) == TMath::Abs(pdgStable[i])) {
     isStable = kTRUE;
     break;
   }
  }
  return isStable;
}


// Resolution
// TObjArray *fArrResoPt = 0x0;
// TObjArray *fArrResoEta = 0x0;
// TObjArray *fArrResoPhi_Pos = 0x0;
// TObjArray *fArrResoPhi_Neg = 0x0;
TString smearingfile;

// void ReadResoFile(TFile *fRes)
// {
//   //
//   // Set resolution arrays
//   //
//
//   TObjArray *ArrResoPt = 0x0;
//   TObjArray *ArrResoEta = 0x0;
//   TObjArray *ArrResoPhi_Pos = 0x0;
//   TObjArray *ArrResoPhi_Neg = 0x0;
//
//   // Set resolutions: default Run 2
//   if(fRes && fRes->IsOpen()){
//     // Oton
//     ArrResoPt = (TObjArray *)fRes->Get("RelPtResArrCocktail");
//     ArrResoEta = (TObjArray *)fRes->Get("EtaResArrVsPt");
//     ArrResoPhi_Pos = (TObjArray *)fRes->Get("PhiPosResArrVsPt");
//     ArrResoPhi_Neg = (TObjArray *)fRes->Get("PhiEleResArrVsPt");
//   }
//   else {
//     printf("no resolution file given!");
//   }
//
//   fArrResoPt = ArrResoPt;
//   fArrResoEta = ArrResoEta;
//   fArrResoPhi_Pos = ArrResoPhi_Pos;
//   fArrResoPhi_Neg = ArrResoPhi_Neg;
//
//   if(!fArrResoPt || !fArrResoEta || !fArrResoPhi_Pos || !fArrResoPhi_Neg) {
//     printf("Did not find the resolution arrays\n");
//   }
//
// }

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
std::string title2d = ";#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c})";
std::string title3d = ";#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});DCA_{ee} (mm)";



void dileptonPairingRec(std::vector<Track *> vec_track_neg,std::vector<Track *> vec_track_pos, bool pairULS, bool MCpidEle, TClonesArray *particles, Int_t iSce);
void dileptonPairingGen(std::vector<GenParticle *> vec_track_neg,std::vector<GenParticle *> vec_track_pos, bool pairULS, TClonesArray *particles, Int_t iSce);

// Pair histograms
TH3F* hMPtDCA_ULS_gen[3];
// TH3F* hMPtDCA_ULS_gen_primary[3];
// TH3F* hMPtDCA_ULS_gen_heavy[3];
// TH3F* hMPtDCA_ULS_gen_charm[3];
// TH3F* hMPtDCA_ULS_gen_beauty[3];
TH3F* hMPtDCA_ULS_rec[3];
// TH3F* hMPtDCA_ULS_rec_primary[3];
// TH3F* hMPtDCA_ULS_rec_heavy[3];
// TH3F* hMPtDCA_ULS_rec_charm[3];
// TH3F* hMPtDCA_ULS_rec_beauty[3];
TH3F* hMPtDCA_ULS_rec_MCpidEle[3];
// TH3F* hMPtDCA_ULS_rec_MCpidEle_primary[3];
// TH3F* hMPtDCA_ULS_rec_MCpidEle_heavy[3];
// TH3F* hMPtDCA_ULS_rec_MCpidEle_charm[3];
// TH3F* hMPtDCA_ULS_rec_MCpidEle_beauty[3];
TH3F* hMPtDCA_LS_gen[3];
// TH3F* hMPtDCA_LS_gen_primary[3];
// TH3F* hMPtDCA_LS_gen_heavy[3];
// TH3F* hMPtDCA_LS_gen_charm[3];
// TH3F* hMPtDCA_LS_gen_beauty[3];
TH3F* hMPtDCA_LS_rec[3];
// TH3F* hMPtDCA_LS_rec_primary[3];
// TH3F* hMPtDCA_LS_rec_heavy[3];
// TH3F* hMPtDCA_LS_rec_charm[3];
// TH3F* hMPtDCA_LS_rec_beauty[3];
TH3F* hMPtDCA_LS_rec_MCpidEle[3];
// TH3F* hMPtDCA_LS_rec_MCpidEle_primary[3];
// TH3F* hMPtDCA_LS_rec_MCpidEle_heavy[3];
// TH3F* hMPtDCA_LS_rec_MCpidEle_charm[3];
// TH3F* hMPtDCA_LS_rec_MCpidEle_beauty[3];
TH3F* hMPtDCA_LS_rec_misIDoneLeg[3];
TH3F* hMPtDCA_LS_rec_misIDtwoLeg[3];
TH3F* hMPtDCA_LS_rec_misIDPion[3];
TH3F* hMPtDCA_LS_rec_misIDhf[3];


void anaEEstudy(
    const char *inputFile = "delphes100k.root", // one fo the LF and charm part. Has 1M charm+beauty events
    const char *outputFile = "anaEEstudyLFcc.root" // merge output files after analysis was run to keep file size moderate
  )
{

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
  for (int iBin = 0; iBin < nEtaBins+1; iBin++) {eta_binning[iBin] = -maxEta + iBin * 2.* maxEta/nEtaBins; /* std::cout << eta_binning[iBin] << " eta bin border " << std::endl;  */}
  for (int iBin = 0; iBin < nPhiBins+1; iBin++) {phi_binning[iBin] = 0. + iBin * maxPhi/nPhiBins;      /* std::cout << phi_binning[iBin] << " phi bin border " << std::endl;  */}

  // ptee
  for(Int_t i=0  ;i<4   ;i++) {
    ptee_bin_c[i] = 0.025 * (i-  0) +  0.0;//from 0 to 0.075 GeV/c, every 0.025 GeV/c
    //printf("bin %f for %d\n",ptee_bin_c[i],i);
  }
  for(Int_t i=4 ;i<=202  ;i++) {
    ptee_bin_c[i] = 0.05  * (i- 4) +  0.1;//from 0.1 to 10 GeV/c, evety 0.05 GeV/c
    //printf("bin %f for %d\n",ptee_bin_c[i],i);
  }
  // mee
  for(Int_t k = 0 ; k < 401; k++){
    mee_bin_c[k] = k*0.01; // 4./400.
  }
  // DCA
  for(int k = 0 ; k < 201; k++){
    dca_bin_c[k] = k*0.1; // 4./400.
  }


  // Get pointers to branches used in this analysis
  auto events = treeReader->UseBranch("Event");
  auto tracks = treeReader->UseBranch("Track");
  auto particles = treeReader->UseBranch("Particle");


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


  // TOF layer
  o2::delphes::TOFLayer toflayer;
  toflayer.setup(tof_radius, tof_length, tof_sigmat,tof_sigma0);

  // RICH detector
  o2::delphes::RICHdetector richdetector;
  richdetector.setup(rich_radius, rich_length);
  richdetector.setIndex(1.03);
  richdetector.setRadiatorLength(2.);
  richdetector.setEfficiency(0.4);
  richdetector.setSigma(7.e-3);

  // smearer
  // look up tables are selected in generateEfficiencies.sh script
  o2::delphes::TrackSmearer smearer;
  smearer.loadTable(11, "./lutCovm.el.dat");
  smearer.loadTable(13, "./lutCovm.mu.dat");
  smearer.loadTable(211, "./lutCovm.pi.dat");
  smearer.loadTable(321, "./lutCovm.ka.dat");
  smearer.loadTable(2212, "./lutCovm.pr.dat");

  smearer.useEfficiency(true);



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
  auto hBeforeSmearing_Pt_Eta_Phi_rec = new TH3F("hBeforeSmearing_Pt_Eta_Phi_rec", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  auto hAfterSmearing_Pt_Eta_Phi_rec = new TH3F("hAfterSmearing_Pt_Eta_Phi_rec", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hAfterKineCuts_Pt_Eta_Phi_rec[nPIDscenarios];
  TH3F* hTrack_ElePos_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hNegTrack_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hPosTrack_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hAllTracks_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Ele_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Pos_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Muon_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Pion_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Kaon_Rec_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Proton_Rec_Pt_Eta_Phi[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_primary_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_primary_Ele_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_primary_Pos_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_hf_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_hf_Ele_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_hf_Pos_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_charm_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_charm_Ele_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_charm_Pos_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_beauty_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_beauty_Ele_rec[nPIDscenarios];
  // TH3F* hPt_Eta_Phi_beauty_Pos_rec[nPIDscenarios];

  for (int j = 0; j < nPIDscenarios; ++j){
    hAfterKineCuts_Pt_Eta_Phi_rec[j] = new TH3F(Form("hAfterKineCuts_Pt_Eta_Phi_rec_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_ElePos_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hNegTrack_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hNegTrack_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hPosTrack_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hPosTrack_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hAllTracks_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hAllTracks_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Ele_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Ele_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Pos_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pos_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Muon_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Muon_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Pion_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pion_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Kaon_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Kaon_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    hTrack_Proton_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Proton_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_primary_rec[j] = new TH3F(Form("hPt_Eta_Phi_primary_rec_sce%i",j+1), ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_primary_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_primary_Ele_rec_sce%i",j+1), ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_primary_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_primary_Pos_rec_sce%i",j+1), ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_hf_rec[j] = new TH3F(Form("hPt_Eta_Phi_hf_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_hf_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_hf_Ele_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_hf_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_hf_Pos_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_charm_rec[j] = new TH3F(Form("hPt_Eta_Phi_charm_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_charm_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_charm_Ele_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_charm_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_charm_Pos_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_beauty_rec[j] = new TH3F(Form("hPt_Eta_Phi_beauty_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_beauty_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_beauty_Ele_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
    // hPt_Eta_Phi_beauty_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_beauty_Pos_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  }


  // Pair histograms
  for (int j = 0; j < nPIDscenarios; ++j){
  // Pair histograms ULS
    hMPtDCA_ULS_gen[j]         = new TH3F(Form("hMPtDCA_ULS_gen_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_gen_primary[j] = new TH3F(Form("hMPtDCA_ULS_gen_primary_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_gen_heavy[j]   = new TH3F(Form("hMPtDCA_ULS_gen_heavy_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_gen_charm[j]   = new TH3F(Form("hMPtDCA_ULS_gen_charm_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_gen_beauty[j]  = new TH3F(Form("hMPtDCA_ULS_gen_beauty_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_ULS_rec[j]         = new TH3F(Form("hMPtDCA_ULS_rec_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_primary[j] = new TH3F(Form("hMPtDCA_ULS_rec_primary_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_heavy[j]   = new TH3F(Form("hMPtDCA_ULS_rec_heavy_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_charm[j]   = new TH3F(Form("hMPtDCA_ULS_rec_charm_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_beauty[j]  = new TH3F(Form("hMPtDCA_ULS_rec_beauty_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_ULS_rec_MCpidEle[j]         = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_MCpidEle_primary[j] = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_primary_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_MCpidEle_heavy[j]   = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_heavy_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_MCpidEle_charm[j]   = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_charm_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_ULS_rec_MCpidEle_beauty[j]  = new TH3F(Form("hMPtDCA_ULS_rec_MCpidEle_beauty_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  // Pair histograms LS
    hMPtDCA_LS_gen[j]         = new TH3F(Form("hMPtDCA_LS_gen_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_gen_primary[j] = new TH3F(Form("hMPtDCA_LS_gen_primary_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_gen_heavy[j]   = new TH3F(Form("hMPtDCA_LS_gen_heavy_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_gen_charm[j]   = new TH3F(Form("hMPtDCA_LS_gen_charm_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_gen_beauty[j]  = new TH3F(Form("hMPtDCA_LS_gen_beauty_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_LS_rec[j]         = new TH3F(Form("hMPtDCA_LS_rec_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_primary[j] = new TH3F(Form("hMPtDCA_LS_rec_primary_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_heavy[j]   = new TH3F(Form("hMPtDCA_LS_rec_heavy_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_charm[j]   = new TH3F(Form("hMPtDCA_LS_rec_charm_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_beauty[j]  = new TH3F(Form("hMPtDCA_LS_rec_beauty_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_LS_rec_MCpidEle[j]         = new TH3F(Form("hMPtDCA_LS_rec_MCpidEle_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_MCpidEle_primary[j] = new TH3F(Form("hMPtDCA_LS_rec_MCpidEle_primary_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_MCpidEle_heavy[j]   = new TH3F(Form("hMPtDCA_LS_rec_MCpidEle_heavy_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_MCpidEle_charm[j]   = new TH3F(Form("hMPtDCA_LS_rec_MCpidEle_charm_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    // hMPtDCA_LS_rec_MCpidEle_beauty[j]  = new TH3F(Form("hMPtDCA_LS_rec_MCpidEle_beauty_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_LS_rec_misIDoneLeg[j]      = new TH3F(Form("hMPtDCA_LS_rec_misIDoneLeg_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_LS_rec_misIDtwoLeg[j]      = new TH3F(Form("hMPtDCA_LS_rec_misIDtwoLeg_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_LS_rec_misIDPion[j]        = new TH3F(Form("hMPtDCA_LS_rec_misIDPion_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
    hMPtDCA_LS_rec_misIDhf[j]          = new TH3F(Form("hMPtDCA_LS_rec_misIDhf_sce%i",j+1),title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  }




  TList* lRecPIDscenario[nPIDscenarios];
  for (int iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    lRecPIDscenario[iPIDscenario] = new TList();
    lRecPIDscenario[iPIDscenario]->SetName(Form("PIDscenario_%i", iPIDscenario+1));
    lRecPIDscenario[iPIDscenario]->SetOwner();
    lRecPIDscenario[iPIDscenario]->Add(hAfterKineCuts_Pt_Eta_Phi_rec[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hTrack_ElePos_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hNegTrack_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hPosTrack_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hAllTracks_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hTrack_Ele_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hTrack_Pos_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hTrack_Muon_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hTrack_Pion_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hTrack_Kaon_Rec_Pt_Eta_Phi[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hTrack_Proton_Rec_Pt_Eta_Phi[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_primary_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_primary_Ele_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_primary_Pos_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_hf_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_hf_Ele_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_hf_Pos_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_charm_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_charm_Ele_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_charm_Pos_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_beauty_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_beauty_Ele_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_beauty_Pos_rec[iPIDscenario]);


    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_primary[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_heavy[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_charm[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_beauty[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_MCpidEle[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_MCpidEle_primary[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_MCpidEle_heavy[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_MCpidEle_charm[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_rec_MCpidEle_beauty[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_primary[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_heavy[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_charm[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_beauty[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_MCpidEle[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_MCpidEle_primary[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_MCpidEle_heavy[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_MCpidEle_charm[iPIDscenario]);
    // lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_MCpidEle_beauty[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_misIDoneLeg[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_misIDtwoLeg[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_misIDPion[iPIDscenario]);
    lRecPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_rec_misIDhf[iPIDscenario]);
    // listRecTracks->Add(lRecPIDscenario);
  }


  TH3F* hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts  = new TH3F("hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts  = new TH3F("hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);


  TH3F* hTrack_All_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Muon_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Pion_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Kaon_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Proton_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Ele_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Pos_Gen_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_ElePos_GenSmeared_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Ele_GenSmeared_Pt_Eta_Phi[nPIDscenarios];
  TH3F* hTrack_Pos_GenSmeared_Pt_Eta_Phi[nPIDscenarios];
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_All_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_All_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Muon_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Muon_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Pion_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pion_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Kaon_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Kaon_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Proton_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Proton_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Ele_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Ele_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Pos_Gen_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pos_Gen_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_GenSmeared_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_GenSmeared_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Ele_GenSmeared_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Ele_GenSmeared_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Pos_GenSmeared_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pos_GenSmeared_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);


  TList* lGenPIDscenario[nPIDscenarios];
  for (int iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    lGenPIDscenario[iPIDscenario] = new TList();
    lGenPIDscenario[iPIDscenario]->SetName(Form("PIDscenario_%i", iPIDscenario+1));
    lGenPIDscenario[iPIDscenario]->SetOwner();
    lGenPIDscenario[iPIDscenario]->Add(hTrack_All_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_ElePos_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Muon_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Pion_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Kaon_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Proton_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Ele_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Pos_Gen_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_ElePos_GenSmeared_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Ele_GenSmeared_Pt_Eta_Phi[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hTrack_Pos_GenSmeared_Pt_Eta_Phi[iPIDscenario]);

    lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_gen[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_gen_primary[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_gen_heavy[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_gen_charm[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_ULS_gen_beauty[iPIDscenario]);
    lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_gen[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_gen_primary[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_gen_heavy[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_gen_charm[iPIDscenario]);
    // lGenPIDscenario[iPIDscenario]->Add(hMPtDCA_LS_gen_beauty[iPIDscenario]);
    // listGenTracks->Add(lGenPIDscenario);
  }

  std::vector<TList*> vecTListGenPIDscenarios;
  for (size_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    vecTListGenPIDscenarios.push_back(lGenPIDscenario[iPIDscenario]);
  }

  // auto hPt_Eta_Phi_primary_gen = new TH3F("hPt_Eta_Phi_primary_gen", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_primary_Ele_gen = new TH3F("hPt_Eta_Phi_primary_Ele_gen", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_primary_Pos_gen = new TH3F("hPt_Eta_Phi_primary_Pos_gen", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_hf_gen = new TH3F("hPt_Eta_Phi_hf_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_hf_Ele_gen = new TH3F("hPt_Eta_Phi_hf_Ele_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_hf_Pos_gen = new TH3F("hPt_Eta_Phi_hf_Pos_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_charm_gen = new TH3F("hPt_Eta_Phi_charm_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_charm_Ele_gen = new TH3F("hPt_Eta_Phi_charm_Ele_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_charm_Pos_gen = new TH3F("hPt_Eta_Phi_charm_Pos_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_beauty_gen = new TH3F("hPt_Eta_Phi_beauty_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_beauty_Ele_gen = new TH3F("hPt_Eta_Phi_beauty_Ele_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // auto hPt_Eta_Phi_beauty_Pos_gen = new TH3F("hPt_Eta_Phi_beauty_Pos_gen", ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);

  //
  // auto hPrimaryM = new TH1F("hPrimaryM",";MotherPDG",1000,0.5,1000.5);
  // auto hPrimaryGM = new TH1F("hPrimaryGM",";MotherPDG",1000,0.5,1000.5);
  // auto hCharmMm = new TH1F("hCharmMm",";MotherPDG",101,399.5,500.5);
  // auto hCharmMb = new TH1F("hCharmMb",";MotherPDG",1001,3099.5,5000.5);
  //
  // auto hCharmGM = new TH1F("hCharmGM",";MotherPDG",1000,0,10000);
  // auto hBeautyM = new TH1F("hBeautyM",";MotherPDG",1000,0,10000);
  // auto hBeautyGM = new TH1F("hBeautyGM",";MotherPDG",1000,0,10000);



  auto hM_Pt_sig_gen = new TH2F("hM_Pt_sig_gen",title2d.c_str(),300.,0,3.,200,0.,20.);

  TH1F* hTime0 = new TH1F("hTime0_sce%i",";t_{0} (ns)", 1000, -1., 1.);

  // TOF histograms
  auto hBetaP_beforeSmearing = new TH2F("hBetaP_beforeSmearing", ";#it{p} GeV/c;#beta", 500, 0., 10., 500, 0.1, 1.1);
  auto hBetaP_afterSmearing = new TH2F("hBetaP_afterSmearing", ";#it{p} GeV/c;#beta", 500, 0., 10., 500, 0.1, 1.1);
  TH2 *hNsigmaP_TOF[5];
  TH2 *hNsigmaP_TOF_trueElec[5];
  TH2 *hNsigmaP_TOF_trueMuon[5];
  TH2 *hNsigmaP_TOF_truePion[5];
  TH2 *hNsigmaP_TOF_trueKaon[5];
  TH2 *hNsigmaP_TOF_trueProton[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueElec[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueMuon[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_TOF_truePion[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueKaon[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueProton[5][nPIDscenarios];

  const char *pname[5] = {"el", "mu", "pi", "ka", "pr"};
  const char *plabel[5] = {"e", "#mu", "#pi", "K", "p"};
  for (int i = 0; i < 5; ++i){  // partilce type loop
    hNsigmaP_TOF[i] = new TH2F(Form("hNsigmaP_%s_TOF", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_TOF_trueElec[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueElec", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_TOF_trueMuon[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueMuon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_TOF_truePion[i] = new TH2F(Form("hNsigmaP_%s_TOF_truePion", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_TOF_trueKaon[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueKaon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_TOF_trueProton[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueProton", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    for (int j = 0; j < nPIDscenarios; ++j){ // PID scenatio loop
      hNsigmaP_afterPIDcuts_TOF[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_TOF_trueElec[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueElec_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_TOF_trueMuon[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueMuon_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_TOF_truePion[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_truePion_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_TOF_trueKaon[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueKaon_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_TOF_trueProton[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueProton_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    }
  }

  // RICH histograms
  TH2 *hNsigmaBeta_RICH[5];
  TH2 *hNsigmaP_RICH[5];
  TH2 *hNsigmaP_RICH_trueElec[5];
  TH2 *hNsigmaP_RICH_trueMuon[5];
  TH2 *hNsigmaP_RICH_truePion[5];
  TH2 *hNsigmaP_RICH_trueKaon[5];
  TH2 *hNsigmaP_RICH_trueProton[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueElec[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueMuon[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_RICH_truePion[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueKaon[5][nPIDscenarios];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueProton[5][nPIDscenarios];
  auto hCherenkovAngleP_beforeSmearing = new TH2F("hCherenkovAngleP_beforeSmearing", ";#it{p} GeV/c; cherenkov angle", 500, 0., 10., 200, 0., 0.4);
  auto hCherenkovAngleP_afterSmearing = new TH2F("hCherenkovAngleP_afterSmearing", ";#it{p} GeV/c; cherenkov angle", 500, 0., 10., 200, 0., 0.4);

  for (int i = 0; i < 5; ++i)  hNsigmaBeta_RICH[i] = new TH2F(Form("hNsigmaBeta_%s_RICH", pname[i]), Form(";#beta ;n#sigma_{%s}", plabel[i]), 200, 0., 1.1, 400, -100., 100.);

  for (int i = 0; i < 5; ++i) { // particle type loop
    hNsigmaP_RICH[i] = new TH2F(Form("hNsigmaP_%s_RICH", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 25.);
    hNsigmaP_RICH_trueElec[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueElec", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_RICH_trueMuon[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueMuon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_RICH_truePion[i] = new TH2F(Form("hNsigmaP_%s_RICH_truePion", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_RICH_trueKaon[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueKaon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    hNsigmaP_RICH_trueProton[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueProton", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    for (int j = 0; j < nPIDscenarios; ++j) { // PID scenatio loop
      hNsigmaP_afterPIDcuts_RICH[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_RICH_trueElec[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueElec_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_RICH_trueMuon[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueMuon_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_RICH_truePion[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_truePion_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_RICH_trueKaon[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueKaon_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
      hNsigmaP_afterPIDcuts_RICH_trueProton[i][j] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueProton_sce%i", pname[i], j+1), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
    }
  }


  for (int iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    for (int i = 0; i < 5; ++i) { // particle type loop
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_TOF[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_TOF_trueElec[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_TOF_trueMuon[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_TOF_truePion[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_TOF_trueKaon[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_TOF_trueProton[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_RICH[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_RICH_trueElec[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_RICH_trueMuon[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_RICH_truePion[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_RICH_trueKaon[i][iPIDscenario]);
      lRecPIDscenario[iPIDscenario]->Add(hNsigmaP_afterPIDcuts_RICH_trueProton[i][iPIDscenario]);
    }
  }

  // // Preshower histograms
  // // logx binning
  // const Int_t nbins = 80;
  // double xmin = 1.e-2;
  // double xmax = 1.e2;
  // double logxmin = std::log10(xmin);
  // double logxmax = std::log10(xmax);
  // double binwidth = (logxmax - logxmin) / nbins;
  // double xbins[nbins + 1];
  // for (Int_t i = 0; i <= nbins; ++i)
  //   xbins[i] = std::pow(10., logxmin + i * binwidth);
  //
  // // Preshower efficiency input histograms
  // TH1* hPreshowerEff[5];
  // auto filePreshowerEff = TFile::Open("preshowerEff.root");
  // hPreshowerEff[0] = (TH1*)filePreshowerEff->Get("preshowerEffElectrons");
  // hPreshowerEff[1] = (TH1*)filePreshowerEff->Get("preshowerEffPions");
  // hPreshowerEff[2] = (TH1*)filePreshowerEff->Get("preshowerEffPions");
  // hPreshowerEff[3] = (TH1*)filePreshowerEff->Get("preshowerEffPions");
  // hPreshowerEff[4] = (TH1*)filePreshowerEff->Get("preshowerEffPions");
  //
  // // Preshower QA output histograms
  // TH1F* hPreShElectronPt[nPIDscenarios];
  // TH1F* hPreShPionPt[nPIDscenarios];
  // for (int j = 0; j < nPIDscenarios; ++j) hPreShElectronPt[j] = new TH1F(Form("PreshowerElectronPt_sce%i",j+1), ";#it{p_{T}} (GeV/#it{c});", nbins, xbins);
  // for (int j = 0; j < nPIDscenarios; ++j) hPreShPionPt[j]     = new TH1F(Form("PreshowerPionPt_sce%i",j+1), ";#it{p_{T}} (GeV/#it{c});", nbins, xbins);
  //
  // for (int iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
  //   lRecPIDscenario[iPIDscenario]->Add(hPreShElectronPt[iPIDscenario]);
  //   lRecPIDscenario[iPIDscenario]->Add(hPreShPionPt[iPIDscenario]);
  // }
  //
  // // PID map
  // std::map<int, int> pidmap = { {11, 0}, {13, 1}, {211, 2}, {321, 3}, {2212, 4} };


  std::vector<TList*> vecTListRecPIDscenarios;
  for (size_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    vecTListRecPIDscenarios.push_back(lRecPIDscenario[iPIDscenario]);
  }

  // Setting Track and Particle vectors
  //-----------------------------------
  std::vector<Track *> vecElectron[nPIDscenarios],vecPositron[nPIDscenarios];
  std::vector<Track *> vecNegTracks[nPIDscenarios],vecPosTracks[nPIDscenarios];
  std::vector<GenParticle *> vecElectronGen[nPIDscenarios],vecPositronGen[nPIDscenarios], vecElectronGenSmeared[nPIDscenarios],vecPositronGenSmeared[nPIDscenarios];
  std::vector<GenParticle *> vecGen, vecNegGen, vecPosGen;
  std::vector<Track *> vecPIDtracks; // , vecTOFtracks

  int nTotalTracks = 0;
  int nTotalParticles = 0;


  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);
    nTracksGen->Fill(tracks->GetEntries());
    nParticles->Fill(particles->GetEntries());
    // nTotalParticles = nTotalParticles + particles->GetEntries();
    // nTotalTracks = nTotalTracks + tracks->GetEntries();
    // cout << "This event has " << particles->GetEntries() << " particles, and " << tracks->GetEntries() << " tracks, the total Number of Particles is now " << nTotalParticles << " and total number of tracks Tracks " << nTotalTracks << endl << endl;

    // std::cout << " NTracks: " << tracks->GetEntries() << std::endl;
    // std::cout << " NParticles: " << particles->GetEntries() << std::endl;



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

      // Check isprimary
      Double_t isphysicalprimary = kTRUE;
      auto pdg = 0;
      pdg = particle->PID;
      if(!IsStable(pdg)) isphysicalprimary = kFALSE;
      // mother
      auto imother  = particle->M1;
      auto mother   = imother != -1 ? (GenParticle *)particles->At(imother) : (GenParticle *)nullptr;
      auto mPid     = 0;
      if(mother) mPid = mother->PID;
      if(IsStable(pdg) && IsStable(mPid)) {
        //printf("pdg and mpdg are %d and %d, not counted\n",pdg,mPid);
        isphysicalprimary = kFALSE;
      }

      if(isphysicalprimary) {
        // increase numer of particle if its in acceptance of  FIT detecotr
        double eta = particle->Eta;
        if( ((2.2<eta&&eta< 5.0) || (-3.4<eta&&eta<-2.3)) && particle->Charge != 0) numParticlesFIT++;
        if( (fabs(eta) < 0.5) && particle->Charge != 0) numParticlesMidRapidity++;
      }
    }






    // centrality selection
    //---------------------
    // accept event if event has more then X particles within the FIT detector accaptance (2.2<eta&&eta< 5.0) || (-3.4<eta&&eta<-2.3)
    nParticlesFIT->Fill(numParticlesFIT);
    nParticlesMidRapidity->Fill(numParticlesMidRapidity);
    if( numParticlesFIT < 5179 ) {
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
      // auto pPID = particle->PID;
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

      double phiGen = TVector2::Phi_0_2pi(particle->Phi);


      // Smearing generated particles
      //-----------------------------
      LV.SetPtEtaPhiM(particle->PT,particle->Eta,phiGen,eMass);
      Double_t charge = -particle->PID;
      LV_smeared = ApplySmearing(fArrResoPt,fArrResoEta,fArrResoPhi_Pos,fArrResoPhi_Neg,LV,charge);

      //Filling suport histograms to see effect of smearing
      if(charge > 0.)  hSmearing_For_Eff_phi_pos->Fill(particle->PT, LV.Phi() - LV_smeared.Phi());
      if(charge < 0.)  hSmearing_For_Eff_phi_neg->Fill(particle->PT, LV.Phi() - LV_smeared.Phi());
                       hSmearing_For_Eff_eta->Fill(particle->PT, LV.Eta() - LV_smeared.Eta());
      if(LV.Pt() > 0.) hSmearing_For_Eff_pt->Fill(particle->PT, (LV.Pt() - LV_smeared.Pt())/LV.Pt());

      // Filling TH3 before kinematic cuts are applied
      if (abs(particle->PID) == 11 )  hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons
      if (abs(particle->PID) == 11 )  hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons


      for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
        // kinematic cuts on particles
        if (kineCuts(particle,iPID_scenario)){
          if (particle->PID == 11 ) vecElectronGen[iPID_scenario].push_back(particle);       // vector filled with generated electrons
          else if (particle->PID == -11 ) vecPositronGen[iPID_scenario].push_back(particle); // vector filled with generated positrons

          // Filling TH3 histograms for generated  Tracks
          //---------------------------------------------
                                            hTrack_All_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of all generated tracks
          if (abs(particle->PID) == 11 )    hTrack_ElePos_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 13 )    hTrack_Muon_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 211 )   hTrack_Pion_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 321 )   hTrack_Kaon_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons
          if (abs(particle->PID) == 2212 )  hTrack_Proton_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons
          if (particle->PID == 11)        hTrack_Ele_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen);    // Pt Eta Phi of generated electrons
          else if (particle->PID == -11)  hTrack_Pos_Gen_Pt_Eta_Phi[iPID_scenario]->Fill(particle->PT,particle->Eta,phiGen);    // Pt Eta Phi of generated positrons
        }

        // applying kinematic cuts on smeared LV vector
        if (!kineCuts(LV_smeared,iPID_scenario)) continue;
        double phiGenSmear = TVector2::Phi_0_2pi(LV_smeared.Phi());
                          // cout << "gen pt eta phi: " << particle->PT << " " << particle->Eta << " " << phiGen << endl;
                          // cout << "gen smeared pt eta phi: " << LV_smeared.Pt() << " " << LV_smeared.Eta() << " " << phiGenSmear << endl;


        // Filling TH3 histograms for generated smeared Tracks
        //----------------------------------------------------
        // look only at electrons
        if (abs(particle->PID) == 11 )  hTrack_ElePos_GenSmeared_Pt_Eta_Phi[iPID_scenario]->Fill(LV_smeared.Pt(),LV_smeared.Eta(),phiGenSmear); // Pt Eta Phi of generated smeared electrons + positrons
        if (particle->PID == 11)        hTrack_Ele_GenSmeared_Pt_Eta_Phi[iPID_scenario]->Fill(LV_smeared.Pt(),LV_smeared.Eta(),phiGenSmear);    // Pt Eta Phi of generated smeared electrons
        else if (particle->PID == -11)  hTrack_Pos_GenSmeared_Pt_Eta_Phi[iPID_scenario]->Fill(LV_smeared.Pt(),LV_smeared.Eta(),phiGenSmear);    // Pt Eta Phi of generated smeared positrons

        if (particle->PID == 11 )       vecElectronGenSmeared[iPID_scenario].push_back(particle);        // vector filled with generated smeared electrons
        else if (particle->PID == -11 ) vecPositronGenSmeared[iPID_scenario].push_back(particle);        // vector filled with generated smeared positrons

        // fill histograms
        // separate electrons & positrons originating from different light flavour (priamary), charm (cc) and buty (bb) decays
        // if (!hasStrangeAncestor(particle, particles) && !hasHeavyAncestor(particle, particles)) {
        //   hPt_Eta_Phi_primary_gen->Fill(particle->PT, particle->Eta, phiGen);
        //   if(particle->PID == 11) hPt_Eta_Phi_primary_Ele_gen->Fill(particle->PT, particle->Eta, phiGen);
        //   else if(particle->PID == -11) hPt_Eta_Phi_primary_Pos_gen->Fill(particle->PT, particle->Eta, phiGen);
        // }
        // else {
        //     if (hasCharmAncestor(particle, particles) || hasBeautyAncestor(particle, particles)){
        //       hPt_Eta_Phi_hf_gen->Fill(particle->PT, particle->Eta, phiGen);
        //       if(particle->PID == 11) hPt_Eta_Phi_hf_Ele_gen->Fill(particle->PT, particle->Eta, phiGen);
        //       else if(particle->PID == -11) hPt_Eta_Phi_hf_Pos_gen->Fill(particle->PT, particle->Eta, phiGen);
        //     }
        //     if(isCharm(mPid) && !hasBeautyAncestor(particle, particles)){
        //       hPt_Eta_Phi_charm_gen->Fill(particle->PT, particle->Eta, phiGen);
        //       if(particle->PID == 11) hPt_Eta_Phi_charm_Ele_gen->Fill(particle->PT, particle->Eta, phiGen);
        //       else if(particle->PID == -11) hPt_Eta_Phi_charm_Pos_gen->Fill(particle->PT, particle->Eta, phiGen);
        //     }
        //     else if (isBeauty(mPid)) {
        //       hPt_Eta_Phi_beauty_gen->Fill(particle->PT, particle->Eta, phiGen);
        //       if(particle->PID == 11) hPt_Eta_Phi_beauty_Ele_gen->Fill(particle->PT, particle->Eta, phiGen);
        //       else if(particle->PID == -11) hPt_Eta_Phi_beauty_Pos_gen->Fill(particle->PT, particle->Eta, phiGen);
        //     }
        //     else if (hasStrangeAncestor(particle, particles)) continue;
        //     // else {
        //     //   std::cout << particle->X << " " << particle->Y << " " << particle->Z << std::endl;
        //     //   std::cout << " --- mother PDG     : " << mother->PID << std::endl;
        //     //   std::cout << " --- grandmother PDG: " << gmother->PID << std::endl;
        //     // }
        // }
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
      // get track and corresponding particle
      auto track = (Track *)tracks->At(itrack);
      auto particle = (GenParticle *)track->Particle.GetObject();

      double phiRec = TVector2::Phi_0_2pi(track->Phi);

      if (bUseTOF){
        // fill beta-p before semaring
        auto p = track->P;
        auto beta = toflayer.getBeta(*track);
        hBetaP_beforeSmearing->Fill(p, beta);
      }

      if (bUseRICH){
        // fill cherenkovAngle-p before semaring
        auto p = track->P;
        auto measAngle = (richdetector.getMeasuredAngle(*track)).first;
        hCherenkovAngleP_beforeSmearing->Fill(p,measAngle);
      }

      // smear track if requested
      if(abs(particle->PID) == 11 ) hBeforeSmearing_Pt_Eta_Phi_rec->Fill(track->PT,track->Eta,phiRec);
      if (bSmear) if (!smearer.smearTrack(*track)) continue; // strange syntax, but works


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

      if( (particle->Eta > -0.5 && particle->Eta<0.5) && particle->Charge != 0 ) hdNdeta_midrap_rec->Fill(particle->Eta);

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


      if(abs(particle->PID) == 11 ) hAfterSmearing_Pt_Eta_Phi_rec->Fill(track->PT,track->Eta,phiRec);


      auto p = track->P;
      auto beta = toflayer.getBeta(*track);
      auto measAngle = (richdetector.getMeasuredAngle(*track)).first;

      // fill beta-p after semaring
      hBetaP_afterSmearing->Fill(p, beta);
      // fill cherenkovAngle-p after semaring
      hCherenkovAngleP_afterSmearing->Fill(p,measAngle);


      //make TOF PID
      std::array<float, 5> PIDdeltatTOF, PIDnsigmaTOF;
      toflayer.makePID(*track, PIDdeltatTOF, PIDnsigmaTOF);
      //make RICH PID
      std::array<float, 5> PIDdeltaangleRICH, PIDnsigmaRICH;
      richdetector.makePID(*track, PIDdeltaangleRICH, PIDnsigmaRICH);

      // if(iPID_scenario == 1){ //Note: 0=Sce1, 1=Sce2 etc... only fill those histograms once, here no PID selection was done yet
        // fill nsigma TOF   (before PID selection is applied)
        hTrackPt->Fill(track->PT);
        hTrackP->Fill(track->P);
        if (toflayer.hasTOF(*track) || richdetector.hasRICH(*track)){
          hTrackPt_hasTOForRICH->Fill(track->PT);
          hTrackP_hasTOForRICH->Fill(track->P);
        }


        if( etaCut(track) && (toflayer.hasTOF(*track) || richdetector.hasRICH(*track))){
          for (int i = 0; i < 5; ++i) {
            if(toflayer.hasTOF(*track)){
              hNsigmaP_TOF[i]->Fill(p, PIDnsigmaTOF[i]);
              if      (abs(particle->PID) == 11)    hNsigmaP_TOF_trueElec[i]->Fill(p, PIDnsigmaTOF[i]);
              else if (abs(particle->PID) == 13)    hNsigmaP_TOF_trueMuon[i]->Fill(p, PIDnsigmaTOF[i]);
              else if (abs(particle->PID) == 211)   hNsigmaP_TOF_truePion[i]->Fill(p, PIDnsigmaTOF[i]);
              else if (abs(particle->PID) == 321)   hNsigmaP_TOF_trueKaon[i]->Fill(p, PIDnsigmaTOF[i]);
              else if (abs(particle->PID) == 2212)  hNsigmaP_TOF_trueProton[i]->Fill(p, PIDnsigmaTOF[i]);
            }

            // fill nsigma RICH  (before PID selection is applied)
            if(richdetector.hasRICH(*track)){
              hNsigmaP_RICH[i]->Fill(p, PIDnsigmaRICH[i]);
              hNsigmaBeta_RICH[i]->Fill(beta, PIDnsigmaRICH[i]);
              if      (abs(particle->PID) == 11)    {hNsigmaP_RICH_trueElec[i]->Fill(p, PIDnsigmaRICH[i]);   }
              else if (abs(particle->PID) == 13)    {hNsigmaP_RICH_trueMuon[i]->Fill(p, PIDnsigmaRICH[i]);   }
              else if (abs(particle->PID) == 211)   {hNsigmaP_RICH_truePion[i]->Fill(p, PIDnsigmaRICH[i]);   }
              else if (abs(particle->PID) == 321)   {hNsigmaP_RICH_trueKaon[i]->Fill(p, PIDnsigmaRICH[i]);   }
              else if (abs(particle->PID) == 2212)  {hNsigmaP_RICH_trueProton[i]->Fill(p, PIDnsigmaRICH[i]); }
            }
          }
          // bool a = fabs(hypot(track->XOuter * 0.1,track->YOuter * 0.1) - 100) < 0.001;
          // bool b = fabs(track->ZOuter * 0.1) < 200;
          // if(p < 0.1) {
          //   cout << " PDG = " << track->PID << ", p = " << p << endl;
          //   cout << " hasTOF = " << toflayer.hasTOF(*track) << ", hasRICH = " << richdetector.hasRICH(*track) << endl;
          //   cout << " TOF NsigmaEle = " << PIDnsigmaTOF[0] << ", NsigmaPi = " << PIDnsigmaTOF[2] << endl;
          //   cout << " RICH NsigmaEle = " << PIDnsigmaRICH[0] << ", NsigmaPi = " << PIDnsigmaRICH[2] << endl;
          //   cout << " TOF hasTOF requirements: (fabs(r - mRadius) < 0.001 ) = " << a << ", (fabs(z) < mLength) = " << b << endl;
          //   cout << endl;
          // }
        }
        // }



      for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {

        // kinatic cuts on tracks
        if (!kineCuts(track, iPID_scenario)) continue;

        if(abs(particle->PID) == 11 ) hAfterKineCuts_Pt_Eta_Phi_rec[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);

        // look only at electrons
        // if(abs(particle->PID) != 11 ) continue;



        //setting variables for i-th PID scenario
        bool   i_useTOFPID = useTOFPID[iPID_scenario];
        bool   i_useRICHPID;
        if (Bz == 0.2)      i_useRICHPID = useRICHPID_B2[iPID_scenario];
        else if (Bz == 0.5) i_useRICHPID = useRICHPID_B5[iPID_scenario];
        else cout << "!!! Error, issue with setting BField. !!!" << endl;
        // bool   i_usePreShPID = usePreShPID[iPID_scenario];
        double i_SigmaTOFEle = nSigmaTOFEle[iPID_scenario];
        double i_SigmaTOFPi = nSigmaTOFPi[iPID_scenario];
        double i_SigmaRICHEle = nSigmaRICHEle[iPID_scenario];
        double i_SigmaRICHPi = nSigmaRICHPi[iPID_scenario];

        // ################## PID selection ##################
        // PID selection, identifying particles by using TOF "OR" RICH.
        // using TOF  electron inclusion and pion rejection
        // using RICH electron inclusion and pion rejection
        bool TOFpid = false;
        bool RICHpid = false;
        // bool PreShpid = false;

        // Preshower PID
        // check if track is identified as an electron by the Preshower
        // auto pdg  = std::abs(track->PID);
        // auto pbin = hPreshowerEff[pidmap[pdg]]->FindBin(p);
        // auto eff  = hPreshowerEff[pidmap[pdg]]->GetBinContent(pbin);
        // if (i_usePreShPID && (gRandom->Uniform() < eff) ){
        //   PreShpid = true;
        //   // and fill some QA histograms for Preshower PID
        //   if(!pidmap[pdg]) hPreShElectronPt[iPID_scenario]->Fill(track->PT);
        //   else hPreShPionPt[iPID_scenario]->Fill(track->PT);
        // }  // select primaries based on 3 sigma DCA cuts
        //


        //TOF PID
        if(i_useTOFPID && (toflayer.hasTOF(*track)) && (p < tof_EleAccep_p_cut)) {
          if (p > tof_PionRej_p_cut && richdetector.hasRICH(*track)) {
            if((fabs(PIDnsigmaTOF[0]) < i_SigmaTOFEle) && (fabs(PIDnsigmaRICH[0]) < nSigmaRICHEle_forTOFPID)) TOFpid = true; // is within 3 sigma of the electron band (TOF)
          }
          else if(p <= tof_PionRej_p_cut){
            if(fabs(PIDnsigmaTOF[0]) < i_SigmaTOFEle) TOFpid = true; // is within 3 sigma of the electron band (TOF)
          }
          else TOFpid = false; // This is rejecting all the heavier partilces which do not have a RICH signal in the pt area of 0.4-0.6 GeV/c


          if(fabs(PIDnsigmaTOF[2]) < i_SigmaTOFPi) TOFpid = false; // is within 3 sigma of the pion band (TOF)
        }
        else TOFpid = false;

        //RICH PID
        if(i_useRICHPID && richdetector.hasRICH(*track)) {
          if(fabs(PIDnsigmaRICH[0]) < i_SigmaRICHEle) RICHpid = true; // is within 3 sigma of the electron band (RICH)
          if((fabs(PIDnsigmaRICH[2]) < i_SigmaRICHPi) && (p > rich_PionRejection_p_cut) ) RICHpid = false; // is within 3 sigma of the pion band (RICH)
        }
        else RICHpid = false;

        // bool a = fabs(hypot(track->XOuter * 0.1,track->YOuter * 0.1) - 100) < 0.001;
        // bool b = fabs(track->ZOuter * 0.1) < 200;
        // if(p < 0.1) {
        //   cout << " TOFpid = " << TOFpid << ", RICHpid = " << RICHpid << endl;
        //   cout << " PDG = " << track->PID << ", p = " << p << endl;
        //   cout << " hasTOF = " << toflayer.hasTOF(*track) << ", hasRICH = " << richdetector.hasRICH(*track) << endl;
        //   cout << " TOF NsigmaEle = " << PIDnsigmaTOF[0] << ", NsigmaPi = " << PIDnsigmaTOF[2] << endl;
        //   cout << " RICH NsigmaEle = " << PIDnsigmaRICH[0] << ", NsigmaPi = " << PIDnsigmaRICH[2] << endl;
        //   cout << " TOF hasTOF requirements: (fabs(r - mRadius) < 0.001 ) = " << a << ", (fabs(z) < mLength) = " << b << endl;
        //   cout << endl;
        // }


        if (!(RICHpid || TOFpid) /*&& !PreShpid*/) continue; // check if TOF or RICH signal is true.
        // ################## end of PID selection ##################

        // fill nsigma after PID cuts histograms    (after PID selection has been applied)
        for (int i = 0; i < 5; ++i) {
          // NSigma plots, separating for true particle ID
          if(toflayer.hasTOF(*track)){
            hNsigmaP_afterPIDcuts_TOF[i][iPID_scenario]->Fill(p, PIDnsigmaTOF[i]);
            if      (abs(particle->PID) == 11)   hNsigmaP_afterPIDcuts_TOF_trueElec[i][iPID_scenario]->Fill(p, PIDnsigmaTOF[i]);
            else if (abs(particle->PID) == 13)   hNsigmaP_afterPIDcuts_TOF_trueMuon[i][iPID_scenario]->Fill(p, PIDnsigmaTOF[i]);
            else if (abs(particle->PID) == 211)  hNsigmaP_afterPIDcuts_TOF_truePion[i][iPID_scenario]->Fill(p, PIDnsigmaTOF[i]);
            else if (abs(particle->PID) == 321)  hNsigmaP_afterPIDcuts_TOF_trueKaon[i][iPID_scenario]->Fill(p, PIDnsigmaTOF[i]);
            else if (abs(particle->PID) == 2212) hNsigmaP_afterPIDcuts_TOF_trueProton[i][iPID_scenario]->Fill(p, PIDnsigmaTOF[i]);
          }
          if(richdetector.hasRICH(*track)){
            hNsigmaP_afterPIDcuts_RICH[i][iPID_scenario]->Fill(p, PIDnsigmaRICH[i]);
            if      (abs(particle->PID) == 11)   hNsigmaP_afterPIDcuts_RICH_trueElec[i][iPID_scenario]->Fill(p, PIDnsigmaRICH[i]);
            else if (abs(particle->PID) == 13)   hNsigmaP_afterPIDcuts_RICH_trueMuon[i][iPID_scenario]->Fill(p, PIDnsigmaRICH[i]);
            else if (abs(particle->PID) == 211)  hNsigmaP_afterPIDcuts_RICH_truePion[i][iPID_scenario]->Fill(p, PIDnsigmaRICH[i]);
            else if (abs(particle->PID) == 321)  hNsigmaP_afterPIDcuts_RICH_trueKaon[i][iPID_scenario]->Fill(p, PIDnsigmaRICH[i]);
            else if (abs(particle->PID) == 2212) hNsigmaP_afterPIDcuts_RICH_trueProton[i][iPID_scenario]->Fill(p, PIDnsigmaRICH[i]);
          }
        }
        // if(!doPID(track, i_useTOFPID, i_useRICHPID, bUsePreSh, tof_EleAccep_p_cut, tof_PionRej_p_cut, rich_PionRejection_p_cut, i_SigmaTOFEle, i_SigmaTOFPi, i_SigmaRICHEle, i_SigmaRICHPi, toflayer, richdetector, PIDnsigmaTOF, PIDnsigmaRICH)) continue;

        // fill histograms
                                    hAllTracks_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec); // Pt Eta Phi of all type of reconstructed particles
        if     (track->Charge < 0)  hNegTrack_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);  // Pt Eta Phi of reconstructed negative tracks
        else if(track->Charge > 0)  hPosTrack_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);  // Pt Eta Phi of reconstructed positive tracks

        // Pt Eta Phi plots for different type of particles
        if      (abs(particle->PID) == 13)    hTrack_Muon_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);
        else if (abs(particle->PID) == 211)   hTrack_Pion_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);
        else if (abs(particle->PID) == 321)   hTrack_Kaon_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);
        else if (abs(particle->PID) == 2212)  hTrack_Proton_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);
        else if (abs(particle->PID) == 11)  sleep(0);
        else std::cout << "Particle not identified!    Particle PID: " << abs(particle->PID) << std::endl;


        // fill track vector
        if (track->Charge < 0) vecNegTracks[iPID_scenario].push_back(track);       //  vector of all reconstructed negative tracks
        else if (track->Charge > 0) vecPosTracks[iPID_scenario].push_back(track);  //  vector of all reconstructed positve tracks

        if (particle->PID == 11 ) vecElectron[iPID_scenario].push_back(track);       // vector of reconstructed electrons
        else if (particle->PID == -11 ) vecPositron[iPID_scenario].push_back(track); // vector of reconstructed positrons

        // look only at electrons
        if(abs(particle->PID) != 11 ) continue;   // considering only electrons & positrons

                                        hTrack_ElePos_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);  // Pt Eta Phi of reconstructed electrons + positrons
        if (particle->PID == 11)        hTrack_Ele_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);     // Pt Eta Phi of reconstructed electrons
        else if (particle->PID == -11)  hTrack_Pos_Rec_Pt_Eta_Phi[iPID_scenario]->Fill(track->PT,track->Eta,phiRec);     // Pt Eta Phi of reconstructed positrons


        // // separate electrons & positrons originating from different light flavour (priamary), charm (cc) and buty (bb) decays
        // if (!hasStrangeAncestor(particle, particles) && !hasHeavyAncestor(particle, particles)) {
        //   hPrimaryM ->Fill(mother->PID);
        //   hPrimaryGM->Fill(gmother->PID);
        //                                  hPt_Eta_Phi_primary_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //   if      (particle->PID == 11)  hPt_Eta_Phi_primary_Ele_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //   else if (particle->PID == -11) hPt_Eta_Phi_primary_Pos_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        // }
        // else {
        //     if (hasCharmAncestor(particle, particles) || hasBeautyAncestor(particle, particles)){
        //                                      hPt_Eta_Phi_hf_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //       if      (particle->PID == 11)  hPt_Eta_Phi_hf_Ele_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //       else if (particle->PID == -11) hPt_Eta_Phi_hf_Pos_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //     }
        //     if(isCharm(mPid) && !hasBeautyAncestor(particle, particles)){
        //       hCharmGM->Fill(gmother->PID);
        //                                      hPt_Eta_Phi_charm_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //       if      (particle->PID == 11)  hPt_Eta_Phi_charm_Ele_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //       else if (particle->PID == -11) hPt_Eta_Phi_charm_Pos_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //     }
        //     else if (isBeauty(mPid)) {
        //     // else if (hasBeautyAncestor(particle, particles) && mPid > 500 && mPid < 600) {
        //       hBeautyM->Fill(mother->PID);
        //       hBeautyGM->Fill(gmother->PID);
        //                                     hPt_Eta_Phi_beauty_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //       if     (particle->PID == 11)  hPt_Eta_Phi_beauty_Ele_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //       else if(particle->PID == -11) hPt_Eta_Phi_beauty_Pos_rec[iPID_scenario]->Fill(track->PT, track->Eta, phiRec);
        //     }
        //     else if (hasStrangeAncestor(particle, particles)) continue;
        //     // else {
        //     //   std::cout << particle->X << " " << particle->Y << " " << particle->Z << std::endl;
        //     //   std::cout << " --- mother PDG     : " << mother->PID << std::endl;
        //     //   std::cout << " --- grandmother PDG: " << gmother->PID << std::endl;
        //     // }
        // }

      }



    }
    vecPIDtracks.clear();

    // cout << endl << " next event " << endl << endl;


    //##################################################
    //############   ULS and LS pairing   ##############
    //############ for Gen and Rec tracks ##############
    //##################################################

    bool pairULS = kTRUE;
    bool MCpidEle = kTRUE;

    for (Int_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
      dileptonPairingRec(vecNegTracks[iPIDscenario], vecPosTracks[iPIDscenario],  pairULS, !MCpidEle, particles, iPIDscenario);
      dileptonPairingRec(vecNegTracks[iPIDscenario], vecNegTracks[iPIDscenario], !pairULS, !MCpidEle, particles, iPIDscenario);
      dileptonPairingRec(vecPosTracks[iPIDscenario], vecPosTracks[iPIDscenario], !pairULS, !MCpidEle, particles, iPIDscenario);
      dileptonPairingRec(vecElectron[iPIDscenario],  vecPositron[iPIDscenario],   pairULS,  MCpidEle, particles, iPIDscenario);
      dileptonPairingRec(vecElectron[iPIDscenario],  vecElectron[iPIDscenario],  !pairULS,  MCpidEle, particles, iPIDscenario);
      dileptonPairingRec(vecPositron[iPIDscenario],  vecPositron[iPIDscenario],  !pairULS,  MCpidEle, particles, iPIDscenario);

      vecNegTracks[iPIDscenario].clear();
      vecPosTracks[iPIDscenario].clear();
      vecElectron[iPIDscenario].clear();
      vecPositron[iPIDscenario].clear();


      dileptonPairingGen(vecElectronGen[iPIDscenario], vecPositronGen[iPIDscenario], pairULS, particles, iPIDscenario);
      dileptonPairingGen(vecElectronGen[iPIDscenario], vecElectronGen[iPIDscenario], !pairULS, particles, iPIDscenario);
      dileptonPairingGen(vecPositronGen[iPIDscenario], vecPositronGen[iPIDscenario], !pairULS, particles, iPIDscenario);

      vecElectronGen[iPIDscenario].clear();
      vecPositronGen[iPIDscenario].clear();
      vecElectronGenSmeared[iPIDscenario].clear();
      vecPositronGenSmeared[iPIDscenario].clear();
    }
    vecNegGen.clear();
    vecPosGen.clear();

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

  for (size_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    fout->mkdir(Form("generated/%s",vecTListGenPIDscenarios.at(iPIDscenario)->GetName()));
    fout->cd(Form("generated/%s",vecTListGenPIDscenarios.at(iPIDscenario)->GetName()));
    vecTListGenPIDscenarios.at(iPIDscenario)->Write();
    vecTListGenPIDscenarios.at(iPIDscenario)->Clear();
  }
  vecTListGenPIDscenarios.clear();

  fout->cd();
  fout->mkdir("reconstructed");
  fout->cd("reconstructed");

  for (size_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    fout->mkdir(Form("reconstructed/%s",vecTListRecPIDscenarios.at(iPIDscenario)->GetName()));
    fout->cd(Form("reconstructed/%s",vecTListRecPIDscenarios.at(iPIDscenario)->GetName()));
    vecTListRecPIDscenarios.at(iPIDscenario)->Write();
    vecTListRecPIDscenarios.at(iPIDscenario)->Clear();
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
  // hAllTracks_Rec_Pt_Eta_Phi->Write();
  // hNegTrack_Rec_Pt_Eta_Phi->Write();
  // hPosTrack_Rec_Pt_Eta_Phi->Write();
  // hTrack_ElePos_Rec_Pt_Eta_Phi->Write();
  // hTrack_Ele_Rec_Pt_Eta_Phi->Write();
  // hTrack_Pos_Rec_Pt_Eta_Phi->Write();
  // hTrack_Muon_Rec_Pt_Eta_Phi->Write();
  // hTrack_Pion_Rec_Pt_Eta_Phi->Write();
  // hTrack_Kaon_Rec_Pt_Eta_Phi->Write();
  // hTrack_Proton_Rec_Pt_Eta_Phi->Write();
  hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts->Write();
  hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts->Write();
  // hTrack_All_Gen_Pt_Eta_Phi->Write();
  // hTrack_ElePos_Gen_Pt_Eta_Phi->Write();
  // hTrack_Muon_Gen_Pt_Eta_Phi->Write();
  // hTrack_Pion_Gen_Pt_Eta_Phi->Write();
  // hTrack_Kaon_Gen_Pt_Eta_Phi->Write();
  // hTrack_Proton_Gen_Pt_Eta_Phi->Write();
  // hTrack_Ele_Gen_Pt_Eta_Phi->Write();
  // hTrack_Pos_Gen_Pt_Eta_Phi->Write();
  // hTrack_ElePos_GenSmeared_Pt_Eta_Phi->Write();
  // hTrack_Ele_GenSmeared_Pt_Eta_Phi->Write();
  // hTrack_Pos_GenSmeared_Pt_Eta_Phi->Write();
  // hPrimaryM->Write();
  // hPrimaryGM->Write();
  // hCharmMm->Write();
  // hCharmMb->Write();
  // hCharmGM->Write();
  // hBeautyM->Write();
  // hBeautyGM->Write();
  hBeforeSmearing_Pt_Eta_Phi_rec->Write();
  hAfterSmearing_Pt_Eta_Phi_rec->Write();
  // hAfterKineCuts_Pt_Eta_Phi_rec->Write();
  // hPt_Eta_Phi_primary_rec->Write();
  // hPt_Eta_Phi_primary_Ele_rec->Write();
  // hPt_Eta_Phi_primary_Pos_rec->Write();
  // hPt_Eta_Phi_hf_rec->Write();
  // hPt_Eta_Phi_hf_Ele_rec->Write();
  // hPt_Eta_Phi_hf_Pos_rec->Write();
  // hPt_Eta_Phi_charm_rec->Write();
  // hPt_Eta_Phi_charm_Ele_rec->Write();
  // hPt_Eta_Phi_charm_Pos_rec->Write();
  // hPt_Eta_Phi_beauty_rec->Write();
  // hPt_Eta_Phi_beauty_Ele_rec->Write();
  // hPt_Eta_Phi_beauty_Pos_rec->Write();
  // hPt_Eta_Phi_primary_gen->Write();
  // hPt_Eta_Phi_primary_Ele_gen->Write();
  // hPt_Eta_Phi_primary_Pos_gen->Write();
  // hPt_Eta_Phi_hf_gen->Write();
  // hPt_Eta_Phi_hf_Ele_gen->Write();
  // hPt_Eta_Phi_hf_Pos_gen->Write();
  // hPt_Eta_Phi_charm_gen->Write();
  // hPt_Eta_Phi_charm_Ele_gen->Write();
  // hPt_Eta_Phi_charm_Pos_gen->Write();
  // hPt_Eta_Phi_beauty_gen->Write();
  // hPt_Eta_Phi_beauty_Ele_gen->Write();
  // hPt_Eta_Phi_beauty_Pos_gen->Write();

  // hMPtDCA_ULS_gen->Write();
  // hMPtDCA_ULS_gen_primary->Write();
  // hMPtDCA_ULS_gen_heavy->Write();
  // hMPtDCA_ULS_gen_charm->Write();
  // hMPtDCA_ULS_gen_beauty->Write();
  // hMPtDCA_LS_gen->Write();
  // hMPtDCA_LS_gen_primary->Write();
  // hMPtDCA_LS_gen_heavy->Write();
  // hMPtDCA_LS_gen_charm->Write();
  // hMPtDCA_LS_gen_beauty->Write();
  // hMPtDCA_ULS_rec->Write();
  // // hMPtDCA_ULS_rec_primary->Write();
  // // hMPtDCA_ULS_rec_heavy->Write();
  // // hMPtDCA_ULS_rec_charm->Write();
  // // hMPtDCA_ULS_rec_beauty->Write();
  // hMPtDCA_ULS_rec_MCpidEle->Write();
  // // hMPtDCA_ULS_rec_MCpidEle_primary->Write();
  // // hMPtDCA_ULS_rec_MCpidEle_heavy->Write();
  // // hMPtDCA_ULS_rec_MCpidEle_charm->Write();
  // // hMPtDCA_ULS_rec_MCpidEle_beauty->Write();
  // hMPtDCA_LS_rec->Write();
  // // hMPtDCA_LS_rec_primary->Write();
  // // hMPtDCA_LS_rec_heavy->Write();
  // // hMPtDCA_LS_rec_charm->Write();
  // // hMPtDCA_LS_rec_beauty->Write();
  // hMPtDCA_LS_rec_MCpidEle->Write();
  // // hMPtDCA_LS_rec_MCpidEle_primary->Write();
  // // hMPtDCA_LS_rec_MCpidEle_heavy->Write();
  // // hMPtDCA_LS_rec_MCpidEle_charm->Write();
  // // hMPtDCA_LS_rec_MCpidEle_beauty->Write();

  // fout->cd("generated/");
  // fout->mkdir("generated/LS");
  // fout->cd("generated/LS");
  // fout->cd();
  // fout->mkdir("reconstructed/");
  // fout->cd("reconstructed/");
  // fout->mkdir("reconstructed/ULS");
  // fout->cd("reconstructed/ULS");


  hM_Pt_sig_gen->Write();

  hBetaP_beforeSmearing->Write();
  hBetaP_afterSmearing->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueElec[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueMuon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_truePion[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueKaon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueProton[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueElec[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_truePion[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueKaon[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueProton[i]->Write();


  hCherenkovAngleP_beforeSmearing->Write();
  hCherenkovAngleP_afterSmearing->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueElec[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueMuon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_truePion[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueKaon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueProton[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueElec[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_truePion[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueKaon[i]->Write();
  // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueProton[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaBeta_RICH[i]->Write();



  // fout->cd("reconstructed/");

  // fout->mkdir("reconstructed/LS");
  // fout->cd("reconstructed/LS");
  fout->Close();



  // stop watch
  watch->Stop();
  watch->Print();

}





// function to ULS or LS pair reconstructed particles
void dileptonPairingRec(std::vector<Track *> vec_track_neg,std::vector<Track *> vec_track_pos, bool pairULS, bool MCpidEle, TClonesArray *particles, Int_t iSce){
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
    // auto imother1 = particle1->M1;
    // auto mother1 = imother1 != -1 ? (GenParticle *)particles->At(imother1) : (GenParticle *)nullptr;
    // auto m1Pid = mother1->PID;
    LV1.SetPtEtaPhiM(track1->PT,track1->Eta,track1->Phi,eMass);

    if(!pairULS) iPos_start = iEle+1;
    for (auto iPos = iPos_start; iPos < vec_track_pos_size; iPos++){
      auto track2 = (Track *) vec_track_pos.at(iPos);
      auto particle2 = (GenParticle *)track2->Particle.GetObject();
      // auto imother2 = particle2->M1;
      // auto mother2 = imother2 != -1 ? (GenParticle *)particles->At(imother2) : (GenParticle *)nullptr;
      // auto m2Pid = mother2->PID;
      LV2.SetPtEtaPhiM(track2->PT,track2->Eta,track2->Phi,eMass);
      LV = LV1 + LV2;

      dca1 = (track1->D0/track1->ErrorD0);
      dca2 = (track2->D0/track2->ErrorD0);
      dca = sqrt((dca1*dca1 + dca2*dca2) / 2);

      if(pairULS && !MCpidEle)  hMPtDCA_ULS_rec[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      if(!pairULS && !MCpidEle) hMPtDCA_LS_rec[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      if(pairULS && MCpidEle)  hMPtDCA_ULS_rec_MCpidEle[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      if(!pairULS && MCpidEle) hMPtDCA_LS_rec_MCpidEle[iSce]->Fill(LV.Mag(),LV.Pt(),dca);

      if((!pairULS && !MCpidEle) && (  (fabs(particle1->PID) != 11)||(fabs(particle2->PID)!=11) ))  hMPtDCA_LS_rec_misIDoneLeg[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      if((!pairULS && !MCpidEle) && (  (fabs(particle1->PID) != 11)&&(fabs(particle2->PID)!=11) ))  hMPtDCA_LS_rec_misIDtwoLeg[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      if((!pairULS && !MCpidEle) && (  (fabs(particle1->PID) == 211)&&(fabs(particle2->PID)==211) ))  hMPtDCA_LS_rec_misIDPion[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      if((!pairULS && !MCpidEle) && ( ((fabs(particle1->PID) != 11) && (hasCharmAncestor(particle1, particles) || hasBeautyAncestor(particle1, particles)) ) || ((fabs(particle2->PID)!=11) && (hasCharmAncestor(particle2, particles) || hasBeautyAncestor(particle2, particles)) )) )  hMPtDCA_LS_rec_misIDhf[iSce]->Fill(LV.Mag(),LV.Pt(),dca);

      // if (mother1 == mother2 && !hasHeavyAncestor(particle1, particles) && !hasStrangeAncestor(particle1, particles)){ // same mother and neutral LF particle, pion or eta
      //   if(pairULS && !MCpidEle)  hMPtDCA_ULS_rec_primary[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && !MCpidEle) hMPtDCA_LS_rec_primary[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(pairULS && MCpidEle)  hMPtDCA_ULS_rec_MCpidEle_primary[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && MCpidEle) hMPtDCA_LS_rec_MCpidEle_primary[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
      // if (mother1 != mother2 && hasHeavyAncestor(particle1, particles) && hasHeavyAncestor(particle2, particles)){
      //   if(pairULS && !MCpidEle)  hMPtDCA_ULS_rec_heavy[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && !MCpidEle) hMPtDCA_LS_rec_heavy[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(pairULS && MCpidEle)  hMPtDCA_ULS_rec_MCpidEle_heavy[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && MCpidEle) hMPtDCA_LS_rec_MCpidEle_heavy[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
      // if ( isCharm(m1Pid) && isCharm(m2Pid) && !hasBeautyAncestor(particle1, particles) && !hasBeautyAncestor(particle2, particles) ){
      // // if(mother1 != mother2  && hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles))
      //   if(pairULS && !MCpidEle)  hMPtDCA_ULS_rec_charm[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && !MCpidEle) hMPtDCA_LS_rec_charm[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(pairULS && MCpidEle)  hMPtDCA_ULS_rec_MCpidEle_charm[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && MCpidEle) hMPtDCA_LS_rec_MCpidEle_charm[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
      // if( isBeauty(m1Pid) && isBeauty(m2Pid) ){
      // // if(mother1 != mother2 && hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles) )
      //   if(pairULS && !MCpidEle)  hMPtDCA_ULS_rec_beauty[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && !MCpidEle) hMPtDCA_LS_rec_beauty[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(pairULS && MCpidEle)  hMPtDCA_ULS_rec_MCpidEle_beauty[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS && MCpidEle) hMPtDCA_LS_rec_MCpidEle_beauty[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
    }
  }
}

// function to ULS or LS pair geneated particles
void dileptonPairingGen(std::vector<GenParticle *> vec_track_neg,std::vector<GenParticle *> vec_track_pos, bool pairULS, TClonesArray *particles, Int_t iSce){
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
    // auto imother1 = particle1->M1;
    // auto mother1 = imother1 != -1 ? (GenParticle *)particles->At(imother1) : (GenParticle *)nullptr;
    // auto m1Pid = mother1->PID;
    LV1.SetPtEtaPhiM(particle1->PT,particle1->Eta,particle1->Phi,eMass);

    if(!pairULS) iPos_start = iEle+1;
    for (int iPos = iPos_start; iPos < vec_track_pos_size; iPos++){
      auto track2 = (Track *) vec_track_pos.at(iPos);
      auto particle2 = (GenParticle *) vec_track_pos.at(iPos);
      // auto imother2 = particle2->M1;
      // auto mother2 = imother2 != -1 ? (GenParticle *)particles->At(imother2) : (GenParticle *)nullptr;
      // auto m2Pid = mother2->PID;
      LV2.SetPtEtaPhiM(particle2->PT,particle2->Eta,particle2->Phi,eMass);
      LV = LV1 + LV2;

      dca1 = (track1->D0/track1->ErrorD0);
      dca2 = (track2->D0/track2->ErrorD0);
      dca = sqrt((dca1*dca1 + dca2*dca2) / 2);

      if(pairULS)  hMPtDCA_ULS_gen[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      if(!pairULS) hMPtDCA_LS_gen[iSce]->Fill(LV.Mag(),LV.Pt(),dca);

      // if (mother1 == mother2 && !hasHeavyAncestor(particle1, particles) && !hasStrangeAncestor(particle1, particles)) // same mother and neutral LF particle, pion or eta
      // {
      //   if(pairULS)  hMPtDCA_ULS_gen_primary[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS) hMPtDCA_LS_gen_primary[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
      // if (mother1 != mother2 && hasHeavyAncestor(particle1, particles) && hasHeavyAncestor(particle2, particles))
      // {
      //   if(pairULS)  hMPtDCA_ULS_gen_heavy[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS) hMPtDCA_LS_gen_heavy[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
      // if ( isCharm(m1Pid) && isCharm(m2Pid) && !hasBeautyAncestor(particle1, particles) && !hasBeautyAncestor(particle2, particles) )
      // // if(mother1 != mother2  && hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles))
      // {
      //   if(pairULS)  hMPtDCA_ULS_gen_charm[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS) hMPtDCA_LS_gen_charm[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
      // if( isBeauty(m1Pid) && isBeauty(m2Pid) )
      // // if(mother1 != mother2 && hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles) )
      // {
      //   if(pairULS)  hMPtDCA_ULS_gen_beauty[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      //   if(!pairULS) hMPtDCA_LS_gen_beauty[iSce]->Fill(LV.Mag(),LV.Pt(),dca);
      // }
    }
  }
}
