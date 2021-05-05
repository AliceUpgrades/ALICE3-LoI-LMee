R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libDelphesO2)

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObjArray.h"


// class GenParticle;
// class Track;
// class TClonesArray;
// #include <GenParticle.h>
// #include <Track.h>
// #include <TClonesArray.h>
// #include <DelphesClasses.h>




bool bSmear    = true;
bool bUsePreSh = true;
bool bUseTOF   = true;
bool bUseRICH  = true;
double Bz = 0.2;            // becomes overwritten by the generateEfficiencies.sh skript
double eMass = 0.000511;

// Cinematic cuts on tracks
// double PtCut = 0.0;     // open cuts
double PtCut = 0.04;       // if B=0.2T becomes overwritten by the generateEfficiencies.sh skript
// double PtCut = 0.08;    // if B=0.5T becomes overwritten by the generateEfficiencies.sh skript
double EtaCut = 1.1;
// double EtaCut = 10.;    // open cuts

Double_t pt_binning[]  = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0};
// Double_t eta_binning[]  = {-5.0,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.0,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0};
Int_t nEtaBins = 10;    Double_t maxEta = 5.;
Int_t nPhiBins = 180;    Double_t maxPhi = 2*TMath::Pi();


// TOF
double tof_radius = 100.; // [cm]
double tof_length = 200.; // [cm]
double tof_sigmat = 0.02; // [ns]
double tof_sigma0 = 0.20; // [ns]

// TOF ele pt acceptance
double tof_EleAccep_pt_cut = 0.6;  // [GeV/c]

double tof_mismatch = 0.01;
std::string tof_mismatch_fname;


// RICH
double rich_radius = 100.;
double rich_length = 200.;

// RICH ele pt acceptance
// double rich_Ele_pt_cut = 1.8; // [GeV/c]
double rich_PionRejection_pt_cut = 1.0; // [GeV/c]



// PID Scenarios
bool useTOFPID[]         = {kTRUE,  kTRUE,  kTRUE};
bool useRICHPID[]        = {kFALSE, kTRUE,  kTRUE};
bool usePreShPID[]       = {kFALSE, kFALSE, kTRUE};
// TOF cuts on tracks
double nSigmaTOFEle[]    = {3.0,    3.0,    3.0};
double nSigmaTOFPi[]     = {3.0,    3.0,    3.0};
// RICH cuts on tracks
double nSigmaRICHEle[]   = {3.0,    3.0,    3.0};
double nSigmaRICHPi[]    = {3.0,    4.0,    4.0};
//Preshower nsigma DCA cuts on tracks
double nSigmaPreShD0[]   = {3.0,    3.0,    3.0};
double nSigmaPreShDZ[]   = {3.0,    3.0,    3.0};


// void makeHistNice(TH1* h, int color){
//   h->SetMarkerColor(color);
//   h->SetLineColor(color);
//   h->SetMarkerStyle(20);
//   h->SetLineWidth(2);
// }

bool hasCommonAncestor(GenParticle* p1,GenParticle* tr2 )
{
  // check if 2 particles have one ancestor
  // fill vector with all ancestors of particle on
  // check if particle2 is included
  // if so, return true if not check for mother of p2
  // maybe do this recursiv???
  return false;
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

// kinematic cuts for tracks
bool kineCuts(Track *tr){
  // check pt and eta for track
  // evaluate as true if criterion is passed
  bool pt = tr->PT > PtCut;
  bool eta = abs(tr->Eta) < EtaCut;
  // all have to be true
  return (pt && eta);
}

// kinematic cuts for generated particles
bool kineCuts(GenParticle *pa){
  // check pt and eta for particle
  // evaluate as true if criterion is passed
  bool pt = pa->PT > PtCut;
  bool eta = abs(pa->Eta) < EtaCut;
  // all have to be true
  return (pt && eta);
}

bool hasHeavyAncestor(GenParticle *particle, TClonesArray *particles)
{
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

bool hasBeautyAncestor(GenParticle *particle, TClonesArray *particles)
{
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
  Double_t smearing = ((TH1D *)(fArrResoPt->At(ptbin)))->GetRandom() * pt;
  const Double_t sPt = pt - smearing;

  // smear eta
  ptbin     = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->FindBin(pt);
  ptbin_max = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->GetNbins();
  if (ptbin < 1)
    ptbin = 1;
  else if (ptbin > ptbin_max)
    ptbin = ptbin_max;
  smearing = ((TH1D *)(fArrResoEta->At(ptbin)))->GetRandom();
  const Double_t sEta = eta - smearing;

  // smear phi
  ptbin     = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->FindBin(pt);
  ptbin_max = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->GetNbins();
  if (ptbin < 1)
    ptbin = 1;
  if (ptbin > ptbin_max)
    ptbin = ptbin_max;
  if (ch > 0) {
    smearing = ((TH1D *)(fArrResoPhi_Pos->At(ptbin)))->GetRandom();
  } else if (ch < 0) {
    smearing = ((TH1D *)(fArrResoPhi_Neg->At(ptbin)))->GetRandom();
  }
  const Double_t sPhi = phi - smearing;

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

auto hMPt_ULS_gen         = new TH2F("hMPt_ULS_gen",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_gen_primary = new TH2F("hMPt_ULS_gen_primary",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_gen_heavy   = new TH2F("hMPt_ULS_gen_heavy",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_gen_charm   = new TH2F("hMPt_ULS_gen_charm",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_gen_beauty  = new TH2F("hMPt_ULS_gen_beauty",title2d.c_str(),300.,0,3.,200,0.,10.);
auto hMPt_ULS_rec         = new TH2F("hMPt_ULS_rec",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_primary = new TH2F("hMPt_ULS_rec_primary",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_heavy   = new TH2F("hMPt_ULS_rec_heavy",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_charm   = new TH2F("hMPt_ULS_rec_charm",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_beauty  = new TH2F("hMPt_ULS_rec_beauty",title2d.c_str(),300.,0,3.,200,0.,10.);
auto hMPt_ULS_rec_MCpidEle         = new TH2F("hMPt_ULS_rec_MCpidEle",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_MCpidEle_primary = new TH2F("hMPt_ULS_rec_MCpidEle_primary",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_MCpidEle_heavy   = new TH2F("hMPt_ULS_rec_MCpidEle_heavy",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_MCpidEle_charm   = new TH2F("hMPt_ULS_rec_MCpidEle_charm",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_ULS_rec_MCpidEle_beauty  = new TH2F("hMPt_ULS_rec_MCpidEle_beauty",title2d.c_str(),300.,0,3.,200,0.,10.);
// Pair histograms LS
auto hMPt_LS_gen         = new TH2F("hMPt_LS_gen",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_gen_primary = new TH2F("hMPt_LS_gen_primary",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_gen_heavy   = new TH2F("hMPt_LS_gen_heavy",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_gen_charm   = new TH2F("hMPt_LS_gen_charm",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_gen_beauty  = new TH2F("hMPt_LS_gen_beauty",title2d.c_str(),300.,0,3.,200,0.,10.);
auto hMPt_LS_rec         = new TH2F("hMPt_LS_rec",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_primary = new TH2F("hMPt_LS_rec_primary",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_heavy   = new TH2F("hMPt_LS_rec_heavy",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_charm   = new TH2F("hMPt_LS_rec_charm",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_beauty  = new TH2F("hMPt_LS_rec_beauty",title2d.c_str(),300.,0,3.,200,0.,10.);
auto hMPt_LS_rec_MCpidEle         = new TH2F("hMPt_LS_rec_MCpidEle",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_MCpidEle_primary = new TH2F("hMPt_LS_rec_MCpidEle_primary",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_MCpidEle_heavy   = new TH2F("hMPt_LS_rec_MCpidEle_heavy",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_MCpidEle_charm   = new TH2F("hMPt_LS_rec_MCpidEle_charm",title2d.c_str(),300.,0,3.,200,0.,10.);
// auto hMPt_LS_rec_MCpidEle_beauty  = new TH2F("hMPt_LS_rec_MCpidEle_beauty",title2d.c_str(),300.,0,3.,200,0.,10.);


// function to ULS or LS pair reconstructed particles
void dileptonPairingRec(std::vector<Track *> vec_track_neg,std::vector<Track *> vec_track_pos, bool pairULS, bool MCpidEle, TClonesArray *particles){
  // pairing reconstructed ULS or LS pairs
  double vec_track_neg_size = vec_track_neg.size();
  double vec_track_pos_size = vec_track_pos.size();
  if(!pairULS && (vec_track_neg_size == 1)) return;
  int    iPos_start = 0;
  int    iNeg_end = vec_track_neg_size;
  if(!pairULS) iNeg_end = vec_track_neg_size-1;
  TLorentzVector LV1,LV2,LV;

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

      if(pairULS && !MCpidEle)  hMPt_ULS_rec->Fill(LV.Mag(),LV.Pt());
      if(!pairULS && !MCpidEle) hMPt_LS_rec->Fill(LV.Mag(),LV.Pt());
      if(pairULS && MCpidEle)  hMPt_ULS_rec_MCpidEle->Fill(LV.Mag(),LV.Pt());
      if(!pairULS && MCpidEle) hMPt_LS_rec_MCpidEle->Fill(LV.Mag(),LV.Pt());

      // if (mother1 == mother2 && !hasHeavyAncestor(particle1, particles) && !hasStrangeAncestor(particle1, particles)){ // same mother and neutral LF particle, pion or eta
      //   if(pairULS && !MCpidEle)  hMPt_ULS_rec_primary->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && !MCpidEle) hMPt_LS_rec_primary->Fill(LV.Mag(),LV.Pt());
      //   if(pairULS && MCpidEle)  hMPt_ULS_rec_MCpidEle_primary->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && MCpidEle) hMPt_LS_rec_MCpidEle_primary->Fill(LV.Mag(),LV.Pt());
      // }
      // if (mother1 != mother2 && hasHeavyAncestor(particle1, particles) && hasHeavyAncestor(particle2, particles)){
      //   if(pairULS && !MCpidEle)  hMPt_ULS_rec_heavy->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && !MCpidEle) hMPt_LS_rec_heavy->Fill(LV.Mag(),LV.Pt());
      //   if(pairULS && MCpidEle)  hMPt_ULS_rec_MCpidEle_heavy->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && MCpidEle) hMPt_LS_rec_MCpidEle_heavy->Fill(LV.Mag(),LV.Pt());
      // }
      // if ( isCharm(m1Pid) && isCharm(m2Pid) && !hasBeautyAncestor(particle1, particles) && !hasBeautyAncestor(particle2, particles) ){
      // // if(mother1 != mother2  && hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles))
      //   if(pairULS && !MCpidEle)  hMPt_ULS_rec_charm->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && !MCpidEle) hMPt_LS_rec_charm->Fill(LV.Mag(),LV.Pt());
      //   if(pairULS && MCpidEle)  hMPt_ULS_rec_MCpidEle_charm->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && MCpidEle) hMPt_LS_rec_MCpidEle_charm->Fill(LV.Mag(),LV.Pt());
      // }
      // if( isBeauty(m1Pid) && isBeauty(m2Pid) ){
      // // if(mother1 != mother2 && hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles) )
      //   if(pairULS && !MCpidEle)  hMPt_ULS_rec_beauty->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && !MCpidEle) hMPt_LS_rec_beauty->Fill(LV.Mag(),LV.Pt());
      //   if(pairULS && MCpidEle)  hMPt_ULS_rec_MCpidEle_beauty->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS && MCpidEle) hMPt_LS_rec_MCpidEle_beauty->Fill(LV.Mag(),LV.Pt());
      // }
    }
  }
}

// function to ULS or LS pair geneated particles
void dileptonPairingGen(std::vector<GenParticle *> vec_track_neg,std::vector<GenParticle *> vec_track_pos, bool pairULS, TClonesArray *particles){
  // pairing reconstructed ULS or LS pairs
  double vec_track_neg_size = vec_track_neg.size();
  double vec_track_pos_size = vec_track_pos.size();
  if(!pairULS && (vec_track_neg_size == 1)) return;
  int    iPos_start = 0;
  int    iNeg_end = vec_track_neg_size;
  if(!pairULS) iNeg_end = vec_track_neg_size-1;
  TLorentzVector LV1,LV2,LV;

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

      if(pairULS)  hMPt_ULS_gen->Fill(LV.Mag(),LV.Pt());
      if(!pairULS) hMPt_LS_gen->Fill(LV.Mag(),LV.Pt());

      // if (mother1 == mother2 && !hasHeavyAncestor(particle1, particles) && !hasStrangeAncestor(particle1, particles)) // same mother and neutral LF particle, pion or eta
      // {
      //   if(pairULS)  hMPt_ULS_gen_primary->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS) hMPt_LS_gen_primary->Fill(LV.Mag(),LV.Pt());
      // }
      // if (mother1 != mother2 && hasHeavyAncestor(particle1, particles) && hasHeavyAncestor(particle2, particles))
      // {
      //   if(pairULS)  hMPt_ULS_gen_heavy->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS) hMPt_LS_gen_heavy->Fill(LV.Mag(),LV.Pt());
      // }
      // if ( isCharm(m1Pid) && isCharm(m2Pid) && !hasBeautyAncestor(particle1, particles) && !hasBeautyAncestor(particle2, particles) )
      // // if(mother1 != mother2  && hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles))
      // {
      //   if(pairULS)  hMPt_ULS_gen_charm->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS) hMPt_LS_gen_charm->Fill(LV.Mag(),LV.Pt());
      // }
      // if( isBeauty(m1Pid) && isBeauty(m2Pid) )
      // // if(mother1 != mother2 && hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles) )
      // {
      //   if(pairULS)  hMPt_ULS_gen_beauty->Fill(LV.Mag(),LV.Pt());
      //   if(!pairULS) hMPt_LS_gen_beauty->Fill(LV.Mag(),LV.Pt());
      // }
    }
  }
}


void anaEEstudy(
    const char *inputFile = "delphes100k.root", // one fo the LF and charm part. Has 1M charm+beauty events
    // const char *inputFile = "delphes_500kBeauty.root", // one fo the LF and charm part. Has 1M charm+beauty events
    const char *outputFile = "anaEEstudyLFcc.root" // merge output files after analysis was run to keep file size moderate
    // const char *outputFile = "dcabb.root" // merge output files after analysis was run to keep file size moderate
    // const Int_t nEvents = -1,              // number of events to analyze (all in case of -1)
    // const Bool_t bStarlight = kFALSE       // special treatment in case of STARLight events
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
  // if(nEvents > -1) numberOfEntries = nEvents;
  // Printf("Process %lld out of %lld entries",numberOfEntries,numberOfEntriesTree);

  Int_t nPtBins  = sizeof(pt_binning)/sizeof(*pt_binning) -1;
  Double_t eta_binning[nEtaBins+1];
  Double_t phi_binning[nPhiBins+1];
  for (int iBin = 0; iBin < nEtaBins+1; iBin++) {eta_binning[iBin] = -maxEta + iBin * 2.* maxEta/nEtaBins; /* std::cout << eta_binning[iBin] << " eta bin border " << std::endl;  */}
  for (int iBin = 0; iBin < nPhiBins+1; iBin++) {phi_binning[iBin] = 0. + iBin * maxPhi/nPhiBins;      /* std::cout << phi_binning[iBin] << " phi bin border " << std::endl;  */}


  // Get pointers to branches used in this analysis
  auto events = treeReader->UseBranch("Event");
  auto tracks = treeReader->UseBranch("Track");
  auto particles = treeReader->UseBranch("Particle");

  // Smearing for generated tracks
  if (Bz==0.5) smearingfile = "../resolutionfiles/resolution_test_5kG.root";
  else if (Bz==0.2) smearingfile = "../resolutionfiles/resolution_test_2kG.root";
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
  auto nTracks = new TH1F("nTracks",";Tracks",500,0,500);
  auto nParticles = new TH1F("nParticles",";Particles",500,0,500);

  //Check smearing of generated tracks
  auto hSmearing_For_Eff_pt = new TH2F("hSmearing_For_Eff_pt","(#it{p}_{T,gen}-#it{p}_{T,rec})/#it{p}_{T,gen};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.5,0.5);
  auto hSmearing_For_Eff_phi_pos = new TH2F("hSmearing_For_Eff_phi_pos","#phi_{gen}-#phi_{rec};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.4,0.1);
  auto hSmearing_For_Eff_phi_neg = new TH2F("hSmearing_For_Eff_phi_neg","#phi_{gen}-#phi_{rec};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.1,0.4);
  auto hSmearing_For_Eff_eta = new TH2F("hSmearing_For_Eff_eta","#eta_{gen}-#eta_{rec};#it{p}_{T} (GeV/#it{c})",400,0.,4.,500,-0.1,0.1);

  // Track histograms
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
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_ElePos_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_ElePos_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hNegTrack_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hNegTrack_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hPosTrack_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hPosTrack_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hAllTracks_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hAllTracks_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Ele_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Ele_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Pos_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pos_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Muon_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Muon_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Pion_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Pion_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Kaon_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Kaon_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  for (int j = 0; j < nPIDscenarios; ++j) hTrack_Proton_Rec_Pt_Eta_Phi[j] = new TH3F(Form("hTrack_Proton_Rec_Pt_Eta_Phi_sce%i",j+1),";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_primary_rec[j] = new TH3F(Form("hPt_Eta_Phi_primary_rec_sce%i",j+1), ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_primary_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_primary_Ele_rec_sce%i",j+1), ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_primary_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_primary_Pos_rec_sce%i",j+1), ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_hf_rec[j] = new TH3F(Form("hPt_Eta_Phi_hf_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_hf_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_hf_Ele_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_hf_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_hf_Pos_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_charm_rec[j] = new TH3F(Form("hPt_Eta_Phi_charm_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_charm_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_charm_Ele_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_charm_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_charm_Pos_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_beauty_rec[j] = new TH3F(Form("hPt_Eta_Phi_beauty_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_beauty_Ele_rec[j] = new TH3F(Form("hPt_Eta_Phi_beauty_Ele_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  // for (int j = 0; j < nPIDscenarios; ++j) hPt_Eta_Phi_beauty_Pos_rec[j] = new TH3F(Form("hPt_Eta_Phi_beauty_Pos_rec_sce%i",j+1), ";p_{T} (GeV/c);#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);


  TList* lPIDscenario[nPIDscenarios];
  for (int iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    lPIDscenario[iPIDscenario] = new TList();
    lPIDscenario[iPIDscenario]->SetName(Form("PIDscenario_%i", iPIDscenario+1));
    lPIDscenario[iPIDscenario]->SetOwner();
    lPIDscenario[iPIDscenario]->Add(hTrack_ElePos_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hNegTrack_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hPosTrack_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hAllTracks_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hTrack_Ele_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hTrack_Pos_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hTrack_Muon_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hTrack_Pion_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hTrack_Kaon_Rec_Pt_Eta_Phi[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hTrack_Proton_Rec_Pt_Eta_Phi[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_primary_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_primary_Ele_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_primary_Pos_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_hf_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_hf_Ele_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_hf_Pos_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_charm_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_charm_Ele_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_charm_Pos_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_beauty_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_beauty_Ele_rec[iPIDscenario]);
    // lPIDscenario[iPIDscenario]->Add(hPt_Eta_Phi_beauty_Pos_rec[iPIDscenario]);
    // listRecTracks->Add(lPIDscenario);
  }

                                                                                std::cout << __LINE__ << " ... COUT LINE ..." << std::endl;


  TH3F* hTrack_ElePos_Gen_Pt_Eta_Phi = new TH3F("hTrack_ElePos_Gen_Pt_Eta_Phi",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_Ele_Gen_Pt_Eta_Phi = new TH3F("hTrack_Ele_Gen_Pt_Eta_Phi",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_Pos_Gen_Pt_Eta_Phi = new TH3F("hTrack_Pos_Gen_Pt_Eta_Phi",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_ElePos_GenSmeared_Pt_Eta_Phi = new TH3F("hTrack_ElePos_GenSmeared_Pt_Eta_Phi",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_Ele_GenSmeared_Pt_Eta_Phi = new TH3F("hTrack_Ele_GenSmeared_Pt_Eta_Phi",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
  TH3F* hTrack_Pos_GenSmeared_Pt_Eta_Phi = new TH3F("hTrack_Pos_GenSmeared_Pt_Eta_Phi",";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", nPtBins, pt_binning, nEtaBins, eta_binning, nPhiBins, phi_binning);
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


  auto hPrimaryM = new TH1F("hPrimaryM",";MotherPDG",1000,0.5,1000.5);
  auto hPrimaryGM = new TH1F("hPrimaryGM",";MotherPDG",1000,0.5,1000.5);
  auto hCharmMm = new TH1F("hCharmMm",";MotherPDG",101,399.5,500.5);
  auto hCharmMb = new TH1F("hCharmMb",";MotherPDG",1001,3099.5,5000.5);

  auto hCharmGM = new TH1F("hCharmGM",";MotherPDG",1000,0,10000);
  auto hBeautyM = new TH1F("hBeautyM",";MotherPDG",1000,0,10000);
  auto hBeautyGM = new TH1F("hBeautyGM",";MotherPDG",1000,0,10000);

  auto hBeforeSmearing_Pt_Eta_Phi_rec = new TH3F("hBeforeSmearing_Pt_Eta_Phi_rec", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 500, 0., 10., 200, -10., 10.,180,0,2*TMath::Pi());
  auto hAfterSmearing_Pt_Eta_Phi_rec = new TH3F("hAfterSmearing_Pt_Eta_Phi_rec", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 500, 0., 10., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_1 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_1", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_2 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_2", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_3 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_3", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_4 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_4", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_5 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_5", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_6 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_6", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_7 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_7", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_8 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_8", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_9 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_9", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  // auto hPt_Eta_Phi_rec_Eta_Cut_10 = new TH3F("hPt_Eta_Phi_rec_Eta_Cut_10", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());



  auto hM_Pt_sig_gen = new TH2F("hM_Pt_sig_gen",title2d.c_str(),300.,0,3.,200,0.,20.);



  // TOF histograms
  auto hBetaP_beforeSmearing = new TH2F("hBetaP_beforeSmearing", ";#it{p} GeV/c;#beta", 500, 0., 10., 500, 0.1, 1.1);
  auto hBetaP_afterSmearing = new TH2F("hBetaP_afterSmearing", ";#it{p} GeV/c;#beta", 500, 0., 10., 500, 0.1, 1.1);
  TH2 *hNsigmaP_TOF[5];
  TH2 *hNsigmaP_TOF_trueElec[5];
  TH2 *hNsigmaP_TOF_trueMuon[5];
  TH2 *hNsigmaP_TOF_truePion[5];
  TH2 *hNsigmaP_TOF_trueKaon[5];
  TH2 *hNsigmaP_TOF_trueProton[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueElec[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueMuon[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_truePion[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueKaon[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueProton[5];
  const char *pname[5] = {"el", "mu", "pi", "ka", "pr"};
  const char *plabel[5] = {"e", "#mu", "#pi", "K", "p"};
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF[i] = new TH2F(Form("hNsigmaP_%s_TOF", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueElec[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueElec", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueMuon[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueMuon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_truePion[i] = new TH2F(Form("hNsigmaP_%s_TOF_truePion", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueKaon[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueKaon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueProton[i] = new TH2F(Form("hNsigmaP_%s_TOF_trueProton", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueElec[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueElec", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueMuon[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueMuon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_truePion[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_truePion", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueKaon[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueKaon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueProton[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_TOF_trueProton", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 100.);


  // RICH histograms
  TH2 *hNsigmaBeta_RICH[5];
  TH2 *hDeltaAngleP_RICH[5];
  TH2 *hDeltaAngleP_RICH_trueElec[5];
  TH2 *hDeltaAngleP_RICH_trueMuon[5];
  TH2 *hDeltaAngleP_RICH_truePion[5];
  TH2 *hDeltaAngleP_RICH_trueKaon[5];
  TH2 *hDeltaAngleP_RICH_trueProton[5];
  TH2 *hNsigmaP_RICH[5];
  TH2 *hNsigmaP_RICH_trueElec[5];
  TH2 *hNsigmaP_RICH_trueMuon[5];
  TH2 *hNsigmaP_RICH_truePion[5];
  TH2 *hNsigmaP_RICH_trueKaon[5];
  TH2 *hNsigmaP_RICH_trueProton[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueElec[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueMuon[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_truePion[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueKaon[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueProton[5];
  auto hCherenkovAngleP_beforeSmearing = new TH2F("hCherenkovAngleP_beforeSmearing", ";#it{p} GeV/c; cherenkov angle", 500, 0., 10., 200, 0., 0.4);
  auto hCherenkovAngleP_afterSmearing = new TH2F("hCherenkovAngleP_afterSmearing", ";#it{p} GeV/c; cherenkov angle", 500, 0., 10., 200, 0., 0.4);
  for (int i = 0; i < 5; ++i) hNsigmaBeta_RICH[i] = new TH2F(Form("hNsigmaBeta_%s_RICH", pname[i]), Form(";#beta ;n#sigma_{%s}", plabel[i]), 200, 0., 1.1, 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH[i] = new TH2F(Form("hDeltaAngleP_%s_RICH", pname[i]), Form(";#it{p} GeV/c;#Delta #alpha_{%s}", plabel[i]), 500, 0., 10., 400, -0.4, 0.4);
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueElec[i] = new TH2F(Form("hDeltaAngleP_%s_RICH_trueElec", pname[i]), Form(";#it{p} GeV/c;#Delta #alpha_{%s}", plabel[i]), 500, 0., 10., 400, -0.4, 0.4);
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueMuon[i] = new TH2F(Form("hDeltaAngleP_%s_RICH_trueMuon", pname[i]), Form(";#it{p} GeV/c;#Delta #alpha_{%s}", plabel[i]), 500, 0., 10., 400, -0.4, 0.4);
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_truePion[i] = new TH2F(Form("hDeltaAngleP_%s_RICH_truePion", pname[i]), Form(";#it{p} GeV/c;#Delta #alpha_{%s}", plabel[i]), 500, 0., 10., 400, -0.4, 0.4);
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueKaon[i] = new TH2F(Form("hDeltaAngleP_%s_RICH_trueKaon", pname[i]), Form(";#it{p} GeV/c;#Delta #alpha_{%s}", plabel[i]), 500, 0., 10., 400, -0.4, 0.4);
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueProton[i] = new TH2F(Form("hDeltaAngleP_%s_RICH_trueProton", pname[i]), Form(";#it{p} GeV/c;#Delta #alpha_{%s}", plabel[i]), 500, 0., 10., 400, -0.4, 0.4);
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH[i] = new TH2F(Form("hNsigmaP_%s_RICH", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueElec[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueElec", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueMuon[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueMuon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -80., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_truePion[i] = new TH2F(Form("hNsigmaP_%s_RICH_truePion", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -65., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueKaon[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueKaon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -40., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueProton[i] = new TH2F(Form("hNsigmaP_%s_RICH_trueProton", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -100., 100.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueElec[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueElec", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueMuon[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueMuon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -80., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_truePion[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_truePion", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -65., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueKaon[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueKaon", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -40., 25.);
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueProton[i] = new TH2F(Form("hNsigmaP_%s_afterPIDcuts_RICH_trueProton", pname[i]), Form(";#it{p} GeV/c;n#sigma_{%s}", plabel[i]), 500, 0., 10., 400, -25., 25.);

  // Preshower histograms
  // logx binning
  const Int_t nbins = 80;
  double xmin = 1.e-2;
  double xmax = 1.e2;
  double logxmin = std::log10(xmin);
  double logxmax = std::log10(xmax);
  double binwidth = (logxmax - logxmin) / nbins;
  double xbins[nbins + 1];
  for (Int_t i = 0; i <= nbins; ++i)
    xbins[i] = std::pow(10., logxmin + i * binwidth);

  TH2 *hHit = new TH2F("hHit", ";#eta;#it{p}_{T} (GeV/#it{c})", 200, -4., 4., nbins, xbins);
  auto hMismatchTemplateOut = new TH1F("hMismatchTemplate", "", 3000., -5., 25.);
  // Preshower efficiency input histograms
  TH1* hPreshowerEff[5];
  auto filePreshowerEff = TFile::Open("preshowerEff.root");
  hPreshowerEff[0] = (TH1*)filePreshowerEff->Get("preshowerEffElectrons");
  hPreshowerEff[1] = (TH1*)filePreshowerEff->Get("preshowerEffPions");
  hPreshowerEff[2] = (TH1*)filePreshowerEff->Get("preshowerEffPions");
  hPreshowerEff[3] = (TH1*)filePreshowerEff->Get("preshowerEffPions");
  hPreshowerEff[4] = (TH1*)filePreshowerEff->Get("preshowerEffPions");

  // Preshower QA output histograms
  TH1F* hPreShElectronPt[nPIDscenarios];
  TH1F* hPreShPionPt[nPIDscenarios];
  for (int j = 0; j < nPIDscenarios; ++j) hPreShElectronPt[j] = new TH1F(Form("PreshowerElectronPt_sce%i",j+1), ";#it{p_{T}} (GeV/#it{c});", nbins, xbins);
  for (int j = 0; j < nPIDscenarios; ++j) hPreShPionPt[j]     = new TH1F(Form("PreshowerPionPt_sce%i",j+1), ";#it{p_{T}} (GeV/#it{c});", nbins, xbins);

  for (int iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    lPIDscenario[iPIDscenario]->Add(hPreShElectronPt[iPIDscenario]);
    lPIDscenario[iPIDscenario]->Add(hPreShPionPt[iPIDscenario]);
  }

  // read mismatch template if requested
  TH1 *hMismatchTemplateIn = nullptr;
  if (!tof_mismatch_fname.empty()) {
    auto fmismatch = TFile::Open(tof_mismatch_fname.c_str());
    hMismatchTemplateIn = (TH1 *)fmismatch->Get("hMismatchTemplate");
    hMismatchTemplateIn->SetDirectory(0);
    fmismatch->Close();
  }
  // PID map
  std::map<int, int> pidmap = { {11, 0}, {13, 1}, {211, 2}, {321, 3}, {2212, 4} };


  std::vector<TList*> vecTListPIDscenarios;
  for (size_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    vecTListPIDscenarios.push_back(lPIDscenario[iPIDscenario]);
  }


  std::vector<Track *> vecElectron[nPIDscenarios],vecPositron[nPIDscenarios];
  std::vector<Track *> vecNegTracks[nPIDscenarios],vecPosTracks[nPIDscenarios];
  std::vector<GenParticle *> vecElectronGen,vecPositronGen;
  std::vector<Track *> vecNegTracks_onlyTOFpid[nPIDscenarios], vecPosTracks_onlyTOFpid[nPIDscenarios];


  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);
    nTracks->Fill(tracks->GetEntries());
    nParticles->Fill(particles->GetEntries());

    // std::cout << " NTracks: " << tracks->GetEntries() << std::endl;
    // std::cout << " NParticles: " << particles->GetEntries() << std::endl;

    //##################################################
    //########## reconstructed Track Loop ##############
    //##################################################

    // loop over tracks and fill vectors for reconstructed electrons and positrons
    for (Int_t itrack = 0; itrack < tracks->GetEntries(); ++itrack) {

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

      // look only at electrons
      // if(abs(particle->PID) != 11 ) continue;


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
      if(abs(particle->PID) == 11 ) hAfterSmearing_Pt_Eta_Phi_rec->Fill(track->PT,track->Eta,phiRec);

      // if(abs(track->Eta) < 10.0) hPt_Eta_Phi_rec_Eta_Cut_10->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 9.0)  hPt_Eta_Phi_rec_Eta_Cut_9->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 8.0)  hPt_Eta_Phi_rec_Eta_Cut_8->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 7.0)  hPt_Eta_Phi_rec_Eta_Cut_7->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 6.0)  hPt_Eta_Phi_rec_Eta_Cut_6->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 5.0)  hPt_Eta_Phi_rec_Eta_Cut_5->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 4.0)  hPt_Eta_Phi_rec_Eta_Cut_4->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 3.0)  hPt_Eta_Phi_rec_Eta_Cut_3->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 2.0)  hPt_Eta_Phi_rec_Eta_Cut_2->Fill(track->PT,track->Eta,phiRec);
      // if(abs(track->Eta) < 1.0)  hPt_Eta_Phi_rec_Eta_Cut_1->Fill(track->PT,track->Eta,phiRec);



      auto p = track->P;
      auto beta = toflayer.getBeta(*track);
      auto measAngle = (richdetector.getMeasuredAngle(*track)).first;

      // fill beta-p after semaring
      hBetaP_afterSmearing->Fill(p, beta);
      // fill cherenkovAngle-p after semaring
      hCherenkovAngleP_afterSmearing->Fill(p,measAngle);

      // fill nsigma TOF   (before PID selection is applied)
      std::array<float, 5> PIDdeltatTOF, PIDnsigmaTOF;
      toflayer.makePID(*track, PIDdeltatTOF, PIDnsigmaTOF);
      for (int i = 0; i < 5; ++i) hNsigmaP_TOF[i]->Fill(p, PIDnsigmaTOF[i]);
      if      (abs(particle->PID) == 11)   for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueElec[i]->Fill(p, PIDnsigmaTOF[i]);
      else if (abs(particle->PID) == 13)   for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueMuon[i]->Fill(p, PIDnsigmaTOF[i]);
      else if (abs(particle->PID) == 211)  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_truePion[i]->Fill(p, PIDnsigmaTOF[i]);
      else if (abs(particle->PID) == 321)  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueKaon[i]->Fill(p, PIDnsigmaTOF[i]);
      else if (abs(particle->PID) == 2212) for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueProton[i]->Fill(p, PIDnsigmaTOF[i]);

      // fill nsigma RICH  (before PID selection is applied)
      std::array<float, 5> PIDdeltaangleRICH, PIDnsigmaRICH;
      richdetector.makePID(*track, PIDdeltaangleRICH, PIDnsigmaRICH);
      for (int i = 0; i < 5; ++i) hNsigmaP_RICH[i]->Fill(p, PIDnsigmaRICH[i]);
      for (int i = 0; i < 5; ++i) hNsigmaBeta_RICH[i]->Fill(beta, PIDnsigmaRICH[i]);
      for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH[i]->Fill(p, PIDdeltaangleRICH[i]);
      if      (abs(particle->PID) == 11)   for (int i = 0; i < 5; ++i) {hNsigmaP_RICH_trueElec[i]->Fill(p, PIDnsigmaRICH[i]);   hDeltaAngleP_RICH_trueElec[i]->Fill(p, PIDdeltaangleRICH[i]);}
      else if (abs(particle->PID) == 13)   for (int i = 0; i < 5; ++i) {hNsigmaP_RICH_trueMuon[i]->Fill(p, PIDnsigmaRICH[i]);   hDeltaAngleP_RICH_trueMuon[i]->Fill(p, PIDdeltaangleRICH[i]);}
      else if (abs(particle->PID) == 211)  for (int i = 0; i < 5; ++i) {hNsigmaP_RICH_truePion[i]->Fill(p, PIDnsigmaRICH[i]);   hDeltaAngleP_RICH_truePion[i]->Fill(p, PIDdeltaangleRICH[i]);}
      else if (abs(particle->PID) == 321)  for (int i = 0; i < 5; ++i) {hNsigmaP_RICH_trueKaon[i]->Fill(p, PIDnsigmaRICH[i]);   hDeltaAngleP_RICH_trueKaon[i]->Fill(p, PIDdeltaangleRICH[i]);}
      else if (abs(particle->PID) == 2212) for (int i = 0; i < 5; ++i) {hNsigmaP_RICH_trueProton[i]->Fill(p, PIDnsigmaRICH[i]); hDeltaAngleP_RICH_trueProton[i]->Fill(p, PIDdeltaangleRICH[i]);}



      // kinatic cuts on tracks
      if (!kineCuts(track)) continue;



      // fill output mismatch template (Preshower)
      auto L = std::sqrt(track->XOuter * track->XOuter +
			                   track->YOuter * track->YOuter +
			                   track->ZOuter * track->ZOuter);
      hMismatchTemplateOut->Fill(track->TOuter * 1.e9 - L / 299.79246);

      // do some random mismatch
      if (hMismatchTemplateIn && gRandom->Uniform() < tof_mismatch) {
        track->TOuter = (hMismatchTemplateIn->GetRandom() + L / 299.79246) * 1.e-9;
      }



      for (size_t iPID_scenario = 0; iPID_scenario < nPIDscenarios; iPID_scenario++) {
        //setting variables for i-th PID scenario
        bool   i_useTOFPID = useTOFPID[iPID_scenario];
        bool   i_useRICHPID = useRICHPID[iPID_scenario];
        bool   i_usePreShPID = usePreShPID[iPID_scenario];
        double i_SigmaTOFEle = nSigmaTOFEle[iPID_scenario];
        double i_SigmaTOFPi = nSigmaTOFPi[iPID_scenario];
        double i_SigmaRICHEle = nSigmaRICHEle[iPID_scenario];
        double i_SigmaRICHPi = nSigmaRICHPi[iPID_scenario];
        double i_SigmaPreShD0 = nSigmaPreShD0[iPID_scenario];
        double i_SigmaPreShDZ = nSigmaPreShDZ[iPID_scenario];


        // std::cout << " *** iScenario = " << iPID_scenario << ", useTOFPID = " << i_useTOFPID << std::endl;
        // std::cout << " *** iScenario = " << iPID_scenario << ", useRICHPID = " << i_useRICHPID << std::endl;
        // std::cout << " *** iScenario = " << iPID_scenario << ", selceted nSigmaTOFEle = " << i_SigmaTOFEle << std::endl;
        // std::cout << " *** iScenario = " << iPID_scenario << ", selceted nSigmaTOFPion = " << i_SigmaTOFPi << std::endl;
        // std::cout << " *** iScenario = " << iPID_scenario << ", selceted nSigmaRICHEle = " << i_SigmaRICHEle << std::endl;
        // std::cout << " *** iScenario = " << iPID_scenario << ", selceted nSigmaRICHPion = " << i_SigmaRICHPi << std::endl;


        // ################## PID selection ##################
        // PID selection, identifying particles by using TOF "OR" RICH.
        // using TOF  electron inclusion and pion rejection
        // using RICH electron inclusion and pion rejection
        bool TOFpid = false;
        bool RICHpid = false;
        // Preshower PID
        if (i_usePreShPID && (fabs(track->D0 / track->ErrorD0) > i_SigmaPreShD0) ) continue; // select primaries based on 3 sigma DCA cuts
        if (i_usePreShPID && (fabs(track->DZ / track->ErrorDZ) > i_SigmaPreShDZ) ) continue; // select primaries based on 3 sigma DCA cuts

        //TOF PID
        if(i_useTOFPID && (toflayer.hasTOF(*track)) && (p < tof_EleAccep_p_cut)) {
          if (p > tof_PionRej_p_cut) {
            if((fabs(PIDnsigmaTOF[0]) < i_SigmaTOFEle) && (fabs(PIDnsigmaRICH[0]) < nSigmaRICHEle_forTOFPID)) TOFpid = true; // is within 3 sigma of the electron band (TOF)
          }
          else {
            if(fabs(PIDnsigmaTOF[0]) < i_SigmaTOFEle) TOFpid = true; // is within 3 sigma of the electron band (TOF)
          }

          if(fabs(PIDnsigmaTOF[2]) < i_SigmaTOFPi) TOFpid = false; // is within 3 sigma of the pion band (TOF)
        }

        //RICH PID
        if(i_useRICHPID && richdetector.hasRICH(*track)) {
          if(fabs(PIDnsigmaRICH[0]) < i_SigmaRICHEle) RICHpid = true; // is within 3 sigma of the electron band (RICH)
          if(fabs( (PIDnsigmaRICH[2]) < i_SigmaRICHPi) && (p > rich_PionRejection_p_cut) ) RICHpid = false; // is within 3 sigma of the pion band (RICH)
        }

        if (!(RICHpid || TOFpid)) continue; // check if TOF or RICH signal is true.
        // ################## end of PID selection ##################
        // std::cout << "particle is accepted -> set acceptedByPID kTURE" << std::endl;
        acceptedByPID[iPID_scenario] = kTRUE;



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


        // fill nsigma after PID cuts     (after PID selection has been applied)
        // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF[i]->Fill(p, PIDnsigmaTOF[i]);
        // for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH[i]->Fill(p, PIDnsigmaRICH[i]);
        // NSigma plots, separating for true particle ID
        // if      (abs(particle->PID) == 11)   for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueElec[i]->Fill(p, PIDnsigmaTOF[i]);
        // else if (abs(particle->PID) == 13)   for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->Fill(p, PIDnsigmaTOF[i]);
        // else if (abs(particle->PID) == 211)  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_truePion[i]->Fill(p, PIDnsigmaTOF[i]);
        // else if (abs(particle->PID) == 321)  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueKaon[i]->Fill(p, PIDnsigmaTOF[i]);
        // else if (abs(particle->PID) == 2212) for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueProton[i]->Fill(p, PIDnsigmaTOF[i]);
        //
        // if      (abs(particle->PID) == 11)   for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueElec[i]->Fill(p, PIDnsigmaRICH[i]);
        // else if (abs(particle->PID) == 13)   for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->Fill(p, PIDnsigmaRICH[i]);
        // else if (abs(particle->PID) == 211)  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_truePion[i]->Fill(p, PIDnsigmaRICH[i]);
        // else if (abs(particle->PID) == 321)  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueKaon[i]->Fill(p, PIDnsigmaRICH[i]);
        // else if (abs(particle->PID) == 2212) for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueProton[i]->Fill(p, PIDnsigmaRICH[i]);


        // check if track is identified as an electron by the Preshower
        auto pdg  = std::abs(track->PID);
        auto pbin = hPreshowerEff[pidmap[pdg]]->FindBin(p);
        auto eff  = hPreshowerEff[pidmap[pdg]]->GetBinContent(pbin);
        if (gRandom->Uniform() < eff){
          // and fill some QA histograms for Preshower PID
          if(!pidmap[pdg]) hPreShElectronPt[iPID_scenario]->Fill(track->PT);
          else hPreShPionPt[iPID_scenario]->Fill(track->PT);
        }


        // fill track vector
        if (track->Charge < 0) vecNegTracks[iPID_scenario].push_back(track);       //  vector of all reconstructed negative tracks
        else if (track->Charge > 0) vecPosTracks[iPID_scenario].push_back(track);  //  vector of all reconstructed positve tracks

        if (particle->PID == 11 ) vecElectron[iPID_scenario].push_back(track);       // vector of reconstructed electrons
        else if (particle->PID == -11 ) vecPositron[iPID_scenario].push_back(track); // vector of reconstructed positrons

        if (track->Charge < 0) vecNegTracks_onlyTOFpid[iPID_scenario].push_back(track);      // vector of all reconstructed negative tracks selcted by TOF only
        else if (track->Charge > 0) vecPosTracks_onlyTOFpid[iPID_scenario].push_back(track); // vector of all reconstructed positve tracks selcted by TOF only

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

      // fill hit histogram with true (eta,pt)
      hHit->Fill(particle->Eta, particle->PT);


    }



    //##################################################
    //############## generated Track Loop ##############
    //##################################################
    TLorentzVector LV; // for smearing of generated particles
    TLorentzVector LV_smeared; // for smearing of generated particles

    // loop over particles and fill vectors for generated electrons and positrons
    for (Int_t iparticle = 0; iparticle < particles->GetEntries(); ++iparticle) {

      // get particle
      auto particle = (GenParticle *)particles->At(iparticle);

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

      double phiGen = TVector2::Phi_0_2pi(particle->Phi);

      // kinatic cuts on particles
      if (!kineCuts(particle)) continue;

        // look only at electrons
      if(abs(particle->PID) != 11 ) continue;

      // Smeared generated particles
     LV.SetPtEtaPhiM(particle->PT,particle->Eta,phiGen,eMass);
     Double_t charge = -particle->PID;
     LV_smeared = ApplySmearing(fArrResoPt,fArrResoEta,fArrResoPhi_Pos,fArrResoPhi_Neg,LV,charge);
     if(charge > 0.) hSmearing_For_Eff_phi_pos->Fill(particle->PT, LV.Phi() - LV_smeared.Phi());
     if(charge < 0.) hSmearing_For_Eff_phi_neg->Fill(particle->PT, LV.Phi() - LV_smeared.Phi());
     hSmearing_For_Eff_eta->Fill(particle->PT, LV.Eta() - LV_smeared.Eta());
     if(LV.Pt() > 0.) hSmearing_For_Eff_pt->Fill(particle->PT, (LV.Pt() - LV_smeared.Pt())/LV.Pt());

      if (particle->PID == 11 ) vecElectronGen.push_back(particle);       // vector filled with generated electrons
      else if (particle->PID == -11 ) vecPositronGen.push_back(particle); // vector filled with generated positrons

                                      hTrack_ElePos_Gen_Pt_Eta_Phi->Fill(particle->PT,particle->Eta,phiGen); // Pt Eta Phi of generated electrons + positrons
      if (particle->PID == 11)        hTrack_Ele_Gen_Pt_Eta_Phi->Fill(particle->PT,particle->Eta,phiGen);    // Pt Eta Phi of generated electrons
      else if (particle->PID == -11)  hTrack_Pos_Gen_Pt_Eta_Phi->Fill(particle->PT,particle->Eta,phiGen);    // Pt Eta Phi of generated positrons

                                      hTrack_ElePos_GenSmeared_Pt_Eta_Phi->Fill(LV_smeared.Pt(),LV_smeared.Eta(),LV_smeared.Phi()); // Pt Eta Phi of generated smeared electrons + positrons
      if (particle->PID == 11)        hTrack_Ele_GenSmeared_Pt_Eta_Phi->Fill(LV_smeared.Pt(),LV_smeared.Eta(),LV_smeared.Phi());    // Pt Eta Phi of generated smeared electrons
      else if (particle->PID == -11)  hTrack_Pos_GenSmeared_Pt_Eta_Phi->Fill(LV_smeared.Pt(),LV_smeared.Eta(),LV_smeared.Phi());    // Pt Eta Phi of generated smeared positrons

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


    //##################################################
    //############   ULS and LS pairing   ##############
    //############ for Gen and Rec tracks ##############
    //##################################################

    bool pairULS = kTRUE;
    bool MCpidEle = kTRUE;
    dileptonPairingRec(vecNegTracks[0], vecPosTracks[0],  pairULS, !MCpidEle, particles);
    dileptonPairingRec(vecNegTracks[0], vecNegTracks[0], !pairULS, !MCpidEle, particles);
    dileptonPairingRec(vecPosTracks[0], vecPosTracks[0], !pairULS, !MCpidEle, particles);
    dileptonPairingRec(vecElectron[0],  vecPositron[0],   pairULS,  MCpidEle, particles);
    dileptonPairingRec(vecElectron[0],  vecElectron[0],  !pairULS,  MCpidEle, particles);
    dileptonPairingRec(vecPositron[0],  vecPositron[0],  !pairULS,  MCpidEle, particles);

    for (size_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
      vecNegTracks[iPIDscenario].clear();
      vecPosTracks[iPIDscenario].clear();
      vecElectron[iPIDscenario].clear();
      vecPositron[iPIDscenario].clear();
    }


    pairULS = kTRUE;
    dileptonPairingGen(vecElectronGen, vecPositronGen, pairULS, particles);
    dileptonPairingGen(vecElectronGen, vecElectronGen, !pairULS, particles);
    dileptonPairingGen(vecPositronGen, vecPositronGen, !pairULS, particles);

    vecElectronGen.clear();
    vecPositronGen.clear();

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
  listGenTracks->Write();

  fout->cd();
  fout->mkdir("reconstructed");
  fout->cd("reconstructed");

  for (size_t iPIDscenario = 0; iPIDscenario < nPIDscenarios; iPIDscenario++) {
    fout->mkdir(Form("reconstructed/%s",vecTListPIDscenarios.at(iPIDscenario)->GetName()));
    fout->cd(Form("reconstructed/%s",vecTListPIDscenarios.at(iPIDscenario)->GetName()));
    vecTListPIDscenarios.at(iPIDscenario)->Write();
    vecTListPIDscenarios.at(iPIDscenario)->Clear();
  }
  vecTListPIDscenarios.clear();


  fout->cd();
  hSmearing_For_Eff_phi_pos->Write();
  hSmearing_For_Eff_phi_neg->Write();
  hSmearing_For_Eff_eta->Write();
  hSmearing_For_Eff_pt->Write();

  nTracks->Write();
  nParticles->Write();
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
  hTrack_ElePos_Gen_Pt_Eta_Phi->Write();
  hTrack_Ele_Gen_Pt_Eta_Phi->Write();
  hTrack_Pos_Gen_Pt_Eta_Phi->Write();
  hTrack_ElePos_GenSmeared_Pt_Eta_Phi->Write();
  hTrack_Ele_GenSmeared_Pt_Eta_Phi->Write();
  hTrack_Pos_GenSmeared_Pt_Eta_Phi->Write();
  hPrimaryM->Write();
  hPrimaryGM->Write();
  hCharmMm->Write();
  hCharmMb->Write();
  hCharmGM->Write();
  hBeautyM->Write();
  hBeautyGM->Write();
  hBeforeSmearing_Pt_Eta_Phi_rec->Write();
  hAfterSmearing_Pt_Eta_Phi_rec->Write();
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

  hMPt_ULS_gen->Write();
  // hMPt_ULS_gen_primary->Write();
  // hMPt_ULS_gen_heavy->Write();
  // hMPt_ULS_gen_charm->Write();
  // hMPt_ULS_gen_beauty->Write();
  hMPt_LS_gen->Write();
  // hMPt_LS_gen_primary->Write();
  // hMPt_LS_gen_heavy->Write();
  // hMPt_LS_gen_charm->Write();
  // hMPt_LS_gen_beauty->Write();
  hMPt_ULS_rec->Write();
  // hMPt_ULS_rec_primary->Write();
  // hMPt_ULS_rec_heavy->Write();
  // hMPt_ULS_rec_charm->Write();
  // hMPt_ULS_rec_beauty->Write();
  hMPt_ULS_rec_MCpidEle->Write();
  // hMPt_ULS_rec_MCpidEle_primary->Write();
  // hMPt_ULS_rec_MCpidEle_heavy->Write();
  // hMPt_ULS_rec_MCpidEle_charm->Write();
  // hMPt_ULS_rec_MCpidEle_beauty->Write();
  hMPt_LS_rec->Write();
  // hMPt_LS_rec_primary->Write();
  // hMPt_LS_rec_heavy->Write();
  // hMPt_LS_rec_charm->Write();
  // hMPt_LS_rec_beauty->Write();
  hMPt_LS_rec_MCpidEle->Write();
  // hMPt_LS_rec_MCpidEle_primary->Write();
  // hMPt_LS_rec_MCpidEle_heavy->Write();
  // hMPt_LS_rec_MCpidEle_charm->Write();
  // hMPt_LS_rec_MCpidEle_beauty->Write();

  // fout->cd("generated/");
  // fout->mkdir("generated/LS");
  // fout->cd("generated/LS");
  // fout->cd();
  // fout->mkdir("reconstructed/");
  // fout->cd("reconstructed/");
  // fout->mkdir("reconstructed/ULS");
  // fout->cd("reconstructed/ULS");

  // hPt_Eta_Phi_rec_Eta_Cut_10->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_9->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_8->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_7->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_6->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_5->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_4->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_3->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_2->Write();
  // hPt_Eta_Phi_rec_Eta_Cut_1->Write();

  hM_Pt_sig_gen->Write();

  hBetaP_beforeSmearing->Write();
  hBetaP_afterSmearing->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueElec[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueMuon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_truePion[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueKaon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_TOF_trueProton[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueElec[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_truePion[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueKaon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_TOF_trueProton[i]->Write();


  hCherenkovAngleP_beforeSmearing->Write();
  hCherenkovAngleP_afterSmearing->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueElec[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueMuon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_truePion[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueKaon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_RICH_trueProton[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueElec[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_truePion[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueKaon[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaP_afterPIDcuts_RICH_trueProton[i]->Write();
  for (int i = 0; i < 5; ++i) hNsigmaBeta_RICH[i]->Write();
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH[i]->Write();
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueElec[i]->Write();
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueMuon[i]->Write();
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_truePion[i]->Write();
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueKaon[i]->Write();
  for (int i = 0; i < 5; ++i) hDeltaAngleP_RICH_trueProton[i]->Write();

  hHit->Write();
  hMismatchTemplateOut->Write();

  // fout->cd("reconstructed/");

  // fout->mkdir("reconstructed/LS");
  // fout->cd("reconstructed/LS");
  fout->Close();



  // stop watch
  watch->Stop();
  watch->Print();

}
