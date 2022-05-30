//******************************************************************
// ALICE 3 hits to AO2D converter tool
//
// This macro imports ALICE 3 tracker hits, runs the ITS tracker
// on them and then saves the tracks and primary vertex in a
// AO2D.root file that is compliant with the O2 framework.
//
// More specifically, it mimics Run 2-converted data so that
// any analysis geared towards that can run on the output of this
// conversion tool.
//
// To be used compiled ('root.exe -q -b ALICE3toAO2D.C+')
//
// Comments, complaints, suggestions? Please write to:
// --- david.dobrigkeit.chinellato@cern.ch
//******************************************************************

#if !defined(__CLING__) || defined(__ROOTCLING__)
//#include "gandiva/projector.h"
#include <string>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TGeoGlobalMagField.h>
#include <vector>
#include <TTimeStamp.h>
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSMFTSimulation/Hit.h"
#include "ITStracking/Configuration.h"
#include "ITStracking/IOUtils.h"
#include "ITStracking/Tracker.h"
#include "ITStracking/TrackerTraitsCPU.h"
#include "ITStracking/Vertexer.h"
#include "ITStracking/VertexerTraits.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ReconstructionDataFormats/TrackLTIntegral.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "SimulationDataFormat/MCTrack.h"
#include "MathUtils/Cartesian.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "ReconstructionDataFormats/TrackParametrizationWithError.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "Framework/DataTypes.h"
#include "Run2LikeAO2DHelper.h"
#include "CommonConstants/PhysicsConstants.h"

#include <Generators/GeneratorFactory.h>
#include "FairPrimaryGenerator.h"
#include "FairGenerator.h"
#include "FairBoxGenerator.h"
#include <FairLogger.h>
#include <SimConfig/SimConfig.h>
#include <Generators/GeneratorFromFile.h>

#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "Framework/DataTypes.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "ITSBase/GeometryTGeo.h"
//#include "AnalysisCore/MC.h"

#include "TrackSmearer.hh"
#include "TrackUtils.hh"

#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "CommonDataFormat/BunchFilling.h"
#endif

const Double_t lIBresolution = 0.00025f;
const Double_t lOBresolution = 0.00100f;
const Double_t liTOFresolution = 20; //ps
const Double_t loTOFresolution = 20; //ps
const Int_t lNStepsLIntegrator = 200;

//Tracking pass configuration
const Float_t lFirstOpening = 1.;
const Float_t lOpening = 3.;
const Float_t lLastOpening = 6.;

using o2::its::MemoryParameters;
using o2::its::TrackingParameters;
using o2::itsmft::Hit;
using std::string;

//using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

constexpr bool kUseSmearing{true};

const Double_t lMagneticField = 5.0;

// Class to hold the track information for the O2 vertexing
class TrackAlice3 : public o2::track::TrackParCov
{
  using TimeEst = o2::dataformats::TimeStampWithError<float, float>;
  
public:
  TrackAlice3() = default;
  ~TrackAlice3() = default;
  TrackAlice3(const TrackAlice3& src) = default;
  TrackAlice3(const o2::track::TrackParCov& src, const float t = 0, const float te = 1, const int label = 0) : o2::track::TrackParCov(src), timeEst{t, te}, mLabel{label} {}
  const TimeEst& getTimeMUS() const { return timeEst; }
  const int mLabel;
  TimeEst timeEst; ///< time estimate in ns
};

class SmearO2KineGenerator : public o2::eventgen::GeneratorFromO2Kine
{
public:
  SmearO2KineGenerator(const char *name) : GeneratorFromO2Kine(name) { };
  bool Init() override {
    auto retval = o2::eventgen::GeneratorFromO2Kine::Init();
    setContinueMode(true);
    return retval; };
  //  bool importParticles() override {
protected:
  
};

float getDetLengthFromEta(const float eta, const float radius)
{
  return 10. * (10. + radius * std::cos(2 * std::atan(std::exp(-eta))));
}

Double_t TrackLength( o2::track::TrackParCov track, Double_t lX0, Double_t lX1 ){
  std::array<float, 3> lPointN;
  std::array<float, 3> lPointNplus;
  Double_t lLength = 0.0;
  track.propagateTo(lX0, lMagneticField);
  for(Int_t iStep=1; iStep<lNStepsLIntegrator; iStep++){
    track.getXYZGlo(lPointN);
    Float_t lPosition = lX0 + (lX1-lX0)*((Float_t)(iStep))/((Float_t)(lNStepsLIntegrator-1));
    track.propagateTo(lPosition, lMagneticField);
    track.getXYZGlo(lPointNplus);
    lLength += std::hypot( lPointNplus[0]-lPointN[0], lPointNplus[1]-lPointN[1], lPointNplus[2]-lPointN[2] );
  }
  return lLength;
}

/// get mass from TParticlePDG
double getMass(int input_pdg){
  double mass = 0;
  if(TDatabasePDG::Instance()){
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(input_pdg);
    if(particle)  mass = particle->Mass();
    else      std::cout << "===> particle mass equal to 0" << std::endl;
  }
  return mass;
}

Double_t Velocity(Double_t lMomentum, Double_t lMass){
  //Momentum p and mass m -> returns speed in centimeters per picosecond
  //Useful for TOF calculations
  Double_t lA = TMath::Power(lMomentum / lMass, 2);
  return 0.0299792458*TMath::Sqrt(lA/(1+lA));
}

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

void MultiALICE3toAO2D(TString lOutputName = "AO2D.root", TString lDataPath = "./", TString lLutPath =  "")
{
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  std::cout << "\e[1;31m      ALICE 3 hits to AO2D converter tool      \e[0;00m" << std::endl;
  std::cout << "\e[1;31m                12-layer version               \e[0;00m" << std::endl;
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;

  cout<<"Main configuration: "<<endl;
  cout<<"Data path for analysis.......: "<<lDataPath.Data()<<endl;
  cout<<"Data path for LUTs...........: "<<lLutPath.Data()<<endl;
  cout<<"Output name..................: "<<lOutputName.Data()<<endl;
  
  //Open Geometry, please
  o2::base::GeometryManager::loadGeometry(lDataPath.Data());
  
  //Open GRP
  const auto grp = o2::parameters::GRPObject::loadFrom(Form("%so2sim_grp.root",lDataPath.Data()));
  if (!grp) {
    LOG(FATAL) << "Cannot run w/o GRP object";
  }

  o2::base::Propagator::initFieldFromGRP(grp);
  auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  if (!field) {
    LOG(FATAL) << "Failed to load ma";
  }
  double origD[3] = {0., 0., 0.};
  
  std::cout << "*- Starting..." << std::endl;
  const string hitsFileName = Form("%so2sim_HitsTRK.root", lDataPath.Data());
  //TChain mcTree("o2sim");
  //mcTree.AddFile("o2sim_Kine.root");
  //mcTree.SetBranchStatus("*", 0); //disable all branches
  //mcTree.SetBranchStatus("MCTrack*", 1);
  //mcTree.SetBranchStatus("MCEventHeader.", 1);
  
  TFile *f = new TFile(Form("%so2sim_Kine.root", lDataPath.Data()),"READ") ;
  
  auto mcTree = (TTree*)f->Get("o2sim");
  o2::dataformats::MCEventHeader* mcHead = nullptr;
  mcTree->SetBranchStatus("*", 0); //disable all branches
  mcTree->SetBranchStatus("MCTrack*", 1);
  mcTree->SetBranchStatus("MCEventHeader.*", 1);
  mcTree->SetBranchAddress("MCEventHeader.", &mcHead);
  
  std::vector<o2::MCTrack>* mcArr = nullptr;
  mcTree->SetBranchAddress("MCTrack", &mcArr);

  //o2::dataformats::MCEventHeader* mcHead;
  //FairMCEventHeader *mcHead;
  //auto mcHead = new o2::dataformats::MCEventHeader;
  //o2::dataformats::MCEventHeader* mcHead = nullptr;
  //mcTree->SetBranchAddress("MCEventHeader.", &mcHead);

  o2::its::Vertexer vertexer(new o2::its::VertexerTraits());

  TChain itsHits("o2sim");

  itsHits.AddFile(hitsFileName.data());

  o2::its::Tracker tracker(new o2::its::TrackerTraitsCPU());
  tracker.setBz(lMagneticField);

  std::uint32_t roFrame;
  std::vector<Hit>* hits = nullptr;
  itsHits.SetBranchAddress("TRKHit", &hits);
  
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  // TRACKING PARAMETERS
  std::vector<float> LayerRadii = {0.5f, 1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f, 100.0f};
  std::vector<float> LayerZ(12);
  for (int i{0}; i < 12; ++i)
  LayerZ[i] = getDetLengthFromEta(1.44, LayerRadii[i]) + 1.;

  // EXPERIMENTAL THREE-PASS VERSION
  std::vector<TrackingParameters> trackParams(3);
  //Tracking parameters for 12 layer setup
  trackParams[0].NLayers = 12;
  trackParams[0].MinTrackLength = 8; //this is the one with fixed params

  //loosely based on run_trac_alice3.C but with extra stuff for the innermost layers
  //FIXME: This may be subject to further tuning and is only a first guess
  std::vector<float> TrackletMaxDeltaZ = {lFirstOpening*0.1f, lFirstOpening*0.1f, lFirstOpening*0.1f, lFirstOpening*0.1f, lFirstOpening*0.3f, lFirstOpening*0.3f, lFirstOpening*0.3f, lFirstOpening*0.3f, lFirstOpening*0.5f, lFirstOpening*0.5f, lFirstOpening*0.5f};
  //std::vector<float> TrackletMaxDeltaZ = {0.5, 0.9, 0.9, 2.2, 3.4, 5.4, 6.7, 9.9, 9.9, 14, 14};;
  std::vector<float> CellMaxDCA = {lFirstOpening*0.05f, lFirstOpening*0.05f, lFirstOpening*0.05f, lFirstOpening*0.04f, lFirstOpening*0.05f, lFirstOpening*0.2f, lFirstOpening*0.4f, lFirstOpening*0.5f, lFirstOpening*0.5f, lFirstOpening*0.5f};
  std::vector<float> CellMaxDeltaZ = {lFirstOpening*0.2f, lFirstOpening*0.2f, lFirstOpening*0.2f,lFirstOpening*0.4f, lFirstOpening*0.5f, lFirstOpening*0.6f, lFirstOpening*3.0f, lFirstOpening*3.0f, lFirstOpening*3.0f,lFirstOpening*3.0f};
  std::vector<float> NeighbourMaxDeltaCurvature = {lFirstOpening*0.012f, lFirstOpening*0.010f, lFirstOpening*0.008f, lFirstOpening*0.0025f, lFirstOpening*0.003f, lFirstOpening*0.0035f, lFirstOpening*0.004f, lFirstOpening*0.004f, lFirstOpening*0.005f};
  std::vector<float> NeighbourMaxDeltaN = {lFirstOpening*0.002f, lFirstOpening*0.002f, lFirstOpening*0.002f, lFirstOpening*0.0090f, lFirstOpening*0.002f, lFirstOpening*0.005f,lFirstOpening* 0.005f,lFirstOpening* 0.005f, lFirstOpening*0.005f};

  //Second pass - looser cuts, shorter tracks
  std::vector<float> SecondTrackletMaxDeltaZ = {lOpening*0.5f, lOpening*0.9f, lOpening*0.9f, lOpening*2.2f, lOpening*3.4f, lOpening*5.4f, lOpening*6.7f, lOpening*9.9f, lOpening*9.9f, lOpening*14.0f, lOpening*14.0f};;
  std::vector<float> SecondCellMaxDCA = {lOpening*0.05f, lOpening*0.05f, lOpening*0.05f, lOpening*0.04f, lOpening*0.05f, lOpening*0.2f, lOpening*0.4f, lOpening*0.5f, lOpening*0.5f, lOpening*0.5f};
  std::vector<float> SecondCellMaxDeltaZ = {lOpening*0.2f, lOpening*0.2f, lOpening*0.2f, lOpening*0.4f, lOpening*0.5f, lOpening*0.6f, lOpening*3.0f, lOpening*3.0f, lOpening*3.0f, lOpening*3.0f};
  std::vector<float> SecondNeighbourMaxDeltaCurvature = {lOpening*0.012f, lOpening*0.010f, lOpening*0.008f, lOpening*0.0025f, lOpening*0.003f, lOpening*0.0035f, lOpening*0.004f, lOpening*0.004f, lOpening*0.005f};
  std::vector<float> SecondNeighbourMaxDeltaN = {lOpening*0.002f, lOpening*0.002f, lOpening*0.002f, lOpening*0.0090f, lOpening*0.002f, lOpening*0.005f, lOpening*0.005f, lOpening*0.005f, lOpening*0.005f};

  //Third pass: loosest, pick up more secondaries
  std::vector<float> ThirdTrackletMaxDeltaZ = {lLastOpening*0.5f, lLastOpening*0.9f, lLastOpening*0.9f, lLastOpening*2.2f, lLastOpening*3.4f, lLastOpening*5.4f, lLastOpening*6.7f, lLastOpening*9.9f, lLastOpening*9.9f, lLastOpening*14.0f, lLastOpening*14.0f};;
  std::vector<float> ThirdCellMaxDCA = {lLastOpening*0.05f, lLastOpening*0.05f, lLastOpening*0.05f, lLastOpening*0.04f, lLastOpening*0.05f, lLastOpening*0.2f, lLastOpening*0.4f, lLastOpening*0.5f, lLastOpening*0.5f, lLastOpening*0.5f};
  std::vector<float> ThirdCellMaxDeltaZ = {lLastOpening*0.2f, lLastOpening*0.2f, lLastOpening*0.2f, lLastOpening*0.4f, lLastOpening*0.5f, lLastOpening*0.6f, lLastOpening*3.0f, lLastOpening*3.0f, lLastOpening*3.0f, lLastOpening*3.0f};
  std::vector<float> ThirdNeighbourMaxDeltaCurvature = {lLastOpening*0.012f, lLastOpening*0.010f, lLastOpening*0.008f, lLastOpening*0.0025f, lLastOpening*0.003f, lLastOpening*0.0035f, lLastOpening*0.004f, lLastOpening*0.004f, lLastOpening*0.005f};
  std::vector<float> ThirdNeighbourMaxDeltaN = {lLastOpening*0.002f, lLastOpening*0.002f, lLastOpening*0.002f, lLastOpening*0.0090f, lLastOpening*0.002f, lLastOpening*0.005f, lLastOpening*0.005f, lLastOpening*0.005f, lLastOpening*0.005f};

  trackParams[0].LayerRadii = LayerRadii;
  trackParams[0].LayerZ = LayerZ;
  trackParams[0].TrackletMaxDeltaPhi = lFirstOpening*0.3;
  trackParams[0].CellMaxDeltaPhi = lFirstOpening*0.15;
  trackParams[0].CellMaxDeltaTanLambda = lFirstOpening*0.03;
  trackParams[0].TrackletMaxDeltaZ = TrackletMaxDeltaZ;
  trackParams[0].CellMaxDCA = CellMaxDCA;
  trackParams[0].CellMaxDeltaZ = CellMaxDeltaZ;
  trackParams[0].NeighbourMaxDeltaCurvature = NeighbourMaxDeltaCurvature;
  trackParams[0].NeighbourMaxDeltaN = NeighbourMaxDeltaN;

  Float_t lFactor = 1000;
  std::vector<MemoryParameters> memParams(3);
  std::vector<float> CellsMemoryCoefficients = {2.3208e-08f * 100, 2.104e-08f * 100, 1.6432e-08f * 100, 1.2412e-08f * 100, 1.3543e-08f * 100, 1.5e-08f * 100, 1.6e-08f * 100, 1.7e-08f * 100};
  std::vector<float> TrackletsMemoryCoefficients = {0.0016353f * lFactor, 0.0013627f * lFactor, 0.000984f * lFactor, 0.00078135f * lFactor, 0.00057934f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor};

  memParams[0].CellsMemoryCoefficients = CellsMemoryCoefficients;
  memParams[0].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
  memParams[0].MemoryOffset = 180000;

  //Pass 2: capture long tracks that are lower in pT, etc
  trackParams[1].NLayers = 12;
  trackParams[1].MinTrackLength = 6; //this is the one with fixed params

  trackParams[1].LayerRadii = LayerRadii;
  trackParams[1].LayerZ = LayerZ;
  trackParams[1].TrackletMaxDeltaPhi = lOpening*0.3;
  trackParams[1].CellMaxDeltaPhi = lOpening*0.15;
  trackParams[1].CellMaxDeltaTanLambda = lOpening*0.03;
  trackParams[1].TrackletMaxDeltaZ = SecondTrackletMaxDeltaZ;
  trackParams[1].CellMaxDCA = SecondCellMaxDCA;
  trackParams[1].CellMaxDeltaZ = SecondCellMaxDeltaZ;
  trackParams[1].NeighbourMaxDeltaCurvature = SecondNeighbourMaxDeltaCurvature;
  trackParams[1].NeighbourMaxDeltaN = SecondNeighbourMaxDeltaN;

  memParams[1].CellsMemoryCoefficients = CellsMemoryCoefficients;
  memParams[1].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
  memParams[1].MemoryOffset = 180000;

  //Pass 3: Secondaries
  trackParams[2].NLayers = 12;
  trackParams[2].MinTrackLength = 6; //this is the one with fixed params

  trackParams[2].LayerRadii = LayerRadii;
  trackParams[2].LayerZ = LayerZ;
  trackParams[2].TrackletMaxDeltaPhi = TMath::Min(lLastOpening*0.3, TMath::Pi());
  trackParams[2].CellMaxDeltaPhi = lLastOpening*0.15;
  trackParams[2].CellMaxDeltaTanLambda = lLastOpening*0.03;
  trackParams[2].TrackletMaxDeltaZ = ThirdTrackletMaxDeltaZ;
  trackParams[2].CellMaxDCA = ThirdCellMaxDCA;
  trackParams[2].CellMaxDeltaZ = ThirdCellMaxDeltaZ;
  trackParams[2].NeighbourMaxDeltaCurvature = ThirdNeighbourMaxDeltaCurvature;
  trackParams[2].NeighbourMaxDeltaN = ThirdNeighbourMaxDeltaN;

  memParams[2].CellsMemoryCoefficients = CellsMemoryCoefficients;
  memParams[2].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
  memParams[2].MemoryOffset = 180000;

  tracker.setParameters(memParams, trackParams);
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

  constexpr int nBins = 100;
  constexpr float minPt = 0.01;
  constexpr float maxPt = 10;
  double newBins[nBins + 1];
  newBins[0] = minPt;
  double factor = pow(maxPt / minPt, 1. / nBins);
  for (int i = 1; i <= nBins; i++) {
    newBins[i] = factor * newBins[i - 1];
  }

  //Histogram output
  //The ALICE 3 hits (cross-check)
  TH2D* hALICE3 = new TH2D("hALICE3", "", 2200, -110, 110, 2200, -110, 110);
  //Track distribution
  TH1D* hNTracks = new TH1D("hNTracks", "", 5000, 0, 5000);
  //Charged particle spectra
  TH1D* hPtSpectraAll = new TH1D("hPtSpectraAll", "", 500,0,50);
  TH1D* hPtSpectraFakeAll = new TH1D("hPtSpectraFakeAll", "", 500,0,50);
  TH1D* hPtSpectra = new TH1D("hPtSpectra", "", 500,0,5);
  TH1D* hPtSpectraFake = new TH1D("hPtSpectraFake", "", 500,0,5);
  
  TString lParticleName[5] = {"el", "mu", "pi", "ka", "pr"};
  int lParticlePDG[5] = {11, 13, 211, 321, 2212};
  
  TH1D* hPtSpectraRecoPID[5];
  TH1D* hPtSpectraGenPID[5];
  
  for(Int_t ii=0; ii<5; ii++){
    hPtSpectraRecoPID[ii]=new TH1D(Form("hPtSpectraRecoPID_%s", lParticleName[ii].Data()), "", 500,0,50);
    hPtSpectraGenPID[ii]=new TH1D(Form("hPtSpectraGenPID_%s", lParticleName[ii].Data()), "", 500,0,50);
  }
  
  TH1D *hHitDensityMidRap[12];
  for(Int_t ii=0; ii<12; ii++){
    hHitDensityMidRap[ii]=new TH1D(Form("hHitDensityMidRapLayer%i", ii), "", 400,0,2000*TMath::Power(LayerRadii[0],1)/TMath::Power(LayerRadii[ii],1));
  }
  
  TH1D *hPartialTransportDebug = new TH1D("hPartialTransportDebug", "", 6,0,6);
  hPartialTransportDebug->GetXaxis()->SetBinLabel(1, "All particles");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(2, "Not tranported");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(3, "Successfully smeared");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(4, "Fully tracked particles");
  
  TH2D *hPIDiTOF = new TH2D("hPIDiTOF", "", 500,0,5,110,0,1.1);
  TH2D *hPIDoTOF = new TH2D("hPIDoTOF", "", 500,0,5,110,0,1.1);
  
  
  o2::delphes::TrackSmearer smearer;
  cout<<"Loading LUTs for primary particles..."<<flush;
  if( lMagneticField > 10 ){
    smearer.loadTable(11,   Form("%slutCovm.2Tesla.20cm.el.dat",lLutPath.Data()));
    smearer.loadTable(13,   Form("%slutCovm.2Tesla.20cm.mu.dat",lLutPath.Data()));
    smearer.loadTable(211,  Form("%slutCovm.2Tesla.20cm.pi.dat",lLutPath.Data()));
    smearer.loadTable(321,  Form("%slutCovm.2Tesla.20cm.ka.dat",lLutPath.Data()));
    smearer.loadTable(2212, Form("%slutCovm.2Tesla.20cm.pr.dat",lLutPath.Data()));
  }else{
    smearer.loadTable(11,   Form("%slutCovm.el.geometry_v1_tof.rmin20cm.5kG.dat",lLutPath.Data()));
    smearer.loadTable(13,   Form("%slutCovm.mu.geometry_v1_tof.rmin20cm.5kG.dat",lLutPath.Data()));
    smearer.loadTable(211,  Form("%slutCovm.pi.geometry_v1_tof.rmin20cm.5kG.dat",lLutPath.Data()));
    smearer.loadTable(321,  Form("%slutCovm.ka.geometry_v1_tof.rmin20cm.5kG.dat",lLutPath.Data()));
    smearer.loadTable(2212, Form("%slutCovm.pr.geometry_v1_tof.rmin20cm.5kG.dat",lLutPath.Data()));
  }
  cout<<" done"<<endl;
  
  cout<<"Initializing track smearer..."<<flush;
  auto gen = new SmearO2KineGenerator( Form("%so2sim_Kine.root", lDataPath.Data()) );
  //gen->SetStartEvent(lEventNumber);
  gen->Init();
  cout<<" done"<<endl;
  
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  std::cout << "*- Setting up output file..." << std::endl;
  // Setup output
  UInt_t fCompress = 101;
  int fBasketSizeEvents = 1000000;                                                // Maximum basket size of the trees for events
  int fBasketSizeTracks = 10000000;                                               // Maximum basket size of the trees for tracks
  TFile* fOutputFile = TFile::Open(lOutputName.Data(), "RECREATE", "O2 AOD", fCompress); // File to store the trees of time frames

  TString lAlienProcId = gSystem -> Getenv( "ALIEN_PROC_ID" );
  cout<<"Alien process ID: "<<lAlienProcId.Data()<<endl;
  Long_t lProcID = lAlienProcId.Atoi();
  Long_t lProcIDLarge = lProcID/1000;
  lProcID = lProcID - 1000*lProcIDLarge;
  cout<<"Base parse: "<<lProcID<<endl;
  UInt_t tfId = lProcID;

  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  Long_t lGoodEvents = 0;
  Long_t fOffsetLabel = 0;
  Float_t lHitCount[12];
  Float_t lRingArea[12];
  
  //Ring with TMath::Abs(z)/r < 0.1 has area: 2*pi*r*(0.2*r)
  for(Int_t ilay=0; ilay<12; ilay++) lRingArea[ilay] = 0.4*TMath::Pi()*LayerRadii[ilay];
  for(Int_t ilay=0; ilay<12; ilay++) cout<<"Ring area "<<lRingArea[ilay]<<" for layer "<<ilay<<endl;
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // one Data Frame per event strategy
  // Create the output directory for the current time frame
  TDirectory* fOutputDir = 0x0; ///! Pointer to the output Root subdirectory
  fOutputDir = fOutputFile->mkdir(Form("DF_%d", tfId));
  fOutputDir->cd();

  //Create output trees in file
  TTree* fTree[kTrees];
  for (Int_t ii = 0; ii < kTrees; ii++) {
    if (gSaveTree[ii]) {
      std::cout << "*- Creating tree " << gTreeName[ii] << "..." << std::endl;
      fTree[ii] = new TTree(gTreeName[ii], gTreeTitle[ii]);
      fTree[ii]->SetAutoFlush(0);
    }
  }
  if (gSaveTree[kEvents]) {
    fTree[kEvents]->Branch("fIndexBCs", &collision.fIndexBCs, "fIndexBCs/I");
    fTree[kEvents]->Branch("fPosX", &collision.fPosX, "fPosX/F");
    fTree[kEvents]->Branch("fPosY", &collision.fPosY, "fPosY/F");
    fTree[kEvents]->Branch("fPosZ", &collision.fPosZ, "fPosZ/F");
    fTree[kEvents]->Branch("fCovXX", &collision.fCovXX, "fCovXX/F");
    fTree[kEvents]->Branch("fCovXY", &collision.fCovXY, "fCovXY/F");
    fTree[kEvents]->Branch("fCovXZ", &collision.fCovXZ, "fCovXZ/F");
    fTree[kEvents]->Branch("fCovYY", &collision.fCovYY, "fCovYY/F");
    fTree[kEvents]->Branch("fCovYZ", &collision.fCovYZ, "fCovYZ/F");
    fTree[kEvents]->Branch("fCovZZ", &collision.fCovZZ, "fCovZZ/F");
    fTree[kEvents]->Branch("fFlags", &collision.fFlags, "fFlags/s");
    fTree[kEvents]->Branch("fChi2", &collision.fChi2, "fChi2/F");
    fTree[kEvents]->Branch("fNumContrib", &collision.fN, "fNumContrib/i");
    fTree[kEvents]->Branch("fCollisionTime", &collision.fCollisionTime, "fCollisionTime/F");
    fTree[kEvents]->Branch("fCollisionTimeRes", &collision.fCollisionTimeRes, "fCollisionTimeRes/F");
    fTree[kEvents]->Branch("fCollisionTimeMask", &collision.fCollisionTimeMask, "fCollisionTimeMask/b");
    fTree[kEvents]->SetBasketSize("*", fBasketSizeEvents);
  }

  if (gSaveTree[kTracks]) {
    fTree[kTracks]->Branch("fIndexCollisions", &tracks.fIndexCollisions, "fIndexCollisions/I");
    fTree[kTracks]->Branch("fTrackType", &tracks.fTrackType, "fTrackType/b");
    //    fTree[kTracks]->Branch("fTOFclsIndex", &tracks.fTOFclsIndex, "fTOFclsIndex/I");
    //    fTree[kTracks]->Branch("fNTOFcls", &tracks.fNTOFcls, "fNTOFcls/I");
    fTree[kTracks]->Branch("fX", &tracks.fX, "fX/F");
    fTree[kTracks]->Branch("fAlpha", &tracks.fAlpha, "fAlpha/F");
    fTree[kTracks]->Branch("fY", &tracks.fY, "fY/F");
    fTree[kTracks]->Branch("fZ", &tracks.fZ, "fZ/F");
    fTree[kTracks]->Branch("fSnp", &tracks.fSnp, "fSnp/F");
    fTree[kTracks]->Branch("fTgl", &tracks.fTgl, "fTgl/F");
    fTree[kTracks]->Branch("fSigned1Pt", &tracks.fSigned1Pt, "fSigned1Pt/F");
  }
    
  if (gSaveTree[kTracksCov]) {
    // Modified covariance matrix
    fTree[kTracksCov]->Branch("fSigmaY", &tracks.fSigmaY, "fSigmaY/F");
    fTree[kTracksCov]->Branch("fSigmaZ", &tracks.fSigmaZ, "fSigmaZ/F");
    fTree[kTracksCov]->Branch("fSigmaSnp", &tracks.fSigmaSnp, "fSigmaSnp/F");
    fTree[kTracksCov]->Branch("fSigmaTgl", &tracks.fSigmaTgl, "fSigmaTgl/F");
    fTree[kTracksCov]->Branch("fSigma1Pt", &tracks.fSigma1Pt, "fSigma1Pt/F");
    fTree[kTracksCov]->Branch("fRhoZY", &tracks.fRhoZY, "fRhoZY/B");
    fTree[kTracksCov]->Branch("fRhoSnpY", &tracks.fRhoSnpY, "fRhoSnpY/B");
    fTree[kTracksCov]->Branch("fRhoSnpZ", &tracks.fRhoSnpZ, "fRhoSnpZ/B");
    fTree[kTracksCov]->Branch("fRhoTglY", &tracks.fRhoTglY, "fRhoTglY/B");
    fTree[kTracksCov]->Branch("fRhoTglZ", &tracks.fRhoTglZ, "fRhoTglZ/B");
    fTree[kTracksCov]->Branch("fRhoTglSnp", &tracks.fRhoTglSnp, "fRhoTglSnp/B");
    fTree[kTracksCov]->Branch("fRho1PtY", &tracks.fRho1PtY, "fRho1PtY/B");
    fTree[kTracksCov]->Branch("fRho1PtZ", &tracks.fRho1PtZ, "fRho1PtZ/B");
    fTree[kTracksCov]->Branch("fRho1PtSnp", &tracks.fRho1PtSnp, "fRho1PtSnp/B");
    fTree[kTracksCov]->Branch("fRho1PtTgl", &tracks.fRho1PtTgl, "fRho1PtTgl/B");
  }
  
  if (gSaveTree[kTracksExtra]) {
    //
    fTree[kTracksExtra]->Branch("fTPCInnerParam", &tracks.fTPCinnerP, "fTPCInnerParam/F");
    fTree[kTracksExtra]->Branch("fFlags", &tracks.fFlags, "fFlags/i");
    fTree[kTracksExtra]->Branch("fITSClusterMap", &tracks.fITSClusterMap, "fITSClusterMap/b");
    fTree[kTracksExtra]->Branch("fTPCNClsFindable", &tracks.fTPCNClsFindable, "fTPCNClsFindable/b");
    fTree[kTracksExtra]->Branch("fTPCNClsFindableMinusFound", &tracks.fTPCNClsFindableMinusFound, "fTPCNClsFindableMinusFound/B");
    fTree[kTracksExtra]->Branch("fTPCNClsFindableMinusCrossedRows", &tracks.fTPCNClsFindableMinusCrossedRows, "fTPCNClsFindableMinusCrossedRows/B");
    fTree[kTracksExtra]->Branch("fTPCNClsShared", &tracks.fTPCNClsShared, "fTPCNClsShared/b");
    fTree[kTracksExtra]->Branch("fTRDPattern", &tracks.fTRDPattern, "fTRDPattern/b");
    fTree[kTracksExtra]->Branch("fITSChi2NCl", &tracks.fITSChi2NCl, "fITSChi2NCl/F");
    fTree[kTracksExtra]->Branch("fTPCChi2NCl", &tracks.fTPCChi2NCl, "fTPCChi2NCl/F");
    fTree[kTracksExtra]->Branch("fTRDChi2", &tracks.fTRDChi2, "fTRDChi2/F");
    fTree[kTracksExtra]->Branch("fTOFChi2", &tracks.fTOFChi2, "fTOFChi2/F");
    fTree[kTracksExtra]->Branch("fTPCSignal", &tracks.fTPCSignal, "fTPCSignal/F");
    fTree[kTracksExtra]->Branch("fTRDSignal", &tracks.fTRDSignal, "fTRDSignal/F");
    fTree[kTracksExtra]->Branch("fTOFSignal", &tracks.fTOFSignal, "fTOFSignal/F");
    fTree[kTracksExtra]->Branch("fLength", &tracks.fLength, "fLength/F");
    fTree[kTracksExtra]->Branch("fTOFExpMom", &tracks.fTOFExpMom, "fTOFExpMom/F");
    fTree[kTracksExtra]->Branch("fTrackEtaEMCAL", &tracks.fTrackEtaEMCAL, "fTrackEtaEMCAL/F");
    fTree[kTracksExtra]->Branch("fTrackPhiEMCAL", &tracks.fTrackPhiEMCAL, "fTrackPhiEMCAL/F");
    fTree[kTracksExtra]->SetBasketSize("*", fBasketSizeTracks);
  }

  
  if (gSaveTree[kMcTrackLabel]) {
    //fTree[kMcTrackLabel]->Branch("fIndexMcParticles", &mctracklabel.fLabel, "fIndexMcParticles/i");
    //fTree[kMcTrackLabel]->Branch("fMcMask", &mctracklabel.fLabelMask, "fMcMask/s");
    fTree[kMcTrackLabel]->Branch("fIndexMcParticles", &mctracklabel.fIndexMcParticles, "fIndexMcParticles/I");
    fTree[kMcTrackLabel]->Branch("fMcMask", &mctracklabel.fMcMask, "fMcMask/s");
    fTree[kMcTrackLabel]->SetBasketSize("*", fBasketSizeTracks);
  }

  if (gSaveTree[kMcCollision]) {
    fTree[kMcCollision]->Branch("fIndexBCs", &mccollision.fIndexBCs, "fIndexBCs/I");
    fTree[kMcCollision]->Branch("fGeneratorsID", &mccollision.fGeneratorsID, "fGeneratorsID/S");
    fTree[kMcCollision]->Branch("fPosX", &mccollision.fPosX, "fPosX/F");
    fTree[kMcCollision]->Branch("fPosY", &mccollision.fPosY, "fPosY/F");
    fTree[kMcCollision]->Branch("fPosZ", &mccollision.fPosZ, "fPosZ/F");
    fTree[kMcCollision]->Branch("fT", &mccollision.fT, "fT/F");
    fTree[kMcCollision]->Branch("fWeight", &mccollision.fWeight, "fWeight/F");
    fTree[kMcCollision]->Branch("fImpactParameter", &mccollision.fImpactParameter, "fImpactParameter/F");
    fTree[kMcCollision]->SetBasketSize("*", fBasketSizeEvents);
  }

  if (gSaveTree[kMcParticle]) {
    fTree[kMcParticle]->Branch("fIndexMcCollisions", &mcparticle.fIndexMcCollisions, "fIndexMcCollisions/I");
    fTree[kMcParticle]->Branch("fPdgCode", &mcparticle.fPdgCode, "fPdgCode/I");
    fTree[kMcParticle]->Branch("fStatusCode", &mcparticle.fStatusCode, "fStatusCode/I");
    fTree[kMcParticle]->Branch("fFlags", &mcparticle.fFlags, "fFlags/b");
    fTree[kMcParticle]->Branch("fIndexMcParticles_Mother0", &mcparticle.fIndexMcParticles_Mother0, "fIndexMcParticles_Mother0/I");
    fTree[kMcParticle]->Branch("fIndexMcParticles_Mother1", &mcparticle.fIndexMcParticles_Mother1, "fIndexMcParticles_Mother1/I");
    fTree[kMcParticle]->Branch("fIndexMcParticles_Daughter0", &mcparticle.fIndexMcParticles_Daughter0, "fIndexMcParticles_Daughter0/I");
    fTree[kMcParticle]->Branch("fIndexMcParticles_Daughter1", &mcparticle.fIndexMcParticles_Daughter1, "fIndexMcParticles_Daughter1/I");
    fTree[kMcParticle]->Branch("fWeight", &mcparticle.fWeight, "fWeight/F");

    fTree[kMcParticle]->Branch("fPx", &mcparticle.fPx, "fPx/F");
    fTree[kMcParticle]->Branch("fPy", &mcparticle.fPy, "fPy/F");
    fTree[kMcParticle]->Branch("fPz", &mcparticle.fPz, "fPz/F");
    fTree[kMcParticle]->Branch("fE", &mcparticle.fE, "fE/F");

    fTree[kMcParticle]->Branch("fVx", &mcparticle.fVx, "fVx/F");
    fTree[kMcParticle]->Branch("fVy", &mcparticle.fVy, "fVy/F");
    fTree[kMcParticle]->Branch("fVz", &mcparticle.fVz, "fVz/F");
    fTree[kMcParticle]->Branch("fVt", &mcparticle.fVt, "fVt/F");
    fTree[kMcParticle]->SetBasketSize("*", fBasketSizeTracks);
  }

  if (gSaveTree[kMcCollisionLabel]) {
    fTree[kMcCollisionLabel]->Branch("fIndexMcCollisions", &mccollisionlabel.fIndexMcCollisions, "fIndexMcCollisions/I");
    fTree[kMcCollisionLabel]->Branch("fMcMask", &mccollisionlabel.fMcMask, "fMcMask/s");
    fTree[kMcCollisionLabel]->SetBasketSize("*", fBasketSizeEvents);
  }

  if (gSaveTree[kBC]) {
    fTree[kBC]->Branch("fRunNumber", &bc.fRunNumber, "fRunNumber/I");
    fTree[kBC]->Branch("fGlobalBC", &bc.fGlobalBC, "fGlobalBC/l");
    fTree[kBC]->Branch("fTriggerMask", &bc.fTriggerMask, "fTriggerMask/l");
    fTree[kBC]->SetBasketSize("*", fBasketSizeEvents);
  }

  
  // Associate branches for Run 2 BC info
  if (gSaveTree[kRun2BCInfo]) {
    fTree[kRun2BCInfo]->Branch("fEventCuts", &run2bcinfo.fEventCuts, "fEventCuts/i");
    fTree[kRun2BCInfo]->Branch("fTriggerMaskNext50", &run2bcinfo.fTriggerMaskNext50, "fTriggerMaskNext50/l");
    fTree[kRun2BCInfo]->Branch("fL0TriggerInputMask", &run2bcinfo.fL0TriggerInputMask, "fL0TriggerInputMask/i");
    fTree[kRun2BCInfo]->Branch("fSPDClustersL0", &run2bcinfo.fSPDClustersL0, "fSPDClustersL0/s");
    fTree[kRun2BCInfo]->Branch("fSPDClustersL1", &run2bcinfo.fSPDClustersL1, "fSPDClustersL1/s");
    fTree[kRun2BCInfo]->Branch("fSPDFiredChipsL0", &run2bcinfo.fSPDFiredChipsL0, "fSPDFiredChipsL0/s");
    fTree[kRun2BCInfo]->Branch("fSPDFiredChipsL1", &run2bcinfo.fSPDFiredChipsL1, "fSPDFiredChipsL1/s");
    fTree[kRun2BCInfo]->Branch("fSPDFiredFastOrL0", &run2bcinfo.fSPDFiredFastOrL0, "fSPDFiredFastOrL0/s");
    fTree[kRun2BCInfo]->Branch("fSPDFiredFastOrL1", &run2bcinfo.fSPDFiredFastOrL1, "fSPDFiredFastOrL1/s");
    fTree[kRun2BCInfo]->Branch("fV0TriggerChargeA", &run2bcinfo.fV0TriggerChargeA, "fV0TriggerChargeA/s");
    fTree[kRun2BCInfo]->Branch("fV0TriggerChargeC", &run2bcinfo.fV0TriggerChargeC, "fV0TriggerChargeC/s");
    fTree[kRun2BCInfo]->SetBasketSize("*", fBasketSizeEvents);
  }
  
  if (gSaveTree[kFDD]) {
    fTree[kFDD]->Branch("fIndexBCs", &fdd.fIndexBCs, "fIndexBCs/I");
    fTree[kFDD]->Branch("fAmplitudeA", fdd.fAmplitudeA, "fAmplitudeA[4]/F");
    fTree[kFDD]->Branch("fAmplitudeC", fdd.fAmplitudeC, "fAmplitudeC[4]/F");
    fTree[kFDD]->Branch("fTimeA", &fdd.fTimeA, "fTimeA/F");
    fTree[kFDD]->Branch("fTimeC", &fdd.fTimeC, "fTimeC/F");
    fTree[kFDD]->Branch("fTriggerMask", &fdd.fTriggerMask, "fTriggerMask/b");
    fTree[kFDD]->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for V0A
  if (gSaveTree[kFV0A]) {
    fTree[kFV0A]->Branch("fIndexBCs", &fv0a.fIndexBCs, "fIndexBCs/I");
    fTree[kFV0A]->Branch("fAmplitude", fv0a.fAmplitude, "fAmplitude[48]/F");
    fTree[kFV0A]->Branch("fTime", &fv0a.fTime, "fTime/F");
    fTree[kFV0A]->Branch("fTriggerMask", &fv0a.fTriggerMask, "fTriggerMask/b");
    fTree[kFV0A]->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for V0C
  if (gSaveTree[kFV0C]) {
    fTree[kFV0C]->Branch("fIndexBCs", &fv0c.fIndexBCs, "fIndexBCs/I");
    fTree[kFV0C]->Branch("fAmplitude", fv0c.fAmplitude, "fAmplitude[32]/F");
    fTree[kFV0C]->Branch("fTime", &fv0c.fTime, "fTime/F");
    fTree[kFV0C]->SetBasketSize("*", fBasketSizeEvents);
  }

  // Associate branches for FT0
  if (gSaveTree[kFT0]) {
    fTree[kFT0]->Branch("fIndexBCs", &ft0.fIndexBCs, "fIndexBCs/I");
    fTree[kFT0]->Branch("fAmplitudeA", ft0.fAmplitudeA, "fAmplitudeA[96]/F");
    fTree[kFT0]->Branch("fAmplitudeC", ft0.fAmplitudeC, "fAmplitudeC[112]/F");
    fTree[kFT0]->Branch("fTimeA", &ft0.fTimeA, "fTimeA/F");
    fTree[kFT0]->Branch("fTimeC", &ft0.fTimeC, "fTimeC/F");
    fTree[kFT0]->Branch("fTriggerMask", &ft0.fTriggerMask, "fTriggerMask/b");
    fTree[kFT0]->SetBasketSize("*", fBasketSizeEvents);
  }

  if (gSaveTree[kZdc]) {
    fTree[kZdc]->Branch("fIndexBCs", &zdc.fIndexBCs, "fIndexBCs/I");
    fTree[kZdc]->Branch("fEnergyZEM1", &zdc.fEnergyZEM1, "fEnergyZEM1/F");
    fTree[kZdc]->Branch("fEnergyZEM2", &zdc.fEnergyZEM2, "fEnergyZEM2/F");
    fTree[kZdc]->Branch("fEnergyCommonZNA", &zdc.fEnergyCommonZNA, "fEnergyCommonZNA/F");
    fTree[kZdc]->Branch("fEnergyCommonZNC", &zdc.fEnergyCommonZNC, "fEnergyCommonZNC/F");
    fTree[kZdc]->Branch("fEnergyCommonZPA", &zdc.fEnergyCommonZPA, "fEnergyCommonZPA/F");
    fTree[kZdc]->Branch("fEnergyCommonZPC", &zdc.fEnergyCommonZPC, "fEnergyCommonZPC/F");
    fTree[kZdc]->Branch("fEnergySectorZNA", &zdc.fEnergySectorZNA, "fEnergySectorZNA[4]/F");
    fTree[kZdc]->Branch("fEnergySectorZNC", &zdc.fEnergySectorZNC, "fEnergySectorZNC[4]/F");
    fTree[kZdc]->Branch("fEnergySectorZPA", &zdc.fEnergySectorZPA, "fEnergySectorZPA[4]/F");
    fTree[kZdc]->Branch("fEnergySectorZPC", &zdc.fEnergySectorZPC, "fEnergySectorZPC[4]/F");
    fTree[kZdc]->Branch("fTimeZEM1", &zdc.fTimeZEM1, "fTimeZEM1/F");
    fTree[kZdc]->Branch("fTimeZEM2", &zdc.fTimeZEM2, "fTimeZEM2/F");
    fTree[kZdc]->Branch("fTimeZNA", &zdc.fTimeZNA, "fTimeZNA/F");
    fTree[kZdc]->Branch("fTimeZNC", &zdc.fTimeZNC, "fTimeZNC/F");
    fTree[kZdc]->Branch("fTimeZPA", &zdc.fTimeZPA, "fTimeZPA/F");
    fTree[kZdc]->Branch("fTimeZPC", &zdc.fTimeZPC, "fTimeZPC/F");
    fTree[kZdc]->SetBasketSize("*", fBasketSizeEvents);
  }
  
  std::cout << "*- Number of events detected: " << itsHits.GetEntries() << std::endl;
  for (int iEvent{0}; iEvent < itsHits.GetEntriesFast(); ++iEvent) {
    itsHits.GetEntry(iEvent);
    mcTree->GetEvent(iEvent);
    o2::its::ROframe event{iEvent, 12};
    std::cout << "*- Processing event " << iEvent << "..." << std::endl;
  
    int id{0};
    for(Int_t ilay=0; ilay<12; ilay++){
      lHitCount[ilay]=0;
    }
    
    Long_t lNchThisEvent = 0;
    Long_t lNchThisEventAll = 0;
    
    for (Long_t iii=0; iii< mcArr->size(); iii++ ){
      auto part = mcArr->at(iii);
      
      //Count primary particles
      if(TMath::Abs(part.GetEta())<0.5 && part.isPrimary() &&
         (
         TMath::Abs( part.GetPdgCode() ) == 211 ||
         TMath::Abs( part.GetPdgCode() ) == 321 ||
         TMath::Abs( part.GetPdgCode() ) == 2212 ||
         TMath::Abs( part.GetPdgCode() ) == 13 ||
         TMath::Abs( part.GetPdgCode() ) == 11
          )
         )
        lNchThisEvent++;
      
      //Count primary particles
      if(TMath::Abs(part.GetEta())<1.5 && part.isPrimary() &&
         (
         TMath::Abs( part.GetPdgCode() ) == 211 ||
         TMath::Abs( part.GetPdgCode() ) == 321 ||
         TMath::Abs( part.GetPdgCode() ) == 2212 ||
         TMath::Abs( part.GetPdgCode() ) == 13 ||
         TMath::Abs( part.GetPdgCode() ) == 11
          )
         )
        lNchThisEventAll++;
    }
    
    for (auto& hit : *hits) {
      int layer{hit.GetDetectorID()};
      float xyz[3]{hit.GetX(), hit.GetY(), hit.GetZ()};
      float r{std::hypot(xyz[0], xyz[1])};
      float phi{std::atan2(-xyz[1], -xyz[0]) + o2::its::constants::math::Pi};

      float lSmearing = lIBresolution;
      if(layer>3) lSmearing = lOBresolution;
      
      if (kUseSmearing) {
        phi = gRandom->Gaus(phi, std::asin(lSmearing / r));
        xyz[0] = r * std::cos(phi);
        xyz[1] = r * std::sin(phi);
        xyz[2] = gRandom->Gaus(xyz[2], lSmearing);
      }
      
      hALICE3->Fill(xyz[0], xyz[1]);

      //if you see radius + epsilon, it's still the N-th layer... likely a bug
      if (r > 99.0 && r < 101) {
        //std::cout << "*- Exception caught at a radius of "<< r << std::endl;
        layer = 11;
      }
      
      if(TMath::Abs(xyz[2])/r<0.1) lHitCount[layer]+=1/lRingArea[layer];
      
      event.addTrackingFrameInfoToLayer(layer, xyz[0], xyz[1], xyz[2], r, phi, std::array<float, 2>{0.f, xyz[2]},
                                        std::array<float, 3>{lSmearing * lSmearing, 0.f, lSmearing * lSmearing});
      event.addClusterToLayer(layer, xyz[0], xyz[1], xyz[2], event.getClustersOnLayer(layer).size());
      event.addClusterLabelToLayer(layer, o2::MCCompLabel(hit.GetTrackID(), iEvent, iEvent, false));
      event.addClusterExternalIndexToLayer(layer, id++);
    }
    for(Int_t ilay=0; ilay<12; ilay++) hHitDensityMidRap[ilay] -> Fill(lHitCount[ilay]);
    roFrame = iEvent;
    std::cout << "*- Event " << iEvent << " finished adding hits." << std::endl;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //LUT APPROACH STEP 1: DO ALL PRIMARY TRACKS
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    o2::math_utils::Point3D<float> zeropos{0.0,0.0,0.0};
    std::array<float, 6> zerocov;
    for(Int_t jj=0; jj<6; jj++) zerocov[jj]=1e-6;
    o2::dataformats::VertexBase zerovtx(zeropos, zerocov);
    o2::dataformats::DCA dca;
    
    //For vertexing
    std::vector<TrackAlice3> smearedTracks;
    std::vector<o2::InteractionRecord> bcData;
    
    //clear old stuff
    gen->clearParticles();
    
    // loop over events
    if(!gen->importParticles()){
      cout<<"PROBLEM RETRIEVING EVENT "<< iEvent<<endl;
    }
    auto particles = gen->getParticles();
    //lWatchSmearing.Start(kFALSE);
    Long_t lParticlesForVertexing = 0;
    // loop over particles
    cout<<"Allocating TOF info arrays..."<<flush;
    TArrayF lInnerTOFSignalLutted(particles.size());
    TArrayF lOuterTOFSignalLutted(particles.size());
    cout<<" done"<<endl;
    
    for (int iparticle = 0; iparticle < particles.size(); ++iparticle){
      //Inner and outer TOF signal initialization
      lInnerTOFSignalLutted[iparticle] = -100;
      lOuterTOFSignalLutted[iparticle] = -100;
      
      auto particle = particles[iparticle];
      
      if( TMath::Abs(particle.Eta())> 1.5) continue; //skip outside barrel
      hPartialTransportDebug->Fill(0.5);
      // only particles to be transported, which are flagged as status code = 1
      // the particles that have been transported already have status code = 0
      if (particle.GetStatusCode() != 1) continue;
      hPartialTransportDebug->Fill(1.5);
      
      // only particles that we know how to smear
      // namely el, mu, pi, ka, pr
      auto pdg = std::abs(particle.GetPdgCode());
      if (pdg != 11 && pdg != 13 && pdg != 211 && pdg != 321 && pdg != 2212) continue;
      
      if( pdg == 3312 || pdg == 3334 ) {
        cout<<"PROBLEM IN EVENT "<< iEvent<<": there's a spurious xi or omega!!! "<<endl;
      }
      
      // convert particle to o2 track and smear it
      O2Track o2track;
      o2::delphes::TrackUtils::convertTParticleToO2Track(particle, o2track);
      
      Float_t lThisTrackLength20 = -100;
      Float_t lThisTrackLength100 = -100;
      
      //Calculate the perfect expected information
      Float_t lX0=-100, lX1=-100;
      if (!o2track.propagateToDCA(zerovtx, tracker.getBz())) {
        //std::cout << "Propagation failed." << std::endl;
      } else {
        lX0 = o2track.getX();
      }
      if (!o2track.getXatLabR(100.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
        lX1 = -100;
      }
      if(lX0>-99.&&lX1>-99.) lThisTrackLength100 = TrackLength(o2track, lX0, lX1);
      o2track.propagateTo(lX0, lMagneticField);
      
      lX0=-100, lX1=-100;
      lThisTrackLength20 = -100;
      if (!o2track.propagateToDCA(zerovtx, tracker.getBz())) {
        //std::cout << "Propagation failed." << std::endl;
      } else {
        lX0 = o2track.getX();
      }
      if (!o2track.getXatLabR(20.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
        lX1 = -100;
      }
      if(lX0>-99.&&lX1>-99.) lThisTrackLength20 = TrackLength(o2track, lX0, lX1);
      o2track.propagateTo(lX0, lMagneticField);
      
      //Calculate expected times with smearing
      Double_t lExpectedTimeInformation20 = lThisTrackLength20/ Velocity(o2track.getP(), getMass(pdg) );
      Double_t lExpectedTimeInformation100 = lThisTrackLength100/ Velocity(o2track.getP(), getMass(pdg) );
      
      if(lThisTrackLength20>-1)
        lInnerTOFSignalLutted[iparticle] = gRandom->Gaus(lExpectedTimeInformation20, loTOFresolution);
      if(lThisTrackLength100>-1)
        lOuterTOFSignalLutted[iparticle] = gRandom->Gaus(lExpectedTimeInformation100, loTOFresolution);
      
      Int_t alabel = iparticle;
      
      float nch = lNchThisEvent;
      if (!smearer.smearTrack(o2track, pdg, nch)) continue;
      hPartialTransportDebug->Fill(2.5);
      
      //Feed vertex finder
      o2::InteractionRecord ir(iEvent, 0);
      const float t = (ir.bc2ns() + gRandom->Gaus(0., 100.)) * 1e-3;
      smearedTracks.push_back(TrackAlice3{o2track, t, 100.f * 1e-3, TMath::Abs(alabel)});
      lParticlesForVertexing++;
    }
    cout<<"Tracks to be used for vertex finding: "<<lParticlesForVertexing<<endl;
    cout<<"PV finding procedure started"<<endl;

//    //True PV finding sequence - using Ruben's PVertexer
//    o2::BunchFilling bcfill;
//    bcfill.setDefault();
//    o2::vertexing::PVertexer vertexer;
//    vertexer.setValidateWithIR(kFALSE);
//    vertexer.setBunchFilling(bcfill);
//    vertexer.init();
//
//    std::vector<o2::MCCompLabel> lblTracks;
//    std::vector<o2::vertexing::PVertex> vertices;
//    std::vector<o2::vertexing::GIndex> vertexTrackIDs;
//    std::vector<o2::vertexing::V2TRef> v2tRefs;
//    std::vector<o2::MCEventLabel> lblVtx;
//    lblVtx.emplace_back(iEvent, 1);
//    std::vector<o2::dataformats::GlobalTrackID> idxVec; // here we will the global IDs of all used tracks
//    idxVec.reserve(smearedTracks.size());
//    for (unsigned i = 0; i < smearedTracks.size(); i++) {
//      lblTracks.emplace_back(smearedTracks[i].mLabel, iEvent, 1, false);
//      idxVec.emplace_back(i, o2::dataformats::GlobalTrackID::ITS);
//    }
//
//    const int n_vertices = vertexer.process(smearedTracks,
//                                            idxVec,
//                                            gsl::span<o2::InteractionRecord>{bcData},
//                                            vertices,
//                                            vertexTrackIDs,
//                                            v2tRefs,
//                                            gsl::span<const o2::MCCompLabel>{lblTracks},
//                                            lblVtx);
//
//    cout<<"PV finding procedure finished, vertices found: "<<n_vertices<<endl;
//
//    if(n_vertices<1) continue;
    
    std::vector<o2::vertexing::PVertex> vertices;

    o2::vertexing::PVertex vertex1;
    vertex1.setX(0.0);
    vertex1.setY(0.0);
    vertex1.setZ(0.0);
    vertex1.setSigmaX2(1e-8);
    vertex1.setSigmaY2(1e-8);
    vertex1.setSigmaZ2(1e-8);
    vertex1.setSigmaXY(1e-8);
    vertex1.setSigmaXZ(1e-8);
    vertex1.setSigmaYZ(1e-8);
    vertices.emplace_back(vertex1);
    
    cout<<"Main PV at: "<<vertices[0].getX()<<", "<<vertices[0].getY()<<", "<<vertices[0].getZ()<<endl;
    
    //hVertexZ->Fill( vertices[0].getZ() );
    
    std::cout << "*- Number of vertices found: " << vertices.size() << endl;
    //hNVertices->Fill(vertices.size());
    //Skip events with no vertex
    if (vertices.size() == 0) {
      std::cout << "*- No primary vertex found, skipping event" << std::endl;
      continue;
    }
    
    //---> MC collision data
    mccollision.fIndexBCs = lGoodEvents;
    if (!mcHead) {
      std::cout << "*- Problem with MC header! " << std::endl;
      return;
    }
    cout<<"*- MC PV at "<<mcHead->GetX()<<", "<<mcHead->GetY()<<", "<<mcHead->GetZ()<<endl;
    
    mccollision.fPosX = mcHead->GetX();
    mccollision.fPosY = mcHead->GetY();
    mccollision.fPosZ = mcHead->GetZ();
    mccollision.fT = mcHead->GetT();
    mccollision.fWeight = 1;
    mccollision.fImpactParameter = mcHead->GetB();

    mccollisionlabel.fIndexMcCollisions = lGoodEvents;
    mccollisionlabel.fMcMask = 0;
    //---> Save fake hits on i-th layer for track

    o2::math_utils::Point3D<float> pos{vertices[0].getX(), vertices[0].getY(), vertices[0].getZ()};
    float_t pvpos[3] = {vertices[0].getX(), vertices[0].getY(), vertices[0].getZ()};
    std::array<float, 6> cov;
    for (Int_t jj = 0; jj < 6; jj++)
      cov[jj] = vertices[0].getCov()[jj];
    o2::dataformats::VertexBase vtx(pos, cov);
    //o2::dataformats::DCA dca;

    //Add primary vertex by hand!
    event.addPrimaryVertex(pvpos[0], pvpos[1], pvpos[2]) ;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    std::cout << "*- Acquiring collision information..." << std::endl;
    //---> Collision data
    Long_t lEventNumber = lGoodEvents;
    collision.fIndexBCs = fTree[kEvents]->GetEntries();
    fdd.fIndexBCs = collision.fIndexBCs;
    ft0.fIndexBCs = collision.fIndexBCs;
    fv0a.fIndexBCs = collision.fIndexBCs;
    fv0c.fIndexBCs = collision.fIndexBCs;
    zdc.fIndexBCs = collision.fIndexBCs;
    collision.fPosX = pvpos[0];
    collision.fPosY = pvpos[1];
    collision.fPosZ = pvpos[2];
    collision.fCovXX = cov[0];
    collision.fCovXY = cov[1];
    collision.fCovXZ = cov[2];
    collision.fCovYY = cov[3];
    collision.fCovYZ = cov[4];
    collision.fCovZZ = cov[5];
    collision.fChi2 = 1; //vertices[0].getChi2();
    collision.fN = -1; //vertices[0].getNContributors();
    collision.fCollisionTime = 10;
    collision.fCollisionTimeRes = 1e-6;
    ft0.fTimeA = 10;
    ft0.fTimeC = 10;

    //---> Dummy trigger mask to ensure nobody rejects this
    bc.fTriggerMask = 0;
    for (Int_t iii = 0; iii < 60; iii++)
      bc.fTriggerMask |= 1ull << iii;
    bc.fRunNumber = 246087; //ah, the good old days

    std::cout << "*- Event " << iEvent << " tracking" << std::endl;
    tracker.clustersToTracks(event);
    auto& lTracks = tracker.getTracks();
    auto& lTracksLabels = tracker.getTrackLabels();
    std::cout << "*- Event " << iEvent << " done tracking, Ntracks = "<<lTracks.size()<<" size of gen "<<particles.size() << std::endl;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // Now merge track lists and label lists
    cout<<"Initial track array size: "<<lTracks.size()<<endl;
    cout<<"Initial label array size: "<<lTracksLabels.size()<<endl;
    hPartialTransportDebug->Fill(3.5, lTracks.size());
    
    for(Int_t iTrack = 0; iTrack<smearedTracks.size(); iTrack++){
      o2::its::TrackITS lConvertedTrack = o2::its::TrackITS( smearedTracks[iTrack] );
      lTracks.emplace_back(lConvertedTrack);
      lTracksLabels.emplace_back(smearedTracks[iTrack].mLabel, iEvent, 1, false);
    }
    cout<<"Merged track array size: "<<lTracks.size()<<endl;
    cout<<"Merged label array size: "<<lTracksLabels.size()<<endl;
    
    hNTracks->Fill(lTracks.size());

    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //---> Track data
    Long_t lNTracks = lTracks.size();
    for (Int_t i = 0; i < lNTracks; i++) {
      //get result from tracker
      auto& lab = lTracksLabels[i];
      auto& track = lTracks[i];
      int trackID = std::abs(lab.getTrackID());
      
      //Propagate to primary vertex as usual
      o2::dataformats::DCA dca1;
      if (!track.propagateToDCA(vtx, tracker.getBz(), &dca1)) {
        std::cout << "Track propagation to primary vertex failed." << std::endl;
      }

      //Fill QA histograms
      hPtSpectraAll->Fill(track.getPt());
      if (lab.isFake())
        hPtSpectraFakeAll->Fill(track.getPt());
      
      auto part = mcArr->at(trackID);
      Float_t lMCEta = part.GetEta();
      Float_t lMCRap = part.GetRapidity();
      Float_t lMCPt = part.GetPt();
      Int_t lMCPDG = part.GetPdgCode();
      if(part.isPrimary()){
      //if (MC::isPhysicalPrimary(mcArr, part)) {
        if(TMath::Abs(lMCEta)<1.0) hPtSpectra->Fill(lMCPt);
        if(TMath::Abs(lMCEta)<1.0 && lab.isFake()) hPtSpectraFake->Fill(lMCPt);
        if(TMath::Abs(lMCRap)<1.0){
          for(Int_t ipdg=0; ipdg<5; ipdg++){
            if(lMCPDG==lParticlePDG[ipdg]) hPtSpectraRecoPID[ipdg]->Fill(lMCPt);
          }
        }
      }

      tracks.fIndexCollisions = lEventNumber;
      tracks.fTrackType = o2::aod::track::TrackTypeEnum::Run2Track; //Make track selection happy, please
      tracks.fFlags = 0x0;
      //Assume it all worked, fool regular selections
      tracks.fFlags |= o2::aod::track::TrackFlagsRun2Enum::ITSrefit;
      tracks.fFlags |= o2::aod::track::TrackFlagsRun2Enum::TPCrefit;
      tracks.fFlags |= o2::aod::track::TrackFlagsRun2Enum::GoldenChi2;

      //Main: X, alpha, track params
      tracks.fX = track.getX();
      tracks.fY = track.getY();
      tracks.fZ = track.getZ();
      tracks.fAlpha = track.getAlpha();
      tracks.fSnp = track.getSnp();
      tracks.fTgl = track.getTgl();
      tracks.fSigned1Pt = track.getQ2Pt();

      // diagonal elements of covariance matrix
      tracks.fSigmaY = TMath::Sqrt(track.getSigmaY2());
      tracks.fSigmaZ = TMath::Sqrt(track.getSigmaZ2());
      tracks.fSigmaSnp = TMath::Sqrt(track.getSigmaSnp2());
      tracks.fSigmaTgl = TMath::Sqrt(track.getSigmaTgl2());
      tracks.fSigma1Pt = TMath::Sqrt(track.getSigma1Pt2());
      // off-diagonal elements of covariance matrix
      tracks.fRhoZY = (Char_t)(128. * track.getSigmaZY() / tracks.fSigmaZ / tracks.fSigmaY);
      tracks.fRhoSnpY = (Char_t)(128. * track.getSigmaSnpY() / tracks.fSigmaSnp / tracks.fSigmaY);
      tracks.fRhoSnpZ = (Char_t)(128. * track.getSigmaSnpZ() / tracks.fSigmaSnp / tracks.fSigmaZ);
      tracks.fRhoTglY = (Char_t)(128. * track.getSigmaTglY() / tracks.fSigmaTgl / tracks.fSigmaY);
      tracks.fRhoTglZ = (Char_t)(128. * track.getSigmaTglZ() / tracks.fSigmaTgl / tracks.fSigmaZ);
      tracks.fRhoTglSnp = (Char_t)(128. * track.getSigmaTglSnp() / tracks.fSigmaTgl / tracks.fSigmaSnp);
      tracks.fRho1PtY = (Char_t)(128. * track.getSigma1PtY() / tracks.fSigma1Pt / tracks.fSigmaY);
      tracks.fRho1PtZ = (Char_t)(128. * track.getSigma1PtZ() / tracks.fSigma1Pt / tracks.fSigmaZ);
      tracks.fRho1PtSnp = (Char_t)(128. * track.getSigma1PtSnp() / tracks.fSigma1Pt / tracks.fSigmaSnp);
      tracks.fRho1PtTgl = (Char_t)(128. * track.getSigma1PtTgl() / tracks.fSigma1Pt / tracks.fSigmaTgl);

      //insist it's good
      //hack: find last layer that's still good
      Float_t lLastLayer = 0;
      for(Int_t ilay=0;ilay<12;ilay++){
        if( track.hasHitOnLayer(ilay)==1) lLastLayer=ilay;
      }
      
      tracks.fITSChi2NCl = 1;
      tracks.fTPCSignal = 1;
      tracks.fTPCChi2NCl = 1.0;
      tracks.fTPCNClsFindable = (UChar_t)(120);
      tracks.fTPCNClsFindableMinusFound = (Char_t)(0);
      tracks.fTPCNClsFindableMinusCrossedRows = (Char_t)(0);
      UChar_t fITSClusterMap = 0u;
      fITSClusterMap |= 0x1 << 0; // flag manually
      fITSClusterMap |= 0x1 << 1; // flag manually
      tracks.fITSClusterMap = fITSClusterMap;

      //MC labels for MC use - negative, yes, but negative with offset
      mctracklabel.fIndexMcParticles = TMath::Abs(lab.getTrackID()) + fOffsetLabel;
      mctracklabel.fMcMask = 0;
      //Tag as fake. Note: used first bit only.
      if (lab.isFake())
        mctracklabel.fMcMask = 1;
      
      //Experimental: TOF information
      Double_t loTOFsignal = -1;
      Double_t liTOFsignal = -1;
      
      //cout<<"Track number "<<i<<" check cluster on layer 6 and 11: "<<track.hasHitOnLayer(6)<<" and "<<track.hasHitOnLayer(11)<<endl;

      Int_t lMCID = lab.getTrackID();

      if( lMCID>=0 ) {
        for (auto& hit : *hits) {
          if(hit.GetTrackID()!=lMCID) continue;

          int layer{hit.GetDetectorID()};
          float xyz[3]{hit.GetX(), hit.GetY(), hit.GetZ()};
          float r{std::hypot(xyz[0], xyz[1])};
          float phi{std::atan2(-xyz[1], -xyz[0]) + o2::its::constants::math::Pi};
          if (r > 99.0 && r < 101) layer = 11;
          //constexpr float LightSpeedCm2S = 299792458.e2;           // C in cm/s
          if( layer == 11 ){
            loTOFsignal = gRandom->Gaus(hit.GetTime(), 1e-12*loTOFresolution);
            //loTOFsignal = hit.GetTime();
          }
          if( layer == 6 )
            liTOFsignal = gRandom->Gaus(hit.GetTime(), 1e-12*liTOFresolution);
            //liTOFsignal = hit.GetTime();
        }
      }
      
      if( loTOFsignal < 0 && lOuterTOFSignalLutted[lMCID] > -99 ) loTOFsignal = lOuterTOFSignalLutted[lMCID]*1e-12;
      if( liTOFsignal < 0 && lInnerTOFSignalLutted[lMCID] > -99 ) liTOFsignal = lInnerTOFSignalLutted[lMCID]*1e-12;
      
      Double_t lFastestTimeoTOF = -1;
      Double_t lFastestTimeiTOF = -1;
      Float_t lX0 = track.getX();
      if(liTOFsignal>0){
        Float_t lX1;
        Bool_t lToOut = track.getXatLabR(20.0,lX1,tracker.getBz(),o2::track::DirOutward);
        std::array<float, 3> lPointN;
        std::array<float, 3> lPointNplus;
        Double_t lLength = 0.0;
        if(lToOut){
          track.propagateParamTo(lX0, tracker.getBz());
          for(Int_t iStep=1; iStep<lNStepsLIntegrator; iStep++){
            track.getXYZGlo(lPointN);
            Float_t lPosition = lX0 + (lX1-lX0)*((Float_t)(iStep))/((Float_t)(lNStepsLIntegrator-1));
            track.propagateParamTo(lPosition, tracker.getBz());
            track.getXYZGlo(lPointNplus);
            lLength += std::hypot( lPointNplus[0]-lPointN[0], lPointNplus[1]-lPointN[1], lPointNplus[2]-lPointN[2] );
          }
        }
        lFastestTimeiTOF = lLength / o2::constants::physics::LightSpeedCm2S; // In seconds!
      }

      if(loTOFsignal>0){
        Float_t lX1;
        Bool_t lToOut = track.getXatLabR(100.0,lX1,tracker.getBz(),o2::track::DirOutward);
        std::array<float, 3> lPointN;
        std::array<float, 3> lPointNplus;
        Double_t lLength = 0.0;
        if(lToOut){
          track.propagateParamTo(lX0, tracker.getBz());
          for(Int_t iStep=1; iStep<lNStepsLIntegrator; iStep++){
            track.getXYZGlo(lPointN);
            Float_t lPosition = lX0 + (lX1-lX0)*((Float_t)(iStep))/((Float_t)(lNStepsLIntegrator-1));
            track.propagateParamTo(lPosition, tracker.getBz());
            track.getXYZGlo(lPointNplus);
            lLength += std::hypot( lPointNplus[0]-lPointN[0], lPointNplus[1]-lPointN[1], lPointNplus[2]-lPointN[2] );
          }
        }
        lFastestTimeoTOF = lLength / o2::constants::physics::LightSpeedCm2S; // In seconds!
      }

      if(lMCID>=0){
        auto part = mcArr->at(lMCID);
        if(part.isPrimary()){
          if( lFastestTimeoTOF > 0) hPIDoTOF -> Fill( part.GetP(), lFastestTimeoTOF/loTOFsignal );
          if( lFastestTimeiTOF > 0) hPIDiTOF -> Fill( part.GetP(), lFastestTimeiTOF/liTOFsignal );
        }
      }
      
      //Outer TOF signal will be the TOF signal
      tracks.fTOFSignal = loTOFsignal;
      //Inner TOF signal will be the TPC signal
      tracks.fTPCSignal = liTOFsignal;
      
      fTree[kTracks]->Fill();
      fTree[kTracksCov]->Fill();
      fTree[kTracksExtra]->Fill();
      fTree[kMcTrackLabel]->Fill();
    }
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //---> MC stack information for de-referencing
    for (Long_t iii = 0; iii < (Long_t)mcArr->size(); iii++) {
      auto part = mcArr->at(iii);
      
      Float_t lMCEta = part.GetEta();
      Float_t lMCRap = part.GetRapidity();
      Float_t lMCPt = part.GetPt();
      Int_t lMCPDG = part.GetPdgCode();
      if(part.isPrimary()){
        if(TMath::Abs(lMCRap)<1.0){
          for(Int_t ipdg=0; ipdg<5; ipdg++){
            if(lMCPDG==lParticlePDG[ipdg]) hPtSpectraGenPID[ipdg]->Fill(lMCPt);
          }
        }
      }

      mcparticle.fIndexMcCollisions = lGoodEvents;

      //Get the kinematic values of the particles
      mcparticle.fPdgCode = part.GetPdgCode();
      mcparticle.fStatusCode = part.isPrimary();

      mcparticle.fFlags = 0;
      if (part.isSecondary())
        mcparticle.fFlags |= MCParticleFlags::ProducedInTransport;

      mcparticle.fIndexMcParticles_Mother0 = part.getMotherTrackId();
      if (mcparticle.fIndexMcParticles_Mother0 > -1)
        mcparticle.fIndexMcParticles_Mother0 += fOffsetLabel;
      mcparticle.fIndexMcParticles_Mother1 = -1;
      mcparticle.fIndexMcParticles_Daughter0 = part.getFirstDaughterTrackId();
      if (mcparticle.fIndexMcParticles_Daughter0 > -1)
        mcparticle.fIndexMcParticles_Daughter0 += fOffsetLabel;
      mcparticle.fIndexMcParticles_Daughter1 = part.getLastDaughterTrackId();
      if (mcparticle.fIndexMcParticles_Daughter1 > -1)
        mcparticle.fIndexMcParticles_Daughter1 += fOffsetLabel;
      mcparticle.fWeight = 1;

      mcparticle.fPx = part.Px();
      mcparticle.fPy = part.Py();
      mcparticle.fPz = part.Pz();
      mcparticle.fE = part.GetEnergy();

      mcparticle.fVx = part.Vx();
      mcparticle.fVy = part.Vy();
      mcparticle.fVz = part.Vz();
      mcparticle.fVt = part.T();
      
      fTree[kMcParticle]->Fill();
    }
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // Go for conversion: save info
    fTree[kEvents]->Fill();
    fTree[kMcCollision]->Fill();
    fTree[kMcCollisionLabel]->Fill();
    fTree[kBC]->Fill();
    fTree[kRun2BCInfo]->Fill();
    fTree[kFDD]->Fill();
    fTree[kFV0A]->Fill();
    fTree[kFV0C]->Fill();
    fTree[kFT0]->Fill();
    fTree[kZdc]->Fill();
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    fOffsetLabel = fTree[kMcParticle]->GetEntries(); //processed total
    //fOffsetLabel = 0;
    lGoodEvents++;
  }
  
  // Write out this time frame
  fOutputDir->cd();
  fTree[kEvents]->Write();
  fTree[kBC]->Write();
  fTree[kRun2BCInfo]->Write();
  fTree[kFDD]->Write();
  fTree[kFV0A]->Write();
  fTree[kFV0C]->Write();
  fTree[kFT0]->Write();
  fTree[kZdc]->Write();
  fTree[kTracks]->Write();
  fTree[kTracksCov]->Write();
  fTree[kTracksExtra]->Write();
  fTree[kMcTrackLabel]->Write();
  fTree[kMcParticle]->Write();
  fTree[kMcCollision]->Write();
  fTree[kMcCollisionLabel]->Write();

  TFile output("conversion-output.root", "recreate");

  //QA of conversion process: the basics
  hALICE3->Write();
  hNTracks->Write();
  hPtSpectraAll->Write();
  hPtSpectraFakeAll->Write();
  hPtSpectra->Write();
  hPtSpectraFake->Write();
  hPIDoTOF -> Write();
  hPIDiTOF -> Write();
  hPartialTransportDebug -> Write();
  
  for(Int_t ii=0; ii<5; ii++) hPtSpectraRecoPID[ii]->Write();
  for(Int_t ii=0; ii<5; ii++) hPtSpectraGenPID[ii]->Write();
  
  for(Int_t ilay=0; ilay<12; ilay++) hHitDensityMidRap[ilay] -> Write();

  std::cout << "*- Saved " << lGoodEvents << " events with a PV (total processed: " << itsHits.GetEntries() << "). Enjoy! \U0001F596" << std::endl;
}
