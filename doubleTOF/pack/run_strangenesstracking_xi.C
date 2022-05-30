#include <string>

#include <TFile.h>
#include <TChain.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TGeoGlobalMagField.h>
#include <vector>

#include <Generators/GeneratorFactory.h>
#include "FairPrimaryGenerator.h"
#include "FairGenerator.h"
#include "FairBoxGenerator.h"
#include <FairLogger.h>
#include <SimConfig/SimConfig.h>
#include <Generators/GeneratorFromFile.h>

#include "Steer/InteractionSampler.h"
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
#include "DataFormatsITSMFT/CompCluster.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "SimulationDataFormat/MCTrack.h"
#include "MathUtils/Cartesian.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "ReconstructionDataFormats/TrackParametrizationWithError.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTruthContainer.h"

#include <Generators/GeneratorFactory.h>
#include "FairPrimaryGenerator.h"
#include "FairGenerator.h"
#include "FairBoxGenerator.h"
#include <FairLogger.h>
#include <SimConfig/SimConfig.h>
#include <Generators/GeneratorFromFile.h>

#include "DetectorsVertexing/DCAFitterN.h"
#include "RecoDecay.h"

#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "Framework/DataTypes.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "ITSBase/GeometryTGeo.h"
#include "TStopwatch.h" //for bechmarking

#include "TrackSmearer.hh"
#include "TrackUtils.hh"

#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "CommonDataFormat/BunchFilling.h"
//

//This switch will enable or disable ALL CUTS BELOW
const Float_t lMasterSwitch = -1; //1: ON, -1: OFF

const Bool_t lAssociateXiDaughters = kTRUE;
const Bool_t lAssociateXic = kTRUE;
const Bool_t lAssociateXicc = kTRUE;
//Ask for hohlweger boolean
const Bool_t lHohlwegerBooleanXic = kFALSE;
const Bool_t lHohlwegerBooleanXicc = kFALSE;

//Tracking pass configuration
const Float_t lFirstOpening = 1;
const Float_t lOpening = 8.;
const Float_t lThirdPass = 30.;
const Float_t lFourthPass = 200.;

//===============================================================
// PRESELECTIONS
//===============================================================

//Remove fundamentally insane combinatorics, please
const Float_t lMinV0Radius = lMasterSwitch*0.475; //in centimeters
const Float_t lMinXiRadius = lMasterSwitch*0.475; //in centimeters
const Float_t lMinDCAxyBachelor = lMasterSwitch*0.0040; //in centimeters
const Float_t lMinDCAzBachelor = lMasterSwitch*0.0040; //in centimeters
const Float_t lMinDCAxyPositive = lMasterSwitch*0.0050; //in centimeters
const Float_t lMinDCAzPositive = lMasterSwitch*0.0040; //in centimeters
const Float_t lMinDCAxyNegative = lMasterSwitch*0.0100; //in centimeters
const Float_t lMinDCAzNegative = lMasterSwitch*0.0050; //in centimeters
const Float_t lMinDCAxyHFPions = lMasterSwitch*0.0010; //in centimeters
const Float_t lMinDCAzHFPions  = lMasterSwitch*0.0010; //in centimeters

const Float_t lMassWindowLambda = 0.012;
const Float_t lMassWindowXi = 0.012;
const Float_t lMassWindowXic = 0.060;
const Float_t lMassWindowXicc = 0.6;

const Float_t lMaxDCAV0Daughters = 1000; //in microns
const Float_t lMaxDCACascadeDaughters = 1200; //in microns

const Float_t lMinXiTransverseMomentum = 0.2;
const Float_t lMaxDCAxyPositive = 5000; //in microns

const Float_t lMaxTimeOffsetInnerTOF = 100; //in picoseconds, for weak decay
const Float_t lMaxTimeOffsetHFPions = 100; //in picoseconds, for HF pions

const Float_t lMaxDCAToPVXiCCCandidate = 50; //in microns

const Float_t lMaxDCAXiCDaughters = 50; //in microns
const Float_t lMaxDCAXiCCDaughters = 50; //in microns
//===============================================================

const Double_t lIBresolution = 0.00025f;
const Double_t lOBresolution = 0.00100f;

const Double_t liTOFresolution = 50; //ps
const Double_t loTOFresolution = 20; //ps
const Int_t lNStepsLIntegrator = 200;

using o2::its::MemoryParameters;
using o2::its::TrackingParameters;
using o2::itsmft::Hit;
using std::string;

using TrackingVertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

constexpr bool kUseSmearing{true};

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

struct particle {
  int pdg = 0;
  int pdgmother = 0;
  int motherID = 0;
  int nLayers = 0;
  float pt;
  float eta;
  float phi;
  float recoPt;
  float recoEta;
  float energyFirst;
  float energyLast;
  float vX;
  float vY;
  float vZ;
  int isReco = 0; // 0 = no, 1 = good, 2 = fake
};

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

float getDetLengthFromEta(const float eta, const float radius)
{
  return 10. * (10. + radius * std::cos(2 * std::atan(std::exp(-eta))));
}

Double_t Velocity(Double_t lMomentum, Double_t lMass){
  //Momentum p and mass m -> returns speed in centimeters per picosecond
  //Useful for TOF calculations
  Double_t lA = TMath::Power(lMomentum / lMass, 2);
  return 0.0299792458*TMath::Sqrt(lA/(1+lA));
}

Double_t TrackLength( o2::track::TrackParCov track, Double_t lX0, Double_t lX1 , Double_t lMagneticField ){
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

void run_strangenesstracking_xi(TString lDataPath = "./", TString lLutPath =  "../")
{
  
  TStopwatch lWatchTracking, lWatchAnalysis, lWatchSmearing;
  
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  std::cout << "\e[1;31m   ALICE 3 Strangeness tracking tool: XiCC      \e[0;00m" << std::endl;
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  
  cout<<"Main configuration: "<<endl;
  cout<<"Data path for analysis.......: "<<lDataPath.Data()<<endl;
  cout<<"Data path for LUTs...........: "<<lLutPath.Data()<<endl;
  
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
  
  //Operational parameters
  const Double_t lMagneticField = field->GetBz(0,0,0);
  cout<<"Magnetic field auto-detected to be "<<lMagneticField<<endl;
  
  Int_t Particle = 3122;
  Int_t lNLayers = 7;
  const string hitsFileName = Form("%so2sim_HitsTRK.root", lDataPath.Data());
  TChain mcTree("o2sim");
  mcTree.AddFile(Form("%so2sim_Kine.root", lDataPath.Data()));
  
  o2::dataformats::MCEventHeader* mcHead = nullptr;
  mcTree.SetBranchStatus("*", 0); //disable all branches
  mcTree.SetBranchStatus("MCTrack*", 1);
  mcTree.SetBranchStatus("MCEventHeader.*", 1);
  mcTree.SetBranchAddress("MCEventHeader.", &mcHead);
  
  std::vector<o2::MCTrack>* mcArr = nullptr;
  mcTree.SetBranchAddress("MCTrack", &mcArr);
  
  TChain itsHits("o2sim");
  
  itsHits.AddFile(hitsFileName.data());
  
  o2::its::Tracker tracker(new o2::its::TrackerTraitsCPU());
  tracker.setBz(lMagneticField);
  tracker.setCorrType(o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrTGeo);
  
  std::uint32_t roFrame;
  std::vector<Hit>* hits = nullptr;
  itsHits.SetBranchAddress("TRKHit", &hits);
  
//  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//  // TRACKING PARAMETERS - CLASSICAL
//  std::vector<TrackingParameters> trackParams(2);
//  //Tracking parameters for 12 layer setup
//  trackParams[0].NLayers = 12;
//  trackParams[0].MinTrackLength = 8; //this is the one with fixed params
//  std::vector<float> LayerRadii = {0.5f, 1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f, 100.0f};
//  std::vector<float> LayerZ(12);
//  for (int i{0}; i < 12; ++i)
//  LayerZ[i] = getDetLengthFromEta(1.44, LayerRadii[i]) + 1.;
//
//  //loosely based on run_trac_alice3.C but with extra stuff for the innermost layers
//  //FIXME: This may be subject to further tuning and is only a first guess
//  std::vector<float> TrackletMaxDeltaZ = {lFirstOpening*0.1f, lFirstOpening*0.1f, lFirstOpening*0.1f, lFirstOpening*0.1f, lFirstOpening*0.3f, lFirstOpening*0.3f, lFirstOpening*0.3f, lFirstOpening*0.3f, lFirstOpening*0.5f, lFirstOpening*0.5f, lFirstOpening*0.5f};
//  //std::vector<float> TrackletMaxDeltaZ = {0.5, 0.9, 0.9, 2.2, 3.4, 5.4, 6.7, 9.9, 9.9, 14, 14};;
//  std::vector<float> CellMaxDCA = {lFirstOpening*0.05f, lFirstOpening*0.05f, lFirstOpening*0.05f, lFirstOpening*0.04f, lFirstOpening*0.05f, lFirstOpening*0.2f, lFirstOpening*0.4f, lFirstOpening*0.5f, lFirstOpening*0.5f, lFirstOpening*0.5f};
//  std::vector<float> CellMaxDeltaZ = {lFirstOpening*0.2f, lFirstOpening*0.2f, lFirstOpening*0.2f,lFirstOpening*0.4f, lFirstOpening*0.5f, lFirstOpening*0.6f, lFirstOpening*3.0f, lFirstOpening*3.0f, lFirstOpening*3.0f,lFirstOpening*3.0f};
//  std::vector<float> NeighbourMaxDeltaCurvature = {lFirstOpening*0.012f, lFirstOpening*0.010f, lFirstOpening*0.008f, lFirstOpening*0.0025f, lFirstOpening*0.003f, lFirstOpening*0.0035f, lFirstOpening*0.004f, lFirstOpening*0.004f, lFirstOpening*0.005f};
//  std::vector<float> NeighbourMaxDeltaN = {lFirstOpening*0.002f, lFirstOpening*0.002f, lFirstOpening*0.002f, lFirstOpening*0.0090f, lFirstOpening*0.002f, lFirstOpening*0.005f,lFirstOpening* 0.005f,lFirstOpening* 0.005f, lFirstOpening*0.005f};
//
//  //Second pass - looser cuts, shorter tracks
//  std::vector<float> SecondTrackletMaxDeltaZ = {lOpening*0.5f, lOpening*0.9f, lOpening*0.9f, lOpening*2.2f, lOpening*3.4f, lOpening*5.4f, lOpening*6.7f, lOpening*9.9f, lOpening*9.9f, lOpening*14.0f, lOpening*14.0f};;
//  std::vector<float> SecondCellMaxDCA = {lOpening*0.05f, lOpening*0.05f, lOpening*0.05f, lOpening*0.04f, lOpening*0.05f, lOpening*0.2f, lOpening*0.4f, lOpening*0.5f, lOpening*0.5f, lOpening*0.5f};
//  std::vector<float> SecondCellMaxDeltaZ = {lOpening*0.2f, lOpening*0.2f, lOpening*0.2f, lOpening*0.4f, lOpening*0.5f, lOpening*0.6f, lOpening*3.0f, lOpening*3.0f, lOpening*3.0f, lOpening*3.0f};
//  std::vector<float> SecondNeighbourMaxDeltaCurvature = {lOpening*0.012f, lOpening*0.010f, lOpening*0.008f, lOpening*0.0025f, lOpening*0.003f, lOpening*0.0035f, lOpening*0.004f, lOpening*0.004f, lOpening*0.005f};
//  std::vector<float> SecondNeighbourMaxDeltaN = {lOpening*0.002f, lOpening*0.002f, lOpening*0.002f, lOpening*0.0090f, lOpening*0.002f, lOpening*0.005f, lOpening*0.005f, lOpening*0.005f, lOpening*0.005f};
//
//  trackParams[0].LayerRadii = LayerRadii;
//  trackParams[0].LayerZ = LayerZ;
//  trackParams[0].TrackletMaxDeltaPhi = lFirstOpening*0.3;
//  trackParams[0].CellMaxDeltaPhi = lFirstOpening*0.15;
//  trackParams[0].CellMaxDeltaTanLambda = lFirstOpening*0.03;
//  trackParams[0].TrackletMaxDeltaZ = TrackletMaxDeltaZ;
//  trackParams[0].CellMaxDCA = CellMaxDCA;
//  trackParams[0].CellMaxDeltaZ = CellMaxDeltaZ;
//  trackParams[0].NeighbourMaxDeltaCurvature = NeighbourMaxDeltaCurvature;
//  trackParams[0].NeighbourMaxDeltaN = NeighbourMaxDeltaN;
//
//  Float_t lFactor = 200;
//  std::vector<MemoryParameters> memParams(2);
//  std::vector<float> CellsMemoryCoefficients = {2.3208e-08f * 30, 2.104e-08f * 30, 1.6432e-08f * 30, 1.2412e-08f * 30, 1.3543e-08f * 30, 1.5e-08f * 30, 1.6e-08f * 30, 1.7e-08f * 30};
//  std::vector<float> TrackletsMemoryCoefficients = {0.0016353f * lFactor, 0.0013627f * lFactor, 0.000984f * lFactor, 0.00078135f * lFactor, 0.00057934f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor};
//  memParams[0].CellsMemoryCoefficients = CellsMemoryCoefficients;
//  memParams[0].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
//  memParams[0].MemoryOffset = 120000;
//
//  //Pass 2: capture long tracks that are lower in pT, etc
//  trackParams[1].NLayers = 12;
//  trackParams[1].MinTrackLength = 6; //this is the one with fixed params
//
//  trackParams[1].LayerRadii = LayerRadii;
//  trackParams[1].LayerZ = LayerZ;
//  trackParams[1].TrackletMaxDeltaPhi = lOpening*0.3;
//  trackParams[1].CellMaxDeltaPhi = lOpening*0.15;
//  trackParams[1].CellMaxDeltaTanLambda = lOpening*0.03;
//  trackParams[1].TrackletMaxDeltaZ = SecondTrackletMaxDeltaZ;
//  trackParams[1].CellMaxDCA = SecondCellMaxDCA;
//  trackParams[1].CellMaxDeltaZ = SecondCellMaxDeltaZ;
//  trackParams[1].NeighbourMaxDeltaCurvature = SecondNeighbourMaxDeltaCurvature;
//  trackParams[1].NeighbourMaxDeltaN = SecondNeighbourMaxDeltaN;
//
//  memParams[1].CellsMemoryCoefficients = CellsMemoryCoefficients;
//  memParams[1].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
//  memParams[1].MemoryOffset = 120000;
//
//  tracker.setParameters(memParams, trackParams);
//  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  
  
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  // TRACKING PARAMETERS
  std::vector<float> LayerRadii = {0.5f, 1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f, 100.0f};
  std::vector<float> LayerZ(12);
  for (int i{0}; i < 12; ++i)
  LayerZ[i] = getDetLengthFromEta(1.44, LayerRadii[i]) + 1.;

  // EXPERIMENTAL THREE-PASS VERSION
  std::vector<TrackingParameters> trackParams(4);
  //Tracking parameters for 12 layer setup
  trackParams[0].NLayers = 12;
  trackParams[0].MinTrackLength = 12; //this is the one with fixed params

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
  std::vector<float> ThirdTrackletMaxDeltaZ = {lThirdPass*0.5f, lThirdPass*0.9f, lThirdPass*0.9f, lThirdPass*2.2f, lThirdPass*3.4f, lThirdPass*5.4f, lThirdPass*6.7f, lThirdPass*9.9f, lThirdPass*9.9f, lThirdPass*14.0f, lThirdPass*14.0f};;
  std::vector<float> ThirdCellMaxDCA = {lThirdPass*0.05f, lThirdPass*0.05f, lThirdPass*0.05f, lThirdPass*0.04f, lThirdPass*0.05f, lThirdPass*0.2f, lThirdPass*0.4f, lThirdPass*0.5f, lThirdPass*0.5f, lThirdPass*0.5f};
  std::vector<float> ThirdCellMaxDeltaZ = {lThirdPass*0.2f, lThirdPass*0.2f, lThirdPass*0.2f, lThirdPass*0.4f, lThirdPass*0.5f, lThirdPass*0.6f, lThirdPass*3.0f, lThirdPass*3.0f, lThirdPass*3.0f, lThirdPass*3.0f};
  std::vector<float> ThirdNeighbourMaxDeltaCurvature = {lThirdPass*0.012f, lThirdPass*0.010f, lThirdPass*0.008f, lThirdPass*0.0025f, lThirdPass*0.003f, lThirdPass*0.0035f, lThirdPass*0.004f, lThirdPass*0.004f, lThirdPass*0.005f};
  std::vector<float> ThirdNeighbourMaxDeltaN = {lThirdPass*0.002f, lThirdPass*0.002f, lThirdPass*0.002f, lThirdPass*0.0090f, lThirdPass*0.002f, lThirdPass*0.005f, lThirdPass*0.005f, lThirdPass*0.005f, lThirdPass*0.005f};
  
  //Fourth pass: bonus
  std::vector<float> FourthTrackletMaxDeltaZ = {lFourthPass*0.5f, lFourthPass*0.9f, lFourthPass*0.9f, lFourthPass*2.2f, lFourthPass*3.4f, lFourthPass*5.4f, lFourthPass*6.7f, lFourthPass*9.9f, lFourthPass*9.9f, lFourthPass*14.0f, lFourthPass*14.0f};;
  std::vector<float> FourthCellMaxDCA = {lFourthPass*0.05f, lFourthPass*0.05f, lFourthPass*0.05f, lFourthPass*0.04f, lFourthPass*0.05f, lFourthPass*0.2f, lFourthPass*0.4f, lFourthPass*0.5f, lFourthPass*0.5f, lFourthPass*0.5f};
  std::vector<float> FourthCellMaxDeltaZ = {lFourthPass*0.2f, lFourthPass*0.2f, lFourthPass*0.2f, lFourthPass*0.4f, lFourthPass*0.5f, lFourthPass*0.6f, lFourthPass*3.0f, lFourthPass*3.0f, lFourthPass*3.0f, lFourthPass*3.0f};
  std::vector<float> FourthNeighbourMaxDeltaCurvature = {lFourthPass*0.012f, lFourthPass*0.010f, lFourthPass*0.008f, lFourthPass*0.0025f, lFourthPass*0.003f, lFourthPass*0.0035f, lFourthPass*0.004f, lFourthPass*0.004f, lFourthPass*0.005f};
  std::vector<float> FourthNeighbourMaxDeltaN = {lFourthPass*0.002f, lFourthPass*0.002f, lFourthPass*0.002f, lFourthPass*0.0090f, lFourthPass*0.002f, lFourthPass*0.005f, lFourthPass*0.005f, lFourthPass*0.005f, lFourthPass*0.005f};

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

  Float_t lFactor = 2000;
  std::vector<MemoryParameters> memParams(4);
  std::vector<float> CellsMemoryCoefficients = {2.3208e-08f * 50, 2.104e-08f * 50, 1.6432e-08f * 50, 1.2412e-08f * 50, 1.3543e-08f * 50, 1.5e-08f * 50, 1.6e-08f * 50, 1.7e-08f * 50};
  std::vector<float> TrackletsMemoryCoefficients = {0.0016353f * lFactor, 0.0013627f * lFactor, 0.000984f * lFactor, 0.00078135f * lFactor, 0.00057934f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor};

  memParams[0].CellsMemoryCoefficients = CellsMemoryCoefficients;
  memParams[0].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
  memParams[0].MemoryOffset = 120000;

  //Pass 2: capture long tracks that are lower in pT, etc
  trackParams[1].NLayers = 12;
  trackParams[1].MinTrackLength = 8; //this is the one with fixed params

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
  trackParams[2].TrackletMaxDeltaPhi = TMath::Min(lThirdPass*0.3, TMath::Pi());
  trackParams[2].CellMaxDeltaPhi = lThirdPass*0.15;
  trackParams[2].CellMaxDeltaTanLambda = lThirdPass*0.03;
  trackParams[2].TrackletMaxDeltaZ = ThirdTrackletMaxDeltaZ;
  trackParams[2].CellMaxDCA = ThirdCellMaxDCA;
  trackParams[2].CellMaxDeltaZ = ThirdCellMaxDeltaZ;
  trackParams[2].NeighbourMaxDeltaCurvature = ThirdNeighbourMaxDeltaCurvature;
  trackParams[2].NeighbourMaxDeltaN = ThirdNeighbourMaxDeltaN;

  memParams[2].CellsMemoryCoefficients = CellsMemoryCoefficients;
  memParams[2].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
  memParams[2].MemoryOffset = 180000;
  
  //Pass 4: Bonus
  trackParams[3].NLayers = 12;
  trackParams[3].MinTrackLength = 6; //this is the one with fixed params

  trackParams[3].LayerRadii = LayerRadii;
  trackParams[3].LayerZ = LayerZ;
  trackParams[3].TrackletMaxDeltaPhi = TMath::Min(lFourthPass*0.3, TMath::Pi());
  trackParams[3].CellMaxDeltaPhi = lFourthPass*0.15;
  trackParams[3].CellMaxDeltaTanLambda = lFourthPass*0.03;
  trackParams[3].TrackletMaxDeltaZ = FourthTrackletMaxDeltaZ;
  trackParams[3].CellMaxDCA = FourthCellMaxDCA;
  trackParams[3].CellMaxDeltaZ = FourthCellMaxDeltaZ;
  trackParams[3].NeighbourMaxDeltaCurvature = FourthNeighbourMaxDeltaCurvature;
  trackParams[3].NeighbourMaxDeltaN = FourthNeighbourMaxDeltaN;

  memParams[3].CellsMemoryCoefficients = CellsMemoryCoefficients;
  memParams[3].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
  memParams[3].MemoryOffset = 180000;

  tracker.setParameters(memParams, trackParams);
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  
  Double_t ptbinlimits[] ={ 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.2,2.4,2.6,2.8,3.0,
    3.3,3.6,3.9,4.2,4.6,5,5.4,5.9, 6.5,7,7.5,8,8.5,9.2,10,11,12,13.5,15,17,20};
  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
  
  TH1D dcav0("dcav0", ";DCA V0 daughters (#mum);Counts", 400, 0.0,400);
  TH1D cosPA("cosPA", ";CosPA;", 400, -1,1);
  TH1D v0radius("v0radius", ";v0radius;", 500, 0,50);
  TH1D hmassk0("hmassk0", ";Mass;", 2000,0,2);
  
  TH1D *hNVertices = new TH1D("hNVertices", "", 10,0,10);
  TH1D *hNTracks = new TH1D("hNTracks", "", 10000,0,10000);
  TH1D *hNLongTracks = new TH1D("hNLongTracks", "", 10000,0,10000);
  TH1D *hNLongPrimaryTracks = new TH1D("hNLongPrimaryTracks", "", 10000,0,10000);
  TH1D *hNPrimaryTracks = new TH1D("hNPrimaryTracks", "", 10000,0,10000);
  
  TH1D *hGenerated = new TH1D("hGenerated", "", ptbinnumb, ptbinlimits);
  TH1D *hPt = new TH1D("hPt", "", ptbinnumb, ptbinlimits);
  TH2D *hMassVsPt = new TH2D("hMassVsPt", "", 150,0,15,2000,0,2);
  TH2D *hDCADauVsPt = new TH2D("hDCADauVsPt", "", 150,0,15,600,0,600);
  TH2D *hPAVsPt = new TH2D("hPAVsPt", "", 150,0,15,200,0,0.5);
  
  TProfile *hDCAProf = new TProfile("hDCAProf","",ptbinnumb, ptbinlimits);
  TProfile *hPAProf = new TProfile("hPAProf","",ptbinnumb, ptbinlimits);
  
  TH1D *hCascFinding = new TH1D("hCascFinding", "", 10,0,10);
  TH1D *hMassXi = new TH1D("hMassXi", "", 200, 1.672-0.050, 1.672+0.050);
  
  TH1D dcacasc("dcacasc", ";DCA V0 daughters (#mum);Counts", 400, 0.0,400);
  TProfile *hCascDCAProf = new TProfile("hCascDCAProf","",ptbinnumb, ptbinlimits);
  TH1D *hPtCascade = new TH1D("hPtCascade", "", 100,0,10);
  TH2D *hCascadeMassVsPt = new TH2D("hCascadeMassVsPt", "", 150,0,15,2000,0,2);
  TH1D* hCascadeCosPA = new TH1D("hCascadeCosPA", "", 400, -1, 1);
  TProfile *hCascadePAProf = new TProfile("hCascadePAProf","",ptbinnumb, ptbinlimits);
  
  TH1D *hTrackCount = new TH1D("hTrackCount", "", 5,0,5);
  TH1D *hCombinatorics = new TH1D("hCombinatorics", "",2,0,2);
  TH1D *hCombinatoricsV0 = new TH1D("hCombinatoricsV0", "",5,0,5);
  
  TH1D *hHitXa = new TH1D("hHitXa", "", 2000,-10,10);
  TH1D *hHitYa = new TH1D("hHitYa", "", 2000,-10,10);
  TH1D *hHitZa = new TH1D("hHitZa", "", 2000,-10,10);
  TH1D *hHitNumber = new TH1D("hHitNumber", "", 10,0,10);
  
  TH1D *hXiR2D = new TH1D("hXiR2D", "", 2000,0,100);
  TH1D *hLambdaR2D = new TH1D("hLambdaR2D", "", 2000,0,100);
  Double_t lRadii[] = {0.0f,0.5f,1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f, 100.0f};
  TH1D *hXiR2Dbinned = new TH1D("hXiR2Dbinned", "", 12, lRadii);
  TH1D *hLambdaR2Dbinned = new TH1D("hLambdaR2Dbinned", "", 12, lRadii);
  
  TH2D *hPtEtaGeneratedXi = new TH2D("hPtEtaGeneratedXi", "", 200,0,20,30,-1.5,1.5);
  TH2D *hPtEtaGeneratedXiC = new TH2D("hPtEtaGeneratedXiC", "", 200,0,20,30,-1.5,1.5);
  TH2D *hPtEtaGeneratedXiCC = new TH2D("hPtEtaGeneratedXiCC", "", 200,0,20,30,-1.5,1.5);
  
  TH2D *hPtYGeneratedXi = new TH2D("hPtYGeneratedXi", "", 200,0,20,30,-1.5,1.5);
  TH2D *hPtYGeneratedXiC = new TH2D("hPtYGeneratedXiC", "", 200,0,20,30,-1.5,1.5);
  TH2D *hPtYGeneratedXiCC = new TH2D("hPtYGeneratedXiCC", "", 200,0,20,30,-1.5,1.5);
  
  Long_t lMaxMult = 10000;
  
  TH2D *hPtMultGeneratedXi_AllEta = new TH2D("hPtMultGeneratedXi_AllEta", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiC_AllEta = new TH2D("hPtMultGeneratedXiC_AllEta", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiCC_AllEta = new TH2D("hPtMultGeneratedXiCC_AllEta", "", 200,0,20,lMaxMult,0,lMaxMult);
  
  TH2D *hPtMultGeneratedXi_AllY = new TH2D("hPtMultGeneratedXi_AllY", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiC_AllY = new TH2D("hPtMultGeneratedXiC_AllY", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiCC_AllY = new TH2D("hPtMultGeneratedXiCC_AllY", "", 200,0,20,lMaxMult,0,lMaxMult);
  
  TH2D *hPtMultGeneratedXi_Eta5 = new TH2D("hPtMultGeneratedXi_Eta5", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiC_Eta5 = new TH2D("hPtMultGeneratedXiC_Eta5", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiCC_Eta5 = new TH2D("hPtMultGeneratedXiCC_Eta5", "", 200,0,20,lMaxMult,0,lMaxMult);
  
  TH2D *hPtMultGeneratedXi_Y5 = new TH2D("hPtMultGeneratedXi_Y5", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiC_Y5 = new TH2D("hPtMultGeneratedXiC_Y5", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiCC_Y5 = new TH2D("hPtMultGeneratedXiCC_Y5", "", 200,0,20,lMaxMult,0,lMaxMult);
  
  TH2D *hPtMultGeneratedXi_Eta10 = new TH2D("hPtMultGeneratedXi_Eta10", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiC_Eta10 = new TH2D("hPtMultGeneratedXiC_Eta10", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiCC_Eta10 = new TH2D("hPtMultGeneratedXiCC_Eta10", "", 200,0,20,lMaxMult,0,lMaxMult);
  
  TH2D *hPtMultGeneratedXi_Y10 = new TH2D("hPtMultGeneratedXi_Y10", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiC_Y10 = new TH2D("hPtMultGeneratedXiC_Y10", "", 200,0,20,lMaxMult,0,lMaxMult);
  TH2D *hPtMultGeneratedXiCC_Y10 = new TH2D("hPtMultGeneratedXiCC_Y10", "", 200,0,20,lMaxMult,0,lMaxMult);
  
  TH1D *hTime = new TH1D("hTime", "Processing time", 6, 0, 6);
  hTime->GetXaxis()->SetBinLabel(1, "Smearing real time");
  hTime->GetXaxis()->SetBinLabel(2, "Tracking real time");
  hTime->GetXaxis()->SetBinLabel(3, "Analysis real time");
  hTime->GetXaxis()->SetBinLabel(4, "Smearing CPU time");
  hTime->GetXaxis()->SetBinLabel(5, "Tracking CPU time");
  hTime->GetXaxis()->SetBinLabel(6, "Analysis CPU time");
  
  TH1D *hPartialTransportDebug = new TH1D("hPartialTransportDebug", "", 6,0,6);
  hPartialTransportDebug->GetXaxis()->SetBinLabel(1, "All particles");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(2, "Not tranported");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(3, "Successfully smeared");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(4, "Fully tracked particles");
  
  TH1D *hQADCAxy = new TH1D("hQADCAxy","", 1000,-500,500);
  TH1D *hQADCAz = new TH1D("hQADCAz","", 1000,-500,500);
  
  TH1D *hQAMCNch = new TH1D("hQAMCNch", "", 6000, 0, 6000);
  
  TH2D *hQANchCorrelation = new TH2D("hQANchCorrelation", "", 600, 0, 12000, 600, 0, 12000);
  
  //Define o2 fitter, 2-prong
  o2::vertexing::DCAFitterN<2> fitterV0, fitterCasc, fitterCascCC;
  o2::vertexing::DCAFitterN<3> retrofitter, fitterCascC;
  
  fitterV0.setBz(lMagneticField);
  fitterV0.setPropagateToPCA(true);
  fitterV0.setMaxR(200.);
  fitterV0.setMinParamChange(1e-5);
  fitterV0.setMinRelChi2Change(0.9);
  fitterV0.setMaxDZIni(1e9);
  fitterV0.setMaxChi2(1e9);
  fitterV0.setUseAbsDCA(true);
  
  fitterCasc.setBz(lMagneticField);
  fitterCasc.setPropagateToPCA(true);
  fitterCasc.setMaxR(200.);
  fitterCasc.setMinParamChange(1e-3);
  fitterCasc.setMinRelChi2Change(0.9);
  fitterCasc.setMaxDZIni(1e9);
  fitterCasc.setMaxChi2(1e9);
  fitterCasc.setUseAbsDCA(true);
  
  fitterCascC.setBz(lMagneticField);
  fitterCascC.setPropagateToPCA(true);
  fitterCascC.setMaxR(200.);
  fitterCascC.setMinParamChange(1e-3);
  fitterCascC.setMinRelChi2Change(0.9);
  fitterCascC.setMaxDZIni(1e9);
  fitterCascC.setMaxChi2(1e9);
  fitterCascC.setUseAbsDCA(true);
  
  fitterCascCC.setBz(lMagneticField);
  fitterCascCC.setPropagateToPCA(true);
  fitterCascCC.setMaxR(200.);
  fitterCascCC.setMinParamChange(1e-3);
  fitterCascCC.setMinRelChi2Change(0.9);
  fitterCascCC.setMaxDZIni(1e9);
  fitterCascCC.setMaxChi2(1e9);
  fitterCascCC.setUseAbsDCA(true);
  
  //Atempt retro-fit: fit 3 tracks to Xi
  retrofitter.setBz(lMagneticField);
  retrofitter.setPropagateToPCA(true);
  retrofitter.setMaxR(200.);
  retrofitter.setMinParamChange(1e-3);
  retrofitter.setMinRelChi2Change(0.9);
  retrofitter.setMaxDZIni(1e9);
  retrofitter.setMaxChi2(1e9);
  retrofitter.setUseAbsDCA(true);
  
  Int_t fParticle, fMinNLayerReco;
  Double_t fDCAV0Dau, fPA, fV0Radius, fMass, fPt, fPtMC, fPtMCMother;
  
  TFile *fileTree = new TFile("xi.treeoutput.root", "RECREATE", "", 109);
  
  TH1D *hEventCounter = new TH1D("hEventCounter", "", 1,0,1);
  TH1D *hEventCounterWithVertex = new TH1D("hEventCounterWithVertex", "", 1,0,1);
  TH1D *hVertexZ = new TH1D("hVertexZ", "", 100, -15,15);
  TH1D *hXiGeneratedPt = new TH1D("hXiGeneratedPt", "", 200,0,20);
  TH1D *hXiCGeneratedPt = new TH1D("hXiCGeneratedPt", "", 200,0,20);
  TH1D *hXiCCGeneratedPt = new TH1D("hXiCCGeneratedPt", "", 200,0,20);
  
  TH1D *hNDirectTrackingLayers = new TH1D("hNDirectTrackingLayers", "",12,0,12);
  TH1D *hNV0DauLayers = new TH1D("hNV0DauLayers", "",12,0,12);
  
  //Store hit densities at each layer
  TH2D *hHitDensityMap[12]; //overall map
  TH1D *hHitDensityCounter[12]; //single event
  TProfile *hHitDensityProfile[12];
  for(Int_t ii=0; ii<12; ii++){
    //Make adaptive binning as required thinking radially outwards
    hHitDensityCounter[ii]=new TH1D(Form("hHitDensityCounter_%i", ii), "", 300,0,150);
    hHitDensityMap[ii]=new TH2D(Form("hHitDensityMap_%i", ii), "", 300,0,150,100,0,20000*TMath::Power(LayerRadii[0],1)/TMath::Power(LayerRadii[ii],1));
    hHitDensityProfile[ii] = new TProfile(Form("hHitDensityProfile_%i", ii), "", 300,0,150);
  }
  
  //Declare
  Float_t fPtXi=0, fPtMCXi=0;
  Float_t fXicPt=0, fXicPtMC=0;
  
  Bool_t fTrueXi = kFALSE;
  Bool_t fTrueXic = kFALSE;
  Bool_t fTrueXicc = kFALSE;
  
  Float_t fXiDCAxyToPV=0;
  Float_t fXiDCAzToPV=0;
  Float_t fXicDCAxyToPV=0;
  Float_t fXicDCAzToPV=0;
  Float_t fXiccDCAxyToPV=0;
  Float_t fXiccDCAzToPV=0;
  
  Int_t fXiHitsAdded=0;;
  Float_t fXicMass=0;
  
  Float_t fPic1DCAxyToPV =0;
  Float_t fPic1DCAzToPV =0;
  Float_t fPic2DCAxyToPV =0;
  Float_t fPic2DCAzToPV =0;
  
  Float_t fV0DauDCA = 0;
  Float_t fXiDauDCA = 0;
  
  Float_t fXiDecayRadius = 0;
  Float_t fV0DecayRadius = 0;
  Float_t fXiDecayRadiusMC = 0;
  Float_t fV0DecayRadiusMC = 0;
  
  Float_t fLambdaMass = 0, fXiMass = 0;
  
  Float_t fXicDecayRadius = 0;
  Float_t fXicDecayDistanceFromPV = 0;
  Float_t fXicDaughterDCA = 0;
  Float_t fXiccDecayRadius = 0;
  Float_t fXiccDecayDistanceFromPV = 0;
  Float_t fXiccDaughterDCA = 0;
  
  Float_t fPrimaryVertexX = 0, fPrimaryVertexY = 0, fPrimaryVertexZ = 0;
  Float_t fXicVertexX = 0, fXicVertexY = 0, fXicVertexZ = 0;
  Float_t fXicVertexXMC = 0, fXicVertexYMC = 0, fXicVertexZMC = 0;
  Float_t fXiccMass=0;
  
  Float_t fPiccDCAxyToPV, fPiccDCAzToPV;
  
  Bool_t fFirstCandidateXiC, fFirstCandidateXiCC;
  
  Int_t fNTracks = 0, fNLongTracks = 0, fNLongPrimaryTracks = 0, fNPrimaryTracks = 0;
  
  TTree *fTreeCandidates = new TTree ( "fTreeCandidates", "Event Characterization Tree" ) ;
  fTreeCandidates->Branch ("fNTracks",    &fNTracks,    "fNTracks/I"    );
  fTreeCandidates->Branch ("fNLongTracks",    &fNLongTracks,    "fNLongTracks/I"    );
  fTreeCandidates->Branch ("fNLongPrimaryTracks",    &fNLongPrimaryTracks,    "fNLongPrimaryTracks/I"    );
  fTreeCandidates->Branch ("fNPrimaryTracks",    &fNPrimaryTracks,    "fNPrimaryTracks/I"    );

  Float_t fPtXiC, fPtMCXiC, fPtXiCC, fPtMCXiCC;
  Float_t fPXiC, fPMCXiC, fPXiCC, fPMCXiCC;
  
  fTreeCandidates->Branch ("fPtXi",    &fPtXi,    "fPtXi/F"    );
  fTreeCandidates->Branch ("fPtXiC",  &fPtXiC,  "fPtXiC/F"  );
  fTreeCandidates->Branch ("fPtXiCC",  &fPtXiCC,  "fPtXiCC/F"  );
  fTreeCandidates->Branch ("fPtMCXi",  &fPtMCXi,  "fPtMCXi/F"  );
  fTreeCandidates->Branch ("fPtMCXiC",  &fPtMCXiC,  "fPtMCXiC/F"  );
  fTreeCandidates->Branch ("fPtMCXiCC",  &fPtMCXiCC,  "fPtMCXiCC/F"  );
  
  fTreeCandidates->Branch ("fPXiC",  &fPXiC,  "fPXiC/F"  );
  fTreeCandidates->Branch ("fPXiCC",  &fPXiCC,  "fPXiCC/F"  );
  fTreeCandidates->Branch ("fPMCXiC",  &fPMCXiC,  "fPMCXiC/F"  );
  fTreeCandidates->Branch ("fPMCXiCC",  &fPMCXiCC,  "fPMCXiCC/F"  );
  
  fTreeCandidates->Branch ("fXiDecayRadiusMC",  &fXiDecayRadiusMC,  "fXiDecayRadiusMC/F"  );
  fTreeCandidates->Branch ("fV0DecayRadiusMC",  &fV0DecayRadiusMC,  "fV0DecayRadiusMC/F"  );
  fTreeCandidates->Branch ("fXiDecayRadius",  &fXiDecayRadius,  "fXiDecayRadius/F"  );
  fTreeCandidates->Branch ("fV0DecayRadius",  &fV0DecayRadius,  "fV0DecayRadius/F"  );
  
  fTreeCandidates->Branch ("fXiDCAxyToPV",  &fXiDCAxyToPV,  "fXiDCAxyToPV/F"  );
  fTreeCandidates->Branch ("fXiDCAzToPV",  &fXiDCAzToPV,  "fXiDCAzToPV/F"  );
  
  fTreeCandidates->Branch ("fV0DauDCA",  &fV0DauDCA,  "fV0DauDCA/F"  );
  fTreeCandidates->Branch ("fXiDauDCA",  &fXiDauDCA,  "fXiDauDCA/F"  );
  fTreeCandidates->Branch ("fTrueXi",  &fTrueXi,  "fTrueXi/O"  );

  fTreeCandidates->Branch ("fXiHitsAdded",  &fXiHitsAdded,  "fXiHitsAdded/I"  );

  fTreeCandidates->Branch ("fLambdaMass",  &fLambdaMass,  "fLambdaMass/F"  );
  fTreeCandidates->Branch ("fXiMass",  &fXiMass,  "fXiMass/F"  );
  
  //PDG Codes
  Int_t fPiccMotherPDG, fPic1MotherPDG, fPic2MotherPDG;
  Float_t fPiccMotherPt, fPic1MotherPt, fPic2MotherPt;
  Int_t fPiccMotherNChain, fPic1MotherNChain, fPic2MotherNChain;
  Int_t fPiccMotherChain[20], fPic1MotherChain[20], fPic2MotherChain[20];
  
  Float_t fXiCCtoXiCLength, fXiCtoXiLength;
  
  Float_t fBachelorDCAxy, fBachelorDCAz, fNegativeDCAxy, fNegativeDCAz, fPositiveDCAxy, fPositiveDCAz;
  
  fTreeCandidates->Branch ("fBachelorDCAxy",  &fBachelorDCAxy,  "fBachelorDCAxy/F"  );
  fTreeCandidates->Branch ("fBachelorDCAz",  &fBachelorDCAz,  "fBachelorDCAz/F"  );
  fTreeCandidates->Branch ("fNegativeDCAxy",  &fNegativeDCAxy,  "fNegativeDCAxy/F"  );
  fTreeCandidates->Branch ("fNegativeDCAz",  &fNegativeDCAz,  "fNegativeDCAz/F"  );
  fTreeCandidates->Branch ("fPositiveDCAxy",  &fPositiveDCAxy,  "fPositiveDCAxy/F"  );
  fTreeCandidates->Branch ("fPositiveDCAz",  &fPositiveDCAz,  "fPositiveDCAz/F"  );
  
  Float_t fV0DCAxyToPV, fV0DCAzToPV, fV0DecayLength, fXiDecayLength;
  fTreeCandidates->Branch ("fV0DCAxyToPV",  &fV0DCAxyToPV,  "fV0DCAxyToPV/F"  );
  fTreeCandidates->Branch ("fV0DCAzToPV",  &fV0DCAzToPV,  "fV0DCAzToPV/F"  );
  fTreeCandidates->Branch ("fV0DecayLength",  &fV0DecayLength,  "fV0DecayLength/F"  );
  fTreeCandidates->Branch ("fXiDecayLength",  &fXiDecayLength,  "fXiDecayLength/F"  );
  
  //Six prong momenta for efficiency composition
  Float_t fPiCCPt, fPiC1Pt, fPiC2Pt, fBachelorPt, fPositivePt, fNegativePt;
  fTreeCandidates->Branch ("fBachelorPt",  &fBachelorPt,  "fBachelorPt/F"  );
  fTreeCandidates->Branch ("fPositivePt",  &fPositivePt,  "fPositivePt/F"  );
  fTreeCandidates->Branch ("fNegativePt",  &fNegativePt,  "fNegativePt/F"  );
  
  Float_t fV0TotalMomentum, fXiTotalMomentum;
  fTreeCandidates->Branch ("fV0TotalMomentum",  &fV0TotalMomentum,  "fV0TotalMomentum/F"  );
  fTreeCandidates->Branch ("fXiTotalMomentum",  &fXiTotalMomentum,  "fXiTotalMomentum/F"  );
  
  Int_t fPiC1Clusters, fPiC2Clusters, fPiCCClusters, fNegativeClusters, fPositiveClusters, fBachelorClusters;
  Float_t fPiC1Chisquare, fPiC2Chisquare, fPiCCChisquare, fNegativeChisquare, fPositiveChisquare, fBachelorChisquare;
  
  fTreeCandidates->Branch ("fNegativeClusters",  &fNegativeClusters,  "fNegativeClusters/I"  );
  fTreeCandidates->Branch ("fPositiveClusters",  &fPositiveClusters,  "fPositiveClusters/I"  );
  fTreeCandidates->Branch ("fBachelorClusters",  &fBachelorClusters,  "fBachelorClusters/I"  );
  
  fTreeCandidates->Branch ("fNegativeChisquare",  &fNegativeChisquare,  "fNegativeChisquare/F"  );
  fTreeCandidates->Branch ("fPositiveChisquare",  &fPositiveChisquare,  "fPositiveChisquare/F"  );
  fTreeCandidates->Branch ("fBachelorChisquare",  &fBachelorChisquare,  "fBachelorChisquare/F"  );
  
  TLorentzVector *lVector = new TLorentzVector();
  Float_t fPiCCEta, fPiC1Eta, fPiC2Eta, fNegativeEta, fPositiveEta, fBachelorEta;

  fTreeCandidates->Branch ("fNegativeEta",  &fNegativeEta,  "fNegativeEta/F"  );
  fTreeCandidates->Branch ("fPositiveEta",  &fPositiveEta,  "fPositiveEta/F"  );
  fTreeCandidates->Branch ("fBachelorEta",  &fBachelorEta,  "fBachelorEta/F"  );
  
  Float_t fXiEta, fXiCEta, fXiCCEta;
  fTreeCandidates->Branch ("fXiEta",  &fXiEta,  "fXiEta/F"  );
  
  Bool_t fPiccUsed;
  Int_t fPicUsed;
  
  Float_t fBachelorTOFSignal, fPositiveTOFSignal, fNegativeTOFSignal, fPic1TOFSignal, fPic2TOFSignal, fPiccTOFSignal;
  Float_t fBachelorTOFSignalMC, fPositiveTOFSignalMC, fNegativeTOFSignalMC, fPic1TOFSignalMC, fPic2TOFSignalMC, fPiccTOFSignalMC;
  Float_t fBachelorExpectedSignal, fPositiveExpectedSignal, fNegativeExpectedSignal, fPic1ExpectedSignal, fPic2ExpectedSignal, fPiccExpectedSignal;
  Float_t fBachelorExpectedSignalFromPV, fPositiveExpectedSignalFromPV, fNegativeExpectedSignalFromPV;
  
  Float_t fBachelorLength, fNegativeLength, fPositiveLength, fCascadeLength, fV0Length;
  Float_t fBachelorLife, fNegativeLife, fPositiveLife, fCascadeLife, fV0Life;
  
  Int_t fPiccPDG, fPic1PDG, fPic2PDG, fNegativePDG, fPositivePDG, fBachelorPDG;
  fTreeCandidates->Branch ("fNegativePDG",  &fNegativePDG,  "fNegativePDG/I"  );
  fTreeCandidates->Branch ("fPositivePDG",  &fPositivePDG,  "fPositivePDG/I"  );
  fTreeCandidates->Branch ("fBachelorPDG",  &fBachelorPDG,  "fBachelorPDG/I"  );
  
  //INNER TOF
  Float_t fBachelorInnerTOF20Signal, fPositiveInnerTOF20Signal, fNegativeInnerTOF20Signal, fPic1InnerTOF20Signal, fPic2InnerTOF20Signal, fPiccInnerTOF20Signal;
  Float_t fBachelorInnerTOF50Signal, fPositiveInnerTOF50Signal, fNegativeInnerTOF50Signal, fPic1InnerTOF50Signal, fPic2InnerTOF50Signal, fPiccInnerTOF50Signal;
  Float_t fBachelorInnerTOFSignalMC, fPositiveInnerTOFSignalMC, fNegativeInnerTOFSignalMC, fPic1InnerTOFSignalMC, fPic2InnerTOFSignalMC, fPiccInnerTOFSignalMC;
  Float_t fBachelorInnerExpectedSignal, fPositiveInnerExpectedSignal, fNegativeInnerExpectedSignal, fPic1InnerExpectedSignal, fPic2InnerExpectedSignal, fPiccInnerExpectedSignal;
  
  fTreeCandidates->Branch ("fBachelorInnerTOF20Signal",  &fBachelorInnerTOF20Signal,  "fBachelorInnerTOF20Signal/F"  );
  fTreeCandidates->Branch ("fPositiveInnerTOF20Signal",  &fPositiveInnerTOF20Signal,  "fPositiveInnerTOF20Signal/F"  );
  fTreeCandidates->Branch ("fNegativeInnerTOF20Signal",  &fNegativeInnerTOF20Signal,  "fNegativeInnerTOF20Signal/F"  );
  
  fTreeCandidates->Branch ("fBachelorInnerTOF50Signal",  &fBachelorInnerTOF50Signal,  "fBachelorInnerTOF50Signal/F"  );
  fTreeCandidates->Branch ("fPositiveInnerTOF50Signal",  &fPositiveInnerTOF50Signal,  "fPositiveInnerTOF50Signal/F"  );
  fTreeCandidates->Branch ("fNegativeInnerTOF50Signal",  &fNegativeInnerTOF50Signal,  "fNegativeInnerTOF50Signal/F"  );
  
  fTreeCandidates->Branch ("fBachelorInnerTOFSignalMC",  &fBachelorInnerTOFSignalMC,  "fBachelorInnerTOFSignalMC/F"  );
  fTreeCandidates->Branch ("fPositiveInnerTOFSignalMC",  &fPositiveInnerTOFSignalMC,  "fPositiveInnerTOFSignalMC/F"  );
  fTreeCandidates->Branch ("fNegativeInnerTOFSignalMC",  &fNegativeInnerTOFSignalMC,  "fNegativeInnerTOFSignalMC/F"  );
  
  fTreeCandidates->Branch ("fBachelorInnerExpectedSignal",  &fBachelorInnerExpectedSignal,  "fBachelorInnerExpectedSignal/F"  );
  fTreeCandidates->Branch ("fPositiveInnerExpectedSignal",  &fPositiveInnerExpectedSignal,  "fPositiveInnerExpectedSignal/F"  );
  fTreeCandidates->Branch ("fNegativeInnerExpectedSignal",  &fNegativeInnerExpectedSignal,  "fNegativeInnerExpectedSignal/F"  );
  
  Bool_t fUsesXiCCProngs, fTrueNegative, fTrueBachelor, fTruePositive, fTruePicc, fTruePic1, fTruePic2;
  fTreeCandidates->Branch ("fTrueNegative",  &fTrueNegative,  "fTrueNegative/O"  );
  fTreeCandidates->Branch ("fTrueBachelor",  &fTrueBachelor,  "fTrueBachelor/O"  );
  fTreeCandidates->Branch ("fTruePositive",  &fTruePositive,  "fTruePositive/O"  );
  
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
  gen->Init();
  cout<<" done"<<endl;
  
  //This is for convenience -> Be careful, please!
  fParticle = Particle;
  fMinNLayerReco = lNLayers;
  
  o2::steer::InteractionSampler irSampler;
  irSampler.setInteractionRate(10000);
  irSampler.init();

  o2::vertexing::PVertexer vertexer;
  vertexer.setValidateWithIR(kFALSE);
  vertexer.setBunchFilling(irSampler.getBunchFilling());
  vertexer.init();
  
  //vertexer QA
  TFile *fQA = new TFile("vertexer.qa.root", "RECREATE");
  TTree *fTreeQA = new TTree("fTreeQA", "");
  
  Float_t fPVx, fPVy, fPVz;
  Float_t fPVxMC=0, fPVyMC=0, fPVzMC=0;
  
  Int_t fNContributors;
  
  fTreeQA->Branch("fPVx", &fPVx, "fPVx/F");
  fTreeQA->Branch("fPVy", &fPVy, "fPVy/F");
  fTreeQA->Branch("fPVz", &fPVz, "fPVz/F");
  fTreeQA->Branch("fPVxMC", &fPVxMC, "fPVxMC/F");
  fTreeQA->Branch("fPVyMC", &fPVyMC, "fPVyMC/F");
  fTreeQA->Branch("fPVzMC", &fPVzMC, "fPVzMC/F");
  
  fTreeQA->Branch("fNContributors", &fNContributors, "fNContributors/I");


  for (int iEvent{0}; iEvent < itsHits.GetEntriesFast(); ++iEvent) {
    fXiDecayRadius = 0;
    fV0DecayRadius = 0;
    fXiDecayRadiusMC = 0;
    fV0DecayRadiusMC = 0;
    std::cout << "*************** Event " << iEvent << " ***************" << std::endl;
    itsHits.GetEntry(iEvent);
    mcTree.GetEvent(iEvent);
    o2::its::ROframe event{iEvent, 12};
    hEventCounter->Fill(0.5);
    
    //Locate Xi decay point
    Int_t lDirectTrackingNLayers = 0;
    Int_t lNLayersLambda = 0;
    cout<<"Checking parenthood - mcArr size "<<mcArr->size()<<endl;
    
    //Check type of MC
    fPtMCXiC = -1;
    fPtMCXiCC = -1;
    
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
      
//      if(part.GetPdgCode()==4232) fPtMCXiC = std::hypot(part.Px(),part.Py());
//      if(part.GetPdgCode()==4422) fPtMCXiCC = std::hypot(part.Px(),part.Py());
//      if(part.GetPdgCode()==4232){
//        fPMCXiC = std::hypot(part.Px(),part.Py(),part.Pz());
//        hXiCGeneratedPt->Fill(fPtMCXiC);
//        hPtEtaGeneratedXiC->Fill(fPtMCXiC, part.GetEta());
//        hPtYGeneratedXiC->Fill(fPtMCXiC, part.GetRapidity());
//      }
//      if(part.GetPdgCode()==4422){
//        fPMCXiCC = std::hypot(part.Px(),part.Py(),part.Pz());
//        hXiCCGeneratedPt->Fill(fPtMCXiCC);
//        hPtEtaGeneratedXiCC->Fill(fPtMCXiCC, part.GetEta());
//        hPtYGeneratedXiCC->Fill(fPtMCXiCC, part.GetRapidity());
//      }
//      if(part.GetPdgCode()==3312){
//        hXiGeneratedPt->Fill(part.GetPt());
//        hPtEtaGeneratedXi->Fill(part.GetPt(), part.GetEta());
//        hPtYGeneratedXi->Fill(part.GetPt(), part.GetRapidity());
//      }
    }
    hQAMCNch->Fill(lNchThisEvent);
    
    int id{0};
    std::map<int, particle> mapPDG;
    Long_t NHitAll=0, NHitInner=0, NHitBachelor=0, NHitV0=0 ;
    std::cout << "*- Filling hits" << std::endl;
    
    //Reset single event counters, please
    for(Int_t ii=0; ii<12; ii++) hHitDensityCounter[ii] -> Reset();
    
    for (auto& hit : *hits) {
      int layer{hit.GetDetectorID()};
      float xyz[3]{hit.GetX(), hit.GetY(), hit.GetZ()};
      float r{std::hypot(xyz[0], xyz[1])};
      float phi{std::atan2(-xyz[1], -xyz[0]) + o2::its::constants::math::Pi};
      
      //if you see radius + epsilon, it's still the N-th layer... likely a bug
      if(r>99.0&&r<101){
        //std::cout << "*- Exception caught at a radius of "<< r << std::endl;
        layer = 11;
      }
      if(layer<0 || layer>11){
        std::cout << "*- Layer "<< layer << ": no such thing! What are you doing, mister?"<< std::endl;
      }
      
      float lSmearing = lIBresolution;
      if(layer>3) lSmearing = lOBresolution;
      
      if (kUseSmearing) {
        phi = gRandom->Gaus(phi, std::asin(lSmearing / r));
        xyz[0] = r * std::cos(phi);
        xyz[1] = r * std::sin(phi);
        xyz[2] = gRandom->Gaus(xyz[2], lSmearing);
      }
      
      hHitDensityCounter[layer] -> Fill(xyz[2]);//Fill with Z coordinate of this hit, please
      
      NHitAll++;
      if(NHitAll%1000==0) std::cout << "*- Hits filled: "<< NHitAll << std::endl;
      event.addTrackingFrameInfoToLayer(layer, xyz[0], xyz[1], xyz[2], r, phi, std::array<float, 2>{0.f, xyz[2]},
                                        std::array<float, 3>{lSmearing * lSmearing, 0.f, lSmearing * lSmearing});
      event.addClusterToLayer(layer, xyz[0], xyz[1], xyz[2], event.getClustersOnLayer(layer).size());
      event.addClusterLabelToLayer(layer, o2::MCCompLabel(hit.GetTrackID(), iEvent, iEvent, false));
      event.addClusterExternalIndexToLayer(layer, id++);
      
      if (mapPDG.find(hit.GetTrackID()) == mapPDG.end() && hit.GetTrackID()>0) {
        mapPDG[hit.GetTrackID()] = particle();
        mapPDG[hit.GetTrackID()].nLayers |= 1 << layer;
        mapPDG[hit.GetTrackID()].pt = std::hypot(hit.GetPx(), hit.GetPy());
        if (hit.GetTrackID() < mcArr->size()) {
          auto part = mcArr->at(hit.GetTrackID());
          mapPDG[hit.GetTrackID()].energyFirst = part.GetEnergy();
          mapPDG[hit.GetTrackID()].pdg = part.GetPdgCode();
          mapPDG[hit.GetTrackID()].pt = part.GetPt();
          mapPDG[hit.GetTrackID()].eta = part.GetEta();
          mapPDG[hit.GetTrackID()].phi = part.GetPhi();
          if(abs(part.getMotherTrackId()) < mcArr->size()){
            auto partmother = mcArr->at( abs(part.getMotherTrackId()) );
            mapPDG[hit.GetTrackID()].pdgmother = partmother.GetPdgCode();
            mapPDG[hit.GetTrackID()].motherID = part.getMotherTrackId();
          }
          mapPDG[hit.GetTrackID()].vX = part.Vx();
          mapPDG[hit.GetTrackID()].vY = part.Vy();
          mapPDG[hit.GetTrackID()].vZ = part.Vz();
        }
      } else {
        mapPDG[hit.GetTrackID()].nLayers |= 1 << layer;
        mapPDG[hit.GetTrackID()].energyLast = hit.GetE();
      }
      //std::cout << "*- Event " << iEvent << " hit.GetTrackID() = " <<hit.GetTrackID() << " endhit"<<std::endl;
    }
    
    for(Int_t layer=0; layer<12; layer++){
      for(Int_t iz=1; iz<hHitDensityCounter[layer]->GetNbinsX()+1; iz++){
        Double_t lHitDensity = hHitDensityCounter[layer]->GetBinContent(iz);
        Double_t lHitPosition = hHitDensityCounter[layer]->GetBinCenter(iz);
        hHitDensityMap[layer] -> Fill(lHitPosition, lHitDensity);
        hHitDensityProfile[layer] -> Fill(lHitPosition, lHitDensity);
      }
    }
    
    roFrame = iEvent;
    //std::cout << "*- Event " << iEvent << " MC loop" << std::endl;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //LUT APPROACH STEP 1: DO ALL PRIMARY TRACKS
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    //Perfect MC vertex from mcHeader
    
    fPVxMC = mcHead->GetX();
    fPVyMC = mcHead->GetY();
    fPVzMC = mcHead->GetZ();

    o2::math_utils::Point3D<float> zeropos{0,0,0};
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
    lWatchSmearing.Start(kFALSE);
    Long_t lParticlesForVertexing = 0;
    // loop over particles
    cout<<"Allocating TOF info arrays..."<<flush;
    TArrayF lInnerTOFSignalLutted(particles.size());
    TArrayF lOuterTOFSignalLutted(particles.size());
    cout<<" done"<<endl;

    o2::InteractionRecord ir = irSampler.generateCollisionTime();
    Long_t lCorruptedTracks=0;
    
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
      if(lX0>-99.&&lX1>-99.) lThisTrackLength100 = TrackLength(o2track, lX0, lX1, lMagneticField);
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
      if(lX0>-99.&&lX1>-99.) lThisTrackLength20 = TrackLength(o2track, lX0, lX1, lMagneticField);
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
      
      if (TMath::IsNaN(o2track.getZ())) {
        lCorruptedTracks++;
        continue;
      }
      hPartialTransportDebug->Fill(2.5);
      
      //Feed vertex finder
      const float t = (ir.bc2ns() + gRandom->Gaus(0., 100.)) * 1e-3;
      smearedTracks.push_back(TrackAlice3{o2track, t, 100.f * 1e-3, TMath::Abs(alabel)});
      lParticlesForVertexing++;
    }
    cout<<"Corrupted tracks: "<<lCorruptedTracks<<endl; 
    cout<<"Tracks to be used for vertex finding: "<<lParticlesForVertexing<<endl;
    cout<<"PV finding procedure started"<<endl;
    
    std::vector<o2::MCCompLabel> lblTracks;
    std::vector<o2::vertexing::PVertex> vertices;
    std::vector<o2::vertexing::GIndex> vertexTrackIDs;
    std::vector<o2::vertexing::V2TRef> v2tRefs;
    std::vector<o2::MCEventLabel> lblVtx;
    lblVtx.emplace_back(iEvent, 1);
    std::vector<o2::dataformats::GlobalTrackID> idxVec; // here we will the global IDs of all used tracks
    idxVec.reserve(smearedTracks.size());
    static int cntEv = 0;
    for (unsigned i = 0; i < smearedTracks.size(); i++) {
      lblTracks.emplace_back(smearedTracks[i].mLabel, iEvent, 1, false);
      idxVec.emplace_back(i, o2::dataformats::GlobalTrackID::ITS);
      if (TMath::IsNaN(smearedTracks[i].getZ())) {
	LOG(INFO) << "Corrupted track: " << smearedTracks[i].asString();
	continue;
      }           
      //LOG(INFO) << "trc " << i << " " << cntEv << " " << smearedTracks[i].getZ() << " " << smearedTracks[i].timeEst.getTimeStamp() << " " << smearedTracks[i].timeEst.getTimeStampError();
    }

    const int n_vertices = vertexer.process(smearedTracks,
                                            idxVec,
                                            gsl::span<o2::InteractionRecord>{bcData},
                                            vertices,
                                            vertexTrackIDs,
                                            v2tRefs,
                                            gsl::span<const o2::MCCompLabel>{lblTracks},
                                            lblVtx);


    cout<<"PV finding procedure finished, vertices found: "<<n_vertices<<endl;
    
    int cntv = 0;
    for ( auto vert : vertices)  {
      LOG(INFO) << "PV: " << cntv++ << " " << cntEv << " " << vert.getX()<<" "<<vert.getY()<<" "<<vert.getZ()<< " " << vert.getNContributors() << " " << vert.getTimeStamp().getTimeStamp() << " " << vert.getTimeStamp().getTimeStampError() << endl;
    }
    cntEv++;
    

    if(n_vertices<1){
      fPVx = -100;
      fPVy = -100;
      fPVz = -100;
      fNContributors=0;
      continue; 
    }else{
      fPVx = vertices[0].getX();
      fPVy = vertices[0].getY();
      fPVz = vertices[0].getZ();
      fNContributors = vertices[0].getNContributors();
      cout<<"Main PV at: "<<vertices[0].getX()<<", "<<vertices[0].getY()<<", "<<vertices[0].getZ()<<endl;
      hVertexZ->Fill( vertices[0].getZ() );
    }
    fTreeQA->Fill();
  
    
    if( TMath::Abs(vertices[0].getZ()) > 10 ) continue; //skip outside center of detector, please
    
    for (Long_t iii=0; iii< mcArr->size(); iii++ ){
      auto part = mcArr->at(iii);
      
      if(part.GetPdgCode()==4232) fPtMCXiC = std::hypot(part.Px(),part.Py());
      if(part.GetPdgCode()==4422) fPtMCXiCC = std::hypot(part.Px(),part.Py());
      if(part.GetPdgCode()==4232){
        fPMCXiC = std::hypot(part.Px(),part.Py(),part.Pz());
        hXiCGeneratedPt->Fill(fPtMCXiC);
        hPtEtaGeneratedXiC->Fill(fPtMCXiC, part.GetEta());
        hPtYGeneratedXiC->Fill(fPtMCXiC, part.GetRapidity());
      }
      if(part.GetPdgCode()==4422){
        fPMCXiCC = std::hypot(part.Px(),part.Py(),part.Pz());
        hXiCCGeneratedPt->Fill(fPtMCXiCC);
        hPtEtaGeneratedXiCC->Fill(fPtMCXiCC, part.GetEta());
        hPtYGeneratedXiCC->Fill(fPtMCXiCC, part.GetRapidity());
      }
      if(part.GetPdgCode()==3312){
        hXiGeneratedPt->Fill(part.GetPt());
        hPtEtaGeneratedXi->Fill(part.GetPt(), part.GetEta());
        hPtYGeneratedXi->Fill(part.GetPt(), part.GetRapidity());
      }
    }
    
    o2::math_utils::Point3D<float> pos{vertices[0].getX(),vertices[0].getY(),vertices[0].getZ()};
    std::array<float, 6> cov;
    for(Int_t jj=0; jj<6; jj++) cov[jj]=vertices[0].getCov()[jj];
    o2::dataformats::VertexBase vtx(pos, cov);
    
//    o2::math_utils::Point3D<float> pos{0,0,0};
//    std::array<float, 6> cov;
//    for(Int_t jj=0; jj<6; jj++) cov[jj]=1e-6;
//    o2::dataformats::VertexBase vtx(pos, cov);
//    o2::dataformats::DCA dca;
//
    hEventCounterWithVertex->Fill(0.5);
//
//    fPrimaryVertexX = vertices[0].getX();
//    fPrimaryVertexY = vertices[0].getY();
//    fPrimaryVertexZ = vertices[0].getZ();
    lWatchSmearing.Stop();
    lWatchTracking.Start(kFALSE);
    
    //event.addPrimaryVertex(pos[0].getX(), pos[0].getY(), pos[0].getZ());
    event.addPrimaryVertex(vertices[0].getX(),vertices[0].getY(),vertices[0].getZ());
    
    std::cout << "*- Event " << iEvent << " tracking" << std::endl;
    tracker.clustersToTracks(event);
    auto& tracks = tracker.getTracks();
    auto& tracksLabels = tracker.getTrackLabels();
    std::cout << "*- Event " << iEvent << " done tracking!" << std::endl;
    
    lWatchTracking.Stop();
    
    lWatchAnalysis.Start(kFALSE);
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // Now merge track lists and label lists
    cout<<"Initial track array size: "<<tracks.size()<<endl;
    cout<<"Initial label array size: "<<tracksLabels.size()<<endl;
    hPartialTransportDebug->Fill(3.5, tracks.size());
    
    for(Int_t iTrack = 0; iTrack<smearedTracks.size(); iTrack++){
      o2::its::TrackITS lConvertedTrack = o2::its::TrackITS( smearedTracks[iTrack] );
      tracks.emplace_back(lConvertedTrack);
      tracksLabels.emplace_back(smearedTracks[iTrack].mLabel, iEvent, 1, false);
    }
    cout<<"Merged track array size: "<<tracks.size()<<endl;
    cout<<"Merged label array size: "<<tracksLabels.size()<<endl;
    
    hQANchCorrelation->Fill(tracks.size(), lNchThisEventAll);
    
    TArrayI negtra(tracks.size());
    TArrayI postra(tracks.size());
    TArrayI bachtra(tracks.size());
    
    TArrayI truenegtra(tracks.size());
    TArrayI truepostra(tracks.size());
    TArrayI truebachtra(tracks.size());
    
    TArrayI pictra(tracks.size());
    TArrayI picctra(tracks.size());
    TArrayI truepic(tracks.size());
    TArrayI truepicfromxicc(tracks.size());
    TArrayI truepicc(tracks.size());
    Long_t nneg=0, npos=0, nbach=0, npic=0, npicc=0;
    
    hNTracks -> Fill(tracks.size());
    fNTracks = tracks.size();
    
    TArrayF DCAxytoPVArray(tracks.size()); //store all DCAs to PV in memory, please
    TArrayF DCAztoPVArray(tracks.size()); //store all DCAs to PV in memory, please
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Loop over tracks, attempt to determine TOF information
    
    //Monte Carlo: perfect signal
    TArrayF lInnerTOFSignalMC(tracks.size());
    TArrayF lOuterTOFSignalMC(tracks.size());
    
    //Reconstructed: smeared signal
    TArrayF lInnerTOFSignal20(tracks.size());
    TArrayF lInnerTOFSignal50(tracks.size());
    TArrayF lOuterTOFSignal(tracks.size());
    TArrayF lTrackLength(tracks.size());
    TArrayF lTrackShortLength(tracks.size());
    
    cout<<"Determining TOF values for inner and outer TOF..."<<endl;
    for (unsigned int i{0}; i < tracks.size(); ++i) {
      auto& lab = tracksLabels[i];
      auto& track = tracks[i];
      int lMCID = lab.getTrackID();
      
      lInnerTOFSignal20[i] = -1;
      lInnerTOFSignal50[i] = -1;
      lOuterTOFSignal[i] = -1;
      lInnerTOFSignalMC[i] = -1;
      lOuterTOFSignalMC[i] = -1;
      
      //Determine length wrt to PV for posterior use
      Float_t lX0=-100, lX1=-100;
      Float_t lThisTrackLength = -100;
      if (!track.propagateToDCA(vtx, tracker.getBz())) {
        //std::cout << "Propagation failed." << std::endl;
      } else {
        lX0 = track.getX();
      }
      if (!track.getXatLabR(100.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
        //std::cout << "Outward propagation failed." << std::endl;
      }
      if(lX0>-99.&&lX1>-99.) lThisTrackLength = TrackLength(track, lX0, lX1, lMagneticField);
      lTrackLength[i] = lThisTrackLength;
      track.propagateTo(lX0, lMagneticField);
      
      lX0=-100, lX1=-100;
      lThisTrackLength = -100;
      if (!track.propagateToDCA(vtx, tracker.getBz())) {
        //std::cout << "Propagation failed." << std::endl;
      } else {
        lX0 = track.getX();
      }
      if (!track.getXatLabR(20.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
        //std::cout << "Outward propagation failed." << std::endl;
      }
      if(lX0>-99.&&lX1>-99.) lThisTrackLength = TrackLength(track, lX0, lX1, lMagneticField);
      lTrackShortLength[i] = lThisTrackLength;
      track.propagateTo(lX0, lMagneticField);
      
      if( lMCID>=0 ) {
        Double_t lLowestTimeOuter = 1e+6;
        Double_t lLowestTimeInner = 1e+6;
        
        for (auto& hit : *hits) {
          if(hit.GetTrackID()!=lMCID) continue;
          
          int layer{hit.GetDetectorID()};
          float xyz[3]{hit.GetX(), hit.GetY(), hit.GetZ()};
          float r{std::hypot(xyz[0], xyz[1])};
          float phi{std::atan2(-xyz[1], -xyz[0]) + o2::its::constants::math::Pi};
          if (r > 99.0 && r < 101) layer = 11;
          //constexpr float LightSpeedCm2S = 299792458.e2;           // C in cm/s
          if( layer == 11 )
            if( 1e+12*hit.GetTime() < lLowestTimeOuter ) lLowestTimeOuter = 1e+12*hit.GetTime();
          if( layer == 6 )
            if( 1e+12*hit.GetTime() < lLowestTimeInner ) lLowestTimeInner = 1e+12*hit.GetTime();
        }
        if(lLowestTimeInner<1e+6-1){
          lInnerTOFSignal20[i] = gRandom->Gaus(lLowestTimeInner, loTOFresolution);
          lInnerTOFSignal50[i] = gRandom->Gaus(lLowestTimeInner, liTOFresolution);
          lInnerTOFSignalMC[i] = lLowestTimeInner;
        }
        if(lLowestTimeOuter<1e+6-1){
          lOuterTOFSignal[i] = gRandom->Gaus(lLowestTimeOuter, loTOFresolution);
          lOuterTOFSignalMC[i] = lLowestTimeOuter;
        }
        
        if(lInnerTOFSignal20[i]<-0.9)
          if(lInnerTOFSignalLutted[lMCID] > -50) lInnerTOFSignal20[i] = lInnerTOFSignalLutted[lMCID];
        if(lOuterTOFSignal[i]<-0.9)
          if(lOuterTOFSignalLutted[lMCID] > -50) lOuterTOFSignal[i] = lOuterTOFSignalLutted[lMCID];
      }
    }
    cout<<"TOF values determined!"<<endl;
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Loop over found tracks, attempt to find V0 candidates
    // simple reset
    
    fPtXi=0;
    fXicPt=0;
    
    fTrueXic = kFALSE;
    
    fXiDCAxyToPV=0;
    fXiDCAzToPV=0;
    fXicDCAxyToPV=0;
    fXicDCAzToPV=0;
    
    fXiHitsAdded=0;
    fXicMass=0;
    
    fNLongTracks = 0;
    fNLongPrimaryTracks = 0;
    fNPrimaryTracks = 0;
    
    Long_t lLongTrackRequirement = -1;
    
    for (unsigned int i{0}; i < tracks.size(); ++i) {
      auto& lab = tracksLabels[i];
      auto& track = tracks[i];
      int trackID = lab.getTrackID();
      
      Int_t lPos, lNeg;
      if( Particle == 3122){
        lPos = 2212; lNeg = -211;
      }
      if( Particle == 310){
        lPos = 211; lNeg= -211;
      }
      Float_t lThisDCAxy=0.0, lThisDCAz=0.0;
      DCAxytoPVArray[i] = 1e+10;
      DCAztoPVArray[i] = 1e+10;
      
      if (!track.propagateToDCA(vtx, tracker.getBz(), &dca)) {
        std::cout << "Propagation failed." << std::endl;
        lThisDCAxy = 1e+10;
        lThisDCAz = 1e+10;
      } else {
        lThisDCAxy = dca.getY();
        lThisDCAz = dca.getZ();
        DCAxytoPVArray[i] = lThisDCAxy;
        DCAztoPVArray[i] = lThisDCAz;
      }
      
      hQADCAxy->Fill(1e+4*lThisDCAxy);
      hQADCAz->Fill(1e+4*lThisDCAz);
      
      fNLongTracks++;
      if( TMath::Abs(1e+4*lThisDCAxy)<25 && TMath::Abs(1e+4*lThisDCAz)<25) fNLongPrimaryTracks++;
      if( TMath::Abs(1e+4*lThisDCAxy)<25 && TMath::Abs(1e+4*lThisDCAz)<25) fNPrimaryTracks++;
      
      Bool_t lIsPos = kFALSE, lIsPosPDG = kFALSE;
      Bool_t lIsNeg = kFALSE, lIsNegPDG = kFALSE;
      
      if(track.getSign()>0) lIsPos = kTRUE;
      if(track.getSign()<0) lIsNeg = kTRUE;
      Bool_t lIsMotherDesired = kFALSE;
      Bool_t lIsGrandMotherDesired = kFALSE;
      if(trackID>=0 && trackID<mcArr->size()){
        auto part = mcArr->at(trackID);
        if( part.GetPdgCode() == lPos ) lIsPosPDG = kTRUE;
        if( part.GetPdgCode() == lNeg ) lIsNegPDG = kTRUE;
        Int_t lMotherTrackID=part.getMotherTrackId();
        if(lMotherTrackID>=0 && lMotherTrackID<mcArr->size()){
          auto partmom = mcArr->at(lMotherTrackID);
          if(partmom.GetPdgCode() == Particle) lIsMotherDesired = kTRUE;
          auto partgrandmomidx = partmom.getMotherTrackId();
          if(partgrandmomidx>=0){
            auto partgrandmom = mcArr->at(partgrandmomidx);
            if(partgrandmom.GetPdgCode() == 3312) lIsGrandMotherDesired = kTRUE;
          }
        }
      }
      
      if( lIsPos && ((lIsMotherDesired&&lIsGrandMotherDesired&&lIsPosPDG)||lAssociateXiDaughters==kFALSE) && ((TMath::Abs(lThisDCAxy) > lMinDCAxyPositive && TMath::Abs(lThisDCAz) > lMinDCAzPositive && TMath::Abs(lThisDCAxy)<lMaxDCAxyPositive) || lMasterSwitch<0 ) ) {
        truepostra[npos] = 0;
        if(lIsMotherDesired) truepostra[npos] = 1;
        postra[npos++]=i;
      }
      if( lIsNeg && ((lIsMotherDesired&&lIsGrandMotherDesired&&lIsNegPDG)||lAssociateXiDaughters==kFALSE) && ((TMath::Abs(lThisDCAxy) > lMinDCAxyNegative && TMath::Abs(lThisDCAz) > lMinDCAzNegative)|| lMasterSwitch<0) )   {
        truenegtra[nneg] = 0;
        if(lIsMotherDesired) truenegtra[nneg] = 1;
        negtra[nneg++]=i;
      }
      
      //Bachelor stuff
      Bool_t lIsBach = kFALSE;
      Bool_t lIsBachPDG = kFALSE;
      if( track.getCharge()<0) lIsBach = kTRUE ; // correct charge!
      lIsMotherDesired = kFALSE;
      if(trackID>=0 && trackID<mcArr->size()){
        auto part = mcArr->at(trackID);
        if( part.GetPdgCode() == -211 ) lIsBachPDG = kTRUE;
        Int_t lMotherTrackID=part.getMotherTrackId();
        if(lMotherTrackID>=0 && lMotherTrackID<mcArr->size()){
          auto partmom = mcArr->at(lMotherTrackID);
          if(partmom.GetPdgCode() == 3312) lIsMotherDesired = kTRUE;
        }
      }
      
      if( lIsBach && ((lIsMotherDesired&&lIsBachPDG)||lAssociateXiDaughters==kFALSE) && ((TMath::Abs(lThisDCAxy) > lMinDCAxyBachelor && TMath::Abs(lThisDCAz) > lMinDCAzBachelor)|| lMasterSwitch<0) ) {
        truebachtra[nbach] = 0;
        if(lIsMotherDesired) truebachtra[nbach] = 1;
        bachtra[nbach++]=i;
      }
      
      //HF pion stuff
      Bool_t lIsPion = kFALSE;
      Bool_t lIsPionPDG = kFALSE;
      Bool_t lIsFromXic = kFALSE;
      Bool_t lIsFromXicFromXicc = kFALSE;
      Bool_t lIsFromXicc = kFALSE;
      if(trackID>=0 && trackID<mcArr->size()){
        auto part = mcArr->at(trackID);
        if( part.GetPdgCode() == 211 ) lIsPionPDG = kTRUE;
        if( track.getCharge() > 0 ) lIsPion = kTRUE; //pion charge
        Int_t lMotherTrackID=part.getMotherTrackId();
        if(lMotherTrackID>=0 && lMotherTrackID<mcArr->size()){
          auto partmom = mcArr->at(lMotherTrackID);
          if(partmom.GetPdgCode() == 4232){
            lIsFromXic = kTRUE;
            //go back once more, if possible
            Int_t lGrandMotherID = partmom.getMotherTrackId();
            if( lGrandMotherID >= 0 ){
              auto partGrandMom = mcArr->at(lGrandMotherID);
              if(partGrandMom.GetPdgCode() == 4422) lIsFromXicFromXicc = kTRUE;
            }
          }
          if(partmom.GetPdgCode() == 4422) lIsFromXicc = kTRUE;
        }
      }
      
      //TOF check
      Double_t lPionExpectedSignal = lTrackShortLength[ i ] / Velocity(track.getP(), 0.139570);
      if( TMath::Abs(lPionExpectedSignal - lInnerTOFSignal20[i] ) > lMaxTimeOffsetHFPions && lMasterSwitch > 0 ) continue;
      
      if(lIsPion==kFALSE) continue; //skip if wrong charge
      
      if( TMath::Abs(lThisDCAxy) > lMinDCAxyHFPions && TMath::Abs(lThisDCAz) > lMinDCAzHFPions ) {
        if(!lAssociateXic||(lIsFromXic&&lIsPionPDG)){
          pictra[npic]=i;
          truepic[npic] = 0;
          truepicfromxicc[npic] = 0;
          if (lIsFromXic) truepic[npic] = 1; //flag correct ones (must be two in this case)
          if (lIsFromXicFromXicc) truepicfromxicc[npic] = 1; //flag correct ones (must be two in this case)
          npic++;
        }
        if(!lAssociateXicc||(lIsFromXicc&&lIsPionPDG)){
          picctra[npicc]=i;
          truepicc[npicc] = 0;
          if (lIsFromXicc) truepicc[npicc] = 1; //flag correct ones (must be two in this case)
          npicc++;
        }
      }
    }
    hNLongTracks -> Fill(fNLongTracks);
    hNLongPrimaryTracks -> Fill(fNLongPrimaryTracks);
    hNPrimaryTracks -> Fill(fNPrimaryTracks);
    
    //Fill generated vs multiplicity
    for (Long_t iii=0; iii< mcArr->size(); iii++ ){
      auto part = mcArr->at(iii);
      Float_t lThisPt = std::hypot(part.Px(),part.Py());
      
      if(part.GetPdgCode()==4232){
        if(TMath::Abs(part.GetEta())<1.5) hPtMultGeneratedXiC_AllEta -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetEta())<1.0) hPtMultGeneratedXiC_Eta10 -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetEta())<0.5) hPtMultGeneratedXiC_Eta5 -> Fill( lThisPt, fNLongTracks);
        
        if(TMath::Abs(part.GetRapidity())<1.5) hPtMultGeneratedXiC_AllY -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetRapidity())<1.0) hPtMultGeneratedXiC_Y10 -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetRapidity())<0.5) hPtMultGeneratedXiC_Y5 -> Fill( lThisPt, fNLongTracks);
      }
      if(part.GetPdgCode()==4422){
        if(TMath::Abs(part.GetEta())<1.5) hPtMultGeneratedXiCC_AllEta -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetEta())<1.0) hPtMultGeneratedXiCC_Eta10 -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetEta())<0.5) hPtMultGeneratedXiCC_Eta5 -> Fill( lThisPt, fNLongTracks);
        
        if(TMath::Abs(part.GetRapidity())<1.5) hPtMultGeneratedXiCC_AllY -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetRapidity())<1.0) hPtMultGeneratedXiCC_Y10 -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetRapidity())<0.5) hPtMultGeneratedXiCC_Y5 -> Fill( lThisPt, fNLongTracks);
      }
      if(part.GetPdgCode()==3312){
        if(TMath::Abs(part.GetEta())<1.5) hPtMultGeneratedXi_AllEta -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetEta())<1.0) hPtMultGeneratedXi_Eta10 -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetEta())<0.5) hPtMultGeneratedXi_Eta5 -> Fill( lThisPt, fNLongTracks);
        
        if(TMath::Abs(part.GetRapidity())<1.5) hPtMultGeneratedXi_AllY -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetRapidity())<1.0) hPtMultGeneratedXi_Y10 -> Fill( lThisPt, fNLongTracks);
        if(TMath::Abs(part.GetRapidity())<0.5) hPtMultGeneratedXi_Y5 -> Fill( lThisPt, fNLongTracks);
      }
    }
  
    hTrackCount->Fill(0.5, npos);
    hTrackCount->Fill(1.5, nneg);
    hTrackCount->Fill(2.5, nbach);
    hTrackCount->Fill(3.5, npic);
    hTrackCount->Fill(4.5, npicc);
    
    cout<<"Valid tracks: "<<npos<<", "<<nneg<<", "<<nbach<<", "<<npic<<", "<<npicc<<endl;
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    for (unsigned int i{0}; i < npos; ++i) {
      for (unsigned int j{0}; j < nneg; ++j) {
        //Try to progate to dca
        hCombinatorics->Fill(0.5);
        auto Track1 = tracks[postra[i]];
        auto Track2 = tracks[negtra[j]];
        
        //Check same mother
        auto& labi = tracksLabels[postra[i]];
        int trackIDi = labi.getTrackID();
        auto& labj = tracksLabels[negtra[j]];
        int trackIDj = labj.getTrackID();
        
        //De-reference from memory, no need to propagate like crazy now
        fPositiveDCAxy = DCAxytoPVArray[postra[i]] * 1.e4;
        fPositiveDCAz = DCAztoPVArray[postra[i]] * 1.e4;
        fNegativeDCAxy = DCAxytoPVArray[negtra[j]] * 1.e4;
        fNegativeDCAz = DCAztoPVArray[negtra[j]] * 1.e4;
        
        fPositiveClusters = Track1.getNumberOfClusters();
        fNegativeClusters = Track2.getNumberOfClusters();
        fPositiveChisquare = Track1.getChi2();
        fNegativeChisquare = Track2.getChi2();
        
        //Not same mother -> discard
        //if( trackIDi >= mcArr->size()) continue;
        //if( trackIDj >= mcArr->size()) continue;
        //auto parti = mcArr->at(trackIDi);
        //auto partj = mcArr->at(trackIDj);
        //if( abs(parti.getMotherTrackId()) != ) continue;
        
        hCombinatoricsV0->Fill(1.5);
        //Get mother momentum from stack
        auto part1 = mcArr->at( trackIDi );
        auto part2 = mcArr->at( trackIDj );
        
        fPositivePDG = part1.GetPdgCode();
        fNegativePDG = part2.GetPdgCode();
        
        Int_t lCascadeID = -1;
        Bool_t lSameMotherPosNeg = kFALSE;
        fPtMCMother = -1;
        fPtMCXi = -1;
        Int_t lCascadeMotherID = -1;
        Int_t lCascadePDGcheck = -1;
        Bool_t lCascadeFromXiCC = kFALSE;
        if( part1.getMotherTrackId()>=0 && part2.getMotherTrackId()>=0 && part1.getMotherTrackId()<mcArr->size() ){
          if( part1.getMotherTrackId() == part2.getMotherTrackId() ) lSameMotherPosNeg = kTRUE;
          auto part3 = mcArr->at( abs(part1.getMotherTrackId()) );
          fPtMCMother = TMath::Sqrt(TMath::Power(part3.Px(),2) + TMath::Power(part3.Py(),2));
          lCascadeID = part3.getMotherTrackId();
          if( lCascadeID >= 0 ){
            auto partCasc = mcArr->at( lCascadeID );
            lCascadePDGcheck = partCasc.GetPdgCode();
            fPtMCXi = partCasc.GetPt();
            lCascadeMotherID = partCasc.getMotherTrackId();
            if( lCascadeMotherID >= 0 ){
              auto partCascMother = mcArr->at( lCascadeMotherID );
              if( partCascMother.GetPdgCode() == 4232 && partCascMother.getMotherTrackId() >= 0){
                auto partCascGrandMother = mcArr->at( partCascMother.getMotherTrackId() );
                if( partCascGrandMother.GetPdgCode() == 4422 ) lCascadeFromXiCC = kTRUE;
              }
            }
          }
        }
        fTrueNegative = kFALSE;
        fTruePositive = kFALSE;
        
        if( lCascadeFromXiCC && truepostra[j] == 1 ) fTruePositive = kTRUE;
        if( lCascadeFromXiCC && truenegtra[j] == 1 ) fTrueNegative = kTRUE;
        
        if(lCascadeID<0&&lAssociateXiDaughters) continue;
        if(lSameMotherPosNeg==kFALSE&&lAssociateXiDaughters) continue;
        if(lCascadePDGcheck!=3312&&lAssociateXiDaughters) continue;
        fV0DecayRadiusMC=-1;
        if( lSameMotherPosNeg == kTRUE ) fV0DecayRadiusMC = std::hypot(part1.Vx(),part1.Vy());
        
        hCombinatoricsV0->Fill(2.5);
        
        fPtMC = TMath::Sqrt(TMath::Power(part1.Px()+part2.Px(),2) + TMath::Power(part1.Py()+part2.Py(),2));
        //pT of V0 via actual V0
        
        int nCand = fitterV0.process(Track1, Track2);
        if (nCand == 0) {
          //std::cout<<"Crap candidate"<<std::endl;
          continue;
        }
        hCombinatoricsV0->Fill(3.5);
        
        float thisdcav0dau = TMath::Sqrt(fitterV0.getChi2AtPCACandidate(0));
        fDCAV0Dau = thisdcav0dau * 1.e4; //in microns!
        if(fDCAV0Dau > lMaxDCAV0Daughters && lMasterSwitch > 0) continue;
        fV0DauDCA = fDCAV0Dau ;
        
        const auto& wvtx = fitterV0.getPCACandidate();
        
        std::array<float, 3> weakdecaypos = {0.};
        std::array<float, 3> pvec0;
        std::array<float, 3> pvec1;
        for (int i = 0; i < 3; i++) {
          weakdecaypos[i] = wvtx[i];
        }
        std::array<float, 3> Track1Pos;
        std::array<float, 3> Track2Pos;
        fitterV0.getTrack(0).getXYZGlo(Track1Pos);
        fitterV0.getTrack(1).getXYZGlo(Track2Pos);
        
        fitterV0.getTrack(0).getPxPyPzGlo(pvec0);
        fitterV0.getTrack(1).getPxPyPzGlo(pvec1);
        
        auto thisv0cospa = RecoDecay::CPA(array{vertices[0].getX(),vertices[0].getY(),vertices[0].getZ()},
                                          array{weakdecaypos[0], weakdecaypos[1], weakdecaypos[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});
        
        //Fiducial: min radius
        auto thisv0radius = std::hypot(weakdecaypos[0], weakdecaypos[1]);
        fV0DecayRadius=thisv0radius;
        
        if( fV0DecayRadius < lMinV0Radius && lMasterSwitch > 0) continue;
        
        Double_t lMassPos, lMassNeg;
        if(Particle==310){
          lMassPos = RecoDecay::getMassPDG(kPiPlus);
          lMassNeg = RecoDecay::getMassPDG(kPiPlus);
        }
        if(Particle==3122){
          lMassPos = RecoDecay::getMassPDG(kProton);
          lMassNeg = RecoDecay::getMassPDG(kPiPlus);
        }
        
        auto massk0 = RecoDecay::M(array{array{pvec0[0], pvec0[1], pvec0[2]}, array{pvec1[0], pvec1[1], pvec1[2]}}, array{lMassPos, lMassNeg});
        fMass = massk0;
        fLambdaMass = fMass;
        
        if( TMath::Abs(fLambdaMass-1.116)> lMassWindowLambda && lMasterSwitch > 0 ) continue; //only look at peak, please
        
        fPt = TMath::Sqrt(TMath::Power(pvec0[0]+pvec1[0],2) + TMath::Power(pvec0[1]+pvec1[1],2) );
        
        Double_t lExpMass=1.115683;
        if(Particle==310) lExpMass = 0.497;
        
        fPA = TMath::ACos(thisv0cospa);
        
        //Now attempt cascade finding, please
        for (unsigned int k{0}; k < nbach; ++k) {
          if(bachtra[k]==negtra[j]) continue; //don't recycle, please
          auto bTrack = tracks[bachtra[k]];
          auto& labBach = tracksLabels[bachtra[k]];
          int BachID = labBach.getTrackID();
          //std::cout<<"---> Cascade finding attempt! "<< std::endl;
          hCascFinding->Fill(0.5);
          
          fBachelorClusters = bTrack.getNumberOfClusters();
          fBachelorChisquare = bTrack.getChi2();
          
          auto lBachelorPart = mcArr->at( BachID );
          fBachelorPDG = lBachelorPart.GetPdgCode();
          fXiDecayRadiusMC = std::hypot(lBachelorPart.Vx(),lBachelorPart.Vy());
          
          Bool_t lSameMotherBachelorV0 = kFALSE;
          if( lSameMotherPosNeg && BachID>=0 && BachID < mcArr->size() ){
            Int_t lCascadeIDviaBach = lBachelorPart.getMotherTrackId();
            if( lCascadeID == lCascadeIDviaBach && lCascadeID>=0 ) lSameMotherBachelorV0 = kTRUE;
          }
          
          fTrueBachelor=kFALSE;
          if( lSameMotherBachelorV0 && truebachtra[k] == 1 ) fTrueBachelor = kTRUE; //bachelor from XiCC decay chain
          
          if (lAssociateXiDaughters==kTRUE && lSameMotherBachelorV0 == kFALSE) continue; //Skip everything not perfect
          
          fBachelorDCAxy = DCAxytoPVArray[bachtra[k]] * 1.e4;
          fBachelorDCAz = DCAztoPVArray[bachtra[k]] * 1.e4;
          
          //Check if associated
          
          std::array<float, 3> pos = {0.};
          std::array<float, 3> posXi = {0.};
          std::array<float, 3> pvecpos = {0.};
          std::array<float, 3> pvecneg = {0.};
          std::array<float, 3> pvecbach = {0.};
          std::array<float, 21> cov0 = {0};
          std::array<float, 21> cov1 = {0};
          std::array<float, 21> covV0 = {0};
          
          //Covariance matrix calculation
          const int momInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
          fitterV0.getTrack(0).getPxPyPzGlo(pvecpos);
          fitterV0.getTrack(1).getPxPyPzGlo(pvecneg);
          fitterV0.getTrack(0).getCovXYZPxPyPzGlo(cov0);
          fitterV0.getTrack(1).getCovXYZPxPyPzGlo(cov1);
          
          for (int i = 0; i < 6; i++) {
            int j = momInd[i];
            covV0[j] = cov0[j]+cov1[j];
          }
          auto covVtxV0 = fitterV0.calcPCACovMatrix();
          
          covV0[0] = covVtxV0(0, 0);
          covV0[1] = covVtxV0(1, 0);
          covV0[2] = covVtxV0(1, 1);
          covV0[3] = covVtxV0(2, 0);
          covV0[4] = covVtxV0(2, 1);
          covV0[5] = covVtxV0(2, 2);
          
          for (int i = 0; i < 21; i++) {
            //std::cout<<" covV0["<<i<<"] = "<< covV0[i]<< std::endl;
          }
          
          const std::array<float, 3> vertex = {(float)wvtx[0], (float)wvtx[1], (float)wvtx[2]};
          const std::array<float, 3> momentum = {pvecpos[0] + pvecneg[0], pvecpos[1] + pvecneg[1], pvecpos[2] + pvecneg[2]};
          
          fV0TotalMomentum = std::hypot( pvecpos[0] + pvecneg[0], pvecpos[1] + pvecneg[1], pvecpos[2] + pvecneg[2]);
          
          auto tV0 = o2::track::TrackParCov(vertex, momentum, covV0, 0);
          
          o2::math_utils::Point3D<float> posV0Point{vertex[0],vertex[1],vertex[2]};
          std::array<float, 6> covV0vertex;
          covV0vertex[0] = covVtxV0(0, 0);
          covV0vertex[1] = covVtxV0(1, 0);
          covV0vertex[2] = covVtxV0(1, 1);
          covV0vertex[3] = covVtxV0(2, 0);
          covV0vertex[4] = covVtxV0(2, 1);
          covV0vertex[5] = covVtxV0(2, 2);
          o2::dataformats::VertexBase vtxV0(posV0Point, covV0vertex);
          
          Float_t lThisDCAxy=0.0, lThisDCAz=0.0;
          if (!tV0.propagateToDCA(vtx, tracker.getBz(), &dca)) {
            std::cout << "Propagation failed." << std::endl;
            lThisDCAxy = 1e+10;
            lThisDCAz = 1e+10;
          } else {
            lThisDCAxy = dca.getY()*1e+4;
            lThisDCAz = dca.getZ()*1e+4;
          }
          fV0DCAxyToPV = lThisDCAxy;
          fV0DCAzToPV = lThisDCAz;
          
          int nCand2 = fitterCasc.process(tV0, bTrack);
          if (nCand2 != 0) {
            fitterCasc.propagateTracksToVertex();
            hCascFinding->Fill(1.5);
            const auto& cascvtx = fitterCasc.getPCACandidate();
            for (int i = 0; i < 3; i++) {
              posXi[i] = cascvtx[i];
            }
          }else{
            continue;
          }
          fXiDauDCA = TMath::Sqrt(fitterCasc.getChi2AtPCACandidate(0))*1e+4 ;
          if(fXiDauDCA>lMaxDCACascadeDaughters && lMasterSwitch > 0) continue;
          float cascradius = std::hypot(posXi[0],posXi[1]);
          fXiDecayRadius = cascradius;
          if( fXiDecayRadius < lMinXiRadius && lMasterSwitch > 0 ) continue;
          
          fV0DecayLength = std::hypot(weakdecaypos[0]-posXi[0],
                                      weakdecaypos[1]-posXi[1],
                                      weakdecaypos[2]-posXi[2]);
          
          fXiDecayLength = std::hypot(posXi[0]-vertices[0].getX(),
                                      posXi[1]-vertices[0].getY(),
                                      posXi[2]-vertices[0].getZ());
          
          //std::cout<<"--------------------------------------------------------------------------------------------------------------------------------------------"<< std::endl;
          //std::cout<<"---> V0 track: "<< std::endl;
          //tV0.print();
          
          //covariance matrix of the cascade
          std::array<float, 21> casccov0 = {0};
          std::array<float, 21> casccov1 = {0};
          std::array<float, 21> casccovcasc = {0};
          std::array<float, 3> pvecv0 = {0.};
          fitterCasc.getTrack(0).getPxPyPzGlo(pvecv0);
          fitterCasc.getTrack(1).getPxPyPzGlo(pvecbach);
          fitterCasc.getTrack(0).getCovXYZPxPyPzGlo(casccov0);
          fitterCasc.getTrack(1).getCovXYZPxPyPzGlo(casccov1);
          for (int i = 0; i < 6; i++) {
            int j = momInd[i];
            casccovcasc[j] = casccov0[j] + casccov1[j];
          }
          auto covVtxCasc = fitterCasc.calcPCACovMatrix();
          casccovcasc[0] = covVtxCasc(0, 0);
          casccovcasc[1] = covVtxCasc(1, 0);
          casccovcasc[2] = covVtxCasc(1, 1);
          casccovcasc[3] = covVtxCasc(2, 0);
          casccovcasc[4] = covVtxCasc(2, 1);
          casccovcasc[5] = covVtxCasc(2, 2);
          
          auto massXi = RecoDecay::M(array{array{pvecpos[0] + pvecneg[0], pvecpos[1] + pvecneg[1], pvecpos[2] + pvecneg[2]}, array{pvecbach[0], pvecbach[1], pvecbach[2]}}, array{1.115683, RecoDecay::getMassPDG(kPiPlus)});
          if( TMath::Abs(massXi-1.322)> lMassWindowXi && lMasterSwitch > 0) continue; //only look at peak, please
          fXiMass = massXi;
          
          const std::array<float, 3> momentumcascade = {pvecpos[0] + pvecneg[0] + pvecbach[0],
            pvecpos[1] + pvecneg[1] + pvecbach[1],
            pvecpos[2] + pvecneg[2] + pvecbach[2]};
          
          fXiTotalMomentum = std::hypot(pvecpos[0] + pvecneg[0] + pvecbach[0],
                                        pvecpos[1] + pvecneg[1] + pvecbach[1],
                                        pvecpos[2] + pvecneg[2] + pvecbach[2]);
          
          fNegativePt = std::hypot(pvecneg[0],pvecneg[1]);
          fPositivePt = std::hypot(pvecpos[0],pvecpos[1]);
          fBachelorPt = std::hypot(pvecbach[0],pvecbach[1]);
          
          lVector->SetPxPyPzE(pvecneg[0],pvecneg[1],pvecneg[2],1e+3); fNegativeEta = lVector->Eta();
          lVector->SetPxPyPzE(pvecpos[0],pvecpos[1],pvecpos[2],1e+3); fPositiveEta = lVector->Eta();
          lVector->SetPxPyPzE(pvecbach[0],pvecbach[1],pvecbach[2],1e+3); fBachelorEta = lVector->Eta();
          
          double pTcasc = TMath::Sqrt(TMath::Power(momentumcascade[0],2)+TMath::Power(momentumcascade[1],2));
          
          float thisdcacascdau = TMath::Sqrt(fitterCasc.getChi2AtPCACandidate(0));
          
          auto thiscasccospa = RecoDecay::CPA(array{vertices[0].getX(),vertices[0].getY(),vertices[0].getZ()},
                                              array{posXi[0], posXi[1], posXi[2]},
                                              array{momentumcascade[0], momentumcascade[0], momentumcascade[0]});
          
          auto tcascade = o2::track::TrackParCov(posXi, momentumcascade, casccovcasc, -1);
          
          o2::math_utils::Point3D<float> posXiPoint{posXi[0],posXi[1],posXi[2]};
          std::array<float, 6> covcascvertex;
          covcascvertex[0] = covVtxCasc(0, 0);
          covcascvertex[1] = covVtxCasc(1, 0);
          covcascvertex[2] = covVtxCasc(1, 1);
          covcascvertex[3] = covVtxCasc(2, 0);
          covcascvertex[4] = covVtxCasc(2, 1);
          covcascvertex[5] = covVtxCasc(2, 2);
          o2::dataformats::VertexBase vtxcascade(posXiPoint, covcascvertex);
          
          o2::dataformats::DCA dcacascade;
          //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          auto tcascadeenhanced = o2::track::TrackParCov(posXi, momentumcascade, casccovcasc, -1);
          //go for strangeness tracking, please
          //tcascadeenhanced.print();
          
          int iHits = 0;
          fTrueXi = kFALSE;
          
          if( truepostra[i] == 1 && truenegtra[j] == 1 && truebachtra[k] == 1 && lSameMotherBachelorV0 ){
            fTrueXi = kTRUE;
            std::cout << " *- Strangeness tracking starting" << std::endl;
            for (auto& hit : *hits) {
              const int layer{hit.GetDetectorID()};
              float xyz[3]{hit.GetX(), hit.GetY(), hit.GetZ()};
              float r{std::hypot(xyz[0], xyz[1])};
              float phi{std::atan2(-xyz[1], -xyz[0]) + o2::its::constants::math::Pi};
              
              //check if before cascade decay point
              if( r > cascradius ) continue;
              
              float lSmearing = lIBresolution;
              if(layer>3) lSmearing = lOBresolution;
              
              if (kUseSmearing) {
                phi = gRandom->Gaus(phi, std::asin(lSmearing / r));
                xyz[0] = r * std::cos(phi);
                xyz[1] = r * std::sin(phi);
                xyz[2] = gRandom->Gaus(xyz[2], lSmearing);
              }
              
              //check if it actually is the 3334
              if(abs(hit.GetTrackID()) >= mcArr->size()) continue;
              auto part = mcArr->at(hit.GetTrackID());
              if(part.GetPdgCode() != 3312) continue;
              if(hit.GetTrackID() != lCascadeID ) continue;
              
              double alpha = tcascadeenhanced.getAlpha();
              
              double xyz1[3]{ TMath::Cos(alpha)*xyz[0]+TMath::Sin(alpha)*xyz[1],
                -TMath::Sin(alpha)*xyz[0]+TMath::Cos(alpha)*xyz[1],
                xyz[2]};
              
              if(!(tcascadeenhanced.propagateTo(xyz1[0],tracker.getBz()))) continue;
              iHits++;
              const o2::track::TrackParametrization<float>::dim2_t hitpoint = {
                static_cast<float>(xyz1[1]),
                static_cast<float>(xyz1[2])
              };
              const o2::track::TrackParametrization<float>::dim3_t hitpointcov = {lSmearing * lSmearing, 0.f, lSmearing * lSmearing};
              std::cout << " *- Inspect after update:" << std::endl;
              tcascadeenhanced.update(hitpoint,hitpointcov);
              tcascadeenhanced.print();
            }
            std::cout << " *- Strangeness tracking found " << iHits<< " relevant hits" << std::endl;
          }
          if( ! fTrueXi && lAssociateXiDaughters ) continue;
          
          o2::dataformats::DCA dcacascade2;
          if (!tcascadeenhanced.propagateToDCA(vtx, tracker.getBz(), &dcacascade2)) {
            std::cout << "Cascade propagation failed." << std::endl;
          } else {
            fPtXi = tcascadeenhanced.getPt();
            fXiHitsAdded=iHits;
            fXiDCAxyToPV = dcacascade2.getY()*1e+4;
            fXiDCAzToPV = dcacascade2.getZ()*1e+4;
            if(iHits>0){
              std::cout << " *- DCAxy: " <<fabs(dcacascade2.getY()*1e+4)<<" microns"<< std::endl;
              std::cout << " *- DCAz: " <<fabs(dcacascade2.getZ()*1e+4)<<" microns"<< std::endl;
            }
          }
          //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          
          //store eta of Xi based on improved track
          fXiEta = tcascadeenhanced.getEta();
          
          if( fPtXi < lMinXiTransverseMomentum && lMasterSwitch > 0) continue;
          //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          
          //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          //Go for Xi_c
          Double_t lXiCMass = 1e+3;
          Double_t lXiCPt = 1e+3;
          Double_t lXiCPiDCAxyTrue = 1e+3;
          Double_t lXiCPiDCAzTrue = 1e+3;
          Double_t ptreco=1e+3, ptrecoSuper=1e+3, massXiC=1e+3, massXiCSuper=1e+3;
          
          Double_t lXicDecayRadius =1e+3;
          Double_t lXicDecayDistanceFromPV =1e+3;
          Double_t lXicDaughterDCA =1e+3;
          
          fXicDecayRadius = 1e+3;
          fXicDecayDistanceFromPV = 1e+3;
          fXicDaughterDCA = 1e+3;
          fXicVertexX = 1e+3;
          fXicVertexY = 1e+3;
          fXicVertexZ = 1e+3;
          
          fFirstCandidateXiC = kTRUE;
          fFirstCandidateXiCC = kTRUE;
          
          fXiccMass=0.0;
          fPic1DCAxyToPV=1e+10;
          fPic2DCAxyToPV=1e+10;
          fPic1DCAzToPV=1e+10;
          fPic2DCAzToPV=1e+10;

          fPiccDCAxyToPV =1e+10;
          fPiccDCAzToPV =1e+10;
          Double_t lXicVertexX = 1e+3, lXicVertexY = 1e+3, lXicVertexZ = 1e+3;
          
          //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          //Calculate expected TOF for daughter prongs, please!
          Double_t lCascadeLength = -1;
          Double_t lBachelorLength = -1;
          Double_t lV0Length = -1;
          Double_t lNegativeLength = -1;
          Double_t lPositiveLength = -1;
          Double_t lCascSpeed, lBachSpeed, lV0Speed, lNegSpeed, lPosSpeed;
          Double_t lBachelorLengthToPV = -1;
          Double_t lNegativeLengthToPV = -1;
          Double_t lPositiveLengthToPV = -1;
          
          float lX0=-100., lX1=-100.;
          
          //Step 1: cascade flight path
          if (!tcascadeenhanced.propagateToDCA(vtx, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Cascade propagation failed." << std::endl;
          } else {
            lX0 = tcascadeenhanced.getX();
          }
          if (!tcascadeenhanced.propagateToDCA(vtxcascade, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Cascade propagation failed." << std::endl;
          } else {
            lX1 = tcascadeenhanced.getX();
          }
          if(lX0>-99.&&lX1>-99.) lCascadeLength = TrackLength(tcascadeenhanced, lX0, lX1, lMagneticField);
          tcascadeenhanced.propagateTo(lX0, lMagneticField);
          
          //Step 1: V0 flight path
          lX0=-100.; lX1=-100.;
          if (!tV0.propagateToDCA(vtxcascade, tracker.getBz(), &dcacascade2)) {
            //std::cout << "V0 propagation failed." << std::endl;
          } else {
            lX0 = tV0.getX();
          }
          if (!tV0.propagateToDCA(vtxV0, tracker.getBz(), &dcacascade2)) {
            //std::cout << "V0 propagation failed." << std::endl;
          } else {
            lX1 = tV0.getX();
          }
          if(lX0>-99.&&lX1>-99.) lV0Length = TrackLength(tV0, lX0, lX1, lMagneticField);
          tV0.propagateTo(lX0, lMagneticField);
          
          //Step 2: Bachelor flight path
          lX0=-100.; lX1=-100.;
          if (!bTrack.propagateToDCA(vtxcascade, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Bachelor propagation failed." << std::endl;
          } else {
            lX0 = bTrack.getX();
          }
          if (!bTrack.getXatLabR(100.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
            //std::cout << "Bachelor propagation failed." << std::endl;
          }
          if(lX0>-99.&&lX1>-99.) lBachelorLength = TrackLength(bTrack, lX0, lX1, lMagneticField);
          //Propagate back to lX0
          bTrack.propagateTo(lX0, lMagneticField);
          
          //Step 3: Positive flight path
          lX0=-100.; lX1=-100.;
          if (!Track1.propagateToDCA(vtxV0, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Positive propagation failed." << std::endl;
          } else {
            lX0 = Track1.getX();
          }
          if (!Track1.getXatLabR(100.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
            //std::cout << "Positive propagation failed." << std::endl;
          }
          if(lX0>-99.&&lX1>-99.) lPositiveLength = TrackLength(Track1, lX0, lX1, lMagneticField);
          Track1.propagateTo(lX0, lMagneticField);
          
          //Step 3: Negative flight path
          lX0=-100.; lX1=-100.;
          if (!Track2.propagateToDCA(vtxV0, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Negative propagation failed." << std::endl;
          } else {
            lX0 = Track2.getX();
          }
          if (!Track2.getXatLabR(100.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
            //std::cout << "Negative propagation failed." << std::endl;
          }
          if(lX0>-99.&&lX1>-99.) lNegativeLength = TrackLength(Track2, lX0, lX1, lMagneticField);
          Track2.propagateTo(lX0, lMagneticField);
          
          //================================================================================
          // Calculate lengths to the inner TOF (at 20 centimeters)
          //================================================================================
          
          //Basic premise: we tracked with at least 6 prongs
          //this restricts the phase space such that we can safely assume that all daughter
          //tracks have been created by the time we reach the inner TOF
          //This may not be the only possibility in the long run but is okay for now!
          
          Double_t lBachelorShortLength = -1;
          Double_t lNegativeShortLength = -1;
          Double_t lPositiveShortLength = -1;
          
          lX0=-100., lX1=-100.;
          
          //Step 2: Bachelor flight path
          lX0=-100.; lX1=-100.;
          if (!bTrack.propagateToDCA(vtxcascade, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Bachelor propagation failed." << std::endl;
          } else {
            lX0 = bTrack.getX();
          }
          if (!bTrack.getXatLabR(20.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
            //std::cout << "Bachelor propagation failed." << std::endl;
          }
          if(lX0>-99.&&lX1>-99.) lBachelorShortLength = TrackLength(bTrack, lX0, lX1, lMagneticField);
          //Propagate back to lX0
          bTrack.propagateTo(lX0, lMagneticField);
          
          //Step 3: Positive flight path
          lX0=-100.; lX1=-100.;
          if (!Track1.propagateToDCA(vtxV0, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Positive propagation failed." << std::endl;
          } else {
            lX0 = Track1.getX();
          }
          if (!Track1.getXatLabR(20.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
            //std::cout << "Positive propagation failed." << std::endl;
          }
          if(lX0>-99.&&lX1>-99.) lPositiveShortLength = TrackLength(Track1, lX0, lX1, lMagneticField);
          Track1.propagateTo(lX0, lMagneticField);
          
          //Step 4: Negative flight path
          lX0=-100.; lX1=-100.;
          if (!Track2.propagateToDCA(vtxV0, tracker.getBz(), &dcacascade2)) {
            //std::cout << "Negative propagation failed." << std::endl;
          } else {
            lX0 = Track2.getX();
          }
          if (!Track2.getXatLabR(20.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
            //std::cout << "Negative propagation failed." << std::endl;
          }
          if(lX0>-99.&&lX1>-99.) lNegativeShortLength = TrackLength(Track2, lX0, lX1, lMagneticField);
          Track2.propagateTo(lX0, lMagneticField);
          
          //================================================================================
          //Retrieve length to PV for pos/neg/bach (pre-calculated)
          //================================================================================
          
          lPositiveLengthToPV = lTrackLength[postra[i]];
          lNegativeLengthToPV = lTrackLength[negtra[j]];
          lBachelorLengthToPV = lTrackLength[bachtra[k]];
          
          //          cout<<"Track lengths: "<<endl;
          //          cout<<"Cascade....: "<<lCascadeLength<<endl;
          //          cout<<"V0.........: "<<lV0Length<<endl;
          //          cout<<"Bachelor...: "<<lBachelorLength<<endl;
          //          cout<<"Positive...: "<<lPositiveLength<<endl;
          //          cout<<"Negative...: "<<lNegativeLength<<endl;
          
          fBachelorLength = lBachelorLength;
          fNegativeLength = lNegativeLength;
          fPositiveLength = lPositiveLength;
          fCascadeLength = lCascadeLength;
          fV0Length = lV0Length;
          
          lCascSpeed = Velocity(fXiTotalMomentum, 1.32171);
          lV0Speed = Velocity(fV0TotalMomentum, 1.115683);
          lBachSpeed = Velocity(bTrack.getP(), 0.139570);
          lPosSpeed = Velocity(Track1.getP(), 0.938272);
          lNegSpeed = Velocity(Track2.getP(), 0.139570);
          
          //          cout<<"Cascade....: "<<lCascSpeed<<endl;
          //          cout<<"V0.........: "<<lV0Speed<<endl;
          //          cout<<"Bachelor...: "<<lBachSpeed<<endl;
          //          cout<<"Positive...: "<<lPosSpeed<<endl;
          //          cout<<"Negative...: "<<lNegSpeed<<endl;
          
          fBachelorLife = (lBachelorLength/lBachSpeed);
          fNegativeLife = (lNegativeLength/lNegSpeed);
          fPositiveLife = (lPositiveLength/lPosSpeed);
          fCascadeLife = (lCascadeLength/lCascSpeed);
          fV0Life = (lV0Length/lV0Speed);
          
          fBachelorExpectedSignal = (lCascadeLength/lCascSpeed) + (lBachelorLength/lBachSpeed);
          fNegativeExpectedSignal = (lCascadeLength/lCascSpeed) + (lV0Length/lV0Speed) + (lNegativeLength/lNegSpeed);
          fPositiveExpectedSignal = (lCascadeLength/lCascSpeed) + (lV0Length/lV0Speed) + (lPositiveLength/lPosSpeed);
          
          fBachelorInnerExpectedSignal = (lCascadeLength/lCascSpeed) + (lBachelorShortLength/lBachSpeed);
          fNegativeInnerExpectedSignal = (lCascadeLength/lCascSpeed) + (lV0Length/lV0Speed) + (lNegativeShortLength/lNegSpeed);
          fPositiveInnerExpectedSignal = (lCascadeLength/lCascSpeed) + (lV0Length/lV0Speed) + (lPositiveShortLength/lPosSpeed);
          
          fBachelorExpectedSignalFromPV = (lBachelorLengthToPV/lBachSpeed);
          fNegativeExpectedSignalFromPV = (lNegativeLengthToPV/lNegSpeed);
          fPositiveExpectedSignalFromPV = (lPositiveLengthToPV/lPosSpeed);
          
          if(lCascadeLength<0||lBachelorLength<0) fBachelorExpectedSignal=-1e+10;
          if(lCascadeLength<0||lV0Length<0||lPositiveLength<0) fPositiveExpectedSignal=-1e+10;
          if(lCascadeLength<0||lV0Length<0||lNegativeLength<0) fNegativeExpectedSignal=-1e+10;
          
          if(lCascadeLength<0||lBachelorLength<0) fBachelorInnerExpectedSignal=-1e+10;
          if(lCascadeLength<0||lV0Length<0||lPositiveLength<0) fPositiveInnerExpectedSignal=-1e+10;
          if(lCascadeLength<0||lV0Length<0||lNegativeLength<0) fNegativeInnerExpectedSignal=-1e+10;
          
          if(lBachelorLengthToPV<0) fBachelorExpectedSignalFromPV = -1e+10;
          if(lPositiveLengthToPV<0) fPositiveExpectedSignalFromPV = -1e+10;
          if(lNegativeLengthToPV<0) fNegativeExpectedSignalFromPV = -1e+10;
          
          fBachelorTOFSignal = lOuterTOFSignal[bachtra[k]];
          fPositiveTOFSignal = lOuterTOFSignal[postra[i]];
          fNegativeTOFSignal = lOuterTOFSignal[negtra[j]];
          fBachelorTOFSignalMC = lOuterTOFSignalMC[bachtra[k]];
          fPositiveTOFSignalMC = lOuterTOFSignalMC[postra[i]];
          fNegativeTOFSignalMC = lOuterTOFSignalMC[negtra[j]];
          fBachelorInnerTOF20Signal = lInnerTOFSignal20[bachtra[k]];
          fPositiveInnerTOF20Signal = lInnerTOFSignal20[postra[i]];
          fNegativeInnerTOF20Signal = lInnerTOFSignal20[negtra[j]];
          fBachelorInnerTOF50Signal = lInnerTOFSignal50[bachtra[k]];
          fPositiveInnerTOF50Signal = lInnerTOFSignal50[postra[i]];
          fNegativeInnerTOF50Signal = lInnerTOFSignal50[negtra[j]];
          fBachelorInnerTOFSignalMC = lInnerTOFSignalMC[bachtra[k]];
          fPositiveInnerTOFSignalMC = lInnerTOFSignalMC[postra[i]];
          fNegativeInnerTOFSignalMC = lInnerTOFSignalMC[negtra[j]];
          
          if( TMath::Abs(fBachelorInnerExpectedSignal - fBachelorInnerTOF20Signal ) > lMaxTimeOffsetInnerTOF && lMasterSwitch > 0 ) continue;
          if( TMath::Abs(fPositiveInnerExpectedSignal - fPositiveInnerTOF20Signal ) > lMaxTimeOffsetInnerTOF && lMasterSwitch > 0 ) continue;
          if( TMath::Abs(fNegativeInnerExpectedSignal - fNegativeInnerTOF20Signal ) > lMaxTimeOffsetInnerTOF && lMasterSwitch > 0 ) continue;
          
          //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
          //Save Xi
          fTreeCandidates->Fill(); //all information should be at hand now
          //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        }//end loop on bachelors
      }
    }//end loop on track pairs
    lWatchAnalysis.Stop();
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    std::cout<<"---> Valid attempts so far = "<<fTreeCandidates->GetEntries()<< std::endl;
  } //end event loop
  
  fQA->Write(); 
  fileTree->cd();
  hEventCounter->Write(); // keep track of event counter, please
  hEventCounterWithVertex->Write();
  hVertexZ->Write();
  hXiGeneratedPt->Write();
  hXiCGeneratedPt->Write();
  hXiCCGeneratedPt->Write();
  fTreeCandidates->Write();
  hNTracks->Write();
  hNLongTracks->Write();
  hNLongPrimaryTracks->Write();
  hNPrimaryTracks->Write();
  hNVertices->Write();
  
  hPtEtaGeneratedXi->Write();
  hPtEtaGeneratedXiC->Write();
  hPtEtaGeneratedXiCC->Write();
  hPtYGeneratedXi->Write();
  hPtYGeneratedXiC->Write();
  hPtYGeneratedXiCC->Write();
  
  hTime->Fill(0.5, lWatchSmearing.RealTime());
  hTime->Fill(1.5, lWatchTracking.RealTime());
  hTime->Fill(2.5, lWatchAnalysis.RealTime());
  hTime->Fill(3.5, lWatchSmearing.CpuTime());
  hTime->Fill(4.5, lWatchTracking.CpuTime());
  hTime->Fill(5.5, lWatchAnalysis.CpuTime());
  
  hTime->Write();
  
  for(Int_t layer=0; layer<12; layer++)
  hHitDensityMap[layer]->Write();
  for(Int_t layer=0; layer<12; layer++)
  hHitDensityProfile[layer]->Write();
  
  hPtMultGeneratedXi_AllEta -> Write();
  hPtMultGeneratedXiC_AllEta -> Write();
  hPtMultGeneratedXiCC_AllEta -> Write();
  
  hPtMultGeneratedXi_AllY -> Write();
  hPtMultGeneratedXiC_AllY -> Write();
  hPtMultGeneratedXiCC_AllY -> Write();
  
  hPtMultGeneratedXi_Eta5 -> Write();
  hPtMultGeneratedXiC_Eta5 -> Write();
  hPtMultGeneratedXiCC_Eta5 -> Write();
  
  hPtMultGeneratedXi_Y5 -> Write();
  hPtMultGeneratedXiC_Y5 -> Write();
  hPtMultGeneratedXiCC_Y5 -> Write();
  
  hPtMultGeneratedXi_Eta10 -> Write();
  hPtMultGeneratedXiC_Eta10 -> Write();
  hPtMultGeneratedXiCC_Eta10 -> Write();
  
  hPtMultGeneratedXi_Y10 -> Write();
  hPtMultGeneratedXiC_Y10 -> Write();
  hPtMultGeneratedXiCC_Y10 -> Write();
  
  hQADCAxy->Write();
  hQADCAz->Write();
  hQAMCNch->Write();
  hQANchCorrelation->Write();
  
  //std::cout<<"---> Xi              = "<<fTreeXi->GetEntries()<< std::endl;
  std::cout<<"---> Xi_c candidates = "<<fTreeCandidates->GetEntries()<< std::endl;
  
  TFile output("xi.histooutput.root", "recreate");
  
  hGenerated->Write();
  hNVertices->Write();
  hNTracks->Write();
  hVertexZ->Write();
  
  hTrackCount->Write();
  hCombinatorics->Write();
  
  hCombinatoricsV0->Write();
  
  cout<<"Performance characteristics: "<<endl;
  cout<<"Smearing time (real).....: "<<lWatchSmearing.RealTime()<<" seconds"<<endl;
  cout<<"Tracking time (real).....: "<<lWatchTracking.RealTime()<<" seconds"<<endl;
  cout<<"Anslysis time (real).....: "<<lWatchAnalysis.RealTime()<<" seconds"<<endl;
  
  cout<<"Smearing time (CPU).....: "<<lWatchSmearing.CpuTime()<<" seconds"<<endl;
  cout<<"Tracking time (CPU).....: "<<lWatchTracking.CpuTime()<<" seconds"<<endl;
  cout<<"Anslysis time (CPU).....: "<<lWatchAnalysis.CpuTime()<<" seconds"<<endl;
  
}

