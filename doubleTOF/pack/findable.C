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
const Float_t lFirstOpening = 1.;
const Float_t lOpening = 2.;
const Float_t lThirdPass = 3.;
const Float_t lFourthPass = 3.;

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

void findable(TString lDataPath = "./", TString lLutPath =  "../")
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

  cout<<"Initializing track smearer..."<<flush;
  
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
  
  auto gen = new SmearO2KineGenerator( Form("%so2sim_Kine.root", lDataPath.Data()) );
  gen->Init();
  cout<<" done"<<endl;
  
  o2::steer::InteractionSampler irSampler;
  irSampler.setInteractionRate(10000);
  irSampler.init();

  o2::vertexing::PVertexer vertexer;
  vertexer.setValidateWithIR(kFALSE);
  vertexer.setBunchFilling(irSampler.getBunchFilling());
  vertexer.init();
  
  //vertexer QA
  TFile *fQA = new TFile("findable.root", "RECREATE");
  
  TH1D *hGenerated = new TH1D("hGenerated", "", 200,0,20);
  
  TH1D *hPartialTransportDebug = new TH1D("hPartialTransportDebug", "", 6,0,6);
  hPartialTransportDebug->GetXaxis()->SetBinLabel(1, "All particles");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(2, "Not tranported");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(3, "Successfully smeared");
  hPartialTransportDebug->GetXaxis()->SetBinLabel(4, "Fully tracked particles");
    
  TTree *fTree = new TTree("fTree", "");
  
  Float_t fPtMC, fEtaMC;
  Int_t fBachHits, fNegHits, fPosHits;
  
  fTree->Branch("fPtMC", &fPtMC, "fPtMC/F");
  fTree->Branch("fEtaMC", &fEtaMC, "fEtaMC/F");
  
  fTree->Branch("fBachHits", &fBachHits, "fBachHits/I");
  fTree->Branch("fPosHits", &fPosHits, "fPosHits/I");
  fTree->Branch("fNegHits", &fNegHits, "fNegHits/I");
  
  Bool_t fBachHitLayer[12];
  Bool_t fPosHitLayer[12];
  Bool_t fNegHitLayer[12];
  
  for(Int_t ih=0; ih<12; ih++) fTree->Branch( Form("fBachHits%i", ih), &(fBachHitLayer[ih]), Form("fBachHits%i/O",ih) );
  for(Int_t ih=0; ih<12; ih++) fTree->Branch( Form("fPosHits%i", ih), &(fPosHitLayer[ih]), Form("fPosHits%i/O",ih) );
  for(Int_t ih=0; ih<12; ih++) fTree->Branch( Form("fNegHits%i", ih), &(fNegHitLayer[ih]), Form("fNegHits%i/O",ih) );
  
  Bool_t fBachTracked, fPosTracked, fNegTracked;
  
  fTree->Branch("fBachTracked", &fBachTracked, "fBachTracked/O");
  fTree->Branch("fPosTracked", &fPosTracked, "fPosTracked/O");
  fTree->Branch("fNegTracked", &fNegTracked, "fNegTracked/O");
  
  Bool_t fCascadeReco;
  
  fTree->Branch("fCascadeReco", &fCascadeReco, "fCascadeReco/O");
  
  //Define o2 fitter, 2-prong
  o2::vertexing::DCAFitterN<2> fitterV0, fitterCasc;
  
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
  
  //Establish indices for the xi
  const Int_t lMax = 100; //absolute max one can store
  
  Int_t lPoIIndex[lMax];
  Float_t lPoIPt[lMax];
  Float_t lPoIEta[lMax];
  Int_t lPoIBachelor[lMax];
  Int_t lPoIPositive[lMax];
  Int_t lPoINegative[lMax];
  Bool_t lPoIBachelorLayer[lMax][12];
  Bool_t lPoIPositiveLayer[lMax][12];
  Bool_t lPoINegativeLayer[lMax][12];
  Int_t lNPoI = 0;
  
  Long_t lNchThisEvent = 0;
  
  for (int iEvent{0}; iEvent < itsHits.GetEntriesFast(); ++iEvent) {
    std::cout << "*************** Event " << iEvent << "/"<<itsHits.GetEntriesFast()<<" ***************" << std::endl;
    itsHits.GetEntry(iEvent);
    mcTree.GetEvent(iEvent);
    o2::its::ROframe event{iEvent, 12};
    
    cout<<"Finding generated Xi in this event..."<<endl;
    lNPoI=0;
    lNchThisEvent=0;
    //Identify Xis in this event: store indices
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
      
      
      if( part.GetPdgCode() == 3312 ){
        hGenerated->Fill(part.GetPt());
        for(Int_t ih=0; ih<12; ih++){
          lPoIBachelorLayer[lNPoI][ih] = kFALSE;
          lPoIPositiveLayer[lNPoI][ih] = kFALSE;
          lPoINegativeLayer[lNPoI][ih] = kFALSE;
        }
        
        lPoIIndex[lNPoI] = iii;
        lPoIPt[lNPoI] = part.GetPt();
        lPoIEta[lNPoI] = part.GetEta();
        
        Int_t lFirstDauID = part.getFirstDaughterTrackId();
        Int_t lLastDauID = part.getLastDaughterTrackId();
        
        //cout<<"To deref: "<<lFirstDauID<<" and "<<lLastDauID<<endl;
        if(lFirstDauID<0||lLastDauID<0) continue;
        
        auto firstdau = mcArr->at(lFirstDauID);
        auto lastdau = mcArr->at(lLastDauID);
        
        lPoIBachelor[lNPoI] = -1;
        lPoINegative[lNPoI] = -1;
        lPoIPositive[lNPoI] = -1;
        
        if(firstdau.GetPdgCode() == -211&&lastdau.GetPdgCode() == 3122){
          lPoIBachelor[lNPoI] = lFirstDauID;
          
          //cout<<"a) Lambda daus to deref: "<<lastdau.getFirstDaughterTrackId()<<" and "<<lastdau.getLastDaughterTrackId()<<endl;
          
          if(lastdau.getFirstDaughterTrackId()<0||lastdau.getLastDaughterTrackId()<0) continue;
          
          auto lambdadau1 = mcArr->at(lastdau.getFirstDaughterTrackId());
          auto lambdadau2 = mcArr->at(lastdau.getLastDaughterTrackId());
          
          if(lambdadau1.GetPdgCode() == -211 && lambdadau2.GetPdgCode() == 2212){
            lPoINegative[lNPoI] = lastdau.getFirstDaughterTrackId();
            lPoIPositive[lNPoI] = lastdau.getLastDaughterTrackId();
          }
          if(lambdadau1.GetPdgCode() == 2212 && lambdadau2.GetPdgCode() == -211){
            lPoINegative[lNPoI] = lastdau.getLastDaughterTrackId();
            lPoIPositive[lNPoI] = lastdau.getFirstDaughterTrackId();
          }
        }
        if(firstdau.GetPdgCode() == 3122&&lastdau.GetPdgCode() == -211){
          lPoIBachelor[lNPoI] = lLastDauID;
          
          //cout<<"b) Lambda daus to deref: "<<firstdau.getFirstDaughterTrackId()<<" and "<<firstdau.getLastDaughterTrackId()<<endl;
          
          if(firstdau.getFirstDaughterTrackId()<0||firstdau.getLastDaughterTrackId()<0) continue;
          
          auto lambdadau1 = mcArr->at(firstdau.getFirstDaughterTrackId());
          auto lambdadau2 = mcArr->at(firstdau.getLastDaughterTrackId());

          if(lambdadau1.GetPdgCode() == -211 && lambdadau2.GetPdgCode() == 2212){
            lPoINegative[lNPoI] = firstdau.getFirstDaughterTrackId();
            lPoIPositive[lNPoI] = firstdau.getLastDaughterTrackId();
          }
          if(lambdadau1.GetPdgCode() == 2212 && lambdadau2.GetPdgCode() == -211){
            lPoINegative[lNPoI] = firstdau.getLastDaughterTrackId();
            lPoIPositive[lNPoI] = firstdau.getFirstDaughterTrackId();
          }
        }
        cout<<"Xi indexed "<<lPoIIndex[lNPoI]<<", daughter indices "<<lPoIBachelor[lNPoI]<<", "<<lPoIPositive[lNPoI]<<", "<<lPoINegative[lNPoI]<<endl;
        lNPoI++;
      }
    }
    cout<<"Found "<<lNPoI<<" Xi in this event. Sweeping for daughter hit information"<<endl;
    
    int id{0};
    
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
      
      //Sweep PoI list
      for(Int_t ip=0; ip<lNPoI; ip++){
        if( lPoIBachelor[ip] == hit.GetTrackID() ) lPoIBachelorLayer[ip][layer] = kTRUE;
        if( lPoIPositive[ip] == hit.GetTrackID() ) lPoIPositiveLayer[ip][layer] = kTRUE;
        if( lPoINegative[ip] == hit.GetTrackID() ) lPoINegativeLayer[ip][layer] = kTRUE;
      }
      
      event.addTrackingFrameInfoToLayer(layer, xyz[0], xyz[1], xyz[2], r, phi, std::array<float, 2>{0.f, xyz[2]},
                                        std::array<float, 3>{lSmearing * lSmearing, 0.f, lSmearing * lSmearing});
      event.addClusterToLayer(layer, xyz[0], xyz[1], xyz[2], event.getClustersOnLayer(layer).size());
      event.addClusterLabelToLayer(layer, o2::MCCompLabel(hit.GetTrackID(), iEvent, iEvent, false));
      event.addClusterExternalIndexToLayer(layer, id++);
    }
    cout<<"Found daughter information. Now summing hit counters and saving all information to tree"<<endl;
        
    roFrame = iEvent;
    
    //Let's try to do the tracking and see if we find the desired tracks
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

    o2::InteractionRecord ir = irSampler.generateCollisionTime();
    Long_t lCorruptedTracks=0;
    
    for (int iparticle = 0; iparticle < particles.size(); ++iparticle){
      //Inner and outer TOF signal initialization
      
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
    
    if(n_vertices<1){
      continue;
    }
    
    int cntv = 0;
    for ( auto vert : vertices)  {
      LOG(INFO) << "PV: " << cntv++ << " " << cntEv << " " << vert.getX()<<" "<<vert.getY()<<" "<<vert.getZ()<< " " << vert.getNContributors() << " " << vert.getTimeStamp().getTimeStamp() << " " << vert.getTimeStamp().getTimeStampError() << endl;
    }
    cntEv++;
    
    if( TMath::Abs(vertices[0].getZ()) > 10 ) continue; //skip outside center of detector, please
    
    o2::math_utils::Point3D<float> pos{vertices[0].getX(),vertices[0].getY(),vertices[0].getZ()};
    std::array<float, 6> cov;
    for(Int_t jj=0; jj<6; jj++) cov[jj]=vertices[0].getCov()[jj];
    o2::dataformats::VertexBase vtx(pos, cov);
    
    //event.addPrimaryVertex(pos[0].getX(), pos[0].getY(), pos[0].getZ());
    event.addPrimaryVertex(vertices[0].getX(),vertices[0].getY(),vertices[0].getZ());
    
    std::cout << "*- Event " << iEvent << " tracking" << std::endl;
    tracker.clustersToTracks(event);
    auto& tracks = tracker.getTracks();
    auto& tracksLabels = tracker.getTrackLabels();
    std::cout << "*- Event " << iEvent << " done tracking!" << std::endl;
    
    std::cout << "*- Determining if xi daughter indices are captured by tracker or not..." << std::endl;
    
    for(Int_t ip=0; ip<lNPoI; ip++){
      fPtMC = lPoIPt[ip];
      fEtaMC = lPoIEta[ip];
      fBachHits=0;
      fNegHits=0;
      fPosHits=0;
      for(Int_t ih=0; ih<12; ih++){
        fBachHitLayer[ih] = kFALSE;
        fPosHitLayer[ih] = kFALSE;
        fNegHitLayer[ih] = kFALSE;
      }
      for(Int_t ih=0; ih<12; ih++){
        if( lPoIBachelorLayer[ip][ih] ) fBachHits++;
        if( lPoIPositiveLayer[ip][ih] ) fPosHits++;
        if( lPoINegativeLayer[ip][ih] ) fNegHits++;
        if( lPoIBachelorLayer[ip][ih] ) fBachHitLayer[ih] = kTRUE;
        if( lPoIPositiveLayer[ip][ih] ) fPosHitLayer[ih] = kTRUE;
        if( lPoINegativeLayer[ip][ih] ) fNegHitLayer[ih] = kTRUE;
      }
      //only fill if in correct channel
      if( lPoIPositive[ip] < 0 || lPoINegative[ip] < 0) continue;
      fTree->Fill();
      
      fBachTracked = kFALSE;
      fPosTracked = kFALSE;
      fNegTracked = kFALSE;
      
      Int_t lBachTrackIndex=-1, lPosTrackIndex=-1, lNegTrackIndex=-1;
      
      for (unsigned int i{0}; i < tracks.size(); ++i) {
        auto& lab = tracksLabels[i];
        auto& track = tracks[i];
        int lMCID = lab.getTrackID();
        
        if( lMCID == lPoIBachelor[ip] ) {
          fBachTracked = kTRUE;
          lBachTrackIndex = i;
        }
        if( lMCID == lPoIPositive[ip] ){
          fPosTracked = kTRUE;
          lPosTrackIndex = i;
        }
        if( lMCID == lPoINegative[ip] ){
          fNegTracked = kTRUE;
          lNegTrackIndex = i;
        }
      }
      fCascadeReco = kFALSE;
      
      //Attempt combination in case tracks correctly recoed
      if(fBachTracked&&fPosTracked&&fNegTracked){
        auto TrackBach = tracks[lBachTrackIndex];
        auto TrackPos = tracks[lPosTrackIndex];
        auto TrackNeg = tracks[lNegTrackIndex];
        
        int nCand = fitterV0.process(TrackPos, TrackNeg);
        if (nCand > 0) {
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
          
          const auto& wvtx = fitterV0.getPCACandidate();
          
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
          
          auto tV0 = o2::track::TrackParCov(vertex, momentum, covV0, 0);
          
          int nCand2 = fitterCasc.process(tV0, TrackBach);
          if (nCand2 > 0) {
            fCascadeReco = kTRUE;
          }
        }
      } //end attempt at cascade reconstruction
    } //end PoI loop
  } //end event loop
  hGenerated->Write();
  hPartialTransportDebug->Write();
  fQA->Write(); 
  
  cout<<"Done!"<<endl;
}

