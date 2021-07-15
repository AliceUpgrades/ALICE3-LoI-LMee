#include <string>

#include <TFile.h>
#include <TChain.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TGeoGlobalMagField.h>
#include <vector>
#include "CommonConstants/PhysicsConstants.h"

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

//Operational parameters
const Bool_t lAssociateXiDaughters = kFALSE;
//Ask for hohlweger boolean
const Bool_t lHohlwegerBooleanXic = kFALSE;
const Bool_t lHohlwegerBooleanXicc = kFALSE;

//Tracking pass configuration
const Float_t lFirstOpening = 1.0;
const Float_t lOpening = 3.0;

//This switch will enable or disable ALL CUTS BELOW
const Float_t lMasterSwitch = +1; //1: ON, -1: OFF

//Remove fundamentally insane combinatorics, please
const Float_t lMinV0Radius = lMasterSwitch*0.5; //in centimeters
const Float_t lMinXiRadius = lMasterSwitch*0.5; //in centimeters
const Float_t lMinDCAxyBachelor = lMasterSwitch*0.0040; //in centimeters
const Float_t lMinDCAzBachelor = lMasterSwitch*0.0040; //in centimeters
const Float_t lMinDCAxyPositive = lMasterSwitch*0.0050; //in centimeters
const Float_t lMinDCAzPositive = lMasterSwitch*0.0040; //in centimeters
const Float_t lMinDCAxyNegative = lMasterSwitch*0.0100; //in centimeters
const Float_t lMinDCAzNegative = lMasterSwitch*0.0050; //in centimeters
const Float_t lMinDCAxyHFPions = lMasterSwitch*0.0010; //in centimeters
const Float_t lMinDCAzHFPions  = lMasterSwitch*0.0010; //in centimeters

const Float_t lMassWindowLambda = 0.012;
const Float_t lMassWindowXi = 0.150;
const Float_t lMaxDCAV0Daughters = 1000; //in microns
const Float_t lMaxDCACascadeDaughters = 1200; //in microns

const Double_t lIBresolution = 0.00025f;
const Double_t lOBresolution = 0.00100f;

const Double_t liTOFresolution = 50; //ps
const Double_t loTOFresolution = 20; //ps
const Int_t lNStepsLIntegrator = 100;

const Double_t lMagneticField = 5.0;

using o2::its::MemoryParameters;
using o2::its::TrackingParameters;
using o2::itsmft::Hit;
using std::string;

using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

constexpr bool kUseSmearing{true};

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

Double_t TrackLength( o2::track::TrackParCov track, Double_t lX0, Double_t lX1 ){
  std::array<float, 3> lPointN;
  std::array<float, 3> lPointNplus;
  Double_t lLength = 0.0;
  track.propagateTo(lX0, lMagneticField);
  track.getXYZGlo(lPointN);
  Double_t lFirstRadius = std::hypot( lPointN[0], lPointN[1] );
  for(Int_t iStep=1; iStep<lNStepsLIntegrator; iStep++){
    track.getXYZGlo(lPointN);
    Float_t lPosition = lX0 + (lX1-lX0)*((Float_t)(iStep))/((Float_t)(lNStepsLIntegrator-1));
    if(!track.propagateTo(lPosition, lMagneticField)) return -100; //error in propagation
    track.getXYZGlo(lPointNplus);
    lLength += std::hypot( lPointNplus[0]-lPointN[0], lPointNplus[1]-lPointN[1], lPointNplus[2]-lPointN[2] );
  }
  Double_t lLastRadius = std::hypot( lPointNplus[0], lPointNplus[1] );
  track.propagateTo(lX0, lMagneticField);
  return lLength;
}
void run_strangenesstracking_01()
{
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  std::cout << "\e[1;31m              ALICE 3 tracking tool      \e[0;00m" << std::endl;
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  
  //Open Geometry, please
  o2::base::GeometryManager::loadGeometry("./");
  
  //Open GRP
  const auto grp = o2::parameters::GRPObject::loadFrom("o2sim_grp.root");
  if (!grp) {
    LOG(FATAL) << "Cannot run w/o GRP object";
  }

  o2::base::Propagator::initFieldFromGRP(grp);
  auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  if (!field) {
    LOG(FATAL) << "Failed to load ma";
  }
  double origD[3] = {0., 0., 0.};
  
  Int_t Particle = 3122;
  Int_t lNLayers = 7;
  const string hitsFileName = "o2sim_HitsTRK.root";
  TChain mcTree("o2sim");
  mcTree.AddFile("o2sim_Kine.root");
  mcTree.SetBranchStatus("*", 0); //disable all branches
  mcTree.SetBranchStatus("MCTrack*", 1);
  
  std::vector<o2::MCTrack>* mcArr = nullptr;
  mcTree.SetBranchAddress("MCTrack", &mcArr);
  
  o2::its::Vertexer vertexer(new o2::its::VertexerTraits());
  
  TChain itsHits("o2sim");
  
  itsHits.AddFile(hitsFileName.data());
  
  o2::its::Tracker tracker(new o2::its::TrackerTraitsCPU());
  tracker.setBz(5.f);
  tracker.setCorrType(o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrTGeo);
  
  std::uint32_t roFrame;
  std::vector<Hit>* hits = nullptr;
  itsHits.SetBranchAddress("TRKHit", &hits);
  
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  // TRACKING PARAMETERS
  std::vector<TrackingParameters> trackParams(2);
  //Tracking parameters for 12 layer setup
  trackParams[0].NLayers = 12;
  trackParams[0].MinTrackLength = 12; //this is the one with fixed params
  std::vector<float> LayerRadii = {0.5f, 1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f, 100.0f};
  std::vector<float> LayerZ(12);
  for (int i{0}; i < 12; ++i)
    LayerZ[i] = getDetLengthFromEta(1.44, LayerRadii[i]) + 1.;

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

  Float_t lFactor = 200;
  std::vector<MemoryParameters> memParams(2);
  std::vector<float> CellsMemoryCoefficients = {2.3208e-08f * 30, 2.104e-08f * 30, 1.6432e-08f * 30, 1.2412e-08f * 30, 1.3543e-08f * 30, 1.5e-08f * 30, 1.6e-08f * 30, 1.7e-08f * 30};
  std::vector<float> TrackletsMemoryCoefficients = {0.0016353f * lFactor, 0.0013627f * lFactor, 0.000984f * lFactor, 0.00078135f * lFactor, 0.00057934f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor, 0.00052217f * lFactor};
  memParams[0].CellsMemoryCoefficients = CellsMemoryCoefficients;
  memParams[0].TrackletsMemoryCoefficients = TrackletsMemoryCoefficients;
  memParams[0].MemoryOffset = 120000;

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
  memParams[1].MemoryOffset = 120000;

  tracker.setParameters(memParams, trackParams);
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  
  Double_t ptbinlimits[] ={ 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.2,2.4,2.6,2.8,3.0,
    3.3,3.6,3.9,4.2,4.6,5,5.4,5.9, 6.5,7,7.5,8,8.5,9.2,10,11,12,13.5,15,17,20};
  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
  
  TH1D *hNVertices = new TH1D("hNVertices", "", 10,0,10);
  TH1D *hNTracks = new TH1D("hNTracks", "", 10000,0,10000);
  
  TH1D *hXiR2D = new TH1D("hXiR2D", "", 2000,0,100);
  TH1D *hLambdaR2D = new TH1D("hLambdaR2D", "", 2000,0,100);
  Double_t lRadii[] = {0.0f,0.5f,1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f, 100.0f};
  TH1D *hXiR2Dbinned = new TH1D("hXiR2Dbinned", "", 12, lRadii);
  TH1D *hLambdaR2Dbinned = new TH1D("hLambdaR2Dbinned", "", 12, lRadii);
  
  Int_t fParticle, fMinNLayerReco;
  
  TFile *fileTree = new TFile("treeoutput.root", "RECREATE", "", 109);
  
  TH1D *hEventCounter = new TH1D("hEventCounter", "", 1,0,1);
  TH1D *hEventCounterWithVertex = new TH1D("hEventCounterWithVertex", "", 1,0,1);

  Float_t fTrackDCAxy, fTrackDCAz;
  Float_t fTrackDCAxyMC, fTrackDCAzMC;
  Float_t fTrackLength, fTrackP, fTOFSignal, fTOFSignalMC, fInnerTOF, fInnerTOF20, fInnerTOFMC, fTrackShortLength;
  Bool_t fIsPrimary;
  Int_t fTrackPDG;
  TTree *fTreeDebug = new TTree ( "fTreeDebug", "Event Characterization Tree" ) ;
  fTreeDebug->Branch ("fIsPrimary",  &fIsPrimary,  "fIsPrimary/O"  );
//  fTreeDebug->Branch ("fTrackDCAxy",  &fTrackDCAxy,  "fTrackDCAxy/F"  );
//  fTreeDebug->Branch ("fTrackDCAz",  &fTrackDCAz,  "fTrackDCAz/F"  );
//
//  //replace PV position for true position (zero)
//  fTreeDebug->Branch ("fTrackDCAxyMC",  &fTrackDCAxyMC,  "fTrackDCAxyMC/F"  );
//  fTreeDebug->Branch ("fTrackDCAzMC",  &fTrackDCAzMC,  "fTrackDCAzMC/F"  );

  fTreeDebug->Branch ("fTrackPDG",  &fTrackPDG,  "fTrackPDG/I"  );

  fTreeDebug->Branch ("fTrackP",  &fTrackP,  "fTrackP/F"  );
  fTreeDebug->Branch ("fTrackLength",  &fTrackLength,  "fTrackLength/F"  );
  fTreeDebug->Branch ("fTrackShortLength",  &fTrackShortLength,  "fTrackShortLength/F"  );
  
  fTreeDebug->Branch ("fTOFSignal",  &fTOFSignal,  "fTOFSignal/F"  );
  fTreeDebug->Branch ("fTOFSignalMC",  &fTOFSignalMC,  "fTOFSignalMC/F"  );
  
  fTreeDebug->Branch ("fInnerTOF",  &fInnerTOF,  "fInnerTOF/F"  );
  fTreeDebug->Branch ("fInnerTOF20",  &fInnerTOF20,  "fInnerTOF20/F"  );
  fTreeDebug->Branch ("fInnerTOFMC",  &fInnerTOFMC,  "fInnerTOFMC/F"  );
  
  Float_t fTrackChisquare;
  Int_t fTrackClusters;
  Float_t fTrackEta;
  fTreeDebug->Branch ("fTrackChisquare",  &fTrackChisquare,  "fTrackChisquare/F"  );
  fTreeDebug->Branch ("fTrackClusters",  &fTrackClusters,  "fTrackClusters/I"  );
  fTreeDebug->Branch ("fTrackEta",  &fTrackEta,  "fTrackEta/F"  );
  
  //This is for convenience -> Be careful, please!
  fParticle = Particle;
  fMinNLayerReco = lNLayers;
  
  for (int iEvent{0}; iEvent < itsHits.GetEntriesFast(); ++iEvent) {
    std::cout << "*************** Event " << iEvent << " ***************" << std::endl;
    itsHits.GetEntry(iEvent);
    mcTree.GetEvent(iEvent);
    o2::its::ROframe event{iEvent, 12};
    hEventCounter->Fill(0.5);
    
    //Locate Xi decay point
    Int_t lDirectTrackingNLayers = 0;
    Int_t lNLayersLambda = 0;
    cout<<"Checking parenthood - mcArr size "<<mcArr->size()<<endl;
    
        
    int id{0};
    std::map<int, particle> mapPDG;
    Long_t NHitAll=0, NHitInner=0, NHitBachelor=0, NHitV0=0 ;
    std::cout << "*- Filling hits" << std::endl;
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
      
      NHitAll++;
      if(NHitAll%1000==0) std::cout << "*- Hits filled: "<< NHitAll << std::endl;
      event.addTrackingFrameInfoToLayer(layer, xyz[0], xyz[1], xyz[2], r, phi, std::array<float, 2>{0.f, xyz[2]},
                                        std::array<float, 3>{lSmearing * lSmearing, 0.f, lSmearing * lSmearing});
      event.addClusterToLayer(layer, xyz[0], xyz[1], xyz[2], event.getClustersOnLayer(layer).size());
      event.addClusterLabelToLayer(layer, o2::MCCompLabel(hit.GetTrackID(), iEvent, iEvent, false));
      event.addClusterExternalIndexToLayer(layer, id++);
      //std::cout << "*- Event " << iEvent << " hit.GetTrackID() = " <<hit.GetTrackID() << " endhit"<<std::endl;
    }
    roFrame = iEvent;
    //std::cout << "*- Event " << iEvent << " MC loop" << std::endl;
    //Loop over generated PoI
    
    vertexer.clustersToVertices(event);
    
    std::vector<Vertex> vertices = vertexer.exportVertices();
    std::cout<<"Number of vertices found: "<<vertices.size()<<endl;
    hNVertices->Fill(vertices.size());
    //Skip events with no vertex
    if(vertices.size()==0){
      std::cout <<"No primary vertex found, skipping event"<<std::endl;
      continue;
    }
    hEventCounterWithVertex->Fill(0.5);
    
    vertices[0].print();
    
    o2::math_utils::Point3D<float> pos{vertices[0].getX(),vertices[0].getY(),vertices[0].getZ()};
    o2::math_utils::Point3D<float> posZero{0.0,0.0,0.0};
    std::array<float, 6> cov;
    for(Int_t jj=0; jj<6; jj++) cov[jj]=vertices[0].getCov()[jj];
    o2::dataformats::VertexBase vtx(pos, cov);
    o2::dataformats::VertexBase vtxZero(posZero, cov);
    o2::dataformats::DCA dca;
    
    std::cout << "*- Event " << iEvent << " tracking" << std::endl;
    tracker.clustersToTracks(event);
    auto& tracks = tracker.getTracks();
    auto& tracksLabels = tracker.getTrackLabels();
    std::cout << "*- Event " << iEvent << " done tracking!" << std::endl;
    
    hNTracks -> Fill(tracks.size());
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //Loop over tracks, attempt to determine TOF information
    
    //Monte Carlo: perfect signal
    TArrayI lInnerTOFSignalMC(tracks.size());
    TArrayI lOuterTOFSignalMC(tracks.size());
    
    //Reconstructed: smeared signal
    TArrayI lInnerTOFSignal(tracks.size());
    TArrayI lInnerTOFSignal20(tracks.size());
    TArrayI lOuterTOFSignal(tracks.size());
    TArrayF lTrackLength(tracks.size());
    TArrayF lTrackShortLength(tracks.size());
    
    cout<<"Determining TOF values for inner and outer TOF..."<<endl;
    for (unsigned int i{0}; i < tracks.size(); ++i) {
      auto& lab = tracksLabels[i];
      auto& track = tracks[i];
      int lMCID = lab.getTrackID();
      
      lInnerTOFSignal[i] = -1;
      lOuterTOFSignal[i] = -1;
      lInnerTOFSignalMC[i] = -1;
      lOuterTOFSignalMC[i] = -1;
      lTrackLength[i] = -1;
      //Determine length wrt to PV for posterior use
      Float_t lX0=-100, lX1=-100;
      Float_t lThisTrackLength = -100;
      if (!track.propagateToDCA(vtx, tracker.getBz())) {
        std::cout << "Propagation failed." << std::endl;
      } else {
        lX0 = track.getX();
      }
      if (!track.getXatLabR(100.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
        std::cout << "Outward propagation failed." << std::endl;
      }
      if(lX0>-99.&&lX1>-99.) lThisTrackLength = TrackLength(track, lX0, lX1);
      lTrackLength[i] = lThisTrackLength;
      track.propagateTo(lX0, lMagneticField);
      
      lThisTrackLength = -100;
      lX0=-100, lX1=-100;
      if (!track.propagateToDCA(vtx, tracker.getBz())) {
        std::cout << "Propagation failed." << std::endl;
      } else {
        lX0 = track.getX();
      }
      if (!track.getXatLabR(20.0,lX1,tracker.getBz(),o2::track::DirOutward)) {
        std::cout << "Outward propagation failed." << std::endl;
      }
      if(lX0>-99.&&lX1>-99.) lThisTrackLength = TrackLength(track, lX0, lX1);
      lTrackShortLength[i] = lThisTrackLength;
      track.propagateTo(lX0, lMagneticField);
      
      Float_t lThisDCAxy=0.0, lThisDCAz=0.0;
      Float_t lThisDCAxyMC=0.0, lThisDCAzMC=0.0;
      if (!track.propagateToDCA(vtxZero, tracker.getBz(), &dca)) {
        std::cout << "Propagation failed." << std::endl;
        lThisDCAxyMC = 1e+10;
        lThisDCAzMC = 1e+10;
      } else {
        lThisDCAxyMC = dca.getY();
        lThisDCAzMC = dca.getZ();
      }
      if (!track.propagateToDCA(vtx, tracker.getBz(), &dca)) {
        std::cout << "Propagation failed." << std::endl;
        lThisDCAxy = 1e+10;
        lThisDCAz = 1e+10;
      } else {
        lThisDCAxy = dca.getY();
        lThisDCAz = dca.getZ();
      }
      
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
          lInnerTOFSignal[i] = gRandom->Gaus(lLowestTimeInner, liTOFresolution);
          lInnerTOFSignal20[i] = gRandom->Gaus(lLowestTimeInner, loTOFresolution);
          lInnerTOFSignalMC[i] = lLowestTimeInner;
        }
        if(lLowestTimeOuter<1e+6-1){
          lOuterTOFSignal[i] = gRandom->Gaus(lLowestTimeOuter, loTOFresolution);
          lOuterTOFSignalMC[i] = lLowestTimeOuter;
        }
      }
      auto part = mcArr->at(lMCID);
      
      fTrackDCAxy = lThisDCAxy;
      fTrackDCAz = lThisDCAz;
      fTrackDCAxyMC = lThisDCAxyMC;
      fTrackDCAzMC = lThisDCAzMC;
      fTrackClusters = track.getNumberOfClusters();
      fTrackChisquare = track.getChi2();
      fTrackEta = track.getEta();
      fTrackP = track.getP();
      fTrackLength = lTrackLength[i];
      fTrackShortLength = lTrackShortLength[i];
      fTrackPDG = part.GetPdgCode();
      fIsPrimary = part.isPrimary();
      fTOFSignal = lOuterTOFSignal[i];
      fTOFSignalMC = lOuterTOFSignalMC[i];
      fInnerTOF = lInnerTOFSignal[i];
      fInnerTOF20 = lInnerTOFSignal20[i];
      fInnerTOFMC = lInnerTOFSignalMC[i];
      fTreeDebug->Fill();
    }
    cout<<"TOF values determined!"<<endl;
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  }
  
  fileTree->cd();
  hEventCounter->Write(); // keep track of event counter, please
  hEventCounterWithVertex->Write();

  fTreeDebug->Write();
  hNTracks->Write();
  hNVertices->Write();
    
  //Small file. if this file was created, job was successful!
  TFile output("histooutput.root", "recreate");
  hNVertices->Write();
  hNTracks->Write();
}

