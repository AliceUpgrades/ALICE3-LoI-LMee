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
// #include "ITStracking/TrackerTraitsCPU.h"
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
Double_t VelocityBeta(Double_t lMomentum, Double_t lMass){
  //Momentum p and mass m -> returns speed in centimeters per picosecond
  //Useful for TOF calculations
  Double_t lA = TMath::Power(lMomentum / lMass, 2);
  return TMath::Sqrt(lA/(1+lA));
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

double makeTOFpid(o2::track::TrackParCov track,int pdg, double lMagneticField = 0.0, double tofPos = 100., double tofRes = 50., bool makeSigma = false)
{
  o2::math_utils::Point3D<float> zeropos{0,0,0};
  std::array<float, 6> zerocov;
  for(Int_t jj=0; jj<6; jj++) zerocov[jj]=1e-6;
  o2::dataformats::VertexBase zerovtx(zeropos, zerocov);
  double lThisTrackLength = -1.;
  if(lMagneticField == 0.0) {
    cout << "Please set magnetic field properly!" << flush;
    return 0;
  }
  Float_t lX0=-100, lX1=-100;
  if (!track.propagateToDCA(zerovtx, lMagneticField)) {
    //std::cout << "Propagation failed." << std::endl;
  } else {
    lX0 = track.getX();
  }
  if (!track.getXatLabR(tofPos,lX1,lMagneticField,o2::track::DirOutward)) {
    lX1 = -100;
  }
  if(lX0>-99.&&lX1>-99.) lThisTrackLength = TrackLength(track, lX0, lX1, lMagneticField);
  track.propagateTo(lX0, lMagneticField);
  //Calculate expected times with smearing
  double lExpectedTimeInformation = lThisTrackLength/ Velocity(track.getP(), getMass(pdg) );
  double lMeasuredTimeInformation = gRandom->Gaus(lExpectedTimeInformation, tofRes);
  if (!makeSigma) return lMeasuredTimeInformation/0.0299792458; // return the beta
  else return ((lExpectedTimeInformation-lMeasuredTimeInformation)/tofRes);
}


void run_doubleTOFele(TString lDataPath = "./", TString lLutPath =  "./")
{
  
  TStopwatch lWatchTracking, lWatchAnalysis, lWatchSmearing;
  
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  std::cout << "\e[1;31m   ALICE 3 Strangeness tracking tool: XiCC      \e[0;00m" << std::endl;
  std::cout << "\e[1;31m***********************************************\e[0;00m" << std::endl;
  
  cout<<"Main configuration: "<<endl;
  cout<<"Data path for analysis.......: "<<lDataPath.Data()<<endl;
  cout<<"Data path for LUTs...........: "<<lLutPath.Data()<<endl;
  
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
  // const Double_t lMagneticField = field->GetBz(0,0,0);
  const Double_t lMagneticField = 20;
  cout << "Magenetic Field: " << lMagneticField << endl;

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
  
  std::uint32_t roFrame;
  std::vector<Hit>* hits = nullptr;
  itsHits.SetBranchAddress("TRKHit", &hits);



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
  // cout<<__LINE__<<endl;
  cout<<"Initializing track smearer..."<<flush;
  auto gen = new SmearO2KineGenerator( Form("%so2sim_Kine.root", lDataPath.Data()) );
  gen->Init();
  cout<<" done"<<endl;

  TH2D hist20{"hist20","hist20",100,0,4,150,0,1.5};
  TH2D hist100{"hist100","hist100",100,0,4,150,0,1.5};
  
  for (int iEvent{0}; iEvent < mcTree.GetEntriesFast(); ++iEvent) { // start event loop
    std::cout << "*************** Event " << iEvent << " ***************" << std::endl;
    mcTree.GetEvent(iEvent);

    auto particles = gen->getParticles();
    lWatchSmearing.Start(kFALSE);
    Long_t lParticlesForVertexing = 0;
  
    cout<<"N particles: "<< mcArr->size() << endl;
    cout<<"N particles: "<< particles.size() << endl;

    // loop over particles
    Long_t lCorruptedTracks=0;
    for (int iparticle = 0; iparticle < particles.size(); ++iparticle){
      cout << "N particles: " << particles.size() << endl;
      auto particle = particles[iparticle];
      
      if( TMath::Abs(particle.Eta())> 1.5) continue; //skip outside barrel
      
      // only particles to be transported, which are flagged as status code = 1
      // the particles that have been transported already have status code = 0
      if (particle.GetStatusCode() != 1) continue;
      
      
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

      double tofbeta20 = makeTOFpid(o2track,pdg, 2.0 , 20., 50.); // add true in the end for a sigma
      cout << "tofbeta20 " << tofbeta20 << endl;
      double tofbeta100 = makeTOFpid(o2track,pdg, 2.0 , 100., 20.); // add true in the end for a sigma
      hist20.Fill(o2track.getP(),tofbeta20);
      hist100.Fill(o2track.getP(),tofbeta100);
      
    }
  }
  TFile *fout = new TFile("out.root","recreate");
  hist20.Write();
  hist100.Write();
  fout->Close();
}

