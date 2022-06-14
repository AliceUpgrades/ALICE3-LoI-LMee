R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libDelphesO2)

// bool nsigma = true;
double Bz = 0.2;
constexpr double eMass = 0.000511;

// TOF geometry
double tof1_radius = 100.; // [cm]
double tof1_length = 200.; // [cm]
double tof1_sigmat = 0.02; // [ns]
double tof1_sigma0 = 0.20; // [ns]


// Cinematic cuts on tracks
double PtCut = 0.04;
double EtaCut = 1.1;


bool kineCuts(Track* tr)
{
  // check pt and eta for track
  // evaluate as true if criterion is passed
  bool pt = tr->PT > PtCut;
  bool eta = fabs(tr->Eta) < EtaCut;
  // all have to be true
  return (pt && eta);
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
  // Use a step by step propagation to estimate the track length
  // start with two space points
  std::array<float, 3> lPointN;
  std::array<float, 3> lPointNplus;
  // track length is 0 in the beginning
  Double_t lLength = 0.0;
  // the number of steps we take to estimate (more is better, but slower, maybe could be checked)
  Int_t lNStepsLIntegrator = 200;
  // propagate the track to the starting point
  track.propagateTo(lX0, lMagneticField);
  for(Int_t iStep=1; iStep<lNStepsLIntegrator; iStep++){
    // get the global point
    track.getXYZGlo(lPointN);
    // calculate next position
    Float_t lPosition = lX0 + (lX1-lX0)*((Float_t)(iStep))/((Float_t)(lNStepsLIntegrator-1));
    // propagate to next position
    track.propagateTo(lPosition, lMagneticField);
    // get global coordinats at new position
    track.getXYZGlo(lPointNplus);
    // calculate distance of the space points and add it to the length
    lLength += std::hypot( lPointNplus[0]-lPointN[0], lPointNplus[1]-lPointN[1], lPointNplus[2]-lPointN[2] );
  }
  return lLength;
}

double makeTOFpid(o2::track::TrackParCov track,int pdg, double lMagneticField = 0.0, double tofPos = 100., double tofRes = 50., bool makeSigma = false, int refPdg = 0)
{
  if((pdg == 0) && makeSigma) {
    cout << "please specify reference particle for n sigma calculation" << endl;
    return 0.;
    }
  // create a vertex point in (0|0|0)
  o2::math_utils::Point3D<float> zeropos{0,0,0};
  std::array<float, 6> zerocov;
  // create a coverance matrix with tiny values
  for(Int_t jj=0; jj<6; jj++) zerocov[jj]=1e-6;
  // create vertex at (0|0|0) with small uncertainties
  o2::dataformats::VertexBase zerovtx(zeropos, zerocov);
  // initiate track length
  double lThisTrackLength = -1.;
  // check that B field has a reasonable value
  if(lMagneticField == 0.0) {
    cout << "Please set magnetic field properly!" << flush;
    return 0;
  }
  // initiate local coordinates
  Float_t lX0=-100, lX1=-100;
  // propagate track to DCA (make sure we start in the beginning)
  if (!track.propagateToDCA(zerovtx, lMagneticField)) {
    //std::cout << "Propagation failed." << std::endl;
  } else {
    lX0 = track.getX();
  }
  // Does the track reach the position we have the TOF at
  if (!track.getXatLabR(tofPos,lX1,lMagneticField,o2::track::DirOutward)) {
    lX1 = -100;
  }
  if(lX0>-99.&&lX1>-99.) lThisTrackLength = TrackLength(track, lX0, lX1, lMagneticField);
  track.propagateTo(lX0, lMagneticField);
  //Calculate expected time
  double lExpectedTimeInformation = lThisTrackLength/ Velocity(track.getP(), getMass(pdg) );
  //Calculate expected times with smearing
  double lMeasuredTimeInformation = gRandom->Gaus(lExpectedTimeInformation, tofRes);
  // Calculate time for mass hypothesis (needed for the n sigma calculation)
  double lTimeInformationHypothesis = lThisTrackLength/ Velocity(track.getP(), getMass(refPdg) ); // this is the one for the nsigma
  // Calculate the beta for measured and expected
  double betaMeasured = lThisTrackLength/lMeasuredTimeInformation/0.0299792458   ; // return the beta l/t/c 
  double betaExpected = lThisTrackLength/lTimeInformationHypothesis/0.0299792458 ;
  //return either beta or nsigma value
  if (!makeSigma) return betaMeasured;
  else return ((lMeasuredTimeInformation - lTimeInformationHypothesis)/tofRes);
}


// TODO: write script to do this...
void doubleTOF(const char* inputFile, const char* outputFile = "output.root")
{
  // if(!inputFile) pritf("No input file specified.\n", );
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  auto treeReader = new ExRootTreeReader(&chain);
  auto numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  auto events = treeReader->UseBranch("Event");
  auto tracks = treeReader->UseBranch("Track");
  auto particles = treeReader->UseBranch("Particle");

  // smearer
  o2::delphes::TrackSmearer smearer;
  smearer.loadTable(11, "./lutCovm.el.dat");
  smearer.loadTable(13, "./lutCovm.mu.dat");
  smearer.loadTable(211, "./lutCovm.pi.dat");
  smearer.loadTable(321, "./lutCovm.ka.dat");
  smearer.loadTable(2212, "./lutCovm.pr.dat");

  // Event histograms
  auto nTracks = new TH1F("nTracks", ";Tracks", 101, -0.5, 100.5);

  // Track histograms
  auto hPt = new TH1F("hPt", ";p_T (GeV/c)", 200, 0, 20);
  auto hBeta1 = new TH2F("hBeta1", ";p (GeV/c);#beta", 200, 0, 10,200,0,1.3);
  auto hBeta20 = new TH2F("hBeta20", ";p (GeV/c);#beta", 200, 0, 10,200,0,1.3);
  auto hBeta50 = new TH2F("hBeta50", ";p (GeV/c);#beta", 200, 0, 10,200,0,1.3);
  
  auto hBetaEleSig = new TH2F("hBetaEleSig", ";p (GeV/c);#sigma_{e}", 200, 0, 10,400,-10,10);
  auto hBetaPiSig = new TH2F("hBetaPiSig", ";p (GeV/c);#sigma_{#pi}", 200, 0, 10,400,-10,10);

  TRandom3 rndm3;
  rndm3.SetSeed(12345);

  int eventCounter = 1;
  TRandom3 rnd;
  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);
    nTracks->Fill(tracks->GetEntries());
    // loop over tracks and fill vectors for electrons and positrons
    for (Int_t itrack = 0; itrack < tracks->GetEntries(); ++itrack) {
      // get track and corresponding particle
      auto track = (Track*)tracks->At(itrack);
      auto particle = (GenParticle*)track->Particle.GetObject();
      int pdg = particle->PID;

      // cut away tracks that are way off.
      if (fabs(track->D0/track->ErrorD0) > 3.) continue; // adopt to just stay in the beampipe?
      if (fabs(track->DZ/track->ErrorDZ) > 3.) continue;

      O2Track o2track;
      o2::delphes::TrackUtils::convertTrackToO2Track(*track, o2track, false);
      // smear track
      if (!smearer.smearTrack(*track)) {cout << "no smearing work. pdg: " << pdg << endl; continue;}
      hPt->Fill(track->PT);


      // lets try the TOF pid
      double tofpid1 = makeTOFpid(o2track,pdg, 5. , 100., 1., false);
      double tofpid20 = makeTOFpid(o2track,pdg, 5. , 100., 20., false);
      double tofpid50 = makeTOFpid(o2track,pdg, 5. , 100., 50., false);
      hBeta1->Fill(track->P,tofpid1);
      hBeta20->Fill(track->P,tofpid20);
      hBeta50->Fill(track->P,tofpid50);

      hBetaEleSig->Fill(track->P,makeTOFpid(o2track,pdg, 5. , 100., 20., true, 11));
      hBetaPiSig->Fill(track->P,makeTOFpid(o2track,pdg, 5. , 100., 20., true, 211));

    }

  }
    
  auto fout = TFile::Open(outputFile, "RECREATE");
  nTracks->Write();
  hPt->Write();
  hBeta1->Write();
  hBeta20->Write();
  hBeta50->Write();
  hBetaEleSig->Write();
  hBetaPiSig->Write();
  fout->Close();
}
