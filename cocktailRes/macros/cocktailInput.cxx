R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libDelphesO2)

// Settings
bool smear = true;
bool nsigma = true;

// Magnetic field
double Bz = 0.2; // [T]

// TOF
double tof_radius = 100.; // [cm]
double tof_length = 200.; // [cm]
double tof_sigmat = 0.02; // [ns]
double tof_sigma0 = 0.20; // [ns]
// electron mass
double eMass = 0.000511;

// Kinematic cuts on tracks/particles
double PtCut  = 0.04;
double EtaCut = 1.10;

// TOF cuts on tracks
double nSigmaTOFEle   = 3.;
double nSigmaTOFOther = -5.;

void makeHistNice(TH1* h, int color){
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetMarkerStyle(20);
  h->SetLineWidth(2);
}

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

bool kineCuts(Track *tr){
  // check pt and eta for track
  // evaluate as true if criterion is passed
  bool pt = tr->PT > PtCut;
  bool eta = abs(tr->Eta) < EtaCut;
  // all have to be true
  return (pt && eta);
}

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


void cocktailInput(
  const char *inputFile = "delphes.root",// input file (need different files for different sources)
  const char *outputFile = "ana.root",   // output files (different sources to be merged externally)
  const Int_t nEvents = -1,              // number of events to analyze (all in case of -1)
  const Bool_t bStarlight = kFALSE       // special treatment in case of STARLight events
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
  if(nEvents > -1)
  numberOfEntries = nEvents;

  Printf("Process %lld out of %lld entries",numberOfEntries,numberOfEntriesTree);

  // Get pointers to branches used in this analysis
  auto events = treeReader->UseBranch("Event");
  auto tracks = treeReader->UseBranch("Track");
  auto particles = treeReader->UseBranch("Particle");

  // TOF layer
  o2::delphes::TOFLayer toflayer;
  toflayer.setup(tof_radius, tof_length, tof_sigmat,tof_sigma0);

  // smearer
  o2::delphes::TrackSmearer smearer;
  // smearer.useEfficiency(false);  // check wich version this is needed in
  if(smear){
      smearer.loadTable(11,   "./lutCovm.el.dat");
      smearer.loadTable(13,   "./lutCovm.mu.dat");
      smearer.loadTable(211,  "./lutCovm.pi.dat");
      smearer.loadTable(321,  "./lutCovm.ka.dat");
      smearer.loadTable(2212, "./lutCovm.pr.dat");
  }

  // histograms

  // Event histograms
  TH1F nTracks ("nTracks",";Tracks",2000,0,2000);

  // Track histograms
  TH1F hist_pt_ele_gen ("hist_pt_ele_gen",";p_{T}^{gen} (GeV/c)",400,0.,4);
  TH1F hist_pt_pos_gen ("hist_pt_pos_gen",";p_{T}^{gen} (GeV/c)",400,0.,4);
  TH1F hist_pt_ele_rec ("hist_pt_ele_rec",";p_{T}^{rec} (GeV/c)",400,0.,4);
  TH1F hist_pt_pos_rec ("hist_pt_pos_rec",";p_{T}^{rec} (GeV/c)",400,0.,4);

  TH1F hist_phi_ele_gen ("hist_phi_ele_gen",";#phi^{gen} (rad)",64,0.0,6.4);
  TH1F hist_phi_ele_rec ("hist_phi_ele_rec",";#phi^{rec} (rad)",64,0.0,6.4);
  TH1F hist_phi_pos_gen ("hist_phi_pos_gen",";#phi^{gen} (rad)",64,0.0,6.4);
  TH1F hist_phi_pos_rec ("hist_phi_pos_rec",";#phi^{rec} (rad)",64,0.0,6.4);

  TH1F hist_eta_ele_gen ("hist_eta_ele_gen",";#eta^{gen}",100,-4.,4);
  TH1F hist_eta_ele_rec ("hist_eta_ele_rec",";#eta^{rec}",100,-4.,4);
  TH1F hist_eta_pos_gen ("hist_eta_pos_gen",";#eta^{gen}",100,-4.,4);
  TH1F hist_eta_pos_rec ("hist_eta_pos_rec",";#eta^{rec}",100,-4.,4);

  TH2F hist_deltaPtRel_ele ("hist_deltaPtRel_ele",";p_{T}^{gen} (GeV/c);(p_{T}^{gen} - p_{T}^{rec}) / p_{T}^{gen}",400,0.,4.,300,-0.5,0.5);
  TH2F hist_deltaEta_ele ("hist_deltaEta_ele",";p_{T}^{gen} (GeV/c;#eta^{gen} - #eta^{rec})",400,0.,4.,400,-0.1,0.1);
  TH2F hist_deltaPhi_ele ("hist_deltaPhi_ele",";p_{T}^{gen} (GeV/c);#phi^{gen} - #phi^{rec}",400,0.,4.,800,-0.1,0.4);

  TH2F hist_deltaPtRel_pos ("hist_deltaPtRel_pos",";p_{T}^{gen} (GeV/c);(p_{T}^{gen} - p_{T}^{rec}) / p_{T}^{gen}",400,0.,4.,300,-0.5,0.5);
  TH2F hist_deltaEta_pos ("hist_deltaEta_pos",";p_{T}^{gen} (GeV/c;#eta^{gen} - #eta^{rec})",400,0.,4.,400,-0.1,0.1);
  TH2F hist_deltaPhi_pos ("hist_deltaPhi_pos",";p_{T}^{gen} (GeV/c);#phi^{gen} - #phi^{rec}",400,0.,4.,400,-0.4,0.1);

  TH3F hist_PtEtaPhi_ele_gen ("hist_PtEtaPhi_ele_gen",";p_{T,e} (GeV/c);#eta;#phi (rad)",400,0.,4.,400,-10.,10.,64,0.,6.4);
  TH3F hist_PtEtaPhi_pos_gen ("hist_PtEtaPhi_pos_gen",";p_{T,e} (GeV/c);#eta;#phi (rad)",400,0.,4.,400,-10.,10.,64,0.,6.4);
  TH3F hist_PtEtaPhi_ele_rec ("hist_PtEtaPhi_ele_rec",";p_{T,e} (GeV/c);#eta;#phi (rad)",400,0.,4.,400,-10.,10.,64,0.,6.4);
  TH3F hist_PtEtaPhi_pos_rec ("hist_PtEtaPhi_pos_rec",";p_{T,e} (GeV/c);#eta;#phi (rad)",400,0.,4.,400,-10.,10.,64,0.,6.4);

  double ptRec,ptGen,etaRec,etaGen,phiRec,phiGen;
  double deltaPtRel,deltaEta,deltaPhi;

  int nTracksCounter = 0;
  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {
    nTracksCounter = 0;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);
    // nTracks.Fill(tracks->GetEntries()); // should be filled after accaptance cuts.

    // loop over tracks and fill vectors for electrons and positrons
    for (Int_t itrack = 0; itrack < tracks->GetEntries(); ++itrack) {

      // get track and corresponding particle
      auto track = (Track *)tracks->At(itrack);
      auto particle = (GenParticle *)track->Particle.GetObject();

      // only use tracks that go to the last layer
      // auto radius = sqrt((track->XOuter)*(track->XOuter)+(track->YOuter)*(track->YOuter));
      // if(radius < 1000.) continue;


      // smear track if requested
      if (smear) if (!smearer.smearTrack(*track)) continue; // strange syntax, but works

      ptRec  = track->PT;
      etaRec = track->Eta;
      phiRec = track->Phi;
      ptGen  = particle->PT;
      etaGen = particle->Eta;
      phiGen = TVector2::Phi_0_2pi(particle->Phi);

      deltaPtRel = (ptGen-ptRec) / ptGen;
      deltaEta   = etaGen-etaRec;
      deltaPhi   = phiGen-phiRec;
      // fill histograms

      // if (!kineCuts(track)) continue;
      nTracksCounter++;
      if (particle->PID == 11 )
      {
        hist_deltaPtRel_ele.Fill(ptGen,deltaPtRel);
        hist_deltaEta_ele.Fill(ptGen,deltaEta);
        hist_deltaPhi_ele.Fill(ptGen,deltaPhi);

        hist_pt_ele_gen.Fill(ptGen);
        hist_eta_ele_gen.Fill(etaGen);
        hist_phi_ele_gen.Fill(phiGen);
        hist_PtEtaPhi_ele_gen.Fill(ptGen,etaGen,phiGen);
        hist_pt_ele_rec.Fill(ptRec);
        hist_eta_ele_rec.Fill(etaRec);
        hist_phi_ele_rec.Fill(phiRec);
        hist_PtEtaPhi_ele_rec.Fill(ptRec,etaRec,phiRec);

      }
      // fill vectors
      else if (particle->PID == -11 )
      {
        hist_deltaPtRel_pos.Fill(ptGen,deltaPtRel);
        hist_deltaEta_pos.Fill(ptGen,deltaEta);
        hist_deltaPhi_pos.Fill(ptGen,deltaPhi);

        hist_pt_pos_gen.Fill(ptGen);
        hist_eta_pos_gen.Fill(etaGen);
        hist_phi_pos_gen.Fill(phiGen);
        hist_PtEtaPhi_pos_gen.Fill(ptGen,etaGen,phiGen);
        hist_pt_pos_rec.Fill(ptRec);
        hist_eta_pos_rec.Fill(etaRec);
        hist_phi_pos_rec.Fill(phiRec);
        hist_PtEtaPhi_pos_rec.Fill(ptRec,etaRec,phiRec);

      }

    }
    nTracks.Fill(nTracksCounter);
  }

  // write histograms to file
  auto fout = TFile::Open(outputFile, "RECREATE");
  nTracks.Write();
  hist_pt_ele_gen.Write();
  hist_pt_pos_gen.Write();
  hist_pt_ele_rec.Write();
  hist_pt_pos_rec.Write();
  hist_phi_ele_gen.Write();
  hist_phi_ele_rec.Write();
  hist_phi_pos_gen.Write();
  hist_phi_pos_rec.Write();
  hist_eta_ele_gen.Write();
  hist_eta_ele_rec.Write();
  hist_eta_pos_gen.Write();
  hist_eta_pos_rec.Write();
  hist_deltaPtRel_ele.Write();
  hist_deltaEta_ele.Write();
  hist_deltaPhi_ele.Write();
  hist_deltaPtRel_pos.Write();
  hist_deltaEta_pos.Write();
  hist_deltaPhi_pos.Write();
  hist_PtEtaPhi_ele_gen.Write();
  hist_PtEtaPhi_pos_gen.Write();
  hist_PtEtaPhi_ele_rec.Write();
  hist_PtEtaPhi_pos_rec.Write();

  fout->Close();



  // stop watch
  watch->Stop();
  watch->Print();

}
