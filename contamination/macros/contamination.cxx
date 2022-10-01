R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libDelphesO2)

bool pairing = true;
bool smear = true;
// bool nsigma = true;
double Bz = 0.2;
constexpr double eMass = 0.000511;

// TOF geometry
double tof_radius = 100.; // [cm]
double tof_length = 200.; // [cm]
double tof_sigmat = 0.02; // [ns]
double tof_sigma0 = 0.20; // [ns]
// RICH params
double rich_radius = 100.; // [cm]
double rich_length = 200.; // [cm]

// Cinematic cuts on tracks
double PtCut = 0.04;
double EtaCut = 1.1;

double pionRej = 1. / 1000.;

void makeHistNice(TH1* h, int color)
{
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetMarkerStyle(20);
  h->SetLineWidth(2);
}

bool kineCuts(Track* tr)
{
  // check pt and eta for track
  // evaluate as true if criterion is passed
  bool pt = tr->PT > PtCut;
  bool eta = fabs(tr->Eta) < EtaCut;
  // all have to be true
  return (pt && eta);
}

bool isSelectedMC(Track* tr)
{
  auto particle = (GenParticle*)tr->Particle.GetObject();
  if (abs(particle->PID) == 11)
    return true;
  else
    return false;
}

// bool isSelectedTOF(Track *tr, TOFLayer *tof)
// {
//
//   return false;
// }
//
// bool isSelectedRICH(Track *tr)
// {
//   return false;
// }

bool hasCommonAncestor(GenParticle* p1, GenParticle* tr2)
{
  // check if 2 particles have one ancestor
  // fill vector with all ancestors of particle on
  // check if particle2 is included
  // if so, return true if not check for mother of p2
  // maybe do this recursiv???
  return false;
}

bool hasStrangeAncestor(GenParticle* particle, TClonesArray* particles)
{
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle*)particles->At(imother);
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

bool hasHeavyAncestor(GenParticle* particle, TClonesArray* particles)
{
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle*)particles->At(imother);
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

bool isBeauty(int pid)
{
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

bool hasBeautyAncestor(GenParticle* particle, TClonesArray* particles)
{
  //check until beauty ancestor is found
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle*)particles->At(imother);
  auto pid = mother->PID;
  if (isBeauty(pid)) return true;
  return hasBeautyAncestor(mother, particles);
}

bool isCharm(int pdg)
{
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

bool hasCharmAncestor(GenParticle* particle, TClonesArray* particles)
{
  auto imother = particle->M1;
  if (imother == -1) return false;
  auto mother = (GenParticle*)particles->At(imother);
  auto pid = mother->PID;
  auto bCharm = false;
  bCharm = isCharm(pid);

  if (hasBeautyAncestor(mother, particles)) {
    bCharm = false;
  }
  return bCharm;
}

// TODO: write script to do this...
void contamination(const char* inputFile, const char* outputFile = "output.root")
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

  // TOF layer
  o2::delphes::TOFLayer toflayer;
  toflayer.setup(tof_radius, tof_length, tof_sigmat, tof_sigma0);
  // RICH detector
  o2::delphes::RICHdetector richdetector;
  richdetector.setup(rich_radius, rich_length);
  richdetector.setIndex(1.03);
  richdetector.setRadiatorLength(2.);
  richdetector.setEfficiency(0.4);
  richdetector.setSigma(7.e-3);

  // smearer
  o2::delphes::TrackSmearer smearer;
  smearer.loadTable(11, "./lutCovm.el.dat");
  smearer.loadTable(13, "./lutCovm.mu.dat");
  smearer.loadTable(211, "./lutCovm.pi.dat");
  smearer.loadTable(321, "./lutCovm.ka.dat");
  smearer.loadTable(2212, "./lutCovm.pr.dat");

  // Event histograms
  auto nTracks = new TH1F("nTracks", ";Tracks", 10000, 0, 10000);
  auto nTracksCent = new TH1F("nTracksCent", ";Tracks", 10000, 0, 10000);
  auto nTracksEle = new TH1F("nTracksEle", ";Tracks", 1000, 0, 1000);
  auto nTracksPos = new TH1F("nTracksPos", ";Tracks", 1000, 0, 1000);

  // Histograms TOF
  auto hTime0 = new TH1F("hTime0", ";t_{0} (ns)", 1000, -1., 1.);
  auto hBetaP = new TH2F("hBetaP", ";#it{p} (GeV/#it{c});#beta", 400, 0., 4., 1000, 0.1, 1.1);

  TH2* hNsigmaP_tof[5];
  TH2* hNsigmaP_tof_after[5];
  const char* pname[5] = {"el", "mu", "pi", "ka", "pr"};
  const char* plabel[5] = {"e", "#mu", "#pi", "K", "p"};
  for (int i = 0; i < 5; ++i) {
    hNsigmaP_tof[i] = new TH2F(Form("hNsigmaP_tof_%s", pname[i]), Form(";#it{p} (GeV/#it{c});n#sigma_{%s}", plabel[i]), 400, 0., 4., 250, -25., 25.);
    hNsigmaP_tof_after[i] = new TH2F(Form("hNsigmaP_tof_after_%s", pname[i]), Form(";#it{p} (GeV/#it{c});n#sigma_{%s}", plabel[i]), 400, 0., 4., 250, -25., 25.);
  }

  // auto hBetaP = new TH2F("hBetaP",";p (GeV/c);#beta",400,0.,4.,1000,0.1,1.1);
  // auto tofNsigma = new TH2F("tofNsigma",";pT (GeV/c); n#sigma_{TOF}",400,0,4,200,-10.,10.);
  auto tofNsigmaEleCut = new TH2F("tofNsigmaEleCut", ";pT (GeV/c); n#sigma_{TOF}", 400, 0, 4, 250, -10., 10.);
  auto tofNsigmaPionCut = new TH2F("tofNsigmaPionCut", ";pT (GeV/c); n#sigma_{TOF}", 400, 0, 4, 250, -10., 10.);

  // histograms RICH
  auto hAngleP = new TH2F("hAngleP", ";#it{p} (GeV/#it{c});#theta (rad)", 200, 0, 20, 250, 0., 0.25);
  TH2* hNsigmaP_rich[5];
  TH2* hNsigmaP_rich_after[5];
  std::map<int, int> pidmap = {{11, 0}, {13, 1}, {211, 2}, {321, 3}, {2212, 4}};
  std::map<int, double> pidmass = {{0, 0.00051099891}, {1, 0.10565800}, {2, 0.13957000}, {3, 0.49367700}, {4, 0.93827200}};
  for (int i = 0; i < 5; ++i) {
    hNsigmaP_rich[i] = new TH2F(Form("hNsigmaP_rich_%s", pname[i]), Form(";#it{p} (GeV/#it{c});n#sigma_{%s}", plabel[i]), 200, 0, 20, 250, -25., 25.);
    hNsigmaP_rich_after[i] = new TH2F(Form("hNsigmaP_rich_after_%s", pname[i]), Form(";#it{p} (GeV/#it{c});n#sigma_{%s}", plabel[i]), 200, 0, 20, 250, -25., 25.);
  }

  // Track histograms
  auto hPt_trackEle = new TH1F("hPt_trackEle", ";p_T (GeV/c)", 200, 0, 20);
  auto hPt_trackEle_trueEle = new TH1F("hPt_trackEle_trueEle", ";p_T (GeV/c)", 200, 0, 20);
  auto hPt_trackPID = new TH1F("hPt_trackPID", ";p_T (GeV/c)", 200, 0, 20);
  // check the pid codes of Mothers
  auto hPdg_mother = new TH1F("hPdg_mother", ";pdg code mother", 601, -0.5, 600.5);

  // Pair histograms
  auto hM_Pt_sameMother = new TH2F("hM_Pt_sameMother", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0., 3., 400, 0., 4.);
  auto hM_Pt_ULS = new TH2F("hM_Pt_ULS", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);
  auto hM_Pt_LSplus = new TH2F("hM_Pt_LSplus", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);
  auto hM_Pt_LSminus = new TH2F("hM_Pt_LSminus", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);
  auto hM_Pt_ULS_trueEle = new TH2F("hM_Pt_ULS_trueEle", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);
  auto hM_Pt_LSplus_trueEle = new TH2F("hM_Pt_LSplus_trueEle", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);
  auto hM_Pt_LSminus_trueEle = new TH2F("hM_Pt_LSminus_trueEle", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);

  auto hM_Pt_ULS_cont = new TH2F("hM_Pt_ULS_cont", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);
  auto hM_Pt_LSplus_cont = new TH2F("hM_Pt_LSplus_cont", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);
  auto hM_Pt_LSminus_cont = new TH2F("hM_Pt_LSminus_cont", ";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)", 300., 0, 3., 400, 0., 4.);

  // auto hM_Pt_sameMother = new TH2F("hM_Pt_sameMother",";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)",300.,0.,3.,400,0.,4.);
  // auto hM_Pt_ULS = new TH2F("hM_Pt_ULS",";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)",300.,0,3.,400,0.,4.);
  // auto hM_Pt_LSplus = new TH2F("hM_Pt_LSplus",";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)",300.,0,3.,400,0.,4.);
  // auto hM_Pt_LSminus = new TH2F("hM_Pt_LSminus",";m_{ee} (Gev/c^2);p_{T,ee} (GeV/c)",300.,0,3.,400,0.,4.);

  std::vector<Track*> vecElectron_mcTruth, vecElectron, vecPositron_mcTruth, vecPositron, vecPIDtracks, vecTOFtracks;

  TRandom3 rndm3;
  rndm3.SetSeed(12345);

  int eventCounter = 1;
  TRandom3 rnd;
  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);

    // loop over tracks and fill vectors for electrons and positrons
    for (Int_t itrack = 0; itrack < tracks->GetEntries(); ++itrack) {
      // get track and corresponding particle
      // get track and corresponding particle
      auto track = (Track*)tracks->At(itrack);
      auto particle = (GenParticle*)track->Particle.GetObject();

      // smear track
      if (!smearer.smearTrack(*track)) continue;

      // cut away tracks that are way off.
      if (fabs(track->D0) > 0.4) continue; // adopt to just stay in the beampipe?
      if (fabs(track->DZ) > 3.) continue;

      // check if has TOF
      if (toflayer.hasTOF(*track)) vecTOFtracks.push_back(track);

      // check if has RICH
      // if (!richdetector.hasRICH(*track)) continue;

      // push all tracks into a vector.
      vecPIDtracks.push_back(track);
    }

    nTracks->Fill(vecPIDtracks.size());
    if (vecPIDtracks.size() < 3750) {
      vecPIDtracks.clear();
      continue;
    } // dirty dirty centrality
    nTracksCent->Fill(vecPIDtracks.size());

    std::array<float, 2> tzero;
    toflayer.eventTime(vecTOFtracks, tzero);
    hTime0->Fill(tzero[0]);
    vecTOFtracks.clear();

    for (auto track : vecPIDtracks) {
      auto particle = (GenParticle*)track->Particle.GetObject();
      auto pid = particle->PID;

      if ((abs(pid) != 11) && (fabs(track->Charge) != 0)) {
        if (rnd.Rndm() > 0.001) continue;
      }

      // TOF PID
      auto p = track->P;
      auto beta = toflayer.getBeta(*track);
      hBetaP->Fill(p, beta);

      // fill nsigma
      std::array<float, 5> deltat, nsigmaTOF;
      toflayer.makePID(*track, deltat, nsigmaTOF);
      for (int i = 0; i < 5; ++i)
        hNsigmaP_tof[i]->Fill(p, nsigmaTOF[i]);

      // RICH PID
      auto measurement = richdetector.getMeasuredAngle(*track);
      auto angle = measurement.first;
      auto anglee = measurement.second;
      if (anglee == 0.) continue;

      hAngleP->Fill(p, angle);

      // make pid
      std::array<float, 5> deltaangle, nsigmaRICH;
      richdetector.makePID(*track, deltaangle, nsigmaRICH);
      for (int i = 0; i < 5; ++i) {
        hNsigmaP_rich[i]->Fill(p, nsigmaRICH[i]);
      }

      // fill pid histos after cuts
      for (int i = 0; i < 5; ++i) {
        hNsigmaP_tof_after[i]->Fill(p, nsigmaTOF[i]);
        hNsigmaP_rich_after[i]->Fill(p, nsigmaRICH[i]);
      }
      hPt_trackPID->Fill(track->PT);
      if (abs(pid) == 11) {
        hPt_trackEle_trueEle->Fill(track->PT);
      }

      if (track->Charge < 0.) {
        if (pairing) {
          vecElectron.push_back(track);
          if (abs(pid) == 11) vecElectron_mcTruth.push_back(track);
        }
      } else if (track->Charge > 0.) {
        if (pairing) {
          vecPositron.push_back(track);
          if (abs(pid) == 11) vecPositron_mcTruth.push_back(track);
        }
      }
    }
    vecPIDtracks.clear();

    nTracksEle->Fill(vecElectron_mcTruth.size());
    nTracksPos->Fill(vecPositron_mcTruth.size());
    // pairing

    if (pairing) {
      // be lazy and write lambda
      auto pairULS = [](const std::vector<Track*>& v1, const std::vector<Track*>& v2, TH2* hist) -> void {
        TLorentzVector LV1, LV2;
        for (auto& t1 : v1) {
          LV1.SetPtEtaPhiM(t1->PT, t1->Eta, t1->Phi, eMass);
          for (auto& t2 : v2) {
            LV2.SetPtEtaPhiM(t2->PT, t2->Eta, t2->Phi, eMass);
            hist->Fill((LV1 + LV2).Mag(), (LV1 + LV2).Pt());
          }
        }
      };
      pairULS(vecElectron_mcTruth, vecPositron_mcTruth, hM_Pt_ULS_trueEle);
      pairULS(vecElectron, vecPositron, hM_Pt_ULS);
      // and another lambda for the LS spectra
      auto pairLS = [](const std::vector<Track*>& v, TH2* hist) -> void {
        TLorentzVector LV1, LV2;
        for (auto t1 = v.begin(); t1 != v.end(); ++t1) {
          for (auto t2 = t1 + 1; t2 != v.end(); ++t2) {
            LV1.SetPtEtaPhiM((*t1)->PT, (*t1)->Eta, (*t1)->Phi, eMass);
            LV2.SetPtEtaPhiM((*t2)->PT, (*t2)->Eta, (*t2)->Phi, eMass);
            hist->Fill((LV1 + LV2).Mag(), (LV1 + LV2).Pt());
          }
        }
      };
      pairLS(vecPositron_mcTruth, hM_Pt_LSplus_trueEle);
      pairLS(vecPositron, hM_Pt_LSplus);
      pairLS(vecElectron_mcTruth, hM_Pt_LSminus_trueEle);
      pairLS(vecElectron, hM_Pt_LSminus);
      // We also want to ask if at least one is from contamination
      auto pairULScont = [](const std::vector<Track*>& v1, const std::vector<Track*>& v2, TH2* hist) -> void {
        TLorentzVector LV1, LV2;
        for (auto& t1 : v1) {
          auto p1 = (GenParticle *)t1->Particle.GetObject();
          auto pid1 = fabs(p1->PID);
          LV1.SetPtEtaPhiM(t1->PT, t1->Eta, t1->Phi, eMass);
          for (auto& t2 : v2) {
            auto p2 = (GenParticle *)t2->Particle.GetObject();
            auto pid2 = fabs(p2->PID);
            if((pid1 == 11) && (pid2 == 11)) continue; // at least one should not be an electron
            LV2.SetPtEtaPhiM(t2->PT, t2->Eta, t2->Phi, eMass);
            hist->Fill((LV1 + LV2).Mag(), (LV1 + LV2).Pt());
          }
        }
      };
      pairULScont(vecElectron, vecPositron, hM_Pt_ULS_cont);
      auto pairLScont = [](const std::vector<Track*>& v, TH2* hist) -> void {
        TLorentzVector LV1, LV2;
        for (auto t1 = v.begin(); t1 != v.end(); ++t1) {
          auto p1 = (GenParticle *)(*t1)->Particle.GetObject();
          auto pid1 = fabs(p1->PID);
          for (auto t2 = t1 + 1; t2 != v.end(); ++t2) {
            auto p2 = (GenParticle *)(*t2)->Particle.GetObject();
            auto pid2 = fabs(p2->PID);
            if((pid1 == 11) && (pid2 == 11)) continue; // at least one should not be an electron
            LV1.SetPtEtaPhiM((*t1)->PT, (*t1)->Eta, (*t1)->Phi, eMass);
            LV2.SetPtEtaPhiM((*t2)->PT, (*t2)->Eta, (*t2)->Phi, eMass);
            hist->Fill((LV1 + LV2).Mag(), (LV1 + LV2).Pt());
          }
        }
      };
      pairLScont(vecPositron, hM_Pt_LSplus_cont);
      pairLScont(vecElectron, hM_Pt_LSminus_cont);



    }
    vecElectron_mcTruth.clear();
    vecElectron.clear();
    vecPositron_mcTruth.clear();
    vecPositron.clear();
    eventCounter++;
  }
  auto fout = TFile::Open(outputFile, "RECREATE");
  nTracks->Write();
  nTracksCent->Write();
  nTracksEle->Write();
  nTracksPos->Write();
  hPt_trackEle_trueEle->Write();
  hPt_trackPID->Write();
  hM_Pt_sameMother->Write();
  hM_Pt_ULS->Write();
  hM_Pt_ULS_trueEle->Write();
  hM_Pt_ULS_cont->Write();
  hM_Pt_LSplus->Write();
  hM_Pt_LSplus_trueEle->Write();
  hM_Pt_LSplus_cont->Write();
  hM_Pt_LSminus->Write();
  hM_Pt_LSminus_trueEle->Write();
  hM_Pt_LSminus_cont->Write();
  hTime0->Write();
  hBetaP->Write();
  hAngleP->Write();
  hPdg_mother->Write();
  for (int i = 0; i < 5; ++i) {
    hNsigmaP_tof[i]->Write();
    hNsigmaP_rich[i]->Write();
    hNsigmaP_tof_after[i]->Write();
    hNsigmaP_rich_after[i]->Write();
  }
  fout->Close();
}