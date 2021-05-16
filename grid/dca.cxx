R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libDelphesO2)

bool smear = true;
bool nsigma = true;
double Bz = 0.2;
double eMass = 0.000511;

// TOF geometry
double tof_radius = 100.;
double tof_length = 200.;
double tof_sigmat = 0.020;
double tof_sigma0 = 0.200;
// RICH params
double rich_radius = 100.; // [cm]
double rich_length = 200.; // [cm]

// Cinematic cuts on tracks
double PtCut = 0.08;
double EtaCut = 1.1;

// charm pair clasification
enum charmPairType {kIsNoCharm = 0, kIsDzeroPair,kIsDplusPair,kIsDmixedPair};

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



int charmPair(int mpid1, int mpid2)
{ if (mpid1 == 431) mpid1 = 421; // set D_s as a D0 since the decay length is very similar 123 vs 150 mum
  if (mpid2 == 431) mpid2 = 421; // set D_s as a D0 since the decay length is very similar 123 vs 150 mum
  if ((mpid1 == 421) && (mpid2 == 421)) return charmPairType::kIsDzeroPair;
  else if ((mpid1 == 411) && (mpid2 == 411)) return charmPairType::kIsDplusPair;
  else if ((mpid1 == 411) && (mpid2 == 421)) return charmPairType::kIsDmixedPair;
  else if ((mpid1 == 421) && (mpid2 == 411)) return charmPairType::kIsDmixedPair;
  else return charmPairType::kIsNoCharm;
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


// TODO: write script to do this...
void dca(
    const char *inputFile = "delphes100k.root", // one fo the LF and charm part. Has 1M charm+beauty events
    // const char *inputFile = "delphes_500kBeauty.root", // one fo the LF and charm part. Has 1M charm+beauty events
    const char *outputFile = "dcaLFcc.root" // merge output files after analysis was run to keep file size moderate
    // const char *outputFile = "dcabb.root" // merge output files after analysis was run to keep file size moderate
  )
{

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
  smearer.loadTable(13, "./lutCovm.pi.dat");
  smearer.loadTable(211, "./lutCovm.ka.dat");
  smearer.loadTable(321, "./lutCovm.pr.dat");
  smearer.loadTable(2212, "./lutCovm.mu.dat");
  // if (Bz == 0.2) {
  //   smearer.loadTable(11, "./lutCovm.el.2kG.20cm.scenario1.dat");
  //   smearer.loadTable(13, "./lutCovm.mu.2kG.20cm.scenario1.dat");
  //   smearer.loadTable(211, "./lutCovm.pi.2kG.20cm.scenario1.dat");
  //   smearer.loadTable(321, "./lutCovm.ka.2kG.20cm.scenario1.dat");
  //   smearer.loadTable(2212, "./lutCovm.pr.2kG.20cm.scenario1.dat");
  // } else if (Bz == 0.5) {
  //   smearer.loadTable(11, "./lutCovm.el.5kG.20cm.scenario1.dat");
  //   smearer.loadTable(13, "./lutCovm.mu.5kG.20cm.scenario1.dat");
  //   smearer.loadTable(211, "./lutCovm.pi.5kG.20cm.scenario1.dat");
  //   smearer.loadTable(321, "./lutCovm.ka.5kG.20cm.scenario1.dat");
  //   smearer.loadTable(2212, "./lutCovm.pr.5kG.20cm.scenario1.dat");
  // } else {
  //   std::cout << " --- invalid Bz field: " << Bz << std::endl;
  //   return;
  // }

  // histograms
  std::string title = ";#it{p}_{T} (GeV/#it{c});d_{0} (mm)";
  std::string title3d = ";#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});d_{0} (mm)";
  std::string title3d2 = ";#it{m}_{ee} (GeV/#it{c}^{2});#Delta #varphi (rad);d_{0} (mm)";
  if (nsigma) title = ";#it{p}_{T} (GeV/#it{c});n#sigma_{d_{0}}";
  if (nsigma) title3d = ";#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});n#sigma_{d_{0}}";
  if (nsigma) title3d2 = ";#it{m}_{ee} (GeV/#it{c}^{2});#Delta #varphi (rad);n#sigma_{d_{0}}";
  // Event histograms
  auto nTracks = new TH1F("nTracks",";Tracks",500,0,500);
  // Track histograms
  auto hPrimaryM = new TH1F("hPrimaryM",";MotherPDG",1000,0.5,1000.5);
  auto hPrimaryGM = new TH1F("hPrimaryGM",";MotherPDG",1000,0.5,1000.5);
  auto hCharmMm = new TH1F("hCharmMm",";MotherPDG",101,399.5,500.5);
  auto hCharmMb = new TH1F("hCharmMb",";MotherPDG",1001,3099.5,5000.5);

  auto hCharmGM = new TH1F("hCharmGM",";MotherPDG",1000,0,10000);
  auto hBeautyM = new TH1F("hBeautyM",";MotherPDG",1000,0,10000);
  auto hBeautyGM = new TH1F("hBeautyGM",";MotherPDG",1000,0,10000);
  auto hDCAxy = new TH2F("hDCAxy", title.c_str(), 400, 0., 20., 1200, -30., 30.);
  auto hPt_Eta_Phi_primary = new TH3F("hPt_Eta_Phi_primary", ";#it{p}_{T} (GeV/#it{c});#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  auto hDCAxy_primary = new TH2F("hDCAxy_primary", title.c_str(), 400, 0., 20., 1200, -30., 30.);
  auto hPt_DCARes_primary = new TH2F("hPt_DCARes_primary",";p_{T} (GeV/c); #Delta DCA (mm)", 400,0,20,600,0,0.3);
  auto hPt_DCAAbs_primary = new TH2F("hPt_DCAAbs_primary",";p_{T} (GeV/c); DCA (mm)", 400,0,20,600,-3.,3);
  auto hPt_Eta_Phi_hf = new TH3F("hPt_Eta_Phi_hf", ";p_{T} (GeV/c);#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  auto hDCAxy_secondary_heavy = new TH2F("hDCAxy_secondary_heavy", title.c_str(), 400, 0., 20., 1200, -30., 30.);
  auto hPt_Eta_Phi_charm = new TH3F("hPt_Eta_Phi_charm", ";p_{T} (GeV/c);#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  auto hDCAxy_secondary_charm = new TH2F("hDCAxy_secondary_charm", title.c_str(), 400, 0., 20., 1200, -30., 30.);

  auto hDCAxy_D0s = new TH2F("hDCAxy_D0", title.c_str(), 400, 0., 20., 1200, -30., 30.);
  auto hDCAxy_Dpm = new TH2F("hDCAxy_Dpm", title.c_str(), 400, 0., 20., 1200, -30., 30.);
  auto hDCAxy_Lc = new TH2F("hDCAxy_Lc", title.c_str(), 400, 0., 20., 1200, -30., 30.);

  auto hPt_DCARes_charm = new TH2F("hPt_DCARes_charm",";p_{T} (GeV/c); #Delta DCA (mm)", 400,0,20,600,0,0.3);
  auto hPt_DCAAbs_charm = new TH2F("hPt_DCAAbs_charm",";p_{T} (GeV/c); DCA (mm)", 400,0,20,600,-3.,3);
  auto hPt_Eta_Phi_beauty = new TH3F("hPt_Eta_Phi_beauty", ";p_{T} (GeV/c);#eta;#varphi (rad)", 200, 0., 20., 200, -10., 10.,180,0,2*TMath::Pi());
  auto hDCAxy_secondary_beauty = new TH2F("hDCAxy_secondary_beauty", title.c_str(), 400, 0., 20., 1200, -30., 30.);
  auto hPt_DCARes_beauty = new TH2F("hPt_DCARes_beauty",";p_{T} (GeV/c); #Delta DCA (mm)", 400,0,20,600,0,0.3);
  auto hPt_DCAAbs_beauty = new TH2F("hPt_DCAAbs_beauty",";p_{T} (GeV/c); DCA (mm)", 400,0,20,600,-3.,3);

  // Pair histograms
  // Binning
   // ptee
   const int n_ptee_bin_c = 202;
   Double_t ptee_bin_c[n_ptee_bin_c+1] = {};
   for(int i=0  ;i<4   ;i++) {
     ptee_bin_c[i] = 0.025 * (i-  0) +  0.0;//from 0 to 0.075 GeV/c, every 0.025 GeV/c
   }
   for(int i=4 ;i<=202  ;i++) {
     ptee_bin_c[i] = 0.05  * (i- 4) +  0.1;//from 0.1 to 10 GeV/c, evety 0.05 GeV/c
   }
   // mee
   Double_t mee_bin_c[401];
   int n_mee_bin_c = 400;
   for(int k = 0 ; k < 401; k++){
     mee_bin_c[k] = k*0.01; // 4./400.
   }
   // DCA
   Double_t dca_bin_c[201];
   int n_dca_bin_c = 200;
   for(int k = 0 ; k < 201; k++){
     dca_bin_c[k] = k*0.1; // 4./400.
   }


  auto hM_Pt_DCAprimary = new TH3F("hM_Pt_DCAprimary",title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  auto hM_Pt_DCAheavy = new TH3F("hM_Pt_DCAheavy",title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  auto hM_Pt_DCAcharm = new TH3F("hM_Pt_DCAcharm",title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  // dedicaded sources
  auto hM_Pt_DCAcharm_Dzero = new TH3F("hM_Pt_DCAcharm_Dzero",title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  auto hM_Pt_DCAcharm_Dplus = new TH3F("hM_Pt_DCAcharm_Dplus",title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  auto hM_Pt_DCAcharm_Dmixed = new TH3F("hM_Pt_DCAcharm_Dmixed",title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);
  auto hM_Pt_DCAbeauty = new TH3F("hM_Pt_DCAbeauty",title3d.c_str(),n_mee_bin_c,mee_bin_c,n_ptee_bin_c,ptee_bin_c,n_dca_bin_c,dca_bin_c);

// mee delta phi dca distributions

  auto hM_dPhi_DCAprimary = new TH3F("hM_dPhi_DCAprimary",title3d2.c_str(),300.,0,3.,20,0.,TMath::Pi(),800,0.,20.);
  auto hM_dPhi_DCAcharm = new TH3F("hM_dPhi_DCAcharm",title3d2.c_str(),300.,0,3.,20,0.,TMath::Pi(),800,0.,20.);
  auto hM_dPhi_DCAbeauty = new TH3F("hM_dPhi_DCAbeauty",title3d2.c_str(),300.,0,3.,20,0.,TMath::Pi(),800,0.,20.);

  std::vector<Track *> vecElectron,vecPositron;

  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);
    nTracks->Fill(tracks->GetEntries());
    // loop over tracks and fill vectors for electrons and positrons
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

      // smear track if requested
      if (smear) if (!smearer.smearTrack(*track)) continue; // strange syntax, but works
      // kinatic cuts on tracks
      if (!kineCuts(track)) continue;

      // look only at electrons
      if(abs(particle->PID) != 11 ) continue;
      if (particle->PID == 11 ) vecElectron.push_back(track);
      else if (particle->PID == -11 ) vecPositron.push_back(track);

      //fill track histograms
      auto D0 = track->D0;
      if (nsigma) D0 /= track->ErrorD0;

      // fill histograms
      hDCAxy->Fill(track->PT, D0);
      if (!hasStrangeAncestor(particle, particles) && !hasHeavyAncestor(particle, particles)) {
        hPrimaryM->Fill(mother->PID);
        hPrimaryGM->Fill(gmother->PID);
        hDCAxy_primary->Fill(track->PT, D0);
        hPt_Eta_Phi_primary->Fill(track->PT, track->Eta, track->Phi);
        hPt_DCARes_primary->Fill(track->PT,track->ErrorD0);
        hPt_DCAAbs_primary->Fill(track->PT,track->D0);
        // if(abs(track->D0) > 1.) cout << "\t>>> LF Mother: "<< mPid << "\t>>> LF GrandMother: " << gmPid << "\n";
      }
      else {
          if (hasCharmAncestor(particle, particles) || hasBeautyAncestor(particle, particles)){
            hDCAxy_secondary_heavy->Fill(track->PT, D0);
            hPt_Eta_Phi_hf->Fill(track->PT, track->Eta, track->Phi);
          }
          if(isCharm(mPid) && !hasBeautyAncestor(particle, particles)){
            // hCharmM->Fill(mother->PID);
            if (mPid < 1000) hCharmMm->Fill(mother->PID);
            if (mPid > 999) hCharmMb->Fill(mother->PID);
            hCharmGM->Fill(gmother->PID);
            hPt_Eta_Phi_charm->Fill(track->PT, track->Eta, track->Phi);
            hDCAxy_secondary_charm->Fill(track->PT, D0);
            hPt_DCARes_charm->Fill(track->PT,track->ErrorD0);
            hPt_DCAAbs_charm->Fill(track->PT,track->D0);
            // Fill track distributions for different species
            switch (abs(mPid)) {
              case 411:  // D+
                hDCAxy_Dpm->Fill(track->PT, D0);
                continue;
              case 421:  // D0
                hDCAxy_D0s->Fill(track->PT, D0);
                continue;
              case 431:  // Ds+
                hDCAxy_D0s->Fill(track->PT, D0);
                continue;
              case 4122: // Lambdac+
                hDCAxy_Lc->Fill(track->PT, D0);
                continue;
            }
          }
          else if (isBeauty(mPid)) {
          // else if (hasBeautyAncestor(particle, particles) && mPid > 500 && mPid < 600) {
            hBeautyM->Fill(mother->PID);
            hBeautyGM->Fill(gmother->PID);
            hPt_Eta_Phi_beauty->Fill(track->PT, track->Eta, track->Phi);
            hDCAxy_secondary_beauty->Fill(track->PT, D0);
            hPt_DCARes_beauty->Fill(track->PT,track->ErrorD0);
            hPt_DCAAbs_beauty->Fill(track->PT,track->D0);
          }
          else if (hasStrangeAncestor(particle, particles)) continue;
          // else {
          //   std::cout << particle->X << " " << particle->Y << " " << particle->Z << std::endl;
          //   std::cout << " --- mother PDG     : " << mother->PID << std::endl;
          //   std::cout << " --- grandmother PDG: " << gmother->PID << std::endl;
          // }
      }
    }

    // pairing
    TLorentzVector LV1,LV2,LV;
    double dca1,dca2,dca;
    for (auto track1 : vecElectron){
      auto particle1 = (GenParticle *)track1->Particle.GetObject();
      auto imother1 = particle1->M1;
      auto mother1 = imother1 != -1 ? (GenParticle *)particles->At(imother1) : (GenParticle *)nullptr;
      auto m1Pid = mother1->PID;
      LV1.SetPtEtaPhiM(track1->PT,track1->Eta,track1->Phi,eMass);
      dca1 = track1->D0/track1->ErrorD0;
      for (auto track2 : vecPositron) {
        auto particle2 = (GenParticle *)track2->Particle.GetObject();
        auto imother2 = particle2->M1;
        auto mother2 = imother2 != -1 ? (GenParticle *)particles->At(imother2) : (GenParticle *)nullptr;
        auto m2Pid = mother2->PID;
        LV2.SetPtEtaPhiM(track2->PT,track2->Eta,track2->Phi,eMass);
        dca2 = track2->D0/track2->ErrorD0;
        LV = LV1 + LV2;
        dca = sqrt( (dca1*dca1 + dca2*dca2)/2 );
        auto deltaPhi = TMath::Abs(track1->Phi - track2->Phi);
        if (mother1 == mother2 && !hasHeavyAncestor(particle1, particles) && !hasStrangeAncestor(particle1, particles)) // same mother and neutral LF particle, pion or eta
          {
          hM_Pt_DCAprimary->Fill(LV.Mag(),LV.Pt(),dca);
          hM_dPhi_DCAprimary->Fill(LV.Mag(),deltaPhi,dca);
          }
        if (mother1 != mother2 && hasHeavyAncestor(particle1, particles) && hasHeavyAncestor(particle2, particles))
          {
          hM_Pt_DCAheavy->Fill(LV.Mag(),LV.Pt(),dca);
          }
        if ( isCharm(m1Pid) && isCharm(m2Pid) && !hasBeautyAncestor(particle1, particles) && !hasBeautyAncestor(particle2, particles) )
        // if(mother1 != mother2  && hasCharmAncestor(particle1, particles) && hasCharmAncestor(particle2, particles))
          {
          hM_Pt_DCAcharm->Fill(LV.Mag(),LV.Pt(),dca);
          printf("mpid1: %d\t mpid2: %d\n",m1Pid,m2Pid);
          if (charmPair(abs(m1Pid),abs(m2Pid)) == charmPairType::kIsDzeroPair)  {hM_Pt_DCAcharm_Dzero->Fill(LV.Mag(),LV.Pt(),dca);}
          if (charmPair(abs(m1Pid),abs(m2Pid)) == charmPairType::kIsDplusPair)  {hM_Pt_DCAcharm_Dplus->Fill(LV.Mag(),LV.Pt(),dca);}
          if (charmPair(abs(m1Pid),abs(m2Pid)) == charmPairType::kIsDmixedPair) {hM_Pt_DCAcharm_Dmixed->Fill(LV.Mag(),LV.Pt(),dca);}
          // hM_dPhi_DCAcharm->Fill(LV.Mag(),deltaPhi,dca);
          }
        if( isBeauty(m1Pid) && isBeauty(m2Pid) )
        // if(mother1 != mother2 && hasBeautyAncestor(particle1, particles) && hasBeautyAncestor(particle2, particles) )
          {
          hM_Pt_DCAbeauty->Fill(LV.Mag(),LV.Pt(),dca);
          // hM_dPhi_DCAbeauty->Fill(LV.Mag(),deltaPhi,dca);
          }
      }
    }
    vecElectron.clear();
    vecPositron.clear();
  }

  // auto c = new TCanvas("c", "c", 800, 800);
  // c->SetLogy();
  // c->SetRightMargin(0.03);
  // c->SetTopMargin(0.03);
  // auto hDCAxy_primary_DCA = hDCAxy_primary->ProjectionY("hDCAxy_primary_DCA");
  // auto hDCAxy_hf_DCA = hDCAxy_secondary_heavy->ProjectionY("hDCAxy_hf_DCA");
  // makeHistNice(hDCAxy_primary_DCA,kBlue+1);
  // hDCAxy_primary_DCA->Rebin(2);
  // makeHistNice(hDCAxy_hf_DCA,kRed);
  // hDCAxy_hf_DCA->Rebin(2);
  // hDCAxy_primary_DCA->SetStats(0);
  // hDCAxy_primary_DCA->SetTitle("");
  // hDCAxy_primary_DCA->GetXaxis()->SetTitle("DCA_{e} (#sigma)");
  // hDCAxy_primary_DCA->GetYaxis()->SetTitle("Counts");
  // hDCAxy_primary_DCA->Draw();
  // hDCAxy_hf_DCA->Draw("same");
  // auto leg = new TLegend(0.6,0.8,0.97,0.97);
  // leg->SetBorderSize(0);
  // leg->SetFillStyle(0);
  // leg->AddEntry(hDCAxy_primary_DCA,"primary electrons","l");
  // leg->AddEntry(hDCAxy_hf_DCA,"heavy flavour electrons","l");
  // leg->Draw("same");
  // c->SaveAs("test.pdf");

  auto fout = TFile::Open(outputFile, "RECREATE");
  nTracks->Write();
  hPrimaryM->Write();
  hPrimaryGM->Write();
  hCharmMm->Write();
  hCharmMb->Write();
  hCharmGM->Write();
  hBeautyM->Write();
  hBeautyGM->Write();
  hPt_Eta_Phi_primary->Write();
  hPt_Eta_Phi_hf->Write();
  hPt_Eta_Phi_charm->Write();
  hPt_Eta_Phi_beauty->Write();
  hDCAxy->Write();
  hDCAxy_primary->Write();
  hDCAxy_secondary_charm->Write();
  hDCAxy_secondary_beauty->Write();
  hPt_DCARes_primary->Write();
  hPt_DCAAbs_primary->Write();
  hPt_DCARes_charm->Write();
  hPt_DCARes_beauty->Write();
  hPt_DCAAbs_charm->Write();
  hPt_DCAAbs_beauty->Write();
  hDCAxy_secondary_heavy->Write();
  hDCAxy_Dpm->Write();
  hDCAxy_D0s->Write();
  hDCAxy_Lc->Write();
  hM_Pt_DCAprimary->Write();
  hM_Pt_DCAheavy->Write();
  hM_Pt_DCAcharm->Write();
  hM_Pt_DCAcharm_Dzero->Write();
  hM_Pt_DCAcharm_Dplus->Write();
  hM_Pt_DCAcharm_Dmixed->Write();
  hM_Pt_DCAbeauty->Write();
  hM_dPhi_DCAprimary->Write();
  hM_dPhi_DCAcharm->Write();
  hM_dPhi_DCAbeauty->Write();
  fout->Close();

}
