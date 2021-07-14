// use bools to select which  histograms are upperIntEdge

bool bPlotEfficiency = kTRUE;
bool bPlotPIDhistograms = kTRUE;
  bool bPlotNSigmaProjections = kTRUE;
bool bPlotTrackContamination = kTRUE;

bool bPlotPairHistograms = kTRUE;
  bool bPlotULS = kTRUE;
  bool bPlotLS = kTRUE;

int ith_PIDscenario = 1;
TString strPIDscenario[] = {"TOF", "TOF+RICH (3#sigma_{#pi}^{RICH} rej)", "TOF+RICH (3.5#sigma_{#pi}^{RICH} rej)", "TOF+RICH (4#sigma_{#pi}^{RICH} rej)"};
// TString strPIDscenario[] = {"TOF only 0.04<pte", "TOF only 0.08<pte"};
// TString strPIDscenario[] = {"TOF+RICH (4#sigma_{#pi} 0.2<pte)", "TOF+RICH (4#sigma_{#pi} 0.08<pte)"};

std::vector<Double_t> vec_proj_bin_p = {0.0, 0.3, 0.5, 0.7, 1.0, 2.0, 4.0, 10.0};
std::vector<Double_t> vec_proj_bin_pt = {0.0, 0.3, 0.5, 0.7, 1.0, 2.0, 4.0, 10.0};
std::vector<Double_t> vec_proj_bin_mass = {0.0, 3.0}; // Intervalls for projection in mass slices

// Double_t pt_bin_proj[]  = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0};
Double_t pt_bin_proj[]  = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,4.0};
// Double_t pt_bin_proj[]  = {0.0,0.1,0.5,1.,4};
Int_t nbinspt_proj  = sizeof(pt_bin_proj)/sizeof(*pt_bin_proj) -1;

void makeHistNice(TH1* h, int color){
  h->SetTitle("");
  h->SetStats(0);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetMarkerStyle(20);
  h->SetLineWidth(2);
}

void make3HistNice(TH1* h, int color){
  h->SetTitle("");
  h->SetStats(0);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetMarkerStyle(20);
  // h->SetMarkerSize(0.6);
  h->SetLineWidth(2);
}

void makeHistNiceTH2(TH2* h, int color){
  h->SetTitle("");
  h->SetStats(0);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->SetMarkerStyle(20);
  h->SetLineWidth(2);
}

void SetTextSettings(TLatex* text, Double_t textSize){
  text->SetNDC();
  text->SetTextColor(1);
  text->SetTextSize(textSize);
  text->SetTextFont(43);
}

void SetStyleAlpha(TH1* h, int color, double alpha)
{
  h->SetStats(0);
  h->SetTitle("");
  h->SetMarkerColor(color);
  h->SetLineColorAlpha(color,0); //make lines transparant
  h->SetFillColorAlpha(color,alpha);
  h->SetFillStyle(1001);
  h->SetMarkerStyle(20);
  h->SetLineWidth(2);
}

void normalizeToBinWidth(TH1* h){
  h->Scale(1.,"width");
  // int nBins = h->GetNbinsX();
  // for (size_t iBin = 1; iBin < nBins+1; iBin++) {
  //   h->SetBinContent(iBin, h->GetBinContent(iBin)/h->GetXaxis()->GetBinWidth(iBin));
  // }
}


void anaPlots(TString inputFile)
{
  TFile *fIn  = TFile::Open(inputFile.Data());

  // select which PID cenario (written in the input file) should be analyzed
  TString pathRecPIDScenario = Form("reconstructed/PIDscenario_%i", ith_PIDscenario);
  TString pathGenPIDScenario = Form("generated/PIDscenario_%i", ith_PIDscenario);

  double BField;
  TString collSystem;
  Bool_t findBField;
  Bool_t findCollSystem;
  findBField = inputFile.Contains("B=0.2");
  if (findBField) BField = 0.2;
  findBField = inputFile.Contains("B=0.5");
  if (findBField) BField = 0.5;
  findCollSystem = inputFile.Contains("PbPb");
  if (findCollSystem) collSystem = "Pb-Pb";
  findCollSystem = inputFile.Contains("pp");
  if (findCollSystem) collSystem = "pp";


  // get number of events after centrality celection
  Int_t nEventsCent = (Int_t) ((TH3F*)fIn->Get("nTracksCent"))->GetEntries();



    TLatex *textBField         = new TLatex(0.7, 0.75 , Form("#it{B} = %g T",BField));
    TLatex *textPIDScenario    = new TLatex(0.7, 0.7 , Form("%s",strPIDscenario[ith_PIDscenario-1].Data()));
    textBField->SetNDC(kTRUE);
    textPIDScenario->SetNDC(kTRUE);
    textBField->SetTextSize(0.04);
    textPIDScenario->SetTextSize(0.02);


    TLatex *textBField_conta         = new TLatex(0.17, 0.68 , Form("#it{B} = %g T",BField));
    TLatex *textPIDScenario_conta    = new TLatex(0.7, 0.9 , Form("%s",strPIDscenario[ith_PIDscenario-1].Data()));
    textBField_conta->SetNDC(kTRUE);
    textPIDScenario_conta->SetNDC(kTRUE);
    textBField_conta->SetTextSize(0.03);
    textPIDScenario_conta->SetTextSize(25);



  // read generated ULS histos
  // fIn->cd("generated/ULS");
  // Track histograms

  // TH3F* hRec_TrackPtEtaPhi_primary = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_primary_rec");
  // TH3F* hRecEle_TrackPtEtaPhi_primary = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_primary_Ele_rec");
  // TH3F* hRecPos_TrackPtEtaPhi_primary = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_primary_Pos_rec");
  // TH3F* hRec_TrackPtEtaPhi_hf = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_hf_rec");
  // TH3F* hRecEle_TrackPtEtaPhi_hf = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_hf_Ele_rec");
  // TH3F* hRecPos_TrackPtEtaPhi_hf = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_hf_Pos_rec");
  // TH3F* hRec_TrackPtEtaPhi_cc = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_charm_rec");
  // TH3F* hRecEle_TrackPtEtaPhi_cc = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_charm_Ele_rec");
  // TH3F* hRecPos_TrackPtEtaPhi_cc = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_charm_Pos_rec");
  // TH3F* hRec_TrackPtEtaPhi_bb = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_beauty_rec");
  // TH3F* hRecEle_TrackPtEtaPhi_bb = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_beauty_Ele_rec");
  // TH3F* hRecPos_TrackPtEtaPhi_bb = (TH3F*) fIn->cd(pathRecPIDScenario)->Get("hPt_Eta_Phi_beauty_Pos_rec");

  // TH3F* hGen_TrackPtEtaPhi_primary = (TH3F*) fIn->Get("hPt_Eta_Phi_primary_gen");
  // TH3F* hGenEle_TrackPtEtaPhi_primary = (TH3F*) fIn->Get("hPt_Eta_Phi_primary_Ele_gen");
  // TH3F* hGenPos_TrackPtEtaPhi_primary = (TH3F*) fIn->Get("hPt_Eta_Phi_primary_Pos_gen");
  // TH3F* hGen_TrackPtEtaPhi_hf = (TH3F*) fIn->Get("hPt_Eta_Phi_hf_gen");
  // TH3F* hGenEle_TrackPtEtaPhi_hf = (TH3F*) fIn->Get("hPt_Eta_Phi_hf_Ele_gen");
  // TH3F* hGenPos_TrackPtEtaPhi_hf = (TH3F*) fIn->Get("hPt_Eta_Phi_hf_Pos_gen");
  // TH3F* hGen_TrackPtEtaPhi_cc = (TH3F*) fIn->Get("hPt_Eta_Phi_charm_gen");
  // TH3F* hGenEle_TrackPtEtaPhi_cc = (TH3F*) fIn->Get("hPt_Eta_Phi_charm_Ele_gen");
  // TH3F* hGenPos_TrackPtEtaPhi_cc = (TH3F*) fIn->Get("hPt_Eta_Phi_charm_Pos_gen");
  // TH3F* hGen_TrackPtEtaPhi_bb = (TH3F*) fIn->Get("hPt_Eta_Phi_beauty_gen");
  // TH3F* hGenEle_TrackPtEtaPhi_bb = (TH3F*) fIn->Get("hPt_Eta_Phi_beauty_Ele_gen");
  // TH3F* hGenPos_TrackPtEtaPhi_bb = (TH3F*) fIn->Get("hPt_Eta_Phi_beauty_Pos_gen");

  TH3F* hRec_Track_Pt_Eta_Phi_BeforeSmearing = (TH3F*) fIn->Get("hBeforeSmearing_Pt_Eta_Phi_rec");
  TH3F* hRec_Track_Pt_Eta_Phi_AfterSmearing = (TH3F*) fIn->Get("hAfterSmearing_Pt_Eta_Phi_rec");

  TH3F* hRec_Track_Pt_Eta_Phi_AfterKineCuts = (TH3F*) fIn->Get(Form("%s/hAfterKineCuts_Pt_Eta_Phi_rec_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_Ele_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Ele_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_Pos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Pos_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_AllTrack_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hAllTracks_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_NegTrack_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hNegTrack_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_PosTrack_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hPosTrack_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_Muon_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Muon_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_Pion_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Pion_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_Kaon_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Kaon_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_Proton_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Proton_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH3F* hGen_Track_ElePos_Pt_Eta_Phi_beforeKineCuts = (TH3F*) fIn->Get("hTrack_ElePos_Gen_Pt_Eta_Phi_beforeKineCuts");
  TH3F* hGen_Track_All_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_All_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_ElePos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_Muon_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Muon_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_Pion_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Pion_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_Kaon_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Kaon_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_Proton_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Proton_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_Ele_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Ele_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_Pos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Pos_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));

  // TH3F* hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts = (TH3F*) fIn->Get("hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts");
  TH3F* hGenSmeared_Track_ElePos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_GenSmeared_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGenSmeared_Track_Ele_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Ele_GenSmeared_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGenSmeared_Track_Pos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Pos_GenSmeared_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));

  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut1 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_1");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut2 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_2");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut3 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_3");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut4 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_4");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut5 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_5");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut6 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_6");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut7 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_7");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut8 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_8");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut9 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_9");
  // TH3F* hRec_Track_Pt_Eta_Phi_EtaCut10 = (TH3F*) fIn->Get("hPt_Eta_Phi_rec_Eta_Cut_10");

  std::vector<TH2F*> vecTOF_PIDplots;
  std::vector<TH2F*> vecRICH_PIDplots;
  TH2F* hRec_TOF_BetaP_beforeSmearing;
  TH2F* hRec_TOF_BetaP_afterSmearing;
  TH2F* hRec_RICH_CherenkovAngleP_beforeSmearing;
  TH2F* hRec_RICH_CherenkovAngleP_afterSmearing;

  TH2 *hNsigmaP_TOF_trueElec[5];
  TH2 *hNsigmaP_TOF_trueMuon[5];
  TH2 *hNsigmaP_TOF_truePion[5];
  TH2 *hNsigmaP_TOF_trueKaon[5];
  TH2 *hNsigmaP_TOF_trueProton[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueElec[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueMuon[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_truePion[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueKaon[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueProton[5];

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

  std::vector<TH2F*> vecHistosNSigma_afterCuts;
  std::vector<TH2F*> vecHistosNSigma;


if (bPlotPIDhistograms) {
  hRec_TOF_BetaP_beforeSmearing = (TH2F*) fIn->Get("hBetaP_beforeSmearing");
  hRec_TOF_BetaP_afterSmearing = (TH2F*) fIn->Get("hBetaP_afterSmearing");
  vecTOF_PIDplots.push_back(hRec_TOF_BetaP_beforeSmearing);
  vecTOF_PIDplots.push_back(hRec_TOF_BetaP_afterSmearing);

  hRec_RICH_CherenkovAngleP_beforeSmearing = (TH2F*) fIn->Get("hCherenkovAngleP_beforeSmearing");
  hRec_RICH_CherenkovAngleP_afterSmearing = (TH2F*) fIn->Get("hCherenkovAngleP_afterSmearing");
  vecRICH_PIDplots.push_back(hRec_RICH_CherenkovAngleP_beforeSmearing);
  vecRICH_PIDplots.push_back(hRec_RICH_CherenkovAngleP_afterSmearing);

  TH2F *hRec_TOF_NSigma[5];
  TH2F *hRec_RICH_NSigma[5];
  const char *pname[5] = {"el", "mu", "pi", "ka", "pr"};
  const char *plabel[5] = {"e", "#mu", "#pi", "K", "p"};
  for (int i = 0; i < 5; ++i) {
    hRec_TOF_NSigma[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_TOF", pname[i]));
    hRec_TOF_NSigma[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hRec_TOF_NSigma[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    vecTOF_PIDplots.push_back(hRec_TOF_NSigma[i]);
    hRec_RICH_NSigma[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_RICH", pname[i]));
    hRec_RICH_NSigma[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hRec_RICH_NSigma[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    vecRICH_PIDplots.push_back(hRec_RICH_NSigma[i]);
    hRec_TOF_NSigma[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_sce%i", pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hRec_TOF_NSigma[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hRec_TOF_NSigma[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    vecTOF_PIDplots.push_back(hRec_TOF_NSigma[i]);
    hRec_RICH_NSigma[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_sce%i", pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hRec_RICH_NSigma[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    hRec_RICH_NSigma[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    vecRICH_PIDplots.push_back(hRec_RICH_NSigma[i]);
  }

  for (int i = 0; i < 5; ++i) {
    // hNsigmaP_TOF[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_TOF", pname[i]));
    hNsigmaP_TOF_trueElec[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_TOF_trueElec", pname[i]));
    hNsigmaP_TOF_trueMuon[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_TOF_trueMuon", pname[i]));
    hNsigmaP_TOF_truePion[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_TOF_truePion", pname[i]));
    hNsigmaP_TOF_trueKaon[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_TOF_trueKaon", pname[i]));
    hNsigmaP_TOF_trueProton[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_TOF_trueProton", pname[i]));
    // hNsigmaP_afterPIDcuts_TOF[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_afterPIDcuts_TOF", pname[i]));
    hNsigmaP_afterPIDcuts_TOF_trueElec[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_trueElec_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
    hNsigmaP_afterPIDcuts_TOF_trueMuon[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_trueMuon_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
    hNsigmaP_afterPIDcuts_TOF_truePion[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_truePion_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
    hNsigmaP_afterPIDcuts_TOF_trueKaon[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_trueKaon_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
    hNsigmaP_afterPIDcuts_TOF_trueProton[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_trueProton_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
    hNsigmaP_TOF_trueElec[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    hNsigmaP_TOF_trueMuon[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    hNsigmaP_TOF_truePion[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    hNsigmaP_afterPIDcuts_TOF_trueElec[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    hNsigmaP_afterPIDcuts_TOF_truePion[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{TOF}", plabel[i]));
    hNsigmaP_TOF_trueElec[i]->SetMarkerColor(kBlue);
    hNsigmaP_TOF_trueMuon[i]->SetMarkerColor(kBlack);
    hNsigmaP_TOF_truePion[i]->SetMarkerColor(kRed);
    hNsigmaP_afterPIDcuts_TOF_trueElec[i]->SetMarkerColor(kBlue);
    hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->SetMarkerColor(kBlack);
    hNsigmaP_afterPIDcuts_TOF_truePion[i]->SetMarkerColor(kRed);
    // hNsigmaP_afterPIDcuts_TOF_trueElec[i]->SetMarkerStyle(20);
    hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->SetMarkerStyle(20);
    hNsigmaP_afterPIDcuts_TOF_truePion[i]->SetMarkerStyle(20);
    // hNsigmaP_afterPIDcuts_TOF_trueElec[i]->SetMarkerSize(0.5);
    hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->SetMarkerSize(0.4);
    hNsigmaP_afterPIDcuts_TOF_truePion[i]->SetMarkerSize(0.4);
    hNsigmaP_TOF_trueElec[i]->SetLineColor(kBlue);
    hNsigmaP_TOF_trueMuon[i]->SetLineColor(kBlack);
    hNsigmaP_TOF_truePion[i]->SetLineColor(kRed);
    hNsigmaP_afterPIDcuts_TOF_trueElec[i]->SetLineColor(kBlue);
    hNsigmaP_afterPIDcuts_TOF_trueMuon[i]->SetLineColor(kBlack);
    hNsigmaP_afterPIDcuts_TOF_truePion[i]->SetLineColor(kRed);
  }


  for (int i = 0; i < 5; ++i) {
    hNsigmaP_RICH[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_RICH", pname[i]));
    hNsigmaP_RICH_trueElec[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_RICH_trueElec", pname[i]));
    hNsigmaP_RICH_trueMuon[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_RICH_trueMuon", pname[i]));
    hNsigmaP_RICH_truePion[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_RICH_truePion", pname[i]));
    hNsigmaP_RICH_trueKaon[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_RICH_trueKaon", pname[i]));
    hNsigmaP_RICH_trueProton[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_RICH_trueProton", pname[i]));
    hNsigmaP_afterPIDcuts_RICH[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_trueElec[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_trueElec_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_trueMuon[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_trueMuon_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_truePion[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_truePion_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_trueKaon[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_trueKaon_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_trueProton[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_trueProton_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_RICH_trueElec[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    hNsigmaP_RICH_trueMuon[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    hNsigmaP_RICH_truePion[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    hNsigmaP_afterPIDcuts_RICH_trueElec[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    hNsigmaP_afterPIDcuts_RICH_truePion[i]->GetYaxis()->SetTitle(Form("N#sigma_{%s}^{RICH}", plabel[i]));
    hNsigmaP_RICH_trueElec[i]->SetMarkerColor(kBlue);
    hNsigmaP_RICH_trueMuon[i]->SetMarkerColor(kBlack);
    hNsigmaP_RICH_truePion[i]->SetMarkerColor(kRed);
    hNsigmaP_afterPIDcuts_RICH_trueElec[i]->SetMarkerColor(kBlue);
    hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->SetMarkerColor(kBlack);
    hNsigmaP_afterPIDcuts_RICH_truePion[i]->SetMarkerColor(kRed);
    // hNsigmaP_afterPIDcuts_RICH_trueElec[i]->SetMarkerStyle(20);
    hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->SetMarkerStyle(20);
    hNsigmaP_afterPIDcuts_RICH_truePion[i]->SetMarkerStyle(20);
    // hNsigmaP_afterPIDcuts_RICH_trueElec[i]->SetMarkerSize(0.5);
    hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->SetMarkerSize(0.4);
    hNsigmaP_afterPIDcuts_RICH_truePion[i]->SetMarkerSize(0.4);
    hNsigmaP_RICH_trueElec[i]->SetLineColor(kBlue);
    hNsigmaP_RICH_trueMuon[i]->SetLineColor(kBlack);
    hNsigmaP_RICH_truePion[i]->SetLineColor(kRed);
    hNsigmaP_afterPIDcuts_RICH_trueElec[i]->SetLineColor(kBlue);
    hNsigmaP_afterPIDcuts_RICH_trueMuon[i]->SetLineColor(kBlack);
    hNsigmaP_afterPIDcuts_RICH_truePion[i]->SetLineColor(kRed);
  }

  // pushback nSigma for one specific particle for every true kind of particle than continue to next nsigma spectra
  // used for nSigma projections in specific momentum intervalls
  vecHistosNSigma_afterCuts.push_back((TH2F*)hNsigmaP_afterPIDcuts_RICH_trueElec[0]);
  vecHistosNSigma_afterCuts.push_back((TH2F*)hNsigmaP_afterPIDcuts_RICH_truePion[0]);
  vecHistosNSigma_afterCuts.push_back((TH2F*)hNsigmaP_afterPIDcuts_RICH_trueElec[2]);
  vecHistosNSigma_afterCuts.push_back((TH2F*)hNsigmaP_afterPIDcuts_RICH_truePion[2]);

  vecHistosNSigma.push_back((TH2F*)hNsigmaP_RICH_trueElec[0]);
  vecHistosNSigma.push_back((TH2F*)hNsigmaP_RICH_truePion[0]);
  vecHistosNSigma.push_back((TH2F*)hNsigmaP_RICH_trueElec[2]);
  vecHistosNSigma.push_back((TH2F*)hNsigmaP_RICH_truePion[2]);
}

  //Pair histograms M,pT,DCA
  TH3F* hMPtDCA_ULS_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_primary_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_primary_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_HF_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_heavy_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_CC_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_charm_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_BB_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_beauty_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_rec_MCpidEle = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH3F* hMPtDCA_LS_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_primary_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_primary_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_HF_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_heavy_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_CC_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_charm_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_BB_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_beauty_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH3F* hMPtDCA_ULS_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_gen_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_primary_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_gen_primary_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_HF_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_gen_heavy_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_CC_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_gen_charm_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_ULS_BB_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_gen_beauty_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));

  TH3F* hMPtDCA_LS_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_gen_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_primary_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_gen_primary_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_HF_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_gen_heavy_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_CC_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_gen_charm_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hMPtDCA_LS_BB_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_gen_beauty_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));


  // create projections and profiles
  // TH1F* ptRecTrackPrim  = (TH1F*) hRec_TrackPtEtaPhi_primary->ProjectionX("rec ptTrackPrim");
  // TH1F* ptRecEleTrackPrim  = (TH1F*) hRecEle_TrackPtEtaPhi_primary->ProjectionX("rec e- ptTrackPrim");
  // TH1F* ptRecPosTrackPrim  = (TH1F*) hRecPos_TrackPtEtaPhi_primary->ProjectionX("rec e+ ptTrackPrim");
  // TH1F* etaRecTrackPrim = (TH1F*) hRec_TrackPtEtaPhi_primary->ProjectionY("rec etaTrackPrim");
  // TH1F* etaRecEleTrackPrim = (TH1F*) hRecEle_TrackPtEtaPhi_primary->ProjectionY("rec e- etaTrackPrim");
  // TH1F* etaRecPosTrackPrim = (TH1F*) hRecPos_TrackPtEtaPhi_primary->ProjectionY("rec e+ etaTrackPrim");
  // TH1F* phiRecTrackPrim = (TH1F*) hRec_TrackPtEtaPhi_primary->ProjectionZ("rec phiTrackPrim");
  // TH1F* phiRecEleTrackPrim = (TH1F*) hRecEle_TrackPtEtaPhi_primary->ProjectionZ("rec e- phiTrackPrim");
  // TH1F* phiRecPosTrackPrim = (TH1F*) hRecPos_TrackPtEtaPhi_primary->ProjectionZ("rec e+ phiTrackPrim");
  // TH1F* ptRecTrackCC  = (TH1F*) hRec_TrackPtEtaPhi_cc->ProjectionX("rec ptTrackCC");
  // TH1F* ptRecEleTrackCC  = (TH1F*) hRecEle_TrackPtEtaPhi_cc->ProjectionX("rec e- ptTrackCC");
  // TH1F* ptRecPosTrackCC  = (TH1F*) hRecPos_TrackPtEtaPhi_cc->ProjectionX("rec e+ ptTrackCC");
  // TH1F* etaRecTrackCC = (TH1F*) hRec_TrackPtEtaPhi_cc->ProjectionY("rec etaTrackCC");
  // TH1F* etaRecEleTrackCC = (TH1F*) hRecEle_TrackPtEtaPhi_cc->ProjectionY("rec e- etaTrackCC");
  // TH1F* etaRecPosTrackCC = (TH1F*) hRecPos_TrackPtEtaPhi_cc->ProjectionY("rec e+ etaTrackCC");
  // TH1F* phiRecTrackCC = (TH1F*) hRec_TrackPtEtaPhi_cc->ProjectionZ("rec phiTrackCC");
  // TH1F* phiRecEleTrackCC = (TH1F*) hRecEle_TrackPtEtaPhi_cc->ProjectionZ("rec e- phiTrackCC");
  // TH1F* phiRecPosTrackCC = (TH1F*) hRecPos_TrackPtEtaPhi_cc->ProjectionZ("rec e+ phiTrackCC");
  // TH1F* ptRecTrackBB  = (TH1F*) hRec_TrackPtEtaPhi_bb->ProjectionX("rec ptTrackBB");
  // TH1F* ptRecEleTrackBB  = (TH1F*) hRecEle_TrackPtEtaPhi_bb->ProjectionX("rec e- ptTrackBB");
  // TH1F* ptRecPosTrackBB  = (TH1F*) hRecPos_TrackPtEtaPhi_bb->ProjectionX("rec e+ ptTrackBB");
  // TH1F* etaRecTrackBB = (TH1F*) hRec_TrackPtEtaPhi_bb->ProjectionY("rec etaTrackBB");
  // TH1F* etaRecEleTrackBB = (TH1F*) hRecEle_TrackPtEtaPhi_bb->ProjectionY("rec e- etaTrackBB");
  // TH1F* etaRecPosTrackBB = (TH1F*) hRecPos_TrackPtEtaPhi_bb->ProjectionY("rec e+ etaTrackBB");
  // TH1F* phiRecTrackBB = (TH1F*) hRec_TrackPtEtaPhi_bb->ProjectionZ("rec phiTrackBB");
  // TH1F* phiRecEleTrackBB = (TH1F*) hRecEle_TrackPtEtaPhi_bb->ProjectionZ("rec e- phiTrackBB");
  // TH1F* phiRecPosTrackBB = (TH1F*) hRecPos_TrackPtEtaPhi_bb->ProjectionZ("rec e+ phiTrackBB");

  // TH1F* ptGenTrackPrim  = (TH1F*) hGen_TrackPtEtaPhi_primary->ProjectionX("gen ptTrackPrim");
  // TH1F* ptGenEleTrackPrim  = (TH1F*) hGenEle_TrackPtEtaPhi_primary->ProjectionX("gen e- ptTrackPrim");
  // TH1F* ptGenPosTrackPrim  = (TH1F*) hGenPos_TrackPtEtaPhi_primary->ProjectionX("gen e+ ptTrackPrim");
  // TH1F* etaGenTrackPrim = (TH1F*) hGen_TrackPtEtaPhi_primary->ProjectionY("gen etaTrackPrim");
  // TH1F* etaGenEleTrackPrim = (TH1F*) hGenEle_TrackPtEtaPhi_primary->ProjectionY("gen e- etaTrackPrim");
  // TH1F* etaGenPosTrackPrim = (TH1F*) hGenPos_TrackPtEtaPhi_primary->ProjectionY("gen e+ etaTrackPrim");
  // TH1F* phiGenTrackPrim = (TH1F*) hGen_TrackPtEtaPhi_primary->ProjectionZ("gen phiTrackPrim");
  // TH1F* phiGenEleTrackPrim = (TH1F*) hGenEle_TrackPtEtaPhi_primary->ProjectionZ("gen e- phiTrackPrim");
  // TH1F* phiGenPosTrackPrim = (TH1F*) hGenPos_TrackPtEtaPhi_primary->ProjectionZ("gen e+ phiTrackPrim");
  // TH1F* ptGenTrackCC  = (TH1F*) hGen_TrackPtEtaPhi_cc->ProjectionX("gen ptTrackCC");
  // TH1F* ptGenEleTrackCC  = (TH1F*) hGenEle_TrackPtEtaPhi_cc->ProjectionX("gen e- ptTrackCC");
  // TH1F* ptGenPosTrackCC  = (TH1F*) hGenPos_TrackPtEtaPhi_cc->ProjectionX("gen e+ ptTrackCC");
  // TH1F* etaGenTrackCC = (TH1F*) hGen_TrackPtEtaPhi_cc->ProjectionY("gen etaTrackCC");
  // TH1F* etaGenEleTrackCC = (TH1F*) hGenEle_TrackPtEtaPhi_cc->ProjectionY("gen e- etaTrackCC");
  // TH1F* etaGenPosTrackCC = (TH1F*) hGenPos_TrackPtEtaPhi_cc->ProjectionY("gen e+ etaTrackCC");
  // TH1F* phiGenTrackCC = (TH1F*) hGen_TrackPtEtaPhi_cc->ProjectionZ("gen phiTrackCC");
  // TH1F* phiGenEleTrackCC = (TH1F*) hGenEle_TrackPtEtaPhi_cc->ProjectionZ("gen e- phiTrackCC");
  // TH1F* phiGenPosTrackCC = (TH1F*) hGenPos_TrackPtEtaPhi_cc->ProjectionZ("gen e+ phiTrackCC");
  // TH1F* ptGenTrackBB  = (TH1F*) hGen_TrackPtEtaPhi_bb->ProjectionX("gen ptTrackBB");
  // TH1F* ptGenEleTrackBB  = (TH1F*) hGenEle_TrackPtEtaPhi_bb->ProjectionX("gen e- ptTrackBB");
  // TH1F* ptGenPosTrackBB  = (TH1F*) hGenPos_TrackPtEtaPhi_bb->ProjectionX("gen e+ ptTrackBB");
  // TH1F* etaGenTrackBB = (TH1F*) hGen_TrackPtEtaPhi_bb->ProjectionY("gen etaTrackBB");
  // TH1F* etaGenEleTrackBB = (TH1F*) hGenEle_TrackPtEtaPhi_bb->ProjectionY("gen e- etaTrackBB");
  // TH1F* etaGenPosTrackBB = (TH1F*) hGenPos_TrackPtEtaPhi_bb->ProjectionY("gen e+ etaTrackBB");
  // TH1F* phiGenTrackBB = (TH1F*) hGen_TrackPtEtaPhi_bb->ProjectionZ("gen phiTrackBB");
  // TH1F* phiGenEleTrackBB = (TH1F*) hGenEle_TrackPtEtaPhi_bb->ProjectionZ("gen e- phiTrackBB");
  // TH1F* phiGenPosTrackBB = (TH1F*) hGenPos_TrackPtEtaPhi_bb->ProjectionZ("gen e+ phiTrackBB");

  TH1F* ptRecTrackBeforeSmearing  = (TH1F*) hRec_Track_Pt_Eta_Phi_BeforeSmearing->ProjectionX("Rec beforeSmearing ptTrack");
  TH1F* etaRecTrackBeforeSmearing = (TH1F*) hRec_Track_Pt_Eta_Phi_BeforeSmearing->ProjectionY("Rec beforeSmearing etaTrack");
  TH1F* phiRecTrackBeforeSmearing = (TH1F*) hRec_Track_Pt_Eta_Phi_BeforeSmearing->ProjectionZ("Rec beforeSmearing phiTrack");

  TH1F* ptRecTrackAfterSmearing  = (TH1F*) hRec_Track_Pt_Eta_Phi_AfterSmearing->ProjectionX("Rec afterSmearing ptTrack");
  TH1F* etaRecTrackAfterSmearing = (TH1F*) hRec_Track_Pt_Eta_Phi_AfterSmearing->ProjectionY("Rec afterSmearing etaTrack");
  TH1F* phiRecTrackAfterSmearing = (TH1F*) hRec_Track_Pt_Eta_Phi_AfterSmearing->ProjectionZ("Rec afterSmearing phiTrack");

  TH1F* ptRecTrackAfterKineCuts  = (TH1F*) hRec_Track_Pt_Eta_Phi_AfterKineCuts->ProjectionX("Rec afterKineCuts ptTrack");
  TH1F* etaRecTrackAfterKineCuts = (TH1F*) hRec_Track_Pt_Eta_Phi_AfterKineCuts->ProjectionY("Rec afterKineCuts etaTrack");
  TH1F* phiRecTrackAfterKineCuts = (TH1F*) hRec_Track_Pt_Eta_Phi_AfterKineCuts->ProjectionZ("Rec afterKineCuts phiTrack");

  TH1F* ptRecTrackElePos  = (TH1F*) hRec_Track_ElePos_Pt_Eta_Phi->ProjectionX("Rec ElePos afterPIDCuts ptTrack");
  TH1F* etaRecTrackElePos = (TH1F*) hRec_Track_ElePos_Pt_Eta_Phi->ProjectionY("Rec ElePos afterPIDCuts etaTrack");
  TH1F* phiRecTrackElePos = (TH1F*) hRec_Track_ElePos_Pt_Eta_Phi->ProjectionZ("Rec ElePos afterPIDCuts phiTrack");

  TH2F* ptEtaRecTrackElePos = (TH2F*) hRec_Track_ElePos_Pt_Eta_Phi->Project3D("yx o");

  TH1F* ptRecTrackEle  = (TH1F*) hRec_Track_Ele_Pt_Eta_Phi->ProjectionX("Rec electrons ptTrack");
  TH1F* etaRecTrackEle = (TH1F*) hRec_Track_Ele_Pt_Eta_Phi->ProjectionY("Rec electrons etaTrack");
  TH1F* phiRecTrackEle = (TH1F*) hRec_Track_Ele_Pt_Eta_Phi->ProjectionZ("Rec electrons phiTrack");

  TH1F* ptRecTrackPos  = (TH1F*) hRec_Track_Pos_Pt_Eta_Phi->ProjectionX("Rec positrons ptTrack");
  TH1F* etaRecTrackPos = (TH1F*) hRec_Track_Pos_Pt_Eta_Phi->ProjectionY("Rec positrons etaTrack");
  TH1F* phiRecTrackPos = (TH1F*) hRec_Track_Pos_Pt_Eta_Phi->ProjectionZ("Rec positrons phiTrack");

  TH1F* ptAllRecTrack  = (TH1F*) hRec_AllTrack_Pt_Eta_Phi->ProjectionX("Rec ptTrack");
  TH1F* etaAllRecTrack = (TH1F*) hRec_AllTrack_Pt_Eta_Phi->ProjectionY("Rec etaTrack");
  TH1F* phiAllRecTrack = (TH1F*) hRec_AllTrack_Pt_Eta_Phi->ProjectionZ("Rec phiTrack");

  TH1F* ptRecNegTrack  = (TH1F*) hRec_NegTrack_Pt_Eta_Phi->ProjectionX("Rec negative ptTrack");
  TH1F* etaRecNegTrack = (TH1F*) hRec_NegTrack_Pt_Eta_Phi->ProjectionY("Rec negative etaTrack");
  TH1F* phiRecNegTrack = (TH1F*) hRec_NegTrack_Pt_Eta_Phi->ProjectionZ("Rec negative phiTrack");

  TH1F* ptRecPosTrack  = (TH1F*) hRec_PosTrack_Pt_Eta_Phi->ProjectionX("Rec positive ptTrack");
  TH1F* etaRecPosTrack = (TH1F*) hRec_PosTrack_Pt_Eta_Phi->ProjectionY("Rec positive etaTrack");
  TH1F* phiRecPosTrack = (TH1F*) hRec_PosTrack_Pt_Eta_Phi->ProjectionZ("Rec positive phiTrack");

  TH1F* ptRecMuonTrack  = (TH1F*) hRec_Track_Muon_Pt_Eta_Phi->ProjectionX("Rec Muon ptTrack");
  TH1F* etaRecMuonTrack = (TH1F*) hRec_Track_Muon_Pt_Eta_Phi->ProjectionY("Rec Muon etaTrack");
  TH1F* phiRecMuonTrack = (TH1F*) hRec_Track_Muon_Pt_Eta_Phi->ProjectionZ("Rec Muon phiTrack");

  TH1F* ptRecPionTrack  = (TH1F*) hRec_Track_Pion_Pt_Eta_Phi->ProjectionX("Rec Pion ptTrack");
  TH1F* etaRecPionTrack = (TH1F*) hRec_Track_Pion_Pt_Eta_Phi->ProjectionY("Rec Pion etaTrack");
  TH1F* phiRecPionTrack = (TH1F*) hRec_Track_Pion_Pt_Eta_Phi->ProjectionZ("Rec Pion phiTrack");

  TH1F* ptRecKaonTrack  = (TH1F*) hRec_Track_Kaon_Pt_Eta_Phi->ProjectionX("Rec Kaon ptTrack");
  TH1F* etaRecKaonTrack = (TH1F*) hRec_Track_Kaon_Pt_Eta_Phi->ProjectionY("Rec Kaon etaTrack");
  TH1F* phiRecKaonTrack = (TH1F*) hRec_Track_Kaon_Pt_Eta_Phi->ProjectionZ("Rec Kaon phiTrack");

  TH1F* ptRecProtonTrack  = (TH1F*) hRec_Track_Proton_Pt_Eta_Phi->ProjectionX("Rec Proton ptTrack");
  TH1F* etaRecProtonTrack = (TH1F*) hRec_Track_Proton_Pt_Eta_Phi->ProjectionY("Rec Proton etaTrack");
  TH1F* phiRecProtonTrack = (TH1F*) hRec_Track_Proton_Pt_Eta_Phi->ProjectionZ("Rec Proton phiTrack");

  TH1F* ptGenTrackElePos_beforeKineCuts  = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionX("Gen ElePos ptTrack beforeKineCuts");
  TH1F* etaGenTrackElePos_beforeKineCuts = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionY("Gen ElePos etaTrack beforeKineCuts");
  TH1F* phiGenTrackElePos_beforeKineCuts = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionZ("Gen ElePos phiTrack beforeKineCuts");

  TH1F* ptGenTrackElePos  = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi->ProjectionX("Gen ElePos ptTrack");
  TH1F* etaGenTrackElePos = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi->ProjectionY("Gen ElePos etaTrack");
  TH1F* phiGenTrackElePos = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi->ProjectionZ("Gen ElePos phiTrack");

  TH2F* ptEtaGenTrackElePos = (TH2F*) hGen_Track_ElePos_Pt_Eta_Phi->Project3D("yx o");

  TH1F* ptGenTrackAll   = (TH1F*) hGen_Track_All_Pt_Eta_Phi->ProjectionX("Gen All ptTrack");
  TH1F* ptGenTrackMuon  = (TH1F*) hGen_Track_Muon_Pt_Eta_Phi->ProjectionX("Gen Muon ptTrack");
  TH1F* ptGenTrackPion  = (TH1F*) hGen_Track_Pion_Pt_Eta_Phi->ProjectionX("Gen Pion ptTrack");
  TH1F* ptGenTrackKaon  = (TH1F*) hGen_Track_Kaon_Pt_Eta_Phi->ProjectionX("Gen Kaon ptTrack");
  TH1F* ptGenTrackProton  = (TH1F*) hGen_Track_Proton_Pt_Eta_Phi->ProjectionX("Gen Proton ptTrack");


  TH1F* ptGenTrackEle  = (TH1F*) hGen_Track_Ele_Pt_Eta_Phi->ProjectionX("Gen electrons ptTrack");
  TH1F* etaGenTrackEle = (TH1F*) hGen_Track_Ele_Pt_Eta_Phi->ProjectionY("Gen electrons etaTrack");
  TH1F* phiGenTrackEle = (TH1F*) hGen_Track_Ele_Pt_Eta_Phi->ProjectionZ("Gen electrons phiTrack");

  TH1F* ptGenTrackPos  = (TH1F*) hGen_Track_Pos_Pt_Eta_Phi->ProjectionX("Gen positrons ptTrack");
  TH1F* etaGenTrackPos = (TH1F*) hGen_Track_Pos_Pt_Eta_Phi->ProjectionY("Gen positrons etaTrack");
  TH1F* phiGenTrackPos = (TH1F*) hGen_Track_Pos_Pt_Eta_Phi->ProjectionZ("Gen positrons phiTrack");

  // TH1F* ptGenSmearedTrackElePos_beforeKineCuts  = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionX("GenSmeared ElePos ptTrack beforeKineCuts");
  // TH1F* etaGenSmearedTrackElePos_beforeKineCuts = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionY("GenSmeared ElePos etaTrack beforeKineCuts");
  // TH1F* phiGenSmearedTrackElePos_beforeKineCuts = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionZ("GenSmeared ElePos phiTrack beforeKineCuts");

  TH1F* ptGenSmearedTrackElePos  = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi->ProjectionX("GenSmeared ElePos ptTrack");
  TH1F* etaGenSmearedTrackElePos = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi->ProjectionY("GenSmeared ElePos etaTrack");
  TH1F* phiGenSmearedTrackElePos = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi->ProjectionZ("GenSmeared ElePos phiTrack");

  TH1F* ptGenSmearedTrackEle  = (TH1F*) hGenSmeared_Track_Ele_Pt_Eta_Phi->ProjectionX("GenSmeared electrons ptTrack");
  TH1F* etaGenSmearedTrackEle = (TH1F*) hGenSmeared_Track_Ele_Pt_Eta_Phi->ProjectionY("GenSmeared electrons etaTrack");
  TH1F* phiGenSmearedTrackEle = (TH1F*) hGenSmeared_Track_Ele_Pt_Eta_Phi->ProjectionZ("GenSmeared electrons phiTrack");

  TH1F* ptGenSmearedTrackPos  = (TH1F*) hGenSmeared_Track_Pos_Pt_Eta_Phi->ProjectionX("GenSmeared positrons ptTrack");
  TH1F* etaGenSmearedTrackPos = (TH1F*) hGenSmeared_Track_Pos_Pt_Eta_Phi->ProjectionY("GenSmeared positrons etaTrack");
  TH1F* phiGenSmearedTrackPos = (TH1F*) hGenSmeared_Track_Pos_Pt_Eta_Phi->ProjectionZ("GenSmeared positrons phiTrack");

  etaGenSmearedTrackPos->GetXaxis()->SetRangeUser(-5.0,5.0);

  // TH1F* ptRecTrackEtaCut_1 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut1->ProjectionX("Rec eta1cut pt Track");
  // TH1F* ptRecTrackEtaCut_2 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut2->ProjectionX("Rec eta2cut pt Track");
  // TH1F* ptRecTrackEtaCut_3 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut3->ProjectionX("Rec eta3cut pt Track");
  // TH1F* ptRecTrackEtaCut_4 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut4->ProjectionX("Rec eta4cut pt Track");
  // TH1F* ptRecTrackEtaCut_5 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut5->ProjectionX("Rec eta5cut pt Track");
  // TH1F* ptRecTrackEtaCut_6 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut6->ProjectionX("Rec eta6cut pt Track");
  // TH1F* ptRecTrackEtaCut_7 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut7->ProjectionX("Rec eta7cut pt Track");
  // TH1F* ptRecTrackEtaCut_8 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut8->ProjectionX("Rec eta8cut pt Track");
  // TH1F* ptRecTrackEtaCut_9 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut9->ProjectionX("Rec eta9cut pt Track");
  // TH1F* ptRecTrackEtaCut_10 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut10->ProjectionX("Rec eta10cut pt Track");
  //
  // TH1F* etaRecTrackEtaCut_1 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut1->ProjectionY("Rec eta1cut eta Track");
  // TH1F* etaRecTrackEtaCut_2 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut2->ProjectionY("Rec eta2cut eta Track");
  // TH1F* etaRecTrackEtaCut_3 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut3->ProjectionY("Rec eta3cut eta Track");
  // TH1F* etaRecTrackEtaCut_4 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut4->ProjectionY("Rec eta4cut eta Track");
  // TH1F* etaRecTrackEtaCut_5 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut5->ProjectionY("Rec eta5cut eta Track");
  // TH1F* etaRecTrackEtaCut_6 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut6->ProjectionY("Rec eta6cut eta Track");
  // TH1F* etaRecTrackEtaCut_7 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut7->ProjectionY("Rec eta7cut eta Track");
  // TH1F* etaRecTrackEtaCut_8 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut8->ProjectionY("Rec eta8cut eta Track");
  // TH1F* etaRecTrackEtaCut_9 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut9->ProjectionY("Rec eta9cut eta Track");
  // TH1F* etaRecTrackEtaCut_10 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut10->ProjectionY("Rec eta10cut eta Track");
  //
  // TH1F* phiRecTrackEtaCut_1 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut1->ProjectionZ("Rec eta1cut phi Track");
  // TH1F* phiRecTrackEtaCut_2 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut2->ProjectionZ("Rec eta2cut phi Track");
  // TH1F* phiRecTrackEtaCut_3 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut3->ProjectionZ("Rec eta3cut phi Track");
  // TH1F* phiRecTrackEtaCut_4 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut4->ProjectionZ("Rec eta4cut phi Track");
  // TH1F* phiRecTrackEtaCut_5 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut5->ProjectionZ("Rec eta5cut phi Track");
  // TH1F* phiRecTrackEtaCut_6 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut6->ProjectionZ("Rec eta6cut phi Track");
  // TH1F* phiRecTrackEtaCut_7 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut7->ProjectionZ("Rec eta7cut phi Track");
  // TH1F* phiRecTrackEtaCut_8 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut8->ProjectionZ("Rec eta8cut phi Track");
  // TH1F* phiRecTrackEtaCut_9 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut9->ProjectionZ("Rec eta9cut phi Track");
  // TH1F* phiRecTrackEtaCut_10 = (TH1F*) hRec_Track_Pt_Eta_Phi_EtaCut10->ProjectionZ("Rec eta10cut phi Track");


    TH1F* hParticlesFIT = (TH1F*) fIn->Get("nParticlesFIT");
    Int_t nentries=0;
    Int_t nbins = hParticlesFIT->GetXaxis()->GetNbins();
    Int_t centralNinetyPercent = hParticlesFIT->GetEntries() * 0.9;
    cout  << nbins << endl;
    for (size_t ibin = 1; ibin < nbins ; ibin++) {
      if(nentries >= centralNinetyPercent) continue; // look for the lowest 90% to find the lower egde of the 10% cental events
      nentries = nentries + hParticlesFIT->GetBinContent(ibin);
      if(nentries >= centralNinetyPercent) cout << "nentries has reached 10% at " << hParticlesFIT->GetXaxis()->GetBinCenter(ibin) << endl;
    }

    TH1F* hdNdeta_midrap = (TH1F*) fIn->Get("hdNdeta_midrap_gen");
    TH1F* hParticlesFITCent = (TH1F*) fIn->Get("nParticlesFITCent");
    TH1F* hParticlesMidRapidity = (TH1F*) fIn->Get("nParticlesMidRapidity");
    TH1F* hParticlesMidRapidityCent = (TH1F*) fIn->Get("nParticlesMidRapidityCent");

    hParticlesFITCent->Sumw2();
    hParticlesMidRapidity->Sumw2();
    hParticlesMidRapidityCent->Sumw2();
    double meandNdeta = hdNdeta_midrap->GetMean(1);
    // cout << " mean of dN/deta = " << meandNdeta << endl;
    hdNdeta_midrap->Scale(1/hParticlesFITCent->GetEntries()/*, "width"*/);
    // double binwidth =hdNdeta_midrap->GetXaxis()->GetBinWidth(hdNdeta_midrap->FindBin(0));
    // cout << "cout binwidth = "<< binwidth << endl;
    // double meandNdeta1 = hdNdeta_midrap->GetMean(1);
    // cout << " mean of dN/deta = " << meandNdeta1 << endl;
    makeHistNice(hdNdeta_midrap,kBlue);
    makeHistNice(hParticlesMidRapidity,kBlue);
    makeHistNice(hParticlesMidRapidityCent,kRed);
    TLatex *textmeandNdeta    = new TLatex(0.7, 0.8    , Form(" mean = %e",meandNdeta));
    TLatex *textEntriesdNdeta    = new TLatex(0.7, 0.77    , Form(" nEntries = %e",hdNdeta_midrap->GetEntries()));
    textmeandNdeta->SetNDC(kTRUE);
    textEntriesdNdeta->SetNDC(kTRUE);
    textmeandNdeta->SetTextSize(0.02);
    textEntriesdNdeta->SetTextSize(0.02);

    auto cdNdeta = new TCanvas("cdNdeta","cdNdeta",800,800);
    cdNdeta->SetTopMargin(0.03);
    cdNdeta->SetRightMargin(0.03);
    cdNdeta->SetLeftMargin(0.13);
    hdNdeta_midrap->GetYaxis()->SetTitle("1/N_{ev}dN_{ch}/d#eta");
    hdNdeta_midrap->GetXaxis()->SetTitle("#eta");
    hdNdeta_midrap->SetMinimum(0.);
    hdNdeta_midrap->Draw("hist p e1");
    textmeandNdeta->Draw("same");
    textEntriesdNdeta->Draw("same");
    textBField_conta->Draw("same");
    // textPIDScenario_conta->Draw("same");
    cdNdeta->SaveAs("./plots/meandNdeta.png");




    double legPosPart[4] = {0.15,0.78,0.45,0.93};
    auto legNpartMidRap = new TLegend(legPosPart[0]+0.05,legPosPart[1],legPosPart[2]+0.05,legPosPart[3]);
    legNpartMidRap->SetBorderSize(0);
    legNpartMidRap->SetFillStyle(0);
    // legNpartMidRap->SetTextSize(0.04);
    legNpartMidRap->AddEntry(hParticlesMidRapidity,"before cent. selection","l");
    legNpartMidRap->AddEntry(hParticlesMidRapidityCent,"after cent. selection","l");

    auto cnParticlesMidRap = new TCanvas("cdNdeta","cdNdeta",800,800);
    cnParticlesMidRap->SetTopMargin(0.03);
    cnParticlesMidRap->SetRightMargin(0.03);
    cnParticlesMidRap->SetLeftMargin(0.13);
    cnParticlesMidRap->SetLogy();
    hParticlesMidRapidity->GetYaxis()->SetTitle("N particles |#eta|<0.5");
    hParticlesMidRapidity->GetXaxis()->SetTitle("#eta");
    hParticlesMidRapidity->GetXaxis()->SetRangeUser(0.,2500);
    // hParticlesMidRapidity->SetMinimum(0.);
    hParticlesMidRapidity->Draw("");
    hParticlesMidRapidityCent->Draw("same");
    legNpartMidRap->Draw("same");
    textBField_conta->Draw("same");
    // textPIDScenario_conta->Draw("same");
    cnParticlesMidRap->SaveAs("./plots/nParticlesMidRapidity.png");



  // Profiles to see M_ee and pT_ee spectras
    TH1F* proj_recULS_Mee = (TH1F*) hMPtDCA_ULS_rec->ProjectionX("proj_recULS_Mee");
    TH1F* proj_recULS_Ptee = (TH1F*) hMPtDCA_ULS_rec->ProjectionY("proj_recULS_Ptee");
    TH1F* proj_recULS_DCA = (TH1F*) hMPtDCA_ULS_rec->ProjectionZ("proj_recULS_DCA");
    // TH1F* proj_recULS_MeePrim = (TH1F*) hMPtDCA_ULS_primary_rec->ProjectionX("proj_recULS_MeePrim");
    // TH1F* proj_recULS_PteePrim = (TH1F*) hMPtDCA_ULS_primary_rec->ProjectionY("proj_recULS_PteePrim");
    // TH1F* proj_recULS_MeeCC = (TH1F*) hMPtDCA_ULS_CC_rec->ProjectionX("proj_recULS_MeeCC");
    // TH1F* proj_recULS_PteeCC = (TH1F*) hMPtDCA_ULS_CC_rec->ProjectionY("proj_recULS_PteeCC");
    // TH1F* proj_recULS_MeeBB = (TH1F*) hMPtDCA_ULS_BB_rec->ProjectionX("proj_recULS_MeeBB");
    // TH1F* proj_recULS_PteeBB = (TH1F*) hMPtDCA_ULS_BB_rec->ProjectionY("proj_recULS_PteeBB");

    // TH1F* proj_recULS_MCpidEle_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle->ProjectionX("proj_recULS_MCpidEle_Mee");
    // TH1F* proj_recULS_MCpidEle_Ptee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle->ProjectionY("proj_recULS_MCpidEle_Ptee");

    TH1F* proj_recLS_Mee = (TH1F*) hMPtDCA_LS_rec->ProjectionX("proj_recLS_Mee");
    TH1F* proj_recLS_Ptee = (TH1F*) hMPtDCA_LS_rec->ProjectionY("proj_recLS_Ptee");
    TH1F* proj_recLS_DCA = (TH1F*) hMPtDCA_LS_rec->ProjectionZ("proj_recLS_DCA");

    TH1F* proj_genULS_Mee = (TH1F*) hMPtDCA_ULS_gen->ProjectionX("proj_genULS_Mee");
    TH1F* proj_genULS_Ptee = (TH1F*) hMPtDCA_ULS_gen->ProjectionY("proj_genULS_Ptee");
    TH1F* proj_genULS_DCA = (TH1F*) hMPtDCA_ULS_gen->ProjectionZ("proj_genULS_DCA");

    TH1F* proj_genLS_Mee = (TH1F*) hMPtDCA_LS_gen->ProjectionX("proj_genLS_Mee");
    TH1F* proj_genLS_Ptee = (TH1F*) hMPtDCA_LS_gen->ProjectionY("proj_genLS_Ptee");
    TH1F* proj_genLS_DCA = (TH1F*) hMPtDCA_LS_gen->ProjectionZ("proj_genLS_DCA");
    // TH1F* proj_genLS_MeePrim = (TH1F*) hMPtDCA_LS_primary_gen->ProjectionX("proj_genLS_MeePrim");
    // TH1F* proj_genLS_MeeCC = (TH1F*) hMPtDCA_LS_CC_gen->ProjectionX("proj_genLS_MeeCC");
    // TH1F* proj_genLS_MeeBB = (TH1F*) hMPtDCA_LS_BB_gen->ProjectionX("proj_genLS_MeeBB");
    // TH1F* proj_genLS_PteePrim = (TH1F*) hMPtDCA_LS_primary_gen->ProjectionY("proj_genLS_PteePrim");
    // TH1F* proj_genLS_PteeCC = (TH1F*) hMPtDCA_LS_CC_gen->ProjectionY("proj_genLS_PteeCC");
    // TH1F* proj_genLS_PteeBB = (TH1F*) hMPtDCA_LS_BB_gen->ProjectionY("proj_genLS_PteeBB");

    proj_recULS_Mee->Sumw2();
    proj_recULS_Ptee->Sumw2();
    proj_recULS_DCA->Sumw2();
    // proj_recULS_MeePrim->Sumw2();
    // proj_recULS_PteePrim->Sumw2();
    // proj_recULS_MeeCC->Sumw2();
    // proj_recULS_PteeCC->Sumw2();
    // proj_recULS_MeeBB->Sumw2();
    // proj_recULS_PteeBB->Sumw2();
    proj_recLS_Mee->Sumw2();
    proj_recLS_Ptee->Sumw2();
    proj_recLS_DCA->Sumw2();
    // proj_recULS_MCpidEle_Mee->Sumw2();
    // proj_recULS_MCpidEle_Ptee->Sumw2();

    proj_genULS_Mee->Sumw2();
    proj_genULS_Ptee->Sumw2();
    proj_genULS_DCA->Sumw2();
    proj_genLS_Mee->Sumw2();
    proj_genLS_Ptee->Sumw2();
    proj_genLS_DCA->Sumw2();
    // proj_genLS_MeePrim->Sumw2();
    // proj_genLS_MeeCC->Sumw2();
    // proj_genLS_MeeBB->Sumw2();
    // proj_genLS_PteePrim->Sumw2();
    // proj_genLS_PteeCC->Sumw2();
    // proj_genLS_PteeBB->Sumw2();

    proj_recULS_Mee->Scale(1./nEventsCent);
    proj_recULS_Ptee->Scale(1./nEventsCent);
    proj_recULS_DCA->Scale(1./nEventsCent);
    // proj_recULS_MeePrim->Scale(1./nEventsCent);
    // proj_recULS_PteePrim->Scale(1./nEventsCent);
    // proj_recULS_MeeCC->Scale(1./nEventsCent);
    // proj_recULS_PteeCC->Scale(1./nEventsCent);
    // proj_recULS_MeeBB->Scale(1./nEventsCent);
    // proj_recULS_PteeBB->Scale(1./nEventsCent);
    proj_recLS_Mee->Scale(1./nEventsCent);
    proj_recLS_Ptee->Scale(1./nEventsCent);
    proj_recLS_DCA->Scale(1./nEventsCent);

    proj_genULS_Mee->Scale(1./nEventsCent);
    proj_genULS_Ptee->Scale(1./nEventsCent);
    proj_genULS_DCA->Scale(1./nEventsCent);
    proj_genLS_Mee->Scale(1./nEventsCent);
    proj_genLS_Ptee->Scale(1./nEventsCent);
    proj_genLS_DCA->Scale(1./nEventsCent);
    // proj_genLS_MeePrim->Scale(1./nEventsCent);
    // proj_genLS_MeeCC->Scale(1./nEventsCent);
    // proj_genLS_MeeBB->Scale(1./nEventsCent);
    // proj_genLS_PteePrim->Scale(1./nEventsCent);
    // proj_genLS_PteeCC->Scale(1./nEventsCent);
    // proj_genLS_PteeBB->Scale(1./nEventsCent);


    normalizeToBinWidth(proj_recULS_Mee);
    normalizeToBinWidth(proj_recULS_Ptee);
    normalizeToBinWidth(proj_recULS_DCA);
    // normalizeToBinWidth(proj_recULS_MeePrim);
    // normalizeToBinWidth(proj_recULS_PteePrim);
    // normalizeToBinWidth(proj_recULS_MeeCC);
    // normalizeToBinWidth(proj_recULS_PteeCC);
    // normalizeToBinWidth(proj_recULS_MeeBB);
    // normalizeToBinWidth(proj_recULS_PteeBB);
    normalizeToBinWidth(proj_recLS_Mee);
    normalizeToBinWidth(proj_recLS_Ptee);
    normalizeToBinWidth(proj_recLS_DCA);

    normalizeToBinWidth(proj_genULS_Mee);
    normalizeToBinWidth(proj_genULS_Ptee);
    normalizeToBinWidth(proj_genULS_DCA);
    normalizeToBinWidth(proj_genLS_Mee);
    normalizeToBinWidth(proj_genLS_Ptee);
    normalizeToBinWidth(proj_genLS_DCA);
    // normalizeToBinWidth(proj_genLS_MeePrim);
    // normalizeToBinWidth(proj_genLS_MeeCC);
    // normalizeToBinWidth(proj_genLS_MeeBB);
    // normalizeToBinWidth(proj_genLS_PteePrim);
    // normalizeToBinWidth(proj_genLS_PteeCC);
    // normalizeToBinWidth(proj_genLS_PteeBB);


  // ptRecTrackPrim->Rebin(2);
  // ptRecEleTrackPrim->Rebin(2);
  // ptRecPosTrackPrim->Rebin(2);
  // etaRecTrackPrim->Rebin(2);
  // etaRecEleTrackPrim->Rebin(2);
  // etaRecPosTrackPrim->Rebin(2);
  // phiRecTrackPrim->Rebin(10);
  // phiRecEleTrackPrim->Rebin(10);
  // phiRecPosTrackPrim->Rebin(10);
  // ptRecTrackCC->Rebin(2);
  // ptRecEleTrackCC->Rebin(2);
  // ptRecPosTrackCC->Rebin(2);
  // etaRecTrackCC->Rebin(2);
  // etaRecEleTrackCC->Rebin(2);
  // etaRecPosTrackCC->Rebin(2);
  // phiRecTrackCC->Rebin(10);
  // phiRecEleTrackCC->Rebin(10);
  // phiRecPosTrackCC->Rebin(10);
  // ptRecTrackBB->Rebin(2);
  // ptRecEleTrackBB->Rebin(2);
  // ptRecPosTrackBB->Rebin(2);
  // etaRecTrackBB->Rebin(2);
  // etaRecEleTrackBB->Rebin(2);
  // etaRecPosTrackBB->Rebin(2);
  // phiRecTrackBB->Rebin(10);
  // phiRecEleTrackBB->Rebin(10);
  // phiRecPosTrackBB->Rebin(10);
  // ptGenTrackPrim->Rebin(2);
  // ptGenEleTrackPrim->Rebin(2);
  // ptGenPosTrackPrim->Rebin(2);
  // etaGenTrackPrim->Rebin(2);
  // etaGenEleTrackPrim->Rebin(2);
  // etaGenPosTrackPrim->Rebin(2);
  // phiGenTrackPrim->Rebin(10);
  // phiGenEleTrackPrim->Rebin(10);
  // phiGenPosTrackPrim->Rebin(10);
  // ptGenTrackCC->Rebin(2);
  // ptGenEleTrackCC->Rebin(2);
  // ptGenPosTrackCC->Rebin(2);
  // etaGenTrackCC->Rebin(2);
  // etaGenEleTrackCC->Rebin(2);
  // etaGenPosTrackCC->Rebin(2);
  // phiGenTrackCC->Rebin(10);
  // phiGenEleTrackCC->Rebin(10);
  // phiGenPosTrackCC->Rebin(10);
  // ptGenTrackBB->Rebin(2);
  // ptGenEleTrackBB->Rebin(2);
  // ptGenPosTrackBB->Rebin(2);
  // etaGenTrackBB->Rebin(2);
  // etaGenEleTrackBB->Rebin(2);
  // etaGenPosTrackBB->Rebin(2);
  // phiGenTrackBB->Rebin(10);
  // phiGenEleTrackBB->Rebin(10);
  // phiGenPosTrackBB->Rebin(10);

  // ptRecTrackBeforeSmearing->Rebin(2);
  // etaRecTrackBeforeSmearing->Rebin(2);
  // phiRecTrackBeforeSmearing->Rebin(10);
  // ptRecTrackAfterSmearing->Rebin(2);
  // etaRecTrackAfterSmearing->Rebin(2);
  // phiRecTrackAfterSmearing->Rebin(10);


 TH1F* ptRecTrackBeforeSmearing_rebin = (TH1F*) ptRecTrackBeforeSmearing->Rebin(nbinspt_proj,"ptRecTrackBeforeSmearing_rebin",&pt_bin_proj[0]);
 TH1F* ptRecTrackAfterSmearing_rebin = (TH1F*) ptRecTrackAfterSmearing->Rebin(nbinspt_proj,"ptRecTrackAfterSmearing_rebin",&pt_bin_proj[0]);
 TH1F* ptRecTrackAfterKineCuts_rebin = (TH1F*) ptRecTrackAfterKineCuts->Rebin(nbinspt_proj,"ptRecTrackAfterKineCuts_rebin",&pt_bin_proj[0]);

 TH1F* ptRecTrackElePos_rebin = (TH1F*) ptRecTrackElePos->Rebin(nbinspt_proj,"ptRecTrackElePos_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackElePos_rebin = (TH1F*) ptGenTrackElePos->Rebin(nbinspt_proj,"ptGenTrackElePos_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackAll_rebin = (TH1F*) ptGenTrackAll->Rebin(nbinspt_proj,"ptGenTrackAll_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackMuon_rebin = (TH1F*) ptGenTrackMuon->Rebin(nbinspt_proj,"ptGenTrackMuon_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackPion_rebin = (TH1F*) ptGenTrackPion->Rebin(nbinspt_proj,"ptGenTrackPion_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackKaon_rebin = (TH1F*) ptGenTrackKaon->Rebin(nbinspt_proj,"ptGenTrackKaon_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackProton_rebin = (TH1F*) ptGenTrackProton->Rebin(nbinspt_proj,"ptGenTrackProton_rebin",&pt_bin_proj[0]);
 TH1F* ptGenSmearedTrackElePos_rebin = (TH1F*) ptGenSmearedTrackElePos->Rebin(nbinspt_proj,"ptGenSmearedTrackElePos_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackElePos_beforeKineCuts_rebin = (TH1F*) ptGenTrackElePos_beforeKineCuts->Rebin(nbinspt_proj,"ptGenTrackElePos_beforeKineCuts_rebin",&pt_bin_proj[0]);
 // TH1F* ptGenSmearedTrackElePos_beforeKineCuts_rebin = (TH1F*) ptGenSmearedTrackElePos_beforeKineCuts->Rebin(nbinspt_proj,"ptGenSmearedTrackElePos_beforeKineCuts_rebin",&pt_bin_proj[0]);


 TH1F* ptAllRecTrack_rebin = (TH1F*) ptAllRecTrack->Rebin(nbinspt_proj,"ptAllRecTrack_rebin",&pt_bin_proj[0]);
 TH1F* ptRecMuonTrack_rebin = (TH1F*) ptRecMuonTrack->Rebin(nbinspt_proj,"ptRecMuonTrack_rebin",&pt_bin_proj[0]);
 TH1F* ptRecPionTrack_rebin = (TH1F*) ptRecPionTrack->Rebin(nbinspt_proj,"ptRecPionTrack_rebin",&pt_bin_proj[0]);
 TH1F* ptRecKaonTrack_rebin = (TH1F*) ptRecKaonTrack->Rebin(nbinspt_proj,"ptRecKaonTrack_rebin",&pt_bin_proj[0]);
 TH1F* ptRecProtonTrack_rebin = (TH1F*) ptRecProtonTrack->Rebin(nbinspt_proj,"ptRecProtonTrack_rebin",&pt_bin_proj[0]);

 TH1F* ptRecNegTrack_rebin = (TH1F*) ptRecNegTrack->Rebin(nbinspt_proj,"ptRecNegTrack_rebin",&pt_bin_proj[0]);
 TH1F* ptRecPosTrack_rebin = (TH1F*) ptRecPosTrack->Rebin(nbinspt_proj,"ptRecPosTrack_rebin",&pt_bin_proj[0]);





 TH1F* ptRecTrackEle_rebin = (TH1F*) ptRecTrackEle->Rebin(nbinspt_proj,"ptRecTrackEle_rebin",&pt_bin_proj[0]);
 TH1F* ptRecTrackPos_rebin = (TH1F*) ptRecTrackPos->Rebin(nbinspt_proj,"ptRecTrackPos_rebin",&pt_bin_proj[0]);
 TH1F* ptGenSmearedTrackEle_rebin = (TH1F*) ptGenSmearedTrackEle->Rebin(nbinspt_proj,"ptGenSmearedTrackEle_rebin",&pt_bin_proj[0]);
 TH1F* ptGenSmearedTrackPos_rebin = (TH1F*) ptGenSmearedTrackPos->Rebin(nbinspt_proj,"ptGenSmearedTrackPos_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackEle_rebin = (TH1F*) ptGenTrackEle->Rebin(nbinspt_proj,"ptGenTrackEle_rebin",&pt_bin_proj[0]);
 TH1F* ptGenTrackPos_rebin = (TH1F*) ptGenTrackPos->Rebin(nbinspt_proj,"ptGenTrackPos_rebin",&pt_bin_proj[0]);



  // ptRecTrackElePos->Rebin(2);
  // etaRecTrackElePos->Rebin(2);
  // phiRecTrackElePos->Rebin(10);

  // ptRecTrackEle->Rebin(2);
  // etaRecTrackEle->Rebin(2);
  // phiRecTrackEle->Rebin(10);
  // ptRecTrackPos->Rebin(2);
  // etaRecTrackPos->Rebin(2);
  // phiRecTrackPos->Rebin(10);
  // ptRecNegTrack->Rebin(2);
  // etaRecNegTrack->Rebin(2);
  // phiRecNegTrack->Rebin(10);
  // ptRecPosTrack->Rebin(2);
  // etaRecPosTrack->Rebin(2);
  // phiRecPosTrack->Rebin(10);

  // ptAllRecTrack->Rebin(2);
  // etaAllRecTrack->Rebin(2);
  // phiAllRecTrack->Rebin(10);
  // ptRecMuonTrack->Rebin(2);
  // etaRecMuonTrack->Rebin(2);
  // phiRecMuonTrack->Rebin(10);
  // ptRecPionTrack->Rebin(2);
  // etaRecPionTrack->Rebin(2);
  // phiRecPionTrack->Rebin(10);
  // ptRecKaonTrack->Rebin(2);
  // etaRecKaonTrack->Rebin(2);
  // phiRecKaonTrack->Rebin(10);
  // ptRecProtonTrack->Rebin(2);
  // etaRecProtonTrack->Rebin(2);
  // phiRecProtonTrack->Rebin(10);

  // ptGenTrackEle->Rebin(2);
  // etaGenTrackEle->Rebin(2);
  // phiGenTrackEle->Rebin(10);
  // ptGenTrackPos->Rebin(2);
  // etaGenTrackPos->Rebin(2);
  // phiGenTrackPos->Rebin(10);

  // ptRecTrackEtaCut_1->Rebin(2);
  // ptRecTrackEtaCut_2->Rebin(2);
  // ptRecTrackEtaCut_3->Rebin(2);
  // ptRecTrackEtaCut_4->Rebin(2);
  // ptRecTrackEtaCut_5->Rebin(2);
  // ptRecTrackEtaCut_6->Rebin(2);
  // ptRecTrackEtaCut_7->Rebin(2);
  // ptRecTrackEtaCut_8->Rebin(2);
  // ptRecTrackEtaCut_9->Rebin(2);
  // ptRecTrackEtaCut_10->Rebin(2);
  //
  // etaRecTrackEtaCut_1->Rebin(2);
  // etaRecTrackEtaCut_2->Rebin(2);
  // etaRecTrackEtaCut_3->Rebin(2);
  // etaRecTrackEtaCut_4->Rebin(2);
  // etaRecTrackEtaCut_5->Rebin(2);
  // etaRecTrackEtaCut_6->Rebin(2);
  // etaRecTrackEtaCut_7->Rebin(2);
  // etaRecTrackEtaCut_8->Rebin(2);
  // etaRecTrackEtaCut_9->Rebin(2);
  // etaRecTrackEtaCut_10->Rebin(2);
  //
  // phiRecTrackEtaCut_1->Rebin(10);
  // phiRecTrackEtaCut_2->Rebin(10);
  // phiRecTrackEtaCut_3->Rebin(10);
  // phiRecTrackEtaCut_4->Rebin(10);
  // phiRecTrackEtaCut_5->Rebin(10);
  // phiRecTrackEtaCut_6->Rebin(10);
  // phiRecTrackEtaCut_7->Rebin(10);
  // phiRecTrackEtaCut_8->Rebin(10);
  // phiRecTrackEtaCut_9->Rebin(10);
  // phiRecTrackEtaCut_10->Rebin(10);

  TH1F* projNSigmas_afterCuts_trueElec;
  TH1F* projNSigmas_afterCuts_truePion;
  TH1F* projNSigmas_trueElec;
  TH1F* projNSigmas_truePion;

  TH1F* ptGenEle;
  TH1F* etaGenEle;
  TH1F* phiGenEle;
  TH1F* ptGenPos;
  TH1F* etaGenPos;
  TH1F* phiGenPos;

  TH1F* ptEffElePosGen;
  TH1F* ptEffElePosGenSmeared;
  TH1F* ptEffElePosGen_woPID;
  TH1F* ptEffEle;
  // TH1F* ptEffElePrim;
  // TH1F* ptEffEleCC;
  // TH1F* ptEffEleBB;
  TH1F* etaEffEle;
  // TH1F* etaEffElePrim;
  // TH1F* etaEffEleCC;
  // TH1F* etaEffEleBB;
  TH1F* phiEffEle;
  // TH1F* phiEffElePrim;
  // TH1F* phiEffEleCC;
  // TH1F* phiEffEleBB;
  TH1F* ptEffPos;
  // TH1F* ptEffPosPrim;
  // TH1F* ptEffPosCC;
  // TH1F* ptEffPosBB;
  TH1F* etaEffPos;
  // TH1F* etaEffPosPrim;
  // TH1F* etaEffPosCC;
  // TH1F* etaEffPosBB;
  TH1F* phiEffPos;
  // TH1F* phiEffPosPrim;
  // TH1F* phiEffPosCC;
  // TH1F* phiEffPosBB;

  TH2F* ptEtaEffElePos;

  TH1F* ptPairEffULS;
  TH1F* massPairEffULS;

  if(bPlotEfficiency){
    //track efficiencies
    // ptEffElePrim = (TH1F*) ptRecEleTrackPrim->Clone();
    // ptEffEleCC = (TH1F*) ptRecEleTrackCC->Clone();
    // ptEffEleBB = (TH1F*) ptRecEleTrackBB->Clone();
    // etaEffElePrim = (TH1F*) etaRecEleTrackPrim->Clone();
    // etaEffEleCC = (TH1F*) etaRecEleTrackCC->Clone();
    // etaEffEleBB = (TH1F*) etaRecEleTrackBB->Clone();
    // phiEffElePrim = (TH1F*) phiRecEleTrackPrim->Clone();
    // phiEffEleCC = (TH1F*) phiRecEleTrackCC->Clone();
    // phiEffEleBB = (TH1F*) phiRecEleTrackBB->Clone();
    // ptEffPosPrim = (TH1F*) ptRecPosTrackPrim->Clone();
    // ptEffPosCC = (TH1F*) ptRecPosTrackCC->Clone();
    // ptEffPosBB = (TH1F*) ptRecPosTrackBB->Clone();
    // etaEffPosPrim = (TH1F*) etaRecPosTrackPrim->Clone();
    // etaEffPosCC = (TH1F*) etaRecPosTrackCC->Clone();
    // etaEffPosBB = (TH1F*) etaRecPosTrackBB->Clone();
    // phiEffPosPrim = (TH1F*) phiRecPosTrackPrim->Clone();
    // phiEffPosCC = (TH1F*) phiRecPosTrackCC->Clone();
    // phiEffPosBB = (TH1F*) phiRecPosTrackBB->Clone();

    // ptEffElePrim->Sumw2();
    // ptEffEleCC->Sumw2();
    // ptEffEleBB->Sumw2();
    // etaEffElePrim->Sumw2();
    // etaEffEleCC->Sumw2();
    // etaEffEleBB->Sumw2();
    // phiEffElePrim->Sumw2();
    // phiEffEleCC->Sumw2();
    // phiEffEleBB->Sumw2();
    // ptEffPosPrim->Sumw2();
    // ptEffPosCC->Sumw2();
    // ptEffPosBB->Sumw2();
    // etaEffPosPrim->Sumw2();
    // etaEffPosCC->Sumw2();
    // etaEffPosBB->Sumw2();
    // phiEffPosPrim->Sumw2();
    // phiEffPosCC->Sumw2();
    // phiEffPosBB->Sumw2();

    // ptEffElePrim->Divide(ptEffElePrim,ptGenEleTrackPrim,1,1,"B");
    // ptEffEleCC->Divide(ptEffEleCC,ptGenEleTrackCC,1,1,"B");
    // ptEffEleBB->Divide(ptEffEleBB,ptGenEleTrackBB,1,1,"B");
    // etaEffElePrim->Divide(etaEffElePrim,etaGenEleTrackPrim,1,1,"B");
    // etaEffEleCC->Divide(etaEffEleCC,etaGenEleTrackCC,1,1,"B");
    // etaEffEleBB->Divide(etaEffEleBB,etaGenEleTrackBB,1,1,"B");
    // phiEffElePrim->Divide(phiEffElePrim,phiGenEleTrackPrim,1,1,"B");
    // phiEffEleCC->Divide(phiEffEleCC,phiGenEleTrackCC,1,1,"B");
    // phiEffEleBB->Divide(phiEffEleBB,phiGenEleTrackBB,1,1,"B");
    // ptEffPosPrim->Divide(ptEffPosPrim,ptGenPosTrackPrim,1,1,"B");
    // ptEffPosCC->Divide(ptEffPosCC,ptGenPosTrackCC,1,1,"B");
    // ptEffPosBB->Divide(ptEffPosBB,ptGenPosTrackBB,1,1,"B");
    // etaEffPosPrim->Divide(etaEffPosPrim,etaGenPosTrackPrim,1,1,"B");
    // etaEffPosCC->Divide(etaEffPosCC,etaGenPosTrackCC,1,1,"B");
    // etaEffPosBB->Divide(etaEffPosBB,etaGenPosTrackBB,1,1,"B");
    // phiEffPosPrim->Divide(phiEffPosPrim,phiGenPosTrackPrim,1,1,"B");
    // phiEffPosCC->Divide(phiEffPosCC,phiGenPosTrackCC,1,1,"B");
    // phiEffPosBB->Divide(phiEffPosBB,phiGenPosTrackBB,1,1,"B");

    // ptGenEle = (TH1F*) ptGenEleTrackPrim->Clone();
    // ptGenEle->Sumw2();
    // ptGenEle->Add(ptGenEleTrackCC, 1);
    // ptGenEle->Add(ptGenEleTrackBB, 1);
    // etaGenEle = (TH1F*) etaGenEleTrackPrim->Clone();
    // etaGenEle->Sumw2();
    // etaGenEle->Add(etaGenEleTrackCC, 1);
    // etaGenEle->Add(etaGenEleTrackBB, 1);
    // phiGenEle = (TH1F*) phiGenEleTrackPrim->Clone();
    // phiGenEle->Sumw2();
    // phiGenEle->Add(phiGenEleTrackCC, 1);
    // phiGenEle->Add(phiGenEleTrackBB, 1);
    // ptGenPos = (TH1F*) ptGenPosTrackPrim->Clone();
    // ptGenPos->Sumw2();
    // ptGenPos->Add(ptGenPosTrackCC, 1);
    // ptGenPos->Add(ptGenPosTrackBB, 1);
    // etaGenPos = (TH1F*) etaGenPosTrackPrim->Clone();
    // etaGenPos->Sumw2();
    // etaGenPos->Add(etaGenPosTrackCC, 1);
    // etaGenPos->Add(etaGenPosTrackBB, 1);
    // phiGenPos = (TH1F*) phiGenPosTrackPrim->Clone();
    // phiGenPos->Sumw2();
    // phiGenPos->Add(phiGenPosTrackCC, 1);
    // phiGenPos->Add(phiGenPosTrackBB, 1);


    ptEffElePosGen = (TH1F*) ptRecTrackElePos_rebin->Clone("eff_pT_singleElePos_Gen");
    ptEffElePosGen->Sumw2();
    ptEffElePosGen->Divide(ptEffElePosGen,ptGenTrackElePos_rebin,1,1,"B");
    ptEffElePosGenSmeared = (TH1F*) ptRecTrackElePos_rebin->Clone("eff_pT_singleElePos_GenSmeared");
    ptEffElePosGenSmeared->Sumw2();
    ptEffElePosGenSmeared->Divide(ptEffElePosGenSmeared,ptGenSmearedTrackElePos_rebin,1,1,"B");

    ptEffElePosGen_woPID = (TH1F*) ptRecTrackAfterKineCuts_rebin->Clone("eff_pT_singleElePos_Gen_woPID");
    ptEffElePosGen_woPID->Sumw2();
    ptEffElePosGen_woPID->Divide(ptEffElePosGen_woPID,ptGenTrackElePos_rebin,1,1,"B");

    ptEffEle = (TH1F*) ptRecTrackEle_rebin->Clone("eff_pT_singleElectrons");
    ptEffEle->Sumw2();
    // ptEffEle->Divide(ptEffEle,ptGenSmearedTrackEle_rebin,1,1,"B");
    ptEffEle->Divide(ptEffEle,ptGenTrackEle_rebin,1,1,"B");
    etaEffEle = (TH1F*) etaRecTrackEle->Clone("eff_eta_singleElectrons");
    etaEffEle->Sumw2();
    // etaEffEle->Divide(etaEffEle,etaGenSmearedTrackEle,1,1,"B");
    etaEffEle->Divide(etaEffEle,etaGenTrackEle,1,1,"B");
    phiEffEle = (TH1F*) phiRecTrackEle->Clone("eff_phi_singleElectrons");
    phiEffEle->Sumw2();
    // phiEffEle->Divide(phiEffEle,phiGenSmearedTrackEle,1,1,"B");
    phiEffEle->Divide(phiEffEle,phiGenTrackEle,1,1,"B");
    // ptEffPos = (TH1F*) ptRecPosTrackPrim->Clone("eff_pT_singlePositrons");
    // ptEffPos->Add(ptRecPosTrackCC, 1);
    // ptEffPos->Add(ptRecPosTrackBB, 1);
    ptEffPos = (TH1F*) ptRecTrackPos_rebin->Clone("eff_pT_singlePositrons");
    ptEffPos->Sumw2();
    // ptEffPos->Divide(ptEffPos,ptGenSmearedTrackPos_rebin,1,1,"B");
    ptEffPos->Divide(ptEffPos,ptGenTrackPos_rebin,1,1,"B");
    etaEffPos = (TH1F*) etaRecTrackPos->Clone("eff_eta_singlePositrons");
    etaEffPos->Sumw2();
    // etaEffPos->Divide(etaEffPos,etaGenSmearedTrackPos,1,1,"B");
    etaEffPos->Divide(etaEffPos,etaGenTrackPos,1,1,"B");
    phiEffPos = (TH1F*) phiRecTrackPos->Clone("eff_phi_singlePositrons");
    phiEffPos->Sumw2();
    // phiEffPos->Divide(phiEffPos,phiGenSmearedTrackPos,1,1,"B");
    phiEffPos->Divide(phiEffPos,phiGenTrackPos,1,1,"B");

    // single track efficiency TH2 as function of pt and eta
    ptEtaEffElePos = (TH2F*) ptEtaRecTrackElePos->Clone("eff_ptEta_singleElePos");
    ptEtaEffElePos->Sumw2();
    ptEtaEffElePos->Divide(ptEtaEffElePos,ptEtaGenTrackElePos,1,1,"B");

   // pair efficiencies
   // // ptPairEffULS = (TH1F*) proj_recULS_Ptee->Clone();
   // ptPairEffULS = (TH1F*) proj_recULS_MCpidEle_Ptee->Clone();
   // // massPairEffULS = (TH1F*) proj_recULS_Mee->Clone();
   // massPairEffULS = (TH1F*) proj_recULS_MCpidEle_Mee->Clone();
   // ptPairEffULS->Sumw2();
   // massPairEffULS->Sumw2();
   // ptPairEffULS->Divide(ptPairEffULS,proj_genULS_Ptee,1,1,"B");
   // massPairEffULS->Divide(massPairEffULS,proj_genULS_Mee,1,1,"B");
  }


  TH1F* hMuonRejectionFactorPt;
  TH1F* hPionRejectionFactorPt;
  TH1F* hKaonRejectionFactorPt;
  TH1F* hProtonRejectionFactorPt;
  TH1F* hTotalRejectionFactorPt;
  TH1F* ptRecContamination;
  ptRecContamination = (TH1F*)  ptAllRecTrack_rebin->Clone("pt_Total_contamination");
  ptRecContamination->Add(ptRecTrackElePos_rebin,-1);


  hMuonRejectionFactorPt = (TH1F*) ptGenTrackMuon_rebin->Clone("pT_Muon_Rejectionfactor");
  hPionRejectionFactorPt = (TH1F*) ptGenTrackPion_rebin->Clone("pT_Pion_Rejectionfactor");
  hKaonRejectionFactorPt = (TH1F*) ptGenTrackKaon_rebin->Clone("pT_Kaon_Rejectionfactor");
  hProtonRejectionFactorPt = (TH1F*) ptGenTrackProton_rebin->Clone("pT_Proton_Rejectionfactor");
  hTotalRejectionFactorPt = (TH1F*) ptGenTrackAll_rebin->Clone("pT_Total_Rejectionfactor");
  hTotalRejectionFactorPt->Add(ptGenTrackElePos_rebin,-1);

  hMuonRejectionFactorPt->Sumw2();
  hPionRejectionFactorPt->Sumw2();
  hKaonRejectionFactorPt->Sumw2();
  hProtonRejectionFactorPt->Sumw2();
  hTotalRejectionFactorPt->Sumw2();

  hMuonRejectionFactorPt->Divide(hMuonRejectionFactorPt, ptRecMuonTrack_rebin,1,1,"B");
  hPionRejectionFactorPt->Divide(hPionRejectionFactorPt, ptRecPionTrack_rebin,1,1,"B");
  hKaonRejectionFactorPt->Divide(hKaonRejectionFactorPt, ptRecKaonTrack_rebin,1,1,"B");
  hProtonRejectionFactorPt->Divide(hProtonRejectionFactorPt, ptRecProtonTrack_rebin,1,1,"B");
  hTotalRejectionFactorPt->Divide(hTotalRejectionFactorPt, ptRecContamination,1,1,"B");


  make3HistNice(hMuonRejectionFactorPt, kBlue+1);
  make3HistNice(hPionRejectionFactorPt, kRed+2);
  make3HistNice(hKaonRejectionFactorPt, kGreen+2);
  make3HistNice(hProtonRejectionFactorPt, kOrange+2);
  make3HistNice(hTotalRejectionFactorPt, kBlack);
  hTotalRejectionFactorPt->SetMarkerStyle(24);
  hTotalRejectionFactorPt->SetMarkerSize(0.9);




  TH1F* hPurityRecPtNeg;
  TH1F* hPurityRecEtaNeg;
  TH1F* hPurityRecPhiNeg;
  TH1F* hPurityRecPtPos;
  TH1F* hPurityRecEtaPos;
  TH1F* hPurityRecPhiPos;
  TH1F* hTotalPureContaminationRecPt;
  TH1F* hTotalPureContaminationRecEta;
  TH1F* hTotalPureContaminationRecPhi;
  TH1F* hPureContaminationRecPtNeg;
  TH1F* hPureContaminationRecEtaNeg;
  TH1F* hPureContaminationRecPhiNeg;
  TH1F* hPureContaminationRecPtPos;
  TH1F* hPureContaminationRecEtaPos;
  TH1F* hPureContaminationRecPhiPos;
  TH1F* hTotalRecTrackPt;
  TH1F* hTotalRecTrackEta;
  TH1F* hTotalRecTrackPhi;
  TH1F* hMuonContaminationRecPt;
  TH1F* hMuonContaminationRecEta;
  TH1F* hMuonContaminationRecPhi;
  TH1F* hPionContaminationRecPt;
  TH1F* hPionContaminationRecEta;
  TH1F* hPionContaminationRecPhi;
  TH1F* hKaonContaminationRecPt;
  TH1F* hKaonContaminationRecEta;
  TH1F* hKaonContaminationRecPhi;
  TH1F* hProtonContaminationRecPt;
  TH1F* hProtonContaminationRecEta;
  TH1F* hProtonContaminationRecPhi;


if (bPlotTrackContamination) {
  hTotalRecTrackPt = (TH1F*)  ptAllRecTrack_rebin->Clone();
  hTotalRecTrackEta = (TH1F*) etaAllRecTrack->Clone();
  hTotalRecTrackPhi = (TH1F*) phiAllRecTrack->Clone();
  hPurityRecPtNeg = (TH1F*)  ptRecNegTrack->Clone();
  hPurityRecEtaNeg = (TH1F*) etaRecNegTrack->Clone();
  hPurityRecPhiNeg = (TH1F*) phiRecNegTrack->Clone();
  hPurityRecPtPos = (TH1F*)  ptRecPosTrack->Clone();
  hPurityRecEtaPos = (TH1F*) etaRecPosTrack->Clone();
  hPurityRecPhiPos = (TH1F*) phiRecPosTrack->Clone();
  hTotalRecTrackPt->Sumw2();
  hTotalRecTrackEta->Sumw2();
  hTotalRecTrackPhi->Sumw2();
  hPurityRecPtNeg->Sumw2();
  hPurityRecEtaNeg->Sumw2();
  hPurityRecPhiNeg->Sumw2();
  hPurityRecPtPos->Sumw2();
  hPurityRecEtaPos->Sumw2();
  hPurityRecPhiPos->Sumw2();
  hPurityRecPtNeg->Divide(ptRecTrackEle,ptRecNegTrack,1,1,"B");
  hPurityRecEtaNeg->Divide(etaRecTrackEle,etaRecNegTrack,1,1,"B");
  hPurityRecPhiNeg->Divide(phiRecTrackEle,phiRecNegTrack,1,1,"B");
  hPurityRecPtPos->Divide(ptRecTrackPos,ptRecPosTrack,1,1,"B");
  hPurityRecEtaPos->Divide(etaRecTrackPos,etaRecPosTrack,1,1,"B");
  hPurityRecPhiPos->Divide(phiRecTrackPos,phiRecPosTrack,1,1,"B");

  hTotalPureContaminationRecPt = (TH1F*)  ptAllRecTrack_rebin->Clone("TotalContaminationPt");
  hTotalPureContaminationRecEta = (TH1F*) etaAllRecTrack->Clone("TotalContaminationEta");
  hTotalPureContaminationRecPhi = (TH1F*) phiAllRecTrack->Clone("TotalContaminationPhi");
  hPureContaminationRecPtNeg = (TH1F*)  ptRecNegTrack_rebin->Clone();
  hPureContaminationRecEtaNeg = (TH1F*) etaRecNegTrack->Clone();
  hPureContaminationRecPhiNeg = (TH1F*) phiRecNegTrack->Clone();
  hPureContaminationRecPtPos = (TH1F*)  ptRecPosTrack_rebin->Clone();
  hPureContaminationRecEtaPos = (TH1F*) etaRecPosTrack->Clone();
  hPureContaminationRecPhiPos = (TH1F*) phiRecPosTrack->Clone();
  hMuonContaminationRecPt = (TH1F*)    ptRecMuonTrack_rebin->Clone();
  hMuonContaminationRecEta = (TH1F*)   etaRecMuonTrack->Clone();
  hMuonContaminationRecPhi = (TH1F*)   phiRecMuonTrack->Clone();
  hPionContaminationRecPt = (TH1F*)    ptRecPionTrack_rebin->Clone();
  hPionContaminationRecEta = (TH1F*)   etaRecPionTrack->Clone();
  hPionContaminationRecPhi = (TH1F*)   phiRecPionTrack->Clone();
  hKaonContaminationRecPt = (TH1F*)    ptRecKaonTrack_rebin->Clone();
  hKaonContaminationRecEta = (TH1F*)   etaRecKaonTrack->Clone();
  hKaonContaminationRecPhi = (TH1F*)   phiRecKaonTrack->Clone();
  hProtonContaminationRecPt = (TH1F*)  ptRecProtonTrack_rebin->Clone();
  hProtonContaminationRecEta = (TH1F*) etaRecProtonTrack->Clone();
  hProtonContaminationRecPhi = (TH1F*) phiRecProtonTrack->Clone();
  hPureContaminationRecPtNeg->Sumw2();
  hPureContaminationRecEtaNeg->Sumw2();
  hPureContaminationRecPhiNeg->Sumw2();
  hPureContaminationRecPtPos->Sumw2();
  hPureContaminationRecEtaPos->Sumw2();
  hPureContaminationRecPhiPos->Sumw2();
  hMuonContaminationRecPt->Sumw2();
  hMuonContaminationRecEta->Sumw2();
  hMuonContaminationRecPhi->Sumw2();
  hPionContaminationRecPt->Sumw2();
  hPionContaminationRecEta->Sumw2();
  hPionContaminationRecPhi->Sumw2();
  hKaonContaminationRecPt->Sumw2();
  hKaonContaminationRecEta->Sumw2();
  hKaonContaminationRecPhi->Sumw2();
  hProtonContaminationRecPt->Sumw2();
  hProtonContaminationRecEta->Sumw2();
  hProtonContaminationRecPhi->Sumw2();
  // hTotalPureContaminationRecPt->Add(ptRecTrackElePos_rebin,-1);
  // hTotalPureContaminationRecEta->Add(etaRecTrackElePos,-1);
  // hTotalPureContaminationRecPhi->Add(phiRecTrackElePos,-1);
  // hTotalPureContaminationRecPt->Sumw2();
  // hTotalPureContaminationRecEta->Sumw2();
  // hTotalPureContaminationRecPhi->Sumw2();

  double binError = 0., binValA =0., binValB=0., errorA = 0., errorB = 0.;
  for (int iBin = 1; iBin < hTotalPureContaminationRecPt->GetNbinsX()+1; iBin++) {
   binValA = ptRecTrackElePos_rebin->GetBinContent(iBin);
   binValB = hTotalPureContaminationRecPt->GetBinContent(iBin);
   errorA = ptRecTrackElePos_rebin->GetBinError(iBin);
   errorB = hTotalPureContaminationRecPt->GetBinError(iBin);
   binError = TMath::Sqrt(TMath::Abs(errorA*errorA - errorB*errorB));
   hTotalPureContaminationRecPt->SetBinContent(iBin,binValB-binValA);
   hTotalPureContaminationRecPt->SetBinError(iBin,binError);
  }
  for (int iBin = 1; iBin < hTotalPureContaminationRecEta->GetNbinsX()+1; iBin++) {
   binValA = etaRecTrackElePos->GetBinContent(iBin);
   binValB = hTotalPureContaminationRecEta->GetBinContent(iBin);
   errorA = etaRecTrackElePos->GetBinError(iBin);
   errorB = hTotalPureContaminationRecEta->GetBinError(iBin);
   binError = TMath::Sqrt(TMath::Abs(errorA*errorA - errorB*errorB));
   hTotalPureContaminationRecEta->SetBinContent(iBin,binValB-binValA);
   hTotalPureContaminationRecEta->SetBinError(iBin,binError);
  }
  for (int iBin = 1; iBin < hTotalPureContaminationRecPhi->GetNbinsX()+1; iBin++) {
   binValA = phiRecTrackElePos->GetBinContent(iBin);
   binValB = hTotalPureContaminationRecPhi->GetBinContent(iBin);
   errorA = phiRecTrackElePos->GetBinError(iBin);
   errorB = hTotalPureContaminationRecPhi->GetBinError(iBin);
   binError = TMath::Sqrt(TMath::Abs(errorA*errorA - errorB*errorB));
   hTotalPureContaminationRecPhi->SetBinContent(iBin,binValB-binValA);
   hTotalPureContaminationRecPhi->SetBinError(iBin,binError);
  }

  hPureContaminationRecPtNeg->Add(ptRecTrackEle_rebin,-1);
  hPureContaminationRecEtaNeg->Add(etaRecTrackEle,-1);
  hPureContaminationRecPhiNeg->Add(phiRecTrackEle,-1);
  hPureContaminationRecPtPos->Add(ptRecTrackPos_rebin,-1);
  hPureContaminationRecEtaPos->Add(etaRecTrackPos,-1);
  hPureContaminationRecPhiPos->Add(phiRecTrackPos,-1);
  hPureContaminationRecPtNeg->Divide(hPureContaminationRecPtNeg,hTotalRecTrackPt,1,1,"B");
  hPureContaminationRecEtaNeg->Divide(hPureContaminationRecEtaNeg,hTotalRecTrackEta,1,1,"B");
  hPureContaminationRecPhiNeg->Divide(hPureContaminationRecPhiNeg,hTotalRecTrackPhi,1,1,"B");
  hPureContaminationRecEtaPos->Divide(hPureContaminationRecEtaPos,hTotalRecTrackEta,1,1,"B");
  hPureContaminationRecPhiPos->Divide(hPureContaminationRecPhiPos,hTotalRecTrackPhi,1,1,"B");
  hPureContaminationRecPtPos->Divide(hPureContaminationRecPtPos,hTotalRecTrackPt,1,1,"B");
  hMuonContaminationRecPt->Divide(hMuonContaminationRecPt, hTotalRecTrackPt,1,1,"B");
  hMuonContaminationRecEta->Divide(hMuonContaminationRecEta, hTotalRecTrackEta,1,1,"B");
  hMuonContaminationRecPhi->Divide(hMuonContaminationRecPhi, hTotalRecTrackPhi,1,1,"B");
  hPionContaminationRecPt->Divide(hPionContaminationRecPt, hTotalRecTrackPt,1,1,"B");
  hPionContaminationRecEta->Divide(hPionContaminationRecEta, hTotalRecTrackEta,1,1,"B");
  hPionContaminationRecPhi->Divide(hPionContaminationRecPhi, hTotalRecTrackPhi,1,1,"B");
  hKaonContaminationRecPt->Divide(hKaonContaminationRecPt, hTotalRecTrackPt,1,1,"B");
  hKaonContaminationRecPhi->Divide(hKaonContaminationRecPhi, hTotalRecTrackPhi,1,1,"B");
  hKaonContaminationRecEta->Divide(hKaonContaminationRecEta, hTotalRecTrackEta,1,1,"B");
  hProtonContaminationRecPt->Divide(hProtonContaminationRecPt, hTotalRecTrackPt,1,1,"B");
  hProtonContaminationRecEta->Divide(hProtonContaminationRecEta, hTotalRecTrackEta,1,1,"B");
  hProtonContaminationRecPhi->Divide(hProtonContaminationRecPhi, hTotalRecTrackPhi,1,1,"B");

  // hTotalPureContaminationRecPt->Divide(hTotalPureContaminationRecPt,hTotalRecTrackPt,1,1,"B");
  // hTotalPureContaminationRecEta->Divide(hTotalPureContaminationRecEta,hTotalRecTrackEta,1,1,"B");
  // hTotalPureContaminationRecPhi->Divide(hTotalPureContaminationRecPhi,hTotalRecTrackPhi,1,1,"B");
  for (int iBin = 1; iBin < hTotalPureContaminationRecPt->GetNbinsX()+1; iBin++) {
    binValA = hTotalPureContaminationRecPt->GetBinContent(iBin);
    binValB = hTotalRecTrackPt->GetBinContent(iBin);
    if (binValB==0. || binValA==0.) continue;
    errorA = hTotalPureContaminationRecPt->GetBinError(iBin);
    errorB = hTotalRecTrackPt->GetBinError(iBin);
    binError = binValA/binValB * TMath::Sqrt((errorA*errorA)/(binValA*binValA) + (errorB*errorB)/(binValB*binValB)- (2*errorA*errorA)/(binValA*binValB));
    hTotalPureContaminationRecPt->SetBinContent(iBin,binValA/binValB);
    hTotalPureContaminationRecPt->SetBinError(iBin,binError);
    // std::cout << " value of contamination in " << iBin <<  "-th Bin: " << binValA/binValB <<  ",  with BinError : " << binError << std::endl;
  }
}




if (bPlotTrackContamination) {
  make3HistNice(hPurityRecPtNeg,kBlue+1);
  make3HistNice(hPurityRecEtaNeg,kBlue+1);
  make3HistNice(hPurityRecPhiNeg,kBlue+1);
  make3HistNice(hPurityRecPtPos,kRed+2);
  make3HistNice(hPurityRecEtaPos,kRed+2);
  make3HistNice(hPurityRecPhiPos,kRed+2);

  make3HistNice(hPureContaminationRecPtNeg,kBlue+1);
  make3HistNice(hPureContaminationRecEtaNeg,kBlue+1);
  make3HistNice(hPureContaminationRecPhiNeg,kBlue+1);
  make3HistNice(hPureContaminationRecPtPos,kRed+2);
  make3HistNice(hPureContaminationRecEtaPos,kRed+2);
  make3HistNice(hPureContaminationRecPhiPos,kRed+2);

  make3HistNice(hMuonContaminationRecPt, kBlue+1);
  make3HistNice(hMuonContaminationRecEta, kBlue+1);
  make3HistNice(hMuonContaminationRecPhi, kBlue+1);
  make3HistNice(hPionContaminationRecPt, kRed+2);
  make3HistNice(hPionContaminationRecEta, kRed+2);
  make3HistNice(hPionContaminationRecPhi, kRed+2);
  make3HistNice(hKaonContaminationRecPt, kGreen+2);
  make3HistNice(hKaonContaminationRecEta, kGreen+2);
  make3HistNice(hKaonContaminationRecPhi, kGreen+2);
  make3HistNice(hProtonContaminationRecPt, kOrange+2);
  make3HistNice(hProtonContaminationRecEta, kOrange+2);
  make3HistNice(hProtonContaminationRecPhi, kOrange+2);
  make3HistNice(hTotalPureContaminationRecPt, kBlack);
  make3HistNice(hTotalPureContaminationRecEta, kBlack);
  make3HistNice(hTotalPureContaminationRecPhi, kBlack);

  hTotalPureContaminationRecPt->SetMarkerStyle(24);
  hTotalPureContaminationRecPt->SetMarkerSize(0.9);
  hTotalPureContaminationRecEta->SetMarkerStyle(24);
  hTotalPureContaminationRecEta->SetMarkerSize(0.9);
  hTotalPureContaminationRecPhi->SetMarkerStyle(24);
  hTotalPureContaminationRecPhi->SetMarkerSize(0.9);




  // hPureContaminationRecPtNeg->SetMarkerStyle(24);
  // hPureContaminationRecEtaNeg->SetMarkerStyle(24);
  // hPureContaminationRecPhiNeg->SetMarkerStyle(24);
  // hPureContaminationRecPtPos->SetMarkerStyle(24);
  // hPureContaminationRecEtaPos->SetMarkerStyle(24);
  // hPureContaminationRecPhiPos->SetMarkerStyle(24);
  // hPureContaminationRecPtNeg->SetMarkerSize(0.9);
  // hPureContaminationRecEtaNeg->SetMarkerSize(0.9);
  // hPureContaminationRecPhiNeg->SetMarkerSize(0.9);
  // hPureContaminationRecPtPos->SetMarkerSize(0.9);
  // hPureContaminationRecEtaPos->SetMarkerSize(0.9);
  // hPureContaminationRecPhiPos->SetMarkerSize(0.9);
}

  makeHistNice(ptGenTrackElePos,kGreen+3);
  makeHistNice(ptGenSmearedTrackElePos,kOrange+1);
  ptGenTrackElePos->SetMarkerStyle(24);
  ptGenSmearedTrackElePos->SetMarkerStyle(26);

  makeHistNice(ptGenTrackElePos_rebin,kGreen+3);
  makeHistNice(ptGenSmearedTrackElePos_rebin,kOrange+1);
  makeHistNice(ptRecTrackElePos_rebin,kBlue+1);
  ptGenTrackElePos_rebin->SetMarkerStyle(24);
  ptGenSmearedTrackElePos_rebin->SetMarkerStyle(25);
  makeHistNice(ptGenTrackElePos,kGreen+3);
  makeHistNice(ptGenSmearedTrackElePos,kOrange+1);
  makeHistNice(ptRecTrackElePos,kBlue+1);
  ptGenTrackElePos->SetMarkerStyle(24);
  ptGenSmearedTrackElePos->SetMarkerStyle(25);
  makeHistNice(ptGenTrackElePos_beforeKineCuts_rebin,kViolet);
  ptGenTrackElePos_beforeKineCuts_rebin->SetMarkerStyle(27);
  // makeHistNice(ptGenSmearedTrackElePos_beforeKineCuts_rebin,kRed+3);
  // ptGenSmearedTrackElePos_beforeKineCuts_rebin->SetMarkerStyle(28);

  // normalizeToBinWidth(ptGenTrackElePos_rebin);
  // normalizeToBinWidth(ptGenSmearedTrackElePos_rebin);
  // normalizeToBinWidth(ptRecTrackElePos_rebin);
  //
  // normalizeToBinWidth(ptGenTrackElePos_beforeKineCuts_rebin);
  // normalizeToBinWidth(ptGenSmearedTrackElePos_beforeKineCuts_rebin);



  makeHistNice(ptRecTrackEle_rebin,kGreen+3);
  makeHistNice(etaRecTrackEle,kGreen+3);
  makeHistNice(phiRecTrackEle,kGreen+3);
  makeHistNice(ptRecTrackPos_rebin,kOrange+1);
  makeHistNice(etaRecTrackPos,kOrange+1);
  makeHistNice(phiRecTrackPos,kOrange+1);
  makeHistNice(ptGenTrackEle_rebin,kBlue+1);
  makeHistNice(etaGenTrackEle,kBlue+1);
  makeHistNice(phiGenTrackEle,kBlue+1);
  makeHistNice(ptGenTrackPos_rebin,kRed+1);
  makeHistNice(etaGenTrackPos,kRed+1);
  makeHistNice(phiGenTrackPos,kRed+1);
  ptGenTrackEle_rebin->SetMarkerStyle(24);
  etaGenTrackEle->SetMarkerStyle(24);
  phiGenTrackEle->SetMarkerStyle(24);
  ptGenTrackPos_rebin->SetMarkerStyle(24);
  etaGenTrackPos->SetMarkerStyle(24);
  phiGenTrackPos->SetMarkerStyle(24);

  // make3HistNice(ptRecTrackPrim,kBlue+1);
  // make3HistNice(etaRecTrackPrim,kBlue+1);
  // make3HistNice(phiRecTrackPrim,kBlue+1);
  // make3HistNice(ptRecTrackCC,kRed+2);
  // make3HistNice(etaRecTrackCC,kRed+2);
  // make3HistNice(phiRecTrackCC,kRed+2);
  // make3HistNice(ptRecTrackBB,kMagenta+1);
  // make3HistNice(etaRecTrackBB,kMagenta+1);
  // make3HistNice(phiRecTrackBB,kMagenta+1);
  //
  // make3HistNice(ptGenTrackPrim,kBlue+1);
  // make3HistNice(etaGenTrackPrim,kBlue+1);
  // make3HistNice(phiGenTrackPrim,kBlue+1);
  // make3HistNice(ptGenTrackCC,kRed+2);
  // make3HistNice(etaGenTrackCC,kRed+2);
  // make3HistNice(phiGenTrackCC,kRed+2);
  // make3HistNice(ptGenTrackBB,kMagenta+1);
  // make3HistNice(etaGenTrackBB,kMagenta+1);
  // make3HistNice(phiGenTrackBB,kMagenta+1);
  //
  // make3HistNice(etaGenEleTrackPrim,kBlue+1);
  // make3HistNice(etaGenEleTrackCC,kRed+2);
  // make3HistNice(etaGenEleTrackBB,kMagenta+1);
  // make3HistNice(etaGenPosTrackPrim,kBlue+1);
  // make3HistNice(etaGenPosTrackCC,kRed+2);
  // make3HistNice(etaGenPosTrackBB,kMagenta+1);
  // make3HistNice(etaRecEleTrackPrim,kBlue+1);
  // make3HistNice(etaRecEleTrackCC,kRed+2);
  // make3HistNice(etaRecEleTrackBB,kMagenta+1);
  // make3HistNice(etaRecPosTrackPrim,kBlue+1);
  // make3HistNice(etaRecPosTrackCC,kRed+2);
  // make3HistNice(etaRecPosTrackBB,kMagenta+1);
  // etaGenEleTrackPrim->SetMarkerStyle(22);
  // etaGenEleTrackCC->SetMarkerStyle(22);
  // etaGenEleTrackBB->SetMarkerStyle(22);
  // etaRecEleTrackPrim->SetMarkerStyle(22);
  // etaRecEleTrackCC->SetMarkerStyle(22);
  // etaRecEleTrackBB->SetMarkerStyle(22);
  // etaGenPosTrackPrim->SetMarkerStyle(24);
  // etaGenPosTrackCC->SetMarkerStyle(24);
  // etaGenPosTrackBB->SetMarkerStyle(24);
  // etaRecPosTrackPrim->SetMarkerStyle(24);
  // etaRecPosTrackCC->SetMarkerStyle(24);
  // etaRecPosTrackBB->SetMarkerStyle(24);
  // etaGenEleTrackPrim->SetMarkerSize(0.4);


  make3HistNice(ptRecTrackBeforeSmearing_rebin,kGreen+1);
  make3HistNice(etaRecTrackBeforeSmearing,kGreen+1);
  make3HistNice(phiRecTrackBeforeSmearing,kGreen+1);
  make3HistNice(ptRecTrackAfterSmearing_rebin,kRed+1);
  make3HistNice(etaRecTrackAfterSmearing,kRed+1);
  make3HistNice(phiRecTrackAfterSmearing,kRed+1);
  make3HistNice(ptRecTrackAfterKineCuts_rebin,kRed+4);
  make3HistNice(etaRecTrackAfterKineCuts,kRed+4);
  make3HistNice(phiRecTrackAfterKineCuts,kRed+4);
  make3HistNice(ptRecTrackElePos_rebin,kBlue+1);
  make3HistNice(etaRecTrackElePos,kBlue+1);
  make3HistNice(phiRecTrackElePos,kBlue+1);

  // normalizeToBinWidth(ptRecTrackBeforeSmearing_rebin);
  // normalizeToBinWidth(ptRecTrackAfterSmearing_rebin);
  // normalizeToBinWidth(ptRecTrackAfterKineCuts_rebin);

  // make3HistNice(ptRecTrackEtaCut_1,kOrange-5);
  // make3HistNice(ptRecTrackEtaCut_2,kOrange-4);
  // make3HistNice(ptRecTrackEtaCut_3,kOrange-3);
  // make3HistNice(ptRecTrackEtaCut_4,kOrange-2);
  // make3HistNice(ptRecTrackEtaCut_5,kOrange-1);
  // make3HistNice(ptRecTrackEtaCut_6,kOrange);
  // make3HistNice(ptRecTrackEtaCut_7,kOrange+1);
  // make3HistNice(ptRecTrackEtaCut_8,kOrange+2);
  // make3HistNice(ptRecTrackEtaCut_9,kOrange+3);
  // make3HistNice(ptRecTrackEtaCut_10,kOrange+4);
  //
  // make3HistNice(etaRecTrackEtaCut_1,kOrange-5);
  // make3HistNice(etaRecTrackEtaCut_2,kOrange-4);
  // make3HistNice(etaRecTrackEtaCut_3,kOrange-3);
  // make3HistNice(etaRecTrackEtaCut_4,kOrange-2);
  // make3HistNice(etaRecTrackEtaCut_5,kOrange-1);
  // make3HistNice(etaRecTrackEtaCut_6,kOrange);
  // make3HistNice(etaRecTrackEtaCut_7,kOrange+1);
  // make3HistNice(etaRecTrackEtaCut_8,kOrange+2);
  // make3HistNice(etaRecTrackEtaCut_9,kOrange+3);
  // make3HistNice(etaRecTrackEtaCut_10,kOrange+4);
  //
  // make3HistNice(phiRecTrackEtaCut_1,kOrange-5);
  // make3HistNice(phiRecTrackEtaCut_2,kOrange-4);
  // make3HistNice(phiRecTrackEtaCut_3,kOrange-3);
  // make3HistNice(phiRecTrackEtaCut_4,kOrange-2);
  // make3HistNice(phiRecTrackEtaCut_5,kOrange-1);
  // make3HistNice(phiRecTrackEtaCut_6,kOrange);
  // make3HistNice(phiRecTrackEtaCut_7,kOrange+1);
  // make3HistNice(phiRecTrackEtaCut_8,kOrange+2);
  // make3HistNice(phiRecTrackEtaCut_9,kOrange+3);
  // make3HistNice(phiRecTrackEtaCut_10,kOrange+4);


  // makeHistNiceTH2(hMPtDCA_ULS_primary_rec,kBlue+1);

  if (bPlotEfficiency) {
    makeHistNice(ptEffEle, kBlue+1);
    makeHistNice(etaEffEle, kBlue+1);
    makeHistNice(phiEffEle, kBlue+1);
    makeHistNice(ptEffPos, kRed+2);
    makeHistNice(etaEffPos, kRed+2);
    makeHistNice(phiEffPos, kRed+2);
    // makeHistNice(ptEffElePrim,kBlue+1);
    // makeHistNice(ptEffEleCC,kRed+2);
    // makeHistNice(ptEffEleBB,kMagenta+1);
    // makeHistNice(etaEffElePrim,kBlue+1);
    // makeHistNice(etaEffEleCC,kRed+2);
    // makeHistNice(etaEffEleBB,kMagenta+1);
    // makeHistNice(phiEffElePrim,kBlue+1);
    // makeHistNice(phiEffEleCC,kRed+2);
    // makeHistNice(phiEffEleBB,kMagenta+1);
    // makeHistNice(ptEffPosPrim,kBlue+1);
    // makeHistNice(ptEffPosCC,kRed+2);
    // makeHistNice(ptEffPosBB,kMagenta+1);
    // makeHistNice(etaEffPosPrim,kBlue+1);
    // makeHistNice(etaEffPosCC,kRed+2);
    // makeHistNice(etaEffPosBB,kMagenta+1);
    // makeHistNice(phiEffPosPrim,kBlue+1);
    // makeHistNice(phiEffPosCC,kRed+2);
    // makeHistNice(phiEffPosBB,kMagenta+1);

    makeHistNice(ptEffElePosGen,kBlue+1);
    makeHistNice(ptEffElePosGenSmeared,kRed+2);
    makeHistNice(ptEffElePosGen_woPID,kViolet+2);
  }

  makeHistNice(proj_recULS_Mee,kBlack);
  makeHistNice(proj_recULS_Ptee,kBlack);
  makeHistNice(proj_recULS_DCA,kBlack);
  // makeHistNice(proj_recULS_MeePrim,kBlue+1);
  // makeHistNice(proj_recULS_PteePrim,kBlue+1);
  // makeHistNice(proj_recULS_MeeCC,kRed+2);
  // makeHistNice(proj_recULS_PteeCC,kRed+2);
  // makeHistNice(proj_recULS_MeeBB,kMagenta+1);
  // makeHistNice(proj_recULS_PteeBB,kMagenta+1);

  makeHistNice(proj_recLS_Mee,kRed+2);
  makeHistNice(proj_recLS_Ptee,kRed+2);
  makeHistNice(proj_recLS_DCA,kRed+2);


  // makeHistNice(proj_genLS_MeePrim,kBlue+1);
  // makeHistNice(proj_genLS_PteePrim,kBlue+1);
  // makeHistNice(proj_genLS_MeeCC,kRed+2);
  // makeHistNice(proj_genLS_PteeCC,kRed+2);
  // makeHistNice(proj_genLS_MeeBB,kMagenta+1);
  // makeHistNice(proj_genLS_PteeBB,kMagenta+1);


  double legPosTrack[4] = {0.55,0.78,0.95,0.93};
  double legPosTrack2[4] = {0.55,0.15,0.95,0.3};
  // double legPosKine[4] = {0.3,0.18,0.85,0.45};
  double legPosKine[4] = {0.3,0.78,0.85,0.93};
  double legPosPair[4] = {0.55,0.78,0.95,0.93};
  double legPosCont[4] = {0.2,0.68,0.55,0.93};
  double legPosBottomLeft[4] = {0.15,0.15,0.55,0.5};
  //make some legends
  // auto legTrack1 = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  // legTrack1->SetBorderSize(0);
  // legTrack1->SetFillStyle(0);
  // legTrack1->AddEntry(ptEffElePrim,"LF #rightarrow e","p");
  // legTrack1->AddEntry(ptEffEleCC,"charm #rightarrow e","p");
  // legTrack1->AddEntry(ptEffEleBB,"beauty #rightarrow e","p");
  //
  // auto legLFCCBB = new TLegend(legPosTrack2[0]+0.05,legPosTrack2[1],legPosTrack2[2]+0.05,legPosTrack2[3]);
  // legLFCCBB->SetBorderSize(0);
  // legLFCCBB->SetFillStyle(0);
  // legLFCCBB->AddEntry(ptEffElePrim,"LF #rightarrow e","p");
  // legLFCCBB->AddEntry(ptEffEleCC,"charm #rightarrow e","p");
  // legLFCCBB->AddEntry(ptEffEleBB,"beauty #rightarrow e","p");

  auto legTrackPosEle_top = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legTrackPosEle_top->SetBorderSize(0);
  legTrackPosEle_top->SetFillStyle(0);
  legTrackPosEle_top->SetTextSize(0.04);
  legTrackPosEle_top->AddEntry(ptEffEle,"electrons","p");
  legTrackPosEle_top->AddEntry(ptEffPos,"positrons","p");

  auto legEff_GenGenSmear = new TLegend(legPosTrack[0]+0.05,legPosTrack[1]+0.05,legPosTrack[2]+0.05,legPosTrack[3]);
  legEff_GenGenSmear->SetBorderSize(0);
  legEff_GenGenSmear->SetFillStyle(0);
  legEff_GenGenSmear->SetTextSize(0.04);
  legEff_GenGenSmear->AddEntry(ptEffElePosGen,"rec/gen","p");
  legEff_GenGenSmear->AddEntry(ptEffElePosGenSmeared,"rec/gen smeared","p");

  auto legEff_Gen = new TLegend(legPosTrack[0]+0.15,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legEff_Gen->SetBorderSize(0);
  legEff_Gen->SetFillStyle(0);
  legEff_Gen->SetTextSize(0.04);
  legEff_Gen->AddEntry(ptEffElePosGen,"rec/gen","p");

  auto legEff_Gen_woPID = new TLegend(legPosTrack[0]+0.15,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legEff_Gen_woPID->SetBorderSize(0);
  legEff_Gen_woPID->SetFillStyle(0);
  legEff_Gen_woPID->SetTextSize(0.04);
  legEff_Gen_woPID->AddEntry(ptEffElePosGen_woPID,"rec/gen wo PID","p");


  auto legTrackEle_top = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legTrackEle_top->SetBorderSize(0);
  legTrackEle_top->SetFillStyle(0);
  legTrackEle_top->SetTextSize(0.04);
  legTrackEle_top->AddEntry(ptEffEle,"electrons","p");

  auto legTrackPos_top = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legTrackPos_top->SetBorderSize(0);
  legTrackPos_top->SetFillStyle(0);
  legTrackPos_top->SetTextSize(0.04);
  legTrackPos_top->AddEntry(ptEffPos,"positrons","p");

  auto legGenRecPos_top = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legGenRecPos_top->SetBorderSize(0);
  legGenRecPos_top->SetFillStyle(0);
  legGenRecPos_top->SetTextSize(0.04);
  legGenRecPos_top->AddEntry(ptGenTrackPos,"gen positrons","p");
  legGenRecPos_top->AddEntry(ptRecTrackPos_rebin,"rec positrons","p");

  auto legGenRecEle_top = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legGenRecEle_top->SetBorderSize(0);
  legGenRecEle_top->SetFillStyle(0);
  legGenRecEle_top->SetTextSize(0.04);
  legGenRecEle_top->AddEntry(ptGenTrackEle_rebin,"gen electrons","p");
  legGenRecEle_top->AddEntry(ptRecTrackEle_rebin,"rec electrons","p");

  auto legGenGenSmearRecElePos_top = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legGenGenSmearRecElePos_top->SetBorderSize(0);
  legGenGenSmearRecElePos_top->SetFillStyle(0);
  legGenGenSmearRecElePos_top->AddEntry(ptGenTrackElePos,"gen tracks","p");
  legGenGenSmearRecElePos_top->AddEntry(ptGenSmearedTrackElePos,"gen smeared tracks","p");
  legGenGenSmearRecElePos_top->AddEntry(ptRecTrackElePos,"rec tracks","p");

  auto legGenGenSmearRecElePos_Plus = new TLegend(legPosBottomLeft[0],legPosBottomLeft[1],legPosBottomLeft[2],legPosBottomLeft[3]);
  legGenGenSmearRecElePos_Plus->SetBorderSize(0);
  legGenGenSmearRecElePos_Plus->SetFillStyle(0);
  legGenGenSmearRecElePos_Plus->AddEntry(ptGenTrackElePos_beforeKineCuts_rebin,"gen","p");
  // legGenGenSmearRecElePos_Plus->AddEntry(ptGenSmearedTrackElePos_beforeKineCuts_rebin,"gen smeared","p");
  legGenGenSmearRecElePos_Plus->AddEntry(ptGenTrackElePos_rebin,"gen + kine cuts","p");
  legGenGenSmearRecElePos_Plus->AddEntry(ptGenSmearedTrackElePos_rebin,"gen smeared + kine cuts","p");
  legGenGenSmearRecElePos_Plus->AddEntry(ptRecTrackBeforeSmearing_rebin,"rec","p");
  legGenGenSmearRecElePos_Plus->AddEntry(ptRecTrackAfterSmearing_rebin,"rec + smeared","p");
  legGenGenSmearRecElePos_Plus->AddEntry(ptRecTrackAfterKineCuts_rebin,"rec + smeared + kine","p");
  legGenGenSmearRecElePos_Plus->AddEntry(ptRecTrackElePos_rebin,"rec + smeared + kine + PID","p");

  // auto legEffULSPair = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  // legEffULSPair->SetBorderSize(0);
  // legEffULSPair->SetFillStyle(0);
  // legEffULSPair->AddEntry(ptPairEffULS,"ULS pair efficiency","p");

  TLegend* legContamination;
  legContamination = new TLegend(0.48,0.73,0.83,0.95);
  legContamination->SetBorderSize(0);
  legContamination->SetFillStyle(0);
  legContamination->SetTextSize(0.04);
  // legContamination->AddEntry(hPurityRecPtNeg,"neg purity","p");
  // legContamination->AddEntry(hPurityRecPtPos,"pos purity","p");
  legContamination->AddEntry(hPureContaminationRecPtNeg,"neg","p");
  legContamination->AddEntry(hPureContaminationRecPtPos,"pos","p");

  TLegend* legPIDContamination;
  // legPIDContamination = new TLegend(legPosCont[0]-0.05,legPosCont[1],legPosCont[2],legPosCont[3]);
  // legPIDContamination = new TLegend(0.48,0.73,0.83,0.95);
  legPIDContamination = new TLegend(0.7,0.76,0.9,0.93);
  legPIDContamination->SetBorderSize(0);
  legPIDContamination->SetFillStyle(0);
  legPIDContamination->SetTextSize(25);
  legPIDContamination->SetTextFont(43);
  legPIDContamination->AddEntry(hTotalPureContaminationRecPt,"Total","p");
  legPIDContamination->AddEntry(hMuonContaminationRecPt,"Muon","p");
  legPIDContamination->AddEntry(hPionContaminationRecPt,"Pion","p");
  // legPIDContamination->AddEntry(hKaonContaminationRecPt,"Kaon","p");
  // legPIDContamination->AddEntry(hProtonContaminationRecPt,"Proton","p");

  TLegend* legRejFactor;
  legRejFactor = new TLegend(0.7,0.8,0.9,0.9);
  legRejFactor->SetBorderSize(0);
  legRejFactor->SetFillStyle(0);
  // legRejFactor->SetTextSize(0.025);
  legRejFactor->SetTextSize(25);
  legRejFactor->SetTextFont(43);
  legRejFactor->AddEntry(hMuonRejectionFactorPt,"Muon","p");
  legRejFactor->AddEntry(hPionRejectionFactorPt,"Pion","p");
  // legRejFactor->AddEntry(hKaonRejectionFactorPt,"Kaon","p");
  // legRejFactor->AddEntry(hProtonRejectionFactorPt,"Proton","p");
  // legRejFactor->AddEntry(hTotalRejectionFactorPt,"Total","p");

  TLegend* legRejFactor_mu;
  legRejFactor_mu = new TLegend(legPosCont[0]-0.05,legPosCont[1],legPosCont[2],legPosCont[3]);
  legRejFactor_mu->SetBorderSize(0);
  legRejFactor_mu->SetFillStyle(0);
  legRejFactor_mu->SetTextSize(0.025);
  legRejFactor_mu->AddEntry(hMuonRejectionFactorPt,"Muon","p");

  TLegend* legRejFactor_pi;
  legRejFactor_pi = new TLegend(0.45,0.85,0.9,0.95);
  legRejFactor_pi->SetBorderSize(0);
  legRejFactor_pi->SetFillStyle(0);
  legRejFactor_pi->SetTextSize(0.025);
  legRejFactor_pi->AddEntry(hPionRejectionFactorPt,"Pions","p");


  TLegend* legTotalContamination;
  legTotalContamination = new TLegend(legPosCont[0]-0.05,legPosCont[1],legPosCont[2],legPosCont[3]);
  legTotalContamination->SetBorderSize(0);
  legTotalContamination->SetFillStyle(0);
  legTotalContamination->SetTextSize(0.04);
  legTotalContamination->AddEntry(hTotalPureContaminationRecPt,"Total contamination","p");

  auto legPIDSeparateColorTOF = new TLegend(0.2, 0.16, 0.35, 0.3);
  legPIDSeparateColorTOF->SetBorderSize(0);
  legPIDSeparateColorTOF->SetFillStyle(0);
  legPIDSeparateColorTOF->SetTextSize(25);
  legPIDSeparateColorTOF->SetTextFont(43);
  TLegendEntry *entryInfo00=legPIDSeparateColorTOF->AddEntry(hNsigmaP_RICH_trueElec[0],"Electrons","");
  TLegendEntry *entryInfo01=legPIDSeparateColorTOF->AddEntry(hNsigmaP_RICH_trueMuon[0],"Muons","");
  TLegendEntry *entryInfo02=legPIDSeparateColorTOF->AddEntry(hNsigmaP_RICH_truePion[0],"Pions","");
  entryInfo00->SetTextColor(kBlue);
  entryInfo01->SetTextColor(kBlack);
  entryInfo02->SetTextColor(kRed);
  // legPIDSeparateColorTOF->AddEntry(hNsigmaP_RICH_trueElec[0],"Electrons","l");
  // legPIDSeparateColorTOF->AddEntry(hNsigmaP_RICH_trueMuon[0],"Muon","l");
  // legPIDSeparateColorTOF->AddEntry(hNsigmaP_RICH_truePion[0],"Pion","l");

  auto legPIDSeparateColor = new TLegend(0.7, 0.74, 0.9, 0.9);
  legPIDSeparateColor->SetBorderSize(0);
  legPIDSeparateColor->SetFillStyle(0);
  legPIDSeparateColor->SetTextSize(25);
  legPIDSeparateColor->SetTextFont(43);
  TLegendEntry *entryInfo03=legPIDSeparateColor->AddEntry(hNsigmaP_RICH_trueElec[0],"Electrons","");
  TLegendEntry *entryInfo04=legPIDSeparateColor->AddEntry(hNsigmaP_RICH_trueMuon[0],"Muons","");
  TLegendEntry *entryInfo05=legPIDSeparateColor->AddEntry(hNsigmaP_RICH_truePion[0],"Pions","");
  entryInfo03->SetTextColor(kBlue);
  entryInfo04->SetTextColor(kBlack);
  entryInfo05->SetTextColor(kRed);
  // legPIDSeparateColor->AddEntry(hNsigmaP_RICH_trueElec[0],"Electrons","l");
  // legPIDSeparateColor->AddEntry(hNsigmaP_RICH_trueMuon[0],"Muon","l");
  // legPIDSeparateColor->AddEntry(hNsigmaP_RICH_truePion[0],"Pion","l");

  auto legSmearLabel = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legSmearLabel->SetBorderSize(0);
  legSmearLabel->SetFillStyle(0);
  legSmearLabel->AddEntry(ptRecTrackBeforeSmearing_rebin,"before Smearing","p");
  legSmearLabel->AddEntry(ptRecTrackAfterSmearing_rebin,"after Smearing","p");
  legSmearLabel->AddEntry(ptRecTrackAfterKineCuts_rebin,"after kine Cuts","p");
  legSmearLabel->AddEntry(ptRecTrackElePos_rebin,"after PID cuts","p");
  // legSmearLabel->AddEntry(etaRecTrackEtaCut_1," |#eta| < 1.0","p");
  // legSmearLabel->AddEntry(etaRecTrackEtaCut_2," |#eta| < 2.0","p");
  // legSmearLabel->AddEntry(etaRecTrackEtaCut_4," |#eta| < 4.0","p");
  // legSmearLabel->AddEntry(etaRecTrackEtaCut_5," |#eta| < 5.0","p");
  // legSmearLabel->AddEntry(etaRecTrackEtaCut_6," |#eta| < 6.0","p");
  // legSmearLabel->AddEntry(etaRecTrackEtaCut_8," |#eta| < 8.0","p");
  // legSmearLabel->AddEntry(etaRecTrackEtaCut_10," |#eta| < 10.0","p");



  // //make some legends
  // auto legProfile = new TLegend(legPosTrack[0],legPosTrack[1],legPosTrack[2],legPosTrack[3]);
  // legProfile->SetBorderSize(0);
  // legProfile->SetFillStyle(0);
  // legProfile->AddEntry("LF","LF #rightarrow e","fp");
  // legProfile->AddEntry("CC","charm #rightarrow e","fp");
  // legProfile->AddEntry("BB","beauty #rightarrow e","fp");

  // auto legKine = new TLegend(legPosKine[0],legPosKine[1],legPosKine[2],legPosKine[3]);
  // legKine->SetBorderSize(0);
  // legKine->SetFillStyle(0);
  // legKine->AddEntry(phiRecTrackPrim,"LF #rightarrow e","l");
  // legKine->AddEntry(phiRecTrackCC,"charm #rightarrow e","l");
  // legKine->AddEntry(phiRecTrackBB,"beauty #rightarrow e","l");

  // auto legKineElePos = new TLegend(legPosKine[0],legPosKine[1],legPosKine[2],legPosKine[3]);
  // legKineElePos->SetBorderSize(0);
  // legKineElePos->SetFillStyle(0);
  // legKineElePos->AddEntry(phiRecTrackPrim,"LF #rightarrow e","p");
  // legKineElePos->AddEntry(phiRecTrackCC,"charm #rightarrow e","p");
  // legKineElePos->AddEntry(phiRecTrackBB,"beauty #rightarrow e","p");
  // // legKineElePos->AddEntry(etaRecEleTrackPrim,"e^{-}","p");
  // // legKineElePos->AddEntry(etaRecPosTrackPrim,"e^{+}","p");

  // auto legPair = new TLegend(legPosPair[0],legPosPair[1],legPosPair[2],legPosPair[3]);
  // legPair->SetBorderSize(0);
  // legPair->SetFillStyle(0);
  // legPair->AddEntry(ptRecTrackPrim,"LF #rightarrow ee","pe");
  // legPair->AddEntry(ptRecTrackCC,"c#bar{c} #rightarrow ee","pe");
  // legPair->AddEntry(ptRecTrackBB,"b#bar{b} #rightarrow ee","pe");


  // auto cKine = new TCanvas("cKine","cKine",900,400);
  // cKine->Divide(3,1);
  // for (size_t i = 1; i < 4; i++) {
  //   /* code */
  //   cKine->cd(i)->SetTopMargin(0.03);
  //   cKine->cd(i)->SetRightMargin(0.03);
  // }
  // cKine->cd(1)->SetLogy();
  // ptRecTrackPrim->GetXaxis()->SetRangeUser(0.,8);
  // ptRecTrackPrim->Draw("pe1");
  // ptRecTrackCC->Draw("pe1 same");
  // ptRecTrackBB->Draw("pe1 same");
  // legKineElePos->Draw("same");
  // cKine->cd(2);
  // // etaRecTrackPrim->SetMaximum(0.03);
  // etaRecTrackPrim->GetXaxis()->SetRangeUser(-2.0,2.0);
  // etaRecTrackPrim->Draw("axis");
  // etaRecTrackCC->Draw("pe1 same");
  // // etaRecEleTrackCC->Draw("pe1 same");
  // // etaRecPosTrackCC->Draw("pe1 same");
  // etaRecTrackBB->Draw("pe1 same");
  // // etaRecEleTrackBB->Draw("pe1 same");
  // // etaRecPosTrackBB->Draw("pe1 same");
  // etaRecTrackPrim->Draw("pe1 same");
  // // etaRecEleTrackPrim->Draw("pe1 same");
  // // etaRecPosTrackPrim->Draw("pe1 same");
  // cKine->cd(3);
  // phiRecTrackPrim->SetMinimum(0.);
  // // phiRecTrackPrim->SetMaximum(0.016);
  // phiRecTrackPrim->Draw("axis");
  // phiRecTrackCC->Draw("pe1 same");
  // phiRecTrackBB->Draw("pe1 same");
  // phiRecTrackPrim->Draw("pe1 same");
  // cKine->SaveAs("./plots/RecPtEtaPhi.png");


  // auto cKineGen = new TCanvas("cKineGen","cKineGen",900,400);
  //
  // cKineGen->Divide(3,1);
  // for (size_t i = 1; i < 4; i++) {
  //   /* code */
  //   cKineGen->cd(i)->SetTopMargin(0.03);
  //   cKineGen->cd(i)->SetRightMargin(0.03);
  // }
  // cKineGen->cd(1)->SetLogy();
  // ptGenTrackPrim->GetXaxis()->SetRangeUser(0.,8);
  // ptGenTrackPrim->Draw("pe1");
  // ptGenTrackCC->Draw("pe1 same");
  // ptGenTrackBB->Draw("pe1 same");
  // legKineElePos->Draw("same");
  // cKineGen->cd(2);
  // // etaGenTrackPrim->SetMaximum(0.03);
  // etaGenTrackPrim->GetXaxis()->SetRangeUser(-2.0,2.0);
  // etaGenTrackPrim->Draw("axis");
  // etaGenTrackCC->Draw("pe1 same");
  // // etaGenEleTrackCC->Draw("pe1 same");
  // // etaGenPosTrackCC->Draw("pe1 same");
  // etaGenTrackBB->Draw("pe1 same");
  // // etaGenEleTrackBB->Draw("pe1 same");
  // // etaGenPosTrackBB->Draw("pe1 same");
  // etaGenTrackPrim->Draw("pe1 same");
  // // etaGenEleTrackPrim->Draw("pe1 same");
  // // etaGenPosTrackPrim->Draw("pe1 same");
  // cKineGen->cd(3);
  // phiGenTrackPrim->SetMinimum(0.);
  // // phiGenTrackPrim->SetMaximum(0.016);
  // phiGenTrackPrim->Draw("axis");
  // phiGenTrackCC->Draw("pe1 same");
  // phiGenTrackBB->Draw("pe1 same");
  // phiGenTrackPrim->Draw("pe1 same");
  // cKineGen->SaveAs("./plots/GenPtEtaPhi.png");

  if (bPlotEfficiency) {
    auto cEffElePt = new TCanvas("cEffElePt","cEffElePt",800,800);
    // cEffElePt->SetLogy();
    cEffElePt->SetTopMargin(0.03);
    cEffElePt->SetRightMargin(0.03);
    cEffElePt->SetLeftMargin(0.13);
    ptEffEle->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffEle->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptEffEle->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEffEle->Draw("hist p e1 ");
    // // ptEffElePrim->SetMaximum(1.);
    // ptEffElePrim->Draw("hist p ");
    // ptEffEleCC->Draw("hist p same");
    // ptEffEleBB->Draw("hist p same");
    legTrackEle_top->Draw("same");
    cEffElePt->SaveAs("./plots/EffElePt.png");

    auto cEffElePosPt = new TCanvas("cEffElePosPt","cEffElePosPt",800,800);
    // cEffElePosPt->SetLogy();
    cEffElePosPt->SetTopMargin(0.03);
    cEffElePosPt->SetRightMargin(0.03);
    cEffElePosPt->SetLeftMargin(0.13);
    ptEffEle->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffEle->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptEffEle->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEffEle->SetMaximum(1.1*std::max({ptEffEle->GetMaximum(),ptEffPos->GetMaximum()}));
    ptEffEle->SetMinimum(0.);
    ptEffEle->Draw("hist p e1 ");
    ptEffPos->Draw("hist p e1 same");
    legTrackPosEle_top->Draw("same");
    textBField->Draw("same");
    textPIDScenario->Draw("same");
    cEffElePosPt->SaveAs("./plots/Eff_ElePos_Pt.png");

    auto cEffAddedElePosPt = new TCanvas("cEffAddedElePosPt","cEffAddedElePosPt",800,800);
    // cEffAddedElePosPt->SetLogy();
    cEffAddedElePosPt->SetTopMargin(0.03);
    cEffAddedElePosPt->SetRightMargin(0.03);
    cEffAddedElePosPt->SetLeftMargin(0.13);
    ptEffElePosGen->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffElePosGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptEffElePosGen->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEffElePosGen->SetMaximum(1.1*std::max({ptEffElePosGen->GetMaximum(),ptEffElePosGenSmeared->GetMaximum()}));
    // ptEffElePosGen->SetMaximum(1.1*std::max({ptEffElePosGen->GetMaximum()}));
    // ptEffElePosGenSmeared->SetMaximum(1.1*std::max({ptEffElePosGenSmeared->GetMaximum()}));
    ptEffElePosGen->SetMinimum(0.);
    ptEffElePosGen->Draw("hist p e1");
    ptEffElePosGenSmeared->Draw("hist p e1 same");
    legEff_GenGenSmear->Draw("same");
    textBField->Draw("same");
    textPIDScenario->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_GenGenSmear.png");
    ptEffElePosGen->Draw("hist p e1 ");
    legEff_Gen->Draw("same");
    textBField->Draw("same");
    textPIDScenario->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_Gen.png");
    ptEffElePosGen_woPID->Draw("hist p e1 ");
    legEff_Gen_woPID->Draw("same");
    textBField->Draw("same");
    textPIDScenario->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_Gen_woPID.png");


    auto cEffEleEta = new TCanvas("cEffEleEta","cEffEleEta",800,800);
    // cEffEleEta->SetLogy();
    cEffEleEta->SetTopMargin(0.03);
    cEffEleEta->SetRightMargin(0.03);
    cEffEleEta->SetLeftMargin(0.13);
    etaEffEle->GetYaxis()->SetTitle("#eta^{rec}/#eta^{gen}");
    etaEffEle->GetXaxis()->SetTitle("#eta");
    etaEffEle->GetXaxis()->SetRangeUser(-2.0,2.0);
    etaEffEle->Draw("hist p e1 ");
    // // etaEffElePrim->SetMaximum(1.);
    // etaEffElePrim->Draw("hist p");
    // etaEffEleCC->Draw("hist p same");
    // etaEffEleBB->Draw("hist p same");
    legTrackEle_top->Draw("same");
    cEffEleEta->SaveAs("./plots/EffEleEta.png");

    auto cEffElePhi = new TCanvas("cEffElePhi","cEffElePhi",800,800);
    // cEffElePhi->SetLogy();
    cEffElePhi->SetTopMargin(0.03);
    cEffElePhi->SetRightMargin(0.03);
    cEffElePhi->SetLeftMargin(0.13);
    phiEffEle->GetYaxis()->SetTitle("#varphi^{rec}/#varphi^{gen}");
    phiEffEle->GetXaxis()->SetTitle("#varphi (rad)");
    phiEffEle->Draw("hist p e1 ");
    // // phiEffElePrim->SetMaximum(1.);
    // phiEffElePrim->GetXaxis()->SetRangeUser(-7.0,7.0);
    // phiEffElePrim->Draw("hist p ");
    // phiEffEleCC->Draw("hist p same");
    // phiEffEleBB->Draw("hist p same");
    legTrackEle_top->Draw("same");
    cEffElePhi->SaveAs("./plots/EffElePhi.png");


    auto cEffPosElePtEtaPhi = new TCanvas("cEffPosElePtEtaPhi","cEffPosElePtEtaPhi",1350,600);
    cEffPosElePtEtaPhi->Divide(3,1);
    for (size_t i = 1; i < 4; i++) {
      cEffPosElePtEtaPhi->cd(i)->SetTopMargin(0.03);
      cEffPosElePtEtaPhi->cd(i)->SetRightMargin(0.03);
    }
    cEffPosElePtEtaPhi->cd(1);
    ptEffEle->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffPos->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffEle->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptEffEle->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEffPos->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEffEle->SetMaximum(1.1*std::max({ptEffEle->GetMaximum(),ptEffPos->GetMaximum()}));
    ptEffEle->SetMinimum(0.);
    ptEffEle->Draw("hist p e1");
    ptEffPos->Draw("same hist p e1");
    legTrackPosEle_top->Draw("same");
    cEffPosElePtEtaPhi->cd(2);
    etaEffEle->GetYaxis()->SetTitle("#eta^{rec}/#eta^{gen}");
    etaEffEle->GetXaxis()->SetTitle("#eta");
    etaEffEle->GetXaxis()->SetRangeUser(-2.0,2.0);
    etaEffEle->SetMaximum(1.1*std::max({etaEffEle->GetMaximum(),etaEffPos->GetMaximum()}));
    etaEffEle->SetMinimum(0.);
    etaEffEle->Draw("hist p e1");
    etaEffPos->Draw("same hist p e1");
    cEffPosElePtEtaPhi->cd(3);
    phiEffEle->GetYaxis()->SetTitle("#varphi^{rec}/#varphi^{gen}");
    phiEffEle->GetXaxis()->SetTitle("#varphi (rad)");
    phiEffEle->SetMaximum(1.1*std::max({phiEffEle->GetMaximum(),phiEffPos->GetMaximum()}));
    phiEffEle->SetMinimum(0.);
    phiEffEle->GetXaxis()->SetRangeUser(-7.0,7.0);
    phiEffEle->Draw("hist p e1");
    phiEffPos->Draw("same hist p e1");
    cEffPosElePtEtaPhi->SaveAs("./plots/Eff_ElePos_PtEtaPhi.png");


    // auto cEffElePtEtaPhi = new TCanvas("cEffElePtEtaPhi","cEffElePtEtaPhi",900,400);
    // cEffElePtEtaPhi->Divide(3,1);
    // for (size_t i = 1; i < 4; i++) {
    //   cEffElePtEtaPhi->cd(i)->SetTopMargin(0.03);
    //   cEffElePtEtaPhi->cd(i)->SetRightMargin(0.03);
    // }
    // cEffElePtEtaPhi->cd(1);
    // ptEffElePrim->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    // ptEffElePrim->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // ptEffElePrim->GetXaxis()->SetRangeUser(0.0,4.0);
    // ptEffElePrim->SetMaximum(1.1*std::max({ptEffElePrim->GetMaximum(),ptEffEleCC->GetMaximum(),ptEffEleBB->GetMaximum()}));
    // ptEffElePrim->SetMinimum(0.);
    // ptEffElePrim->Draw("hist p ");
    // ptEffEleCC->Draw("hist p same");
    // ptEffEleBB->Draw("hist p same");
    // legTrack1->Draw("same");
    // cEffElePtEtaPhi->cd(2);
    // etaEffElePrim->GetYaxis()->SetTitle("#eta^{rec}/#eta^{gen}");
    // etaEffElePrim->GetXaxis()->SetTitle("#eta");
    // etaEffElePrim->GetXaxis()->SetRangeUser(-2.0,2.0);
    // etaEffElePrim->SetMaximum(1.1*std::max({etaEffElePrim->GetMaximum(),etaEffEleCC->GetMaximum(),etaEffEleBB->GetMaximum()}));
    // etaEffElePrim->SetMinimum(0.);
    // etaEffElePrim->Draw("hist p");
    // etaEffEleCC->Draw("hist p same");
    // etaEffEleBB->Draw("hist p same");
    // cEffElePtEtaPhi->cd(3);
    // phiEffElePrim->GetYaxis()->SetTitle("#varphi^{rec}/#varphi^{gen}");
    // phiEffElePrim->GetXaxis()->SetTitle("#varphi (rad)");
    // phiEffElePrim->SetMaximum(1.1*std::max({phiEffElePrim->GetMaximum(),phiEffEleCC->GetMaximum(),phiEffEleBB->GetMaximum()}));
    // phiEffElePrim->SetMinimum(0.);
    // phiEffElePrim->GetXaxis()->SetRangeUser(-7.0,7.0);
    // phiEffElePrim->Draw("hist p ");
    // phiEffEleCC->Draw("hist p same");
    // phiEffEleBB->Draw("hist p same");
    // cEffElePtEtaPhi->SaveAs("./plots/EffElePtEtaPhi.png");


    // auto cEffPosPtEtaPhi = new TCanvas("cEffPosPtEtaPhi","cEffPosPtEtaPhi",900,400);
    // cEffPosPtEtaPhi->Divide(3,1);
    // for (size_t i = 1; i < 4; i++) {
    //   cEffPosPtEtaPhi->cd(i)->SetTopMargin(0.03);
    //   cEffPosPtEtaPhi->cd(i)->SetRightMargin(0.03);
    // }
    // cEffPosPtEtaPhi->cd(1);
    // ptEffPosPrim->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    // ptEffPosPrim->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // ptEffPosPrim->GetXaxis()->SetRangeUser(0.0,4.0);
    // ptEffPosPrim->SetMaximum(1.1*std::max({ptEffPosPrim->GetMaximum(),ptEffPosCC->GetMaximum(),ptEffPosBB->GetMaximum()}));
    // ptEffPosPrim->SetMinimum(0.);
    // ptEffPosPrim->Draw("hist p ");
    // ptEffPosCC->Draw("hist p same");
    // ptEffPosBB->Draw("hist p same");
    // legTrack1->Draw("same");
    // cEffPosPtEtaPhi->cd(2);
    // etaEffPosPrim->GetYaxis()->SetTitle("#eta^{rec}/#eta^{gen}");
    // etaEffPosPrim->GetXaxis()->SetTitle("#eta");
    // etaEffPosPrim->GetXaxis()->SetRangeUser(-2.0,2.0);
    // etaEffPosPrim->SetMaximum(1.1*std::max({etaEffPosPrim->GetMaximum(),etaEffPosCC->GetMaximum(),etaEffPosBB->GetMaximum()}));
    // etaEffPosPrim->SetMinimum(0.);
    // etaEffPosPrim->Draw("hist p");
    // etaEffPosCC->Draw("hist p same");
    // etaEffPosBB->Draw("hist p same");
    // cEffPosPtEtaPhi->cd(3);
    // phiEffPosPrim->GetYaxis()->SetTitle("#varphi^{rec}/#varphi^{gen}");
    // phiEffPosPrim->GetXaxis()->SetTitle("#varphi (rad)");
    // phiEffPosPrim->SetMaximum(1.1*std::max({phiEffPosPrim->GetMaximum(),phiEffPosCC->GetMaximum(),phiEffPosBB->GetMaximum()}));
    // phiEffPosPrim->SetMinimum(0.);
    // phiEffPosPrim->GetXaxis()->SetRangeUser(-7.0,7.0);
    // phiEffPosPrim->Draw("hist p ");
    // phiEffPosCC->Draw("hist p same");
    // phiEffPosBB->Draw("hist p same");
    // cEffPosPtEtaPhi->SaveAs("./plots/EffPosPtEtaPhi.png");



    auto cEffPosPt = new TCanvas("cEffPosPt","cEffPosPt",800,800);
    // cEffPosPt->SetLogy();
    cEffPosPt->SetTopMargin(0.03);
    cEffPosPt->SetRightMargin(0.03);
    cEffPosPt->SetLeftMargin(0.13);
    ptEffPos->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffPos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptEffPos->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEffPos->Draw("hist p e1");
    // ptEffPosPrim->SetMaximum(1.);
    // ptEffPosPrim->Draw("hist p ");
    // ptEffPosCC->Draw("hist p same");
    // ptEffPosBB->Draw("hist p same");
    legTrackPos_top->Draw("same");
    cEffPosPt->SaveAs("./plots/EffPosPt.png");

    auto cEffPosEta = new TCanvas("cEffPosEta","cEffPosEta",800,800);
    // cEffPosEta->SetLogy();
    cEffPosEta->SetTopMargin(0.03);
    cEffPosEta->SetRightMargin(0.03);
    cEffPosEta->SetLeftMargin(0.13);
    etaEffPos->GetYaxis()->SetTitle("#eta^{rec}/#eta^{gen}");
    etaEffPos->GetXaxis()->SetTitle("#eta");
    etaEffPos->GetXaxis()->SetRangeUser(-2.0,2.0);
    etaEffPos->Draw("hist p e1");
    // etaEffPosPrim->SetMaximum(1.);
    // etaEffPosPrim->Draw("hist p ");
    // etaEffPosCC->Draw("hist p same");
    // etaEffPosBB->Draw("hist p same");
    legTrackPos_top->Draw("same");
    cEffPosEta->SaveAs("./plots/EffPosEta.png");

    auto cEffPosPhi = new TCanvas("cEffPosPhi","cEffPosPhi",800,800);
    // cEffPosPhi->SetLogy();
    cEffPosPhi->SetTopMargin(0.03);
    cEffPosPhi->SetRightMargin(0.03);
    cEffPosPhi->SetLeftMargin(0.13);
    phiEffPos->GetYaxis()->SetTitle("#varphi^{rec}/#varphi^{gen}");
    phiEffPos->GetXaxis()->SetTitle("#varphi (rad)");
    phiEffPos->Draw("hist p e1");
    // phiEffPosPrim->SetMaximum(1.);
    // phiEffPosPrim->GetXaxis()->SetRangeUser(-7.0,7.0);
    // phiEffPosPrim->Draw("hist p ");
    // phiEffPosCC->Draw("hist p same");
    // phiEffPosBB->Draw("hist p same");
    legTrackPos_top->Draw("same");
    cEffPosPhi->SaveAs("./plots/EffPosPhi.png");

    auto cGenGenSmearRecElePosPt = new TCanvas("cGenGenSmearRecElePosPt","cGenGenSmearRecElePosPt",800,800);
    cGenGenSmearRecElePosPt->SetLogy();
    cGenGenSmearRecElePosPt->SetTopMargin(0.03);
    cGenGenSmearRecElePosPt->SetRightMargin(0.03);
    cGenGenSmearRecElePosPt->SetLeftMargin(0.13);
    ptGenTrackElePos_rebin->GetYaxis()->SetTitle("N^{gen} Tracks");
    ptGenTrackElePos_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptGenTrackElePos_rebin->GetXaxis()->SetRangeUser(0.0,4.0);
    // ptGenTrackElePos_rebin->SetMaximum(10*std::max({ptGenTrackElePos_rebin->GetMaximum(),ptGenSmearedTrackElePos_rebin->GetMaximum(),ptRecTrackElePos_rebin->GetMaximum()}));
    ptGenTrackElePos_rebin->Draw("hist p e1");
    ptGenSmearedTrackElePos_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_rebin->Draw("hist p e1 same");
    legGenGenSmearRecElePos_top->Draw("same");
    textBField->Draw("same");
    textPIDScenario->Draw("same");
    cGenGenSmearRecElePosPt->SaveAs("./plots/GenGenSmearRecElePos_Pt.png");
    ptGenTrackElePos_rebin->GetXaxis()->SetRangeUser(0.,0.5);
    ptGenTrackElePos_rebin->Draw("hist p e1");
    ptGenSmearedTrackElePos_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_rebin->Draw("hist p e1 same");
    textBField->Draw("same");
    textPIDScenario->Draw("same");
    legGenGenSmearRecElePos_Plus->Draw("same");
    ptGenTrackElePos_beforeKineCuts_rebin->Draw("hist p e1 same");
    // ptGenSmearedTrackElePos_beforeKineCuts_rebin->Draw("hist p e1 same");
    ptRecTrackBeforeSmearing_rebin->Draw("hist p e1 same");
    ptRecTrackAfterSmearing_rebin->Draw("hist p e1 same");
    ptRecTrackAfterKineCuts_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_rebin->Draw("hist p e1 same");
    cGenGenSmearRecElePosPt->SaveAs("./plots/GenGenSmearRecElePos_Pt_0-0.5+woPIDcuts.png");

    auto cGenRecElePt = new TCanvas("cGenRecElePt","cGenRecElePt",800,800);
    // cGenRecElePt->SetLogy();
    cGenRecElePt->SetTopMargin(0.03);
    cGenRecElePt->SetRightMargin(0.03);
    cGenRecElePt->SetLeftMargin(0.13);
    ptGenTrackEle_rebin->GetYaxis()->SetTitle("N^{gen} Tracks");
    ptGenTrackEle_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptGenTrackEle_rebin->GetXaxis()->SetRangeUser(0.0,4.0);
    ptGenTrackEle_rebin->SetMaximum(1.1*std::max({ptGenTrackEle->GetMaximum(),ptRecTrackEle->GetMaximum()}));
    ptGenTrackEle_rebin->Draw("hist p e1");
    ptRecTrackEle_rebin->Draw("hist p e1 same");
    legGenRecEle_top->Draw("same");
    cGenRecElePt->SaveAs("./plots/GenRecEle_Pt.png");


    auto cGenRecPosPt = new TCanvas("cGenRecPosPt","cGenRecPosPt",800,800);
    // cGenRecPosPt->SetLogy();
    cGenRecPosPt->SetTopMargin(0.03);
    cGenRecPosPt->SetRightMargin(0.03);
    cGenRecPosPt->SetLeftMargin(0.13);
    ptGenTrackPos_rebin->GetYaxis()->SetTitle("N^{gen} Tracks");
    ptGenTrackPos_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptGenTrackPos_rebin->GetXaxis()->SetRangeUser(0.0,4.0);
    ptGenTrackPos_rebin->SetMaximum(1.1*std::max({ptGenTrackPos->GetMaximum(),ptRecTrackPos->GetMaximum()}));
    ptGenTrackPos_rebin->Draw("hist p e1");
    ptRecTrackPos_rebin->Draw("hist p e1 same");
    legGenRecPos_top->Draw("same");
    cGenRecPosPt->SaveAs("./plots/GenRecPos_Pt.png");



    auto cEffElePosPtEta = new TCanvas("cEffElePosPtEta","cEffElePosPtEta",800,800);
    gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
    cEffElePosPtEta->SetLogz();
    cEffElePosPtEta->SetLogx();
    cEffElePosPtEta->SetTopMargin(0.03);
    cEffElePosPtEta->SetRightMargin(0.13);
    cEffElePosPtEta->SetLeftMargin(0.13);
    ptEtaEffElePos->SetTitle("");
    ptEtaEffElePos->GetYaxis()->SetTitle("#eta");
    ptEtaEffElePos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptEtaEffElePos->GetZaxis()->SetTitle("track efficiency");
    // ptEtaEffElePos->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEtaEffElePos->GetXaxis()->SetRangeUser(0.0,0.4);
    ptEtaEffElePos->GetYaxis()->SetRangeUser(-1.2,1.2);
    ptEtaEffElePos->GetZaxis()->SetRangeUser(0,1.0);
    ptEtaEffElePos->Draw("Colz");
    // legTrackEle_top->Draw("same");
    cEffElePosPtEta->SaveAs("./plots/EffElePosPtEta.png");



    // // pair efficiencies
    // auto cPairEffULSPt = new TCanvas("cPairEffULSPt","cPairEffULSPt",800,800);
    // // cPairEffULSPt->SetLogy();
    // cPairEffULSPt->SetTopMargin(0.03);
    // cPairEffULSPt->SetRightMargin(0.03);
    // cPairEffULSPt->SetLeftMargin(0.13);
    // ptPairEffULS->GetYaxis()->SetTitle("p_{T,ee}^{rec}/p_{T,ee}^{gen}");
    // ptPairEffULS->GetXaxis()->SetTitle("p_{T,ee} (GeV/#it{c})");
    // ptPairEffULS->GetXaxis()->SetRangeUser(0.0,4.0);
    // ptPairEffULS->Draw("hist p e1 ");
    // legEffULSPair->Draw("same");
    // cPairEffULSPt->SaveAs("./plots/Eff_ULSPair_Pt.png");



  }

  auto cBeforeAfterSmearing = new TCanvas("cBeforeAfterSmearing","cBeforeAfterSmearing",1800,800);

  cBeforeAfterSmearing->Divide(3,1);
  for (size_t i = 1; i < 4; i++) {
    /* code */
    cBeforeAfterSmearing->cd(i)->SetTopMargin(0.03);
    cBeforeAfterSmearing->cd(i)->SetRightMargin(0.03);
  }
  cBeforeAfterSmearing->cd(1)->SetLogy();
  ptRecTrackBeforeSmearing->GetXaxis()->SetRangeUser(0.,8);
  // ptRecTrackBeforeSmearing->SetMinimum(0.01);
  ptRecTrackBeforeSmearing->SetMaximum(1.1*std::max({ptRecTrackBeforeSmearing_rebin->GetMaximum(),ptRecTrackAfterSmearing_rebin->GetMaximum(),ptRecTrackElePos_rebin->GetMaximum()}));
  ptRecTrackBeforeSmearing_rebin->Draw("pe1");
  ptRecTrackAfterSmearing_rebin->Draw("pe1 same");
  ptRecTrackAfterKineCuts_rebin->Draw("pe1 same");
  // ptRecTrackEtaCut_10->Draw("pe1 same");
  // ptRecTrackEtaCut_9->Draw("pe1 same");
  // ptRecTrackEtaCut_8->Draw("pe1 same");
  // ptRecTrackEtaCut_7->Draw("pe1 same");
  // ptRecTrackEtaCut_6->Draw("pe1 same");
  // ptRecTrackEtaCut_5->Draw("pe1 same");
  // ptRecTrackEtaCut_4->Draw("pe1 same");
  // ptRecTrackEtaCut_3->Draw("pe1 same");
  // ptRecTrackEtaCut_2->Draw("pe1 same");
  // ptRecTrackEtaCut_1->Draw("pe1 same");
  ptRecTrackElePos_rebin->Draw("pe1 same");
  legSmearLabel->Draw("same");
  cBeforeAfterSmearing->cd(2);
  etaRecTrackBeforeSmearing->GetXaxis()->SetRangeUser(-2.0,2.0);
  etaRecTrackBeforeSmearing->SetMinimum(0.);
  etaRecTrackBeforeSmearing->SetMaximum(1.1*std::max({etaRecTrackBeforeSmearing->GetMaximum(),etaRecTrackAfterSmearing->GetMaximum(),etaRecTrackElePos->GetMaximum()}));
  etaRecTrackBeforeSmearing->Draw("axis");
  etaRecTrackBeforeSmearing->Draw("pe1 same");
  etaRecTrackAfterSmearing->Draw("pe1 same");
  etaRecTrackAfterKineCuts->Draw("pe1 same");
  // etaRecTrackEtaCut_10->Draw("pe1 same");
  // etaRecTrackEtaCut_8->Draw("pe1 same");
  // etaRecTrackEtaCut_9->Draw("pe1 same");
  // etaRecTrackEtaCut_7->Draw("pe1 same");
  // etaRecTrackEtaCut_6->Draw("pe1 same");
  // etaRecTrackEtaCut_5->Draw("pe1 same");
  // etaRecTrackEtaCut_4->Draw("pe1 same");
  // etaRecTrackEtaCut_3->Draw("pe1 same");
  // etaRecTrackEtaCut_2->Draw("pe1 same");
  // etaRecTrackEtaCut_1->Draw("pe1 same");
  etaRecTrackElePos->Draw("pe1 same");
  cBeforeAfterSmearing->cd(3);
  phiRecTrackBeforeSmearing->SetMinimum(0.);
  phiRecTrackBeforeSmearing->SetMaximum(1.1*std::max({phiRecTrackBeforeSmearing->GetMaximum(),phiRecTrackAfterSmearing->GetMaximum(),phiRecTrackElePos->GetMaximum()}));
  phiRecTrackBeforeSmearing->Draw("axis");
  phiRecTrackBeforeSmearing->Draw("pe1 same");
  phiRecTrackAfterSmearing->Draw("pe1 same");
  phiRecTrackAfterKineCuts->Draw("pe1 same");
  // phiRecTrackEtaCut_10->Draw("pe1 same");
  // phiRecTrackEtaCut_9->Draw("pe1 same");
  // phiRecTrackEtaCut_8->Draw("pe1 same");
  // phiRecTrackEtaCut_7->Draw("pe1 same");
  // phiRecTrackEtaCut_6->Draw("pe1 same");
  // phiRecTrackEtaCut_5->Draw("pe1 same");
  // phiRecTrackEtaCut_4->Draw("pe1 same");
  // phiRecTrackEtaCut_3->Draw("pe1 same");
  // phiRecTrackEtaCut_2->Draw("pe1 same");
  // phiRecTrackEtaCut_1->Draw("pe1 same");
  phiRecTrackElePos->Draw("pe1 same");
  cBeforeAfterSmearing->SaveAs("./plots/RecPtEtaPhiBeforeAfterSmearing.png");

  Double_t lx_low; Double_t lx_high; Double_t ly_low; Double_t ly_high;
  lx_low = 0.1; ly_low = 0.74; lx_high = 0.35; ly_high = 0.95;
  auto legendInfo = new TLegend(lx_low,ly_low,lx_high,ly_high);
  TLegendEntry *entryInfo0=legendInfo->AddEntry("ALICE 3","ALICE 3 Study","");
  TLegendEntry *entryInfo1=legendInfo->AddEntry("collisionSystem",Form("0-10%s %s, #sqrt{#it{s}_{NN}} = 5.02 TeV","%",collSystem.Data()),"");
  // TLegendEntry *entryInfo4=legendInfo->AddEntry("Generator",Form("Phythia8 Angantyr, #it{B} = %g T",BField),"");
  // TLegendEntry *entryInfo2;
  // // TLegendEntry *entryInfo3;
  // if(BField == 0.2){
  //   if(ith_PIDscenario == 1){
  //     // entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.04 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //     entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.03 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //   }
  //   if(ith_PIDscenario == 2){
  //     // entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.08 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //     entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.03 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //   }
  // }
  // else if (BField == 0.5){
  //   if(ith_PIDscenario == 1){
  //     // entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.2 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //     entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.03 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //   }
  //   if(ith_PIDscenario >= 2){
  //     // entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.08 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //     entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.03 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
  //   }
  // }

  // TLegendEntry *entryInfo5=legendInfo->AddEntry("Layout",Form("Layout v1, |#it{#eta}_{e}| < 1.1, #it{B} = %g T",BField),"");
  TLegendEntry *entryInfo5=legendInfo->AddEntry("Layout",Form("TOF 20cm, |#it{#eta}_{e}| < 1.1, #it{B} = %g T",BField),"");

  TLegendEntry *entryInfo3=legendInfo->AddEntry("PID",Form("%s PID",strPIDscenario[ith_PIDscenario-1].Data()),"");
  legendInfo->SetBorderSize(0);
  legendInfo->SetFillColorAlpha(0, 0.0);
  // legendInfo->SetTextSize(0.025);
  legendInfo->SetTextFont(43);
  legendInfo->SetTextSize(25);

  auto legendInfoTOF = new TLegend(0.4, 0.12, 0.87, 0.3);
  legendInfoTOF->AddEntry(entryInfo0,entryInfo0->GetLabel(),"");
  legendInfoTOF->AddEntry(entryInfo1,entryInfo1->GetLabel(),"");
  legendInfoTOF->AddEntry(entryInfo5,entryInfo5->GetLabel(),"");
  legendInfoTOF->AddEntry(entryInfo3,entryInfo3->GetLabel(),"");
  legendInfoTOF->SetBorderSize(0);
  legendInfoTOF->SetFillColorAlpha(0, 0.0);
  // legendInfoTOF->SetTextSize(0.025);
  legendInfoTOF->SetTextFont(43);
  legendInfoTOF->SetTextSize(25);

  auto legendInfoNoPID = new TLegend(lx_low, ly_low, lx_high, ly_high);
  legendInfoNoPID->AddEntry(entryInfo0,entryInfo0->GetLabel(),"");
  legendInfoNoPID->AddEntry(entryInfo1,entryInfo1->GetLabel(),"");
  legendInfoNoPID->AddEntry(entryInfo5,entryInfo5->GetLabel(),"");
  legendInfoNoPID->AddEntry("noPID","no PID cuts","");
  legendInfoNoPID->SetBorderSize(0);
  legendInfoNoPID->SetFillColorAlpha(0, 0.0);
  legendInfoNoPID->SetTextFont(43);
  legendInfoNoPID->SetTextSize(25);

  auto legendInfoNoPIDtof = new TLegend(0.4, 0.12, 0.87, 0.3);
  legendInfoNoPIDtof->AddEntry(entryInfo0,entryInfo0->GetLabel(),"");
  legendInfoNoPIDtof->AddEntry(entryInfo1,entryInfo1->GetLabel(),"");
  legendInfoNoPIDtof->AddEntry(entryInfo5,entryInfo5->GetLabel(),"");
  legendInfoNoPIDtof->AddEntry("noPID","no PID cuts","");
  legendInfoNoPIDtof->SetBorderSize(0);
  legendInfoNoPIDtof->SetFillColorAlpha(0, 0.0);
  legendInfoNoPIDtof->SetTextFont(43);
  legendInfoNoPIDtof->SetTextSize(25);

  auto legendInfoCont = (TLegend*)legendInfo->Clone();
  legendInfoCont->AddEntry("Generator","Phythia8 Angantyr","");



  if (bPlotPIDhistograms) {
    // Plotting PID plots for TOF and RICH, + NSigma PID
    TString DocumentPath = "./plots/PID_histograms";
    gSystem->Exec(Form("mkdir -p %s",DocumentPath.Data()));
    auto cPID = new TCanvas("cPID","cPID",800,800);
    auto cPID_logx = new TCanvas("cPID_logx","cPID_logx",800,800);

    for (size_t i = 0; i < vecTOF_PIDplots.size(); i++) {
      gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
      // cPID->cd();
      // cPID->SetLogz();
      // cPID->SetTopMargin(0.03);
      // cPID->SetRightMargin(0.07);
      // cPID->SetLeftMargin(0.13);
      // vecTOF_PIDplots.at(i)->Draw("COL ");
      // legendInfo->Draw("same");
      // if(i % 2 == 0) legendInfoNoPIDtof->Draw("same");
      // else legendInfoTOF->Draw("same");
      // cPID->SetTicks();
      // cPID->SaveAs(Form("./plots/PID_histograms/%s.png",vecTOF_PIDplots.at(i)->GetName()));

      cPID_logx->cd();
      cPID_logx->SetLogx();
      cPID_logx->SetLogz();
      cPID_logx->SetTopMargin(0.03);
      cPID_logx->SetRightMargin(0.07);
      cPID_logx->SetLeftMargin(0.13);
      vecTOF_PIDplots.at(i)->GetYaxis()->SetRangeUser(-15,25);
      vecTOF_PIDplots.at(i)->GetXaxis()->SetRangeUser(0.01,10);
      vecTOF_PIDplots.at(i)->Draw("COL ");

      if(i % 2 == 0) legendInfoNoPIDtof->Draw("same");
      else legendInfoTOF->Draw("same");
      cPID_logx->SetTicks();
      cPID_logx->SaveAs(Form("./plots/PID_histograms/%s_logx.pdf",vecTOF_PIDplots.at(i)->GetName()));
    }
    for (size_t i = 0; i < vecRICH_PIDplots.size(); i++) {
      gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
      // cPID->cd();
      // cPID->SetLogz();
      // // cPID->SetTopMargin(0.03);
      // // cPID->SetRightMargin(0.03);
      // // cPID->SetLeftMargin(0.13);
      // vecRICH_PIDplots.at(i)->Draw("COL ");
      // legendInfo->Draw("same");
      // if(i % 2 == 0) legendInfoNoPID->Draw("same");
      // else legendInfo->Draw("same");
      // cPID->SetTicks();
      // cPID->SaveAs(Form("./plots/PID_histograms/%s.png",vecRICH_PIDplots.at(i)->GetName()));

      cPID_logx->cd();
      cPID_logx->SetLogx();
      cPID_logx->SetLogz();
      // cPID_logx->SetTopMargin(0.03);
      // cPID_logx->SetRightMargin(0.03);
      // cPID_logx->SetLeftMargin(0.13);
      // vecRICH_PIDplots.at(i)->GetYaxis()->SetRangeUser(-25,15);
      vecRICH_PIDplots.at(i)->Draw("COL");
      if(i % 2 == 0) {
        legendInfoNoPID->Draw("same");
        vecRICH_PIDplots.at(i)->GetYaxis()->SetRangeUser(-25,15);
      }
      else {
        legendInfo->Draw("same");
        vecRICH_PIDplots.at(i)->GetYaxis()->SetRangeUser(-6,10);
      }
      cPID_logx->SetTicks();
      cPID_logx->SaveAs(Form("./plots/PID_histograms/%s_logx.pdf",vecRICH_PIDplots.at(i)->GetName()));
    }
  }

  if (bPlotTrackContamination) {
    auto cPosNegContaminationPtEtaPhi = new TCanvas("cPosNegContaminationPtEtaPhi","cPosNegContaminationPtEtaPhi",1350,600);
    cPosNegContaminationPtEtaPhi->Divide(3,1);
    for (size_t i = 1; i < 4; i++) {
      cPosNegContaminationPtEtaPhi->cd(i)->SetTopMargin(0.03);
      cPosNegContaminationPtEtaPhi->cd(i)->SetRightMargin(0.03);
    }
    cPosNegContaminationPtEtaPhi->cd(1);
    hPureContaminationRecPtNeg->GetYaxis()->SetTitle("contamination");
    hPureContaminationRecPtNeg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPureContaminationRecPtNeg->GetXaxis()->SetRangeUser(0.0,4.0);
    hPureContaminationRecPtNeg->SetMaximum(1.1);
    // hPureContaminationRecPtNeg->SetMaximum(1.1*std::max({hPureContaminationRecPtNeg->GetMaximum(),hPureContaminationRecPtPos->GetMaximum()}));
    hPureContaminationRecPtNeg->SetMinimum(0.);
    // hPurityRecPtNeg->Draw("hist p e1");
    // hPurityRecPtPos->Draw("same hist p e1");
    hPureContaminationRecPtNeg->Draw("hist p e1");
    hPureContaminationRecPtPos->Draw("same hist p e1");
    legContamination->Draw("same");
    cPosNegContaminationPtEtaPhi->cd(2);
    hPureContaminationRecEtaNeg->GetYaxis()->SetTitle("contamination");
    hPureContaminationRecEtaNeg->GetXaxis()->SetTitle("#eta");
    hPureContaminationRecEtaNeg->GetXaxis()->SetRangeUser(-2.0,2.0);
    hPureContaminationRecEtaNeg->SetMaximum(0.1);
    // hPureContaminationRecEtaNeg->SetMaximum(1.1*std::max({hPureContaminationRecEtaNeg->GetMaximum(),hPureContaminationRecEtaPos->GetMaximum()}));
    hPureContaminationRecEtaNeg->SetMinimum(0.);
    // hPurityRecEtaNeg->Draw("hist p e1");
    // hPurityRecEtaPos->Draw("same hist p e1");
    hPureContaminationRecEtaNeg->Draw("hist p e1");
    hPureContaminationRecEtaPos->Draw("same hist p e1");
    cPosNegContaminationPtEtaPhi->cd(3);
    hPureContaminationRecPhiNeg->GetYaxis()->SetTitle("contamination");
    hPureContaminationRecPhiNeg->GetXaxis()->SetTitle("#varphi (rad)");
    hPureContaminationRecPhiNeg->SetMaximum(0.1);
    // hPureContaminationRecPhiNeg->SetMaximum(1.1*std::max({hPureContaminationRecPhiNeg->GetMaximum(),hPureContaminationRecPhiPos->GetMaximum()}));
    hPureContaminationRecPhiNeg->SetMinimum(0.);
    hPureContaminationRecPhiNeg->GetXaxis()->SetRangeUser(-7.0,7.0);
    // hPurityRecPhiNeg->Draw("hist p e1");
    // hPurityRecPhiPos->Draw("same hist p e1");
    hPureContaminationRecPhiNeg->Draw("hist p e1");
    hPureContaminationRecPhiPos->Draw("same hist p e1");
    cPosNegContaminationPtEtaPhi->SaveAs("./plots/PID_PosNeg_Contamiantion.png");


    auto cContaminationPtEtaPhi = new TCanvas("cContaminationPtEtaPhi","cContaminationPtEtaPhi",1350,600);
    cContaminationPtEtaPhi->Divide(3,1);
    for (size_t i = 1; i < 4; i++) {
      cContaminationPtEtaPhi->cd(i)->SetTopMargin(0.03);
      cContaminationPtEtaPhi->cd(i)->SetRightMargin(0.03);
    }
    cContaminationPtEtaPhi->cd(1);
    hMuonContaminationRecPt->GetYaxis()->SetTitle("contamination");
    hMuonContaminationRecPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMuonContaminationRecPt->GetXaxis()->SetRangeUser(0.0,4.0);
    hMuonContaminationRecPt->SetMaximum(1.1);
    hMuonContaminationRecPt->SetMinimum(0.);
    hMuonContaminationRecPt->Draw("hist p e1");
    hPionContaminationRecPt->Draw("same hist p e1");
    hKaonContaminationRecPt->Draw("same hist p e1");
    hProtonContaminationRecPt->Draw("same hist p e1");
    hTotalPureContaminationRecPt->Draw("same hist p e1");
    legPIDContamination->Draw("same");
    cContaminationPtEtaPhi->cd(2);
    hMuonContaminationRecEta->GetYaxis()->SetTitle("contamination");
    hMuonContaminationRecEta->GetXaxis()->SetTitle("#eta");
    hMuonContaminationRecEta->GetXaxis()->SetRangeUser(-2.0,2.0);
    hMuonContaminationRecEta->SetMaximum(0.1);
    hMuonContaminationRecEta->SetMinimum(0.);
    hMuonContaminationRecEta->Draw("hist p e1");
    hPionContaminationRecEta->Draw("same hist p e1");
    hKaonContaminationRecEta->Draw("same hist p e1");
    hProtonContaminationRecEta->Draw("same hist p e1");
    hTotalPureContaminationRecEta->Draw("same hist p e1");
    cContaminationPtEtaPhi->cd(3);
    hMuonContaminationRecPhi->GetYaxis()->SetTitle("contamination");
    hMuonContaminationRecPhi->GetXaxis()->SetTitle("#varphi (rad)");
    hMuonContaminationRecPhi->SetMaximum(0.1);
    hMuonContaminationRecPhi->SetMinimum(0.);
    hMuonContaminationRecPhi->GetXaxis()->SetRangeUser(-7.0,7.0);
    hMuonContaminationRecPhi->Draw("hist p e1");
    hPionContaminationRecPhi->Draw("same hist p e1");
    hKaonContaminationRecPhi->Draw("same hist p e1");
    hProtonContaminationRecPhi->Draw("same hist p e1");
    hTotalPureContaminationRecPhi->Draw("same hist p e1");
    cContaminationPtEtaPhi->SaveAs("./plots/PID_Contamiantion.png");



    auto cContaminationPt = new TCanvas("cContaminationPt","cContaminationPt",800,800);
    cContaminationPt->SetTopMargin(0.03);
    cContaminationPt->SetRightMargin(0.03);
    cContaminationPt->SetLeftMargin(0.13);
    hMuonContaminationRecPt->GetYaxis()->SetTitle("contamination");
    hMuonContaminationRecPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hMuonContaminationRecPt->GetXaxis()->SetRangeUser(0.0,0.4);
    // hMuonContaminationRecPt->GetXaxis()->SetRangeUser(0.0,4.0);
    hMuonContaminationRecPt->SetMaximum(0.12);
    hMuonContaminationRecPt->SetMinimum(0.);
    hMuonContaminationRecPt->Draw("hist p e1");
    hPionContaminationRecPt->Draw("same hist p e1");
    hKaonContaminationRecPt->Draw("same hist p e1");
    hProtonContaminationRecPt->Draw("same hist p e1");
    hTotalPureContaminationRecPt->Draw("same hist p e1");
    legPIDContamination->Draw("same");
    // textBField_conta->Draw("same");
    // textPIDScenario_conta->Draw("same");
    legendInfoCont->Draw("same");
    gPad->SetTicks();
    cContaminationPt->SaveAs("./plots/PID_Pt_Contamiantion.pdf");
    // hTotalPureContaminationRecPt->Draw("hist p e1");
    // legTotalContamination->Draw("same");
    // textBField_conta->Draw("same");
    // textPIDScenario_conta->Draw("same");
    // cContaminationPt->SaveAs("./plots/PID_Pt_TotalContamiantion.png");
  }



  auto cRejectionFactorPt = new TCanvas("cRejectionFactorPt","cRejectionFactorPt",800,800);
  cRejectionFactorPt->SetTopMargin(0.03);
  cRejectionFactorPt->SetRightMargin(0.03);
  cRejectionFactorPt->SetLeftMargin(0.13);
  cRejectionFactorPt->SetLogy();
  hPionRejectionFactorPt->GetYaxis()->SetTitle("rejection factor");
  hPionRejectionFactorPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  // hPionRejectionFactorPt->GetXaxis()->SetRangeUser(0.0,4.0);
  // hPionRejectionFactorPt->SetMaximum(10*std::max({hMuonRejectionFactorPt->GetMaximum(),hPionRejectionFactorPt->GetMaximum(),hKaonRejectionFactorPt->GetMaximum(),hProtonRejectionFactorPt->GetMaximum()/*,hTotalRejectionFactorPt->GetMaximum()*/}));
  // hPionRejectionFactorPt->SetMaximum(1.1);
  // hPionRejectionFactorPt->SetMinimum(0.1*std::min({hMuonRejectionFactorPt->GetMinimum(),hPionRejectionFactorPt->GetMinimum(),hKaonRejectionFactorPt->GetMinimum(),hProtonRejectionFactorPt->GetMinimum()/*,hTotalRejectionFactorPt->GetMaximum()*/}));
  // hPionRejectionFactorPt->SetMaximum(100000000000);
  hPionRejectionFactorPt->SetMaximum(100000000);
  hPionRejectionFactorPt->SetMinimum(0.1);
  hPionRejectionFactorPt->GetXaxis()->SetRangeUser(0.0,0.4);
  hMuonRejectionFactorPt->SetMinimum(0.1);
  hPionRejectionFactorPt->Draw("hist p e1");
  hMuonRejectionFactorPt->Draw("same hist p e1");
  hKaonRejectionFactorPt->Draw("same hist p e1");
  hProtonRejectionFactorPt->Draw("same hist p e1");
  // hTotalRejectionFactorPt->Draw("same hist p e1");
  legRejFactor->Draw("same");
  legendInfo->Draw("same");
  // textBField_conta->Draw("same");
  // textPIDScenario_conta->Draw("same");
  cRejectionFactorPt->SetTicks();
  cRejectionFactorPt->SaveAs("./plots/PID_Pt_RejectionFactor.pdf");
  hPionRejectionFactorPt->SetMinimum(1000);
  hPionRejectionFactorPt->SetMaximum(100000000);
  // hPionRejectionFactorPt->GetXaxis()->SetRangeUser(0.0,3.5);
  hPionRejectionFactorPt->Draw("hist p e1");
  TLine* line = new TLine(0.,100000.,4.,100000.);
  line->SetLineStyle(2);
  line->Draw("same");
  legRejFactor_pi->Draw("same");
  // textBField_conta->Draw("same");
  textPIDScenario_conta->Draw("same");
  legendInfo->Draw("same");
  cRejectionFactorPt->SetTicks();
  cRejectionFactorPt->SaveAs("./plots/PID_Pt_RejectionFactor_pion.png");
  // hMuonRejectionFactorPt->Draw("hist p e1");
  // legRejFactor_mu->Draw("same");
  // textBField_conta->Draw("same");
  // textPIDScenario_conta->Draw("same");
  // cRejectionFactorPt->SaveAs("./plots/PID_Pt_RejectionFactor_muon.png");


  if (bPlotPIDhistograms) {
    const char *namesEleMuPi[3] = {"Ele", "Mu", "Pi"};

    for (size_t iNSig = 0; iNSig < (sizeof(namesEleMuPi)/sizeof(namesEleMuPi[0])); iNSig++) {
      auto cNSigmaSeparatePID_TOF = new TCanvas("cNSigmaPionSeparatePID_TOF","cNSigmaPionSeparatePID_TOF",800,800);
      cNSigmaSeparatePID_TOF->SetTopMargin(0.03);
      cNSigmaSeparatePID_TOF->SetRightMargin(0.03);
      cNSigmaSeparatePID_TOF->SetLeftMargin(0.13);
      cNSigmaSeparatePID_TOF->SetLogx();
      // hNsigmaP_TOF_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,4.0);
      // hNsigmaP_TOF_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,0.4);
      hNsigmaP_TOF_truePion[iNSig]->GetYaxis()->SetRangeUser(-15.,25.);
      hNsigmaP_TOF_truePion[iNSig]->Draw("");
      hNsigmaP_TOF_trueElec[iNSig]->Draw("same");
      hNsigmaP_TOF_trueMuon[iNSig]->Draw("same");
      legPIDSeparateColorTOF->Draw("same");
      legendInfoNoPIDtof->Draw("same");
      cNSigmaSeparatePID_TOF->SetTicks();
      cNSigmaSeparatePID_TOF->SaveAs(Form("./plots/PID_histograms/NSigma%s_TOF_SeparatePiEleMuPID_noCuts.pdf",namesEleMuPi[iNSig]));

      auto cNSigmaSeparatePID_TOF_afterCuts = new TCanvas("cNSigmaPionSeparatePID_TOF_afterCuts","cNSigmaPionSeparatePID_TOF_afterCuts",800,800);
      cNSigmaSeparatePID_TOF_afterCuts->SetTopMargin(0.03);
      cNSigmaSeparatePID_TOF_afterCuts->SetRightMargin(0.03);
      cNSigmaSeparatePID_TOF_afterCuts->SetLeftMargin(0.13);
      cNSigmaSeparatePID_TOF_afterCuts->SetLogx();
      // hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,4.0);
      hNsigmaP_afterPIDcuts_TOF_trueElec[iNSig]->GetXaxis()->SetRangeUser(0.0,0.4);
      hNsigmaP_afterPIDcuts_TOF_trueElec[iNSig]->GetYaxis()->SetRangeUser(-15.,25.);
      // hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->SetMaximum(1.1);
      // hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->SetMinimum(0.);
      hNsigmaP_afterPIDcuts_TOF_trueElec[iNSig]->Draw("");
      hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->Draw("same");
      hNsigmaP_afterPIDcuts_TOF_trueMuon[iNSig]->Draw("same");
      legPIDSeparateColorTOF->Draw("same");
      legendInfoTOF->Draw("same");
      cNSigmaSeparatePID_TOF_afterCuts->SetTicks();
      cNSigmaSeparatePID_TOF_afterCuts->SaveAs(Form("./plots/PID_histograms/NSigma%s_TOF_SeparatePiEleMuPID_afterCuts.pdf",namesEleMuPi[iNSig]));



      auto cNSigmaSeparatePID_RICH = new TCanvas("cNSigmaPionSeparatePID_RICH","cNSigmaPionSeparatePID_RICH",800,800);
      cNSigmaSeparatePID_RICH->SetTopMargin(0.03);
      cNSigmaSeparatePID_RICH->SetRightMargin(0.03);
      cNSigmaSeparatePID_RICH->SetLeftMargin(0.13);
      cNSigmaSeparatePID_RICH->SetLogx();
      // hNsigmaP_RICH_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,4.0);
      hNsigmaP_RICH_truePion[iNSig]->GetYaxis()->SetRangeUser(-25.,15.);
      // hNsigmaP_RICH_truePion[iNSig]->SetMaximum(1.1);
      // hNsigmaP_RICH_truePion[iNSig]->SetMinimum(0.);
      hNsigmaP_RICH_truePion[iNSig]->Draw("");
      hNsigmaP_RICH_trueElec[iNSig]->Draw("same");
      hNsigmaP_RICH_trueMuon[iNSig]->Draw("same");
      legPIDSeparateColor->Draw("same");
      legendInfoNoPID->Draw("same");
      cNSigmaSeparatePID_RICH->SetTicks();
      cNSigmaSeparatePID_RICH->SaveAs(Form("./plots/PID_histograms/NSigma%s_RICH_SeparatePiEleMuPID_noCuts.pdf",namesEleMuPi[iNSig]));


      auto cNSigmaSeparatePID_RICH_afterCuts = new TCanvas("cNSigmaPionSeparatePID_RICH_afterCuts","cNSigmaPionSeparatePID_RICH_afterCuts",800,800);
      cNSigmaSeparatePID_RICH_afterCuts->SetTopMargin(0.03);
      cNSigmaSeparatePID_RICH_afterCuts->SetRightMargin(0.03);
      cNSigmaSeparatePID_RICH_afterCuts->SetLeftMargin(0.13);
      cNSigmaSeparatePID_RICH_afterCuts->SetLogx();
      // hNsigmaP_afterPIDcuts_RICH_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,4.0);
      hNsigmaP_afterPIDcuts_RICH_trueElec[iNSig]->GetYaxis()->SetRangeUser(-6.,10.);
      // hNsigmaP_afterPIDcuts_RICH_truePion[iNSig]->SetMaximum(1.1);
      // hNsigmaP_afterPIDcuts_RICH_truePion[iNSig]->SetMinimum(0.);
      hNsigmaP_afterPIDcuts_RICH_trueElec[iNSig]->Draw("");
      hNsigmaP_afterPIDcuts_RICH_truePion[iNSig]->Draw("same");
      hNsigmaP_afterPIDcuts_RICH_trueMuon[iNSig]->Draw("same");
      legPIDSeparateColor->Draw("same");
      legendInfo->Draw("same");
      cNSigmaSeparatePID_RICH_afterCuts->SetTicks();
      cNSigmaSeparatePID_RICH_afterCuts->SaveAs(Form("./plots/PID_histograms/NSigma%s_RICH_SeparatePiEleMuPID_afterCuts.pdf",namesEleMuPi[iNSig]));



      // //projections in RICH
      // auto cPprojectionsPID_RICH = new TCanvas("cPprojectionsPID_RICH","cPprojectionsPID_RICHcPprojectionsPID_RICH",800,800);
      // cPprojectionsPID_RICH->SetTopMargin(0.03);
      // cPprojectionsPID_RICH->SetRightMargin(0.03);
      // cPprojectionsPID_RICH->SetLeftMargin(0.13);
      // cPprojectionsPID_RICH->SetLogx();
      // cPprojectionsPID_RICH->SetLogy();
      // // hNsigmaP_RICH_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,4.0);
      // hNsigmaP_RICH_truePion[iNSig]->GetYaxis()->SetRangeUser(-25.,25.);
      // // hNsigmaP_RICH_truePion[iNSig]->SetMaximum(1.1);
      // // hNsigmaP_RICH_truePion[iNSig]->SetMinimum(0.);
      // auto h1 = hNsigmaP_RICH_truePion[iNSig]->ProjectionX();
      // auto h2 = hNsigmaP_RICH_trueElec[iNSig]->ProjectionX();
      // auto h3 = hNsigmaP_RICH_trueMuon[iNSig]->ProjectionX();
      // h1->SetMarkerColor(kRed);
      // h2->SetMarkerColor(kBlue);
      // h3->SetMarkerColor(kBlack);
      // h1->Draw("");
      // h2->Draw("same");
      // h3->Draw("same");
      // legPIDSeparateColor->Draw("same");
      // legendInfo->Draw("same");
      // textBeforePID->Draw("same");
      // textRICH->Draw("same");
      // cPprojectionsPID_RICH->SetTicks();
      // cPprojectionsPID_RICH->SaveAs(Form("./plots/PID_histograms/NSigma%s_projX_RICH_SeparatePiEleMuPID_noCuts.png",namesEleMuPi[iNSig]));
      //
      // auto cSigmaProjectionsPID_RICH = new TCanvas("cSigmaProjectionsPID_RICH","cSigmaProjectionsPID_RICH",800,800);
      // cSigmaProjectionsPID_RICH->SetTopMargin(0.03);
      // cSigmaProjectionsPID_RICH->SetRightMargin(0.03);
      // cSigmaProjectionsPID_RICH->SetLeftMargin(0.13);
      // cSigmaProjectionsPID_RICH->SetLogx();
      // cSigmaProjectionsPID_RICH->SetLogy();
      // // hNsigmaP_RICH_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,4.0);
      // hNsigmaP_RICH_truePion[iNSig]->GetYaxis()->SetRangeUser(-25.,25.);
      // // hNsigmaP_RICH_truePion[iNSig]->SetMaximum(1.1);
      // // hNsigmaP_RICH_truePion[iNSig]->SetMinimum(0.);
      // auto h4 = hNsigmaP_RICH_truePion[iNSig]->ProjectionY();
      // auto h5 = hNsigmaP_RICH_trueElec[iNSig]->ProjectionY();
      // auto h6 = hNsigmaP_RICH_trueMuon[iNSig]->ProjectionY();
      // h4->SetMarkerColor(kRed);
      // h5->SetMarkerColor(kBlue);
      // h6->SetMarkerColor(kBlack);
      // h4->Draw("");
      // h5->Draw("same");
      // h6->Draw("same");
      // legPIDSeparateColor->Draw("same");
      // legendInfo->Draw("same");
      // textBeforePID->Draw("same");
      // textRICH->Draw("same");
      // cSigmaProjectionsPID_RICH->SetTicks();
      // cSigmaProjectionsPID_RICH->SaveAs(Form("./plots/PID_histograms/NSigma%s_projY_RICH_SeparatePiEleMuPID_noCuts.png",namesEleMuPi[iNSig]));

    }


    Double_t startTextX         = .16;
    Double_t startTextY         = .8;
    Double_t textSize           = 20;
    auto cProjNsigma = new TCanvas("cProjNsigma","cProjNsigma",800,800);
    auto cProjNsigma_afterCuts = new TCanvas("cProjNsigma_afterCuts","cProjNsigma_afterCuts",800,800);
    cProjNsigma->SetTopMargin(0.03);
    cProjNsigma->SetRightMargin(0.03);
    cProjNsigma->SetLeftMargin(0.13);
    cProjNsigma_afterCuts->SetTopMargin(0.03);
    cProjNsigma_afterCuts->SetRightMargin(0.03);
    cProjNsigma_afterCuts->SetLeftMargin(0.13);
    const char *pnameElePi[2] = {"el", "pi"};
    if(bPlotNSigmaProjections){
      for (size_t j = 0; j < vec_proj_bin_p.size()-1;  j++) {      // loop over all p projection intervalls
        TLatex *p_Intervall    = new TLatex(startTextX, startTextY    , Form("%g < #font[12]{p} < %g  #frac{GeV}{#font[12]{c}}",vec_proj_bin_p.at(j) ,vec_proj_bin_p.at(j+1)));
        SetTextSettings(p_Intervall,textSize);
        Int_t startbinY = vecHistosNSigma_afterCuts.at(0)->GetYaxis()->FindBin(vec_proj_bin_p[j]);      // select start bin of pt projection
        Int_t endbinY   = vecHistosNSigma_afterCuts.at(0)->GetYaxis()->FindBin(vec_proj_bin_p[j+1]);    // select end bin of pt projection
        int i = 0;
        for (size_t iHist = 0; iHist < vecHistosNSigma_afterCuts.size()-1; iHist+=2) {
          projNSigmas_afterCuts_trueElec = (TH1F*) vecHistosNSigma_afterCuts.at(iHist)->ProjectionY(Form("hProjNsigmaP_%s_afterPIDcuts_RICH_Proj_p%g:%g", pnameElePi[i],vec_proj_bin_p.at(j),vec_proj_bin_p.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram
          projNSigmas_afterCuts_truePion = (TH1F*) vecHistosNSigma_afterCuts.at(iHist+1)->ProjectionY(Form("hProjNsigmaP_%s_afterPIDCuts_RICH_Proj_p%g:%g", pnameElePi[i],vec_proj_bin_p.at(j),vec_proj_bin_p.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram
          projNSigmas_trueElec = (TH1F*) vecHistosNSigma.at(iHist)->ProjectionY(Form("hProjNsigmaP_%s_RICH_Proj_p%g:%g", pnameElePi[i],vec_proj_bin_p.at(j),vec_proj_bin_p.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram
          projNSigmas_truePion = (TH1F*) vecHistosNSigma.at(iHist+1)->ProjectionY(Form("hProjNsigmaP_%s_RICH_Proj_p%g:%g", pnameElePi[i],vec_proj_bin_p.at(j),vec_proj_bin_p.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram
          projNSigmas_afterCuts_trueElec->GetXaxis()->SetRangeUser(-10.,10.);
          projNSigmas_afterCuts_truePion->GetXaxis()->SetRangeUser(-10.,10.);
          projNSigmas_trueElec->GetXaxis()->SetRangeUser(-10.,10.);
          projNSigmas_truePion->GetXaxis()->SetRangeUser(-10.,10.);

          TLatex *elecScale_afterCuts    = new TLatex(startTextX, startTextY-0.03    , Form("e Scale x%g"  ,projNSigmas_afterCuts_trueElec->Integral()));
          TLatex *pionScale_afterCuts    = new TLatex(startTextX, startTextY-0.06    , Form("#pi Scale x%g",projNSigmas_afterCuts_truePion->Integral()));
          TLatex *elecScale    = new TLatex(startTextX, startTextY-0.04    , Form("e Scale x%g"  ,projNSigmas_trueElec->Integral()));
          TLatex *pionScale    = new TLatex(startTextX, startTextY-0.07    , Form("#pi Scale x%g",projNSigmas_truePion->Integral()));
          SetTextSettings(elecScale_afterCuts,textSize);
          SetTextSettings(pionScale_afterCuts,textSize);
          SetTextSettings(elecScale,textSize);
          SetTextSettings(pionScale,textSize);

          makeHistNice(projNSigmas_afterCuts_trueElec,kBlue+1);
          makeHistNice(projNSigmas_afterCuts_truePion,kRed+2);
          makeHistNice(projNSigmas_trueElec,kBlue+1);
          makeHistNice(projNSigmas_truePion,kRed+2);
          projNSigmas_trueElec->SetMarkerStyle(24);
          projNSigmas_truePion->SetMarkerStyle(24);

          projNSigmas_afterCuts_trueElec->Scale(1./projNSigmas_afterCuts_trueElec->Integral());
          projNSigmas_afterCuts_truePion->Scale(1./projNSigmas_afterCuts_truePion->Integral());
          projNSigmas_trueElec->Scale(1./projNSigmas_trueElec->Integral());
          projNSigmas_truePion->Scale(1./projNSigmas_truePion->Integral());

          // // Normalise Histograms to Number of Events
          // TF1* fTF1dummy = new TF1("dummy","1",-100,100);
          // projY->Divide(fTF1dummy,nEvents);
          double legPosElePi[4] = {0.16,0.83,0.35,0.93};

          auto legTrackElePi_top = new TLegend(legPosElePi[0],legPosElePi[1],legPosElePi[2],legPosElePi[3]);
          legTrackElePi_top->SetBorderSize(0);
          legTrackElePi_top->SetFillStyle(0);
          legTrackElePi_top->SetTextSize(0.025);
          legTrackElePi_top->AddEntry(projNSigmas_trueElec,"electrons","p");
          legTrackElePi_top->AddEntry(projNSigmas_truePion,"pions","p");

          cProjNsigma->cd();
          projNSigmas_trueElec->SetMaximum(1.1*std::max({projNSigmas_trueElec->GetMaximum(),projNSigmas_truePion->GetMaximum()}));
          projNSigmas_trueElec->Draw("hist p e1");
          projNSigmas_truePion->Draw("same hist p e1");
          p_Intervall->Draw("same");
          elecScale->Draw("same");
          pionScale->Draw("same");
          legTrackElePi_top->Draw("same");
          cProjNsigma->SaveAs(Form("./plots/PID_histograms/%s.png",projNSigmas_trueElec->GetName()));

          cProjNsigma_afterCuts->cd();
          projNSigmas_afterCuts_trueElec->SetMaximum(1.1*std::max({projNSigmas_afterCuts_trueElec->GetMaximum(),projNSigmas_afterCuts_truePion->GetMaximum()}));
          projNSigmas_afterCuts_trueElec->Draw("hist p e1");
          projNSigmas_afterCuts_truePion->Draw("same hist p e1");
          p_Intervall->Draw("same");
          elecScale_afterCuts->Draw("same");
          pionScale_afterCuts->Draw("same");
          legTrackElePi_top->Draw("same");
          cProjNsigma_afterCuts->SaveAs(Form("./plots/PID_histograms/%s.png",projNSigmas_afterCuts_trueElec->GetName()));
          i++;
        }
      }
    }
  }


  if (bPlotPairHistograms) {
    if (bPlotULS) {
      //ULS as TH2 and projections of mee and ptee
      // auto cMee_primCCBB = new TCanvas("cMee","cMee",800,800);
      // cMee_primCCBB->SetLogy();
      // cMee_primCCBB->SetTopMargin(0.03);
      // cMee_primCCBB->SetRightMargin(0.03);
      // cMee_primCCBB->SetLeftMargin(0.13);
      // proj_recULS_MeePrim->GetYaxis()->SetTitle("counts");
      // proj_recULS_MeePrim->Draw("hist p e1");
      // proj_recULS_MeeCC->Draw("same hist p e1");
      // proj_recULS_MeeBB->Draw("same hist p e1");
      // legPair->Draw("same");
      // cMee_primCCBB->SaveAs("./plots/Mee_primCCBB.png");

      // auto cPtee_primCCBB = new TCanvas("cPtee","cPtee",800,800);
      // cPtee_primCCBB->SetLogy();
      // cPtee_primCCBB->SetTopMargin(0.03);
      // cPtee_primCCBB->SetRightMargin(0.03);
      // cPtee_primCCBB->SetLeftMargin(0.13);
      // proj_recULS_PteePrim->GetYaxis()->SetTitle("counts");
      // proj_recULS_PteePrim->SetMaximum(0.8);
      // proj_recULS_PteePrim->Draw("hist p e1");
      // proj_recULS_PteeCC->Draw("same hist p e1");
      // proj_recULS_PteeBB->Draw("same hist p e1");
      // legPair->Draw("same");
      // cPtee_primCCBB->SaveAs("./plots/Ptee_primCCBB.png");

      auto cprojMeePteeDCA = new TCanvas("cprojMeePteeDCA","cprojMeePteeDCA",1800,800);
      cprojMeePteeDCA->Divide(3,1);
      for (size_t i = 1; i < 4; i++) {
        cprojMeePteeDCA->cd(i)->SetTopMargin(0.03);
        cprojMeePteeDCA->cd(i)->SetRightMargin(0.03);
        gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
      }
      Double_t lx_low; Double_t lx_high; Double_t ly_low; Double_t ly_high;
      lx_low = 0.45; ly_low = 0.7; lx_high = 0.9; ly_high = 0.89;
      auto legendInfo = new TLegend(lx_low,ly_low,lx_high,ly_high);
      TLegendEntry *entryInfo1=legendInfo->AddEntry("collisionSystem",Form("0-10%s %s, #sqrt{s} = 5.5 TeV","%",collSystem.Data()),"");
      TLegendEntry *entryInfo2;
      // TLegendEntry *entryInfo3;
      if(BField == 0.2){
        if(ith_PIDscenario == 1){
          entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.04 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
          // entryInfo3=legendInfo->AddEntry("PairPt" ,"#font[12]{p}_{T,ee} > 0.08 GeV/#font[12]{c}","");
        }
        if(ith_PIDscenario == 2){
          entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.08 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
          // entryInfo3=legendInfo->AddEntry("PairPt" ,"#font[12]{p}_{T,ee} > 0.16 GeV/#font[12]{c}","");
        }
      }
      else if (BField == 0.5){
        if(ith_PIDscenario == 1){
          entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.2 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
          // entryInfo3=legendInfo->AddEntry("PairPt" ,"#font[12]{p}_{T,ee} > 0.4 GeV/#font[12]{c}","");
        }
        if(ith_PIDscenario == 2){
          entryInfo2=legendInfo->AddEntry("SinglePt","#font[12]{p}_{T,e} > 0.08 GeV/#font[12]{c}, |#it{#eta}_{e}| < 1.1","");
          // entryInfo3=legendInfo->AddEntry("PairPt" ,"#font[12]{p}_{T,ee} > 0.16 GeV/#font[12]{c}","");
        }
      }
      auto legendPtee = new TLegend(lx_low,ly_low,lx_high,ly_high);
      TLegendEntry *entryInfo4=legendPtee->AddEntry("PairMass" ,"0 < #font[12]{m}_{ee} < 3.0 GeV/#font[12]{c^{2}}","");
      TLegendEntry *entryULS=legendPtee->AddEntry(proj_recULS_Ptee ,"ULS","p");
      TLegendEntry *entryLS=legendPtee->AddEntry(proj_recLS_Ptee ,"LS","p");
      legendInfo->SetBorderSize(0);
      legendPtee->SetBorderSize(0);
      legendInfo->SetFillColorAlpha(0, 0.0);
      legendPtee->SetFillColorAlpha(0, 0.0);
      legendInfo->SetTextSize(0.035);
      legendPtee->SetTextSize(0.035);
      cprojMeePteeDCA->cd(1);
      cprojMeePteeDCA->cd(1)->SetLogy();
      proj_recULS_Mee->GetYaxis()->SetTitle("Yield");
      proj_recULS_Mee->GetXaxis()->SetTitle("m_{ee} (GeV/#it{c}^{2})");
      // proj_recULS_Mee->RebinX(2);
      // proj_recLS_Mee->RebinX(2);
      proj_recULS_Mee->GetXaxis()->SetRangeUser(0.0,3.0);
      proj_recULS_Mee->Draw("hist p e1");
      proj_recLS_Mee->Draw("same hist p e1");
      legendInfo->Draw("same");
      cprojMeePteeDCA->cd(2);
      cprojMeePteeDCA->cd(2)->SetLogy();
      proj_recULS_Ptee->GetYaxis()->SetTitle("Yield");
      proj_recULS_Ptee->GetXaxis()->SetTitle("p_{T,ee} (GeV/#it{c})");
      // proj_recULS_Ptee->RebinX(2);
      proj_recULS_Ptee->GetXaxis()->SetRangeUser(0.0,4.0);
      proj_recULS_Ptee->Draw("hist p e1");
      proj_recLS_Ptee->Draw("same hist p e1");
      legendPtee->Draw("same");

      cprojMeePteeDCA->cd(3);
      cprojMeePteeDCA->cd(3)->SetLogy();
      proj_recULS_DCA->GetYaxis()->SetTitle("Yield");
      proj_recULS_DCA->GetXaxis()->SetTitle("DCA_{ee} (#sigma)");
      proj_recULS_DCA->GetXaxis()->SetRangeUser(0.0,20.0);
      proj_recULS_DCA->Draw("hist p e1");
      proj_recLS_DCA->Draw("same hist p e1");
      cprojMeePteeDCA->SaveAs(Form("./plots/projMeePteeDCAee_sce%i.png",ith_PIDscenario));

      // auto cMeePteeTH2 = new TCanvas("cMeePteeTH2","cMeePteeTH2",800,800);
      // cMeePteeTH2->SetTopMargin(0.03);
      // cMeePteeTH2->SetRightMargin(0.03);
      // cMeePteeTH2->SetLeftMargin(0.13);
      // hMPtDCA_ULS_rec->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      // hMPtDCA_ULS_rec->GetXaxis()->SetTitle("m_{ee} (GeV/#it{c}^{2})");
      // hMPtDCA_ULS_rec->SetMaximum(0.8);
      // hMPtDCA_ULS_rec->Draw("");
      // cMeePteeTH2->SaveAs("./plots/MeePtee.png");

      // auto cprojMee = new TCanvas("cprojMee","cprojMee",800,800);
      //   cprojMee->SetTopMargin(0.03);
      //   cprojMee->SetRightMargin(0.03);
      //   gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
      // auto legendULSLS = new TLegend(lx_low+0.06,ly_low-0.1,lx_high+0.06,ly_high-0.2);
      // TLegendEntry *entry5=legendULSLS->AddEntry(proj_recULS_Ptee ,"ULS","p");
      // TLegendEntry *entry6=legendULSLS->AddEntry(proj_recLS_Ptee ,"LS","p");
      //   legendULSLS->SetBorderSize(0);
      //   legendULSLS->SetFillColorAlpha(0, 0.0);
      //   legendULSLS->SetTextSize(0.035);
      //   // cprojMee->SetLogy();
      // TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
      //   pad1->SetBottomMargin(0);
      //   pad1->Draw();
      // cprojMee->cd();
      // TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
      //   pad2->SetTopMargin(0);
      //   pad2->Draw();
      //   pad2->cd();
      //   pad1->SetTopMargin(pad1->GetTopMargin()*0.5);
      //   pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
      //   pad2->SetBottomMargin(pad2->GetBottomMargin()*2.0);
      // pad1->cd();
      // pad1->SetLogy();
      // proj_recULS_Mee->GetYaxis()->SetTitle("Yield");
      // proj_recULS_Mee->GetXaxis()->SetTitle("m_{ee} (GeV/#it{c}^{2})");
      // proj_recULS_Mee->RebinX(2);
      // proj_recLS_Mee->RebinX(2);
      // proj_recULS_Mee->GetXaxis()->SetRangeUser(0.0,3.0);
      // proj_recULS_Mee->Draw("hist p e1");
      // proj_recLS_Mee->Draw("same hist p e1");
      // legendInfo->Draw("same");
      // legendULSLS->Draw("same");
      // pad2->cd();
      // TLine *line = new TLine(0,1,3,1);
      // line->SetLineStyle(2);
      // TH1F* ratioULSLS = (TH1F*) proj_recULS_Mee->Clone();
      // ratioULSLS->Divide(proj_recLS_Mee);
      // ratioULSLS->GetYaxis()->SetTitle("Ratio ULS/LS");
      // ratioULSLS->GetYaxis()->SetRangeUser(0.,3.);
      // ratioULSLS->GetXaxis()->SetLabelSize(0.08);
      // ratioULSLS->GetYaxis()->SetLabelSize(0.08);
      // ratioULSLS->GetXaxis()->SetTitleSize(0.08);
      // ratioULSLS->GetYaxis()->SetTitleSize(0.08);
      // ratioULSLS->GetXaxis()->SetTitleOffset(1.);
      // ratioULSLS->GetYaxis()->SetTitleOffset(0.5);
      // ratioULSLS->Draw("ep E1");
      // line->Draw("same");
      // cprojMee->SaveAs("./plots/projMee_ULSLS.png");
      //
      //
      // TString DocumentPathULS = "./plots/ULSee/Projections";
      // gSystem->Exec(Form("mkdir -p %s",DocumentPathULS.Data()));
      // auto cprojMee_PteeInterval = new TCanvas("cprojMee_PteeInterval","cprojMee_PteeInterval",800,800);
      // auto cprojPtee_MeeInterval = new TCanvas("cprojPtee_MeeInterval","cprojPtee_MeeInterval",800,800);
      // cprojMee_PteeInterval->SetTopMargin(0.03);      cprojPtee_MeeInterval->SetTopMargin(0.03);
      // cprojMee_PteeInterval->SetRightMargin(0.03);    cprojPtee_MeeInterval->SetRightMargin(0.03);
      // cprojMee_PteeInterval->SetLeftMargin(0.13);     cprojPtee_MeeInterval->SetLeftMargin(0.13);
      // Double_t posTextX         = .65;
      // Double_t posTextY         = .85;
      // TH1F* projX;
      // TH1F* projY;
      // for (size_t j = 0; j < vec_proj_bin_pt.size()-1;  j++) {      // loop over all pt projection intervalls
      //   for (size_t k = 0; k < vec_proj_bin_mass.size()-1; k++) {     // loop over all mass projection intervalls
      //     TLatex *pT_Intervall    = new TLatex(posTextX, posTextY    , Form("%g < #font[12]{p}_{T,ee} < %g  #frac{GeV}{#font[12]{c}}",vec_proj_bin_pt.at(j) ,vec_proj_bin_pt.at(j+1)));
      //     TLatex *mass_Intervall  = new TLatex(posTextX-0.05, posTextY-.05, Form("#font[12]{m}_{ee}  Intervall = %g - %g  #frac{GeV}{#font[12]{c^{2}}}",vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1)));
      //     SetTextSettings(pT_Intervall,textSize);
      //     SetTextSettings(mass_Intervall,textSize);
      //     Int_t startbinX = hMPtDCA_ULS_rec->GetXaxis()->FindBin(vec_proj_bin_mass[k]);    // select start bin of mass projection
      //     Int_t endbinX   = hMPtDCA_ULS_rec->GetXaxis()->FindBin(vec_proj_bin_mass[k+1]);  // select end bin of mass projection
      //     Int_t startbinY = hMPtDCA_ULS_rec->GetYaxis()->FindBin(vec_proj_bin_pt[j]);      // select start bin of pt projection
      //     Int_t endbinY   = hMPtDCA_ULS_rec->GetYaxis()->FindBin(vec_proj_bin_pt[j+1]);    // select end bin of pt projection
      //     projX = (TH1F*) hMPtDCA_ULS_rec->ProjectionX(Form("%s_ProjMass%g:%g_pt%g:%g",hMPtDCA_ULS_rec->GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram
      //     // projY = (TH1F*) hMPtDCA_ULS_rec->ProjectionY(Form("%s_ProjPt%g:%g_mass%g:%g",hMPtDCA_ULS_rec->GetName(),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1)),startbinX,endbinX)->Clone();  // Pt projection Histogram
      //     projY = (TH1F*) hMPtDCA_ULS_rec->ProjectionY(Form("%s_ProjPt_mass%g:%g",hMPtDCA_ULS_rec->GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1)),startbinX,endbinX)->Clone();  // Pt projection Histogram
      //     projX->RebinX(10);
      //     projY->RebinX(2);
      //     makeHistNice(projX,kBlue+1);
      //     makeHistNice(projY,kBlue+1);
      //     cprojMee_PteeInterval->cd();
      //     cprojMee_PteeInterval->cd()->SetLogy();
      //     projX->GetYaxis()->SetTitle("Yield");
      //     projX->Draw("hist p e1");
      //     pT_Intervall->Draw("same");
      //     cprojPtee_MeeInterval->cd();
      //     cprojPtee_MeeInterval->cd()->SetLogy();
      //     projY->GetYaxis()->SetTitle("Yield");
      //     projY->GetXaxis()->SetRangeUser(0,10);
      //     projY->Draw("hist p e1");
      //     mass_Intervall->Draw("same");
      //     cprojMee_PteeInterval->SaveAs(Form("./plots/ULSee/Projections/%s.png",projX->GetName()));
      //     cprojPtee_MeeInterval->SaveAs(Form("./plots/ULSee/Projections/%s.png",projY->GetName()));
      //   }
      // }
    }



    if (bPlotLS) {
      //LS++ + LS--
      // auto cLS_Mee_primCCBB = new TCanvas("cLS_Mee_primCCBB","cLS_Mee_primCCBB",800,800);
      // cLS_Mee_primCCBB->SetLogy();
      // cLS_Mee_primCCBB->SetTopMargin(0.03);
      // cLS_Mee_primCCBB->SetRightMargin(0.03);
      // cLS_Mee_primCCBB->SetLeftMargin(0.13);
      // proj_genLS_MeePrim->GetYaxis()->SetTitle("counts");
      // proj_genLS_MeePrim->Draw("hist p e1");
      // proj_genLS_MeeCC->Draw("same hist p e1");
      // proj_genLS_MeeBB->Draw("same hist p e1");
      // legPair->Draw("same");
      // cLS_Mee_primCCBB->SaveAs("./plots/bgkLS_Mee_primCCBB.png");

      // auto cLS_Ptee_primCCBB = new TCanvas("cLS_Ptee_primCCBB","cLS_Ptee_primCCBB",800,800);
      // cLS_Ptee_primCCBB->SetLogy();
      // cLS_Ptee_primCCBB->SetTopMargin(0.03);
      // cLS_Ptee_primCCBB->SetRightMargin(0.03);
      // cLS_Ptee_primCCBB->SetLeftMargin(0.13);
      // proj_genLS_PteePrim->GetYaxis()->SetTitle("counts");
      // proj_genLS_PteePrim->SetMaximum(0.8);
      // proj_genLS_PteePrim->Draw("hist p e1");
      // proj_genLS_PteeCC->Draw("same hist p e1");
      // proj_genLS_PteeBB->Draw("same hist p e1");
      // legPair->Draw("same");
      // cLS_Ptee_primCCBB->SaveAs("./plots/bgkLS_Ptee_primCCBB.png");
    }
  }


//   TString nameEffRootFile = inputFile.Data();
//   TString endDirName = "5.6.21_largeStat";
//   if (inputFile.Contains("data")) nameEffRootFile.ReplaceAll("./data/prod/anaEEstudy.", "");
//   if (inputFile.Contains("B2_100k")) nameEffRootFile.ReplaceAll(Form("../grid/output/B2_100k_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
//   if (inputFile.Contains("B5_100k")) nameEffRootFile.ReplaceAll(Form("../grid/output/B5_100k_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
//   if (inputFile.Contains("B2_2M_502TeV")) nameEffRootFile.ReplaceAll(Form("../grid/output/B2_2M_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
//   if (inputFile.Contains("B5_2M_502TeV")) nameEffRootFile.ReplaceAll(Form("../grid/output/B5_2M_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
//   nameEffRootFile.ReplaceAll(".root", "");
//   TFile *fOut = TFile::Open(Form("./data/TrackEff_%s_PIDscenario%i.root",nameEffRootFile.Data(),ith_PIDscenario),"RECREATE");
//
//   if (bPlotEfficiency) {
//     ptEffEle->SetTitle("eff_singleElectrons_Rec/GenSmeared");
//     ptEffPos->SetTitle("eff_singlePositrons_Rec/GenSmeared");
//     ptEffEle->GetXaxis()->SetRangeUser(0.0,10.0);
//     ptEffPos->GetXaxis()->SetRangeUser(0.0,10.0);
//     ptEffEle->Write();
//     ptEffPos->Write();
//     ptEffElePosGen->SetTitle("eff_Track_Rec/Generated");
//     ptEffElePosGenSmeared->SetTitle("eff_Track_Rec/GeneratedSmeared");
//     ptEffElePosGen->GetXaxis()->SetRangeUser(0.0,10.0);
//     ptEffElePosGenSmeared->GetXaxis()->SetRangeUser(0.0,10.0);
//     ptEffElePosGen->Write();
//     ptEffElePosGenSmeared->Write();
//   }
//
//   if (bPlotTrackContamination) {
//     hTotalPureContaminationRecPt->SetTitle("Total Contamination");
//     hTotalPureContaminationRecPt->GetXaxis()->SetRangeUser(0.0,4.0);
//     hTotalPureContaminationRecPt->Write();
//   }
//   fOut->Close();


}
