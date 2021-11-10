// use bools to select which  histograms are upperIntEdge

bool bPlotEfficiency = kTRUE;
bool bPlotPIDhistograms = kFALSE;
  bool bPlotNSigmaProjections = kTRUE;
bool bPlotTrackContamination = kTRUE;
bool bPlotLFHFcontributions = kTRUE;
bool bElePionRatios = kTRUE;

bool bPlotPairHistograms = kTRUE;
  bool bPlotULS = kTRUE;
  bool bPlotLS = kTRUE;

bool bUseAdditionalFiles = kTRUE;
  bool bReadPionPt = kTRUE;
  bool bCocktailFile = kTRUE;
  bool bReadSingleElectronSepc = kTRUE;

bool bwriteFileLFHFweights = kFALSE;
bool bGenerateBackground = kTRUE;

int ith_PIDscenario = 1;
int nPIDscenarios = 3;
bool bSemiCentral = true;
Bool_t DCAcut = true;

// TString strPIDscenario[] = {"TOF", "TOF+RICH (3#sigma_{#pi}^{RICH} rej)", "TOF+RICH (3.5#sigma_{#pi}^{RICH} rej)", "TOF+RICH (4#sigma_{#pi}^{RICH} rej)"};
// TString strPIDscenario[] = {"0.08 < #it{p}_{T,e}, TOF only", "0.0 < #it{p}_{T,e}, TOF only"};
// TString strPIDscenario[] = {"0.04 < #it{p}_{T,e}, iTOF only", "0.0 < #it{p}_{T,e}, iTOF only"};
// TString strPIDscenario[] = {"0.2 < #it{p}_{T,e}, TOF", "0.2 < #it{p}_{T,e}, RICH", "0.2 < #it{p}_{T,e}, TOF+RICH"};
TString strPIDscenario[] = {"0.2 < #it{p}_{T,e}, oTOF+RICH woPF", "0.2 < #it{p}_{T,e}, oTOF+RICH wPF", "0.2 < #it{p}_{T,e}, TOF+RICH wPF_iTOF","0.2 < #it{p}_{T,e}, TOF+RICH wPF_iTOF","0.2 < #it{p}_{T,e}, TOF+RICH wPF_iTOF"};
// TString strPIDscenario[] = {"0.2 < #it{p}_{T,e}, TOF+RICH", "0.2 < #it{p}_{T,e}, TOF+PS", "0.2 < #it{p}_{T,e}, TOF+RICH+PS","0.2 < #it{p}_{T,e}, TOF+RICH", "0.2 < #it{p}_{T,e}, TOF+PS", "0.2 < #it{p}_{T,e}, TOF+RICH+PS"};
// TString strPIDscenario[] = {"0.2 < #it{p}_{T,e}, TOF", "0.2 < #it{p}_{T,e}, PS", "0.2 < #it{p}_{T,e}, TOF+PS", "0.2 < #it{p}_{T,e}, TOF+PS"};
// TString strPIDscenario[] = {"0.2 < #it{p}_{T,e}, PS"};
// TString strPIDscenario[] = {"0.08 < #it{p}_{T,e}, iTOF only"};
// TString strPIDscenario[] = {"TOF+RICH (4#sigma_{#pi} 0.2<pte)", "TOF+RICH (4#sigma_{#pi} 0.08<pte)"};

Double_t etaCutsVec08[] = {0.8,0.8,0.8,0.8,0.8,0.8};
// Double_t etaCutsVec[] = {0.8,1.25,1.75,2.5,3.0,4.0};
Double_t etaCutsVec175[] = {1.75,1.75,1.75,1.75,1.75};

std::vector<Double_t> vec_proj_bin_p = {0.0, 0.3, 0.5, 0.7, 1.0, 2.0, 4.0, 10.0};
std::vector<Double_t> vec_proj_bin_pt = {0.0, 0.3, 0.5, 0.7, 1.0, 2.0, 4.0, 10.0};
std::vector<Double_t> vec_proj_bin_mass = {0.0, 3.0}; // Intervalls for projection in mass slices

// Double_t pt_bin_proj[]  = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0};
// Double_t pt_bin_proj[]  = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,4.0};
// Double_t pt_bin_proj[]     = {0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0};
// Double_t pt_bin_proj_10[]  = {0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.};
// Double_t pte_bin_2D[]      = {0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,5.,6.,7.,8.,9.,10.};
// Double_t pte_bin_2D[]      = {0.0,0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,5.,6.,7.,8.,9.,10.};
Double_t pte_bin_2D[]         = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.2,3.4,3.6,3.8,4.,5.,6.,7.,8.,9.,10.};
Double_t pte_bin_2D_01_10[]   = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.2,3.4,3.6,3.8,4.,5.,6.,7.,8.,9.,10.};
// Double_t pte_bin_2D[]      = {0.0,0.1,0.2,0.3,0.4,0.5,0.7,0.9,1.1,1.3,1.5,2.,2.5,3.,4.,5.,6.,7.,8.,10.};
// Double_t pt_bin_proj[]  = {0.0,0.1,0.5,1.,4};
// Int_t nbinspt_proj  = sizeof(pt_bin_proj)/sizeof(*pt_bin_proj) -1;
// Int_t nbinspt_proj_10  = sizeof(pt_bin_proj_10)/sizeof(*pt_bin_proj_10) -1;
Int_t nbinspt_proj_2D  = sizeof(pte_bin_2D)/sizeof(*pte_bin_2D) -1;
Int_t nbinspt_proj_2D_01_10  = sizeof(pte_bin_2D_01_10)/sizeof(*pte_bin_2D_01_10) -1;


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

void makeRatioNice(TH1 *hist, TString xlabel, TString ylabel){
  // TGaxis::SetMaxDigits(2);
  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetXaxis()->SetLabelSize(0.1);
  hist->GetYaxis()->SetLabelSize(0.1);
  hist->GetXaxis()->SetTitleSize(0.1);
  hist->GetYaxis()->SetTitleSize(0.1);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(0.4);
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

void ScaleBinWidth(TH1* h, Bool_t div = kTRUE){
  if(!h) return;
  for(int ibin=0; ibin<h->GetNbinsX();++ibin){
    Double_t const _width = h->GetBinWidth(ibin);
    Double_t const _val = h->GetBinContent(ibin);
    Double_t const _error = h->GetBinError(ibin);
    if(_width <= 0.) continue;
    if ( div ) {
      h->SetBinContent(ibin, _val / _width);
      h->SetBinError(ibin, _error / _width);
    } else {
      h->SetBinContent(ibin, _val * _width);
      h->SetBinError(ibin, _error * _width);
    }
  }
}

void invariantYield(TH1* h){
  double pt = 0;
  for (int i = 1;i <= h->GetNbinsX(); ++i)
  {
    pt = h->GetBinCenter(i);
    h->SetBinContent(i, h->GetBinContent(i)/pt);
  }
  double fac = 1/(2*TMath::Pi());
  h->Scale(fac);
}

void reverseInvariantYield(TH1* h){
  double pt = 0;
  for (int i = 1;i <= h->GetNbinsX(); ++i)
  {
    pt = h->GetBinCenter(i);
    h->SetBinContent(i, h->GetBinContent(i)*pt);
  }
  double fac = (2*TMath::Pi());
  h->Scale(fac);
}

//__________________________________________________
Double_t ComputeIntegral(TH1 *g,Double_t low,Double_t high) {
  Double_t step = (high-low)/1000;
  Double_t sum = 0.;

  for(Int_t k=0; k < 1000; k++){
    Double_t xb = low + step/2. + k*step;
    Double_t evalvalue = g->Interpolate(xb);
    sum += evalvalue*step;
  }
  return sum;
}

TH1F *RatioToTheoryy(TH1 *ffit, TH1 *ginput) {
  //
  // Ratio of TGraphAsymmErrors to fit function
  //


  TH1F *graphh = (TH1F *) ginput->Clone(Form("Ratio_%s",ginput->GetName()));

  Int_t npoints = graphh->GetNbinsX();
  TAxis *ptaxis = ginput->GetXaxis();

  for(Int_t k=0; k < npoints; k++) {

    Double_t ex_low = ptaxis->GetBinLowEdge(k+1);
    Double_t ex_high = ptaxis->GetBinUpEdge(k+1);
    Double_t error = ginput->GetBinError(k+1);

    Double_t content = ginput->GetBinContent(k+1);
    Double_t y = ginput->GetBinContent(k+1);

    //Double_t evalfit = ffit->Eval(x);
    Double_t evalfit = ComputeIntegral(ffit,ex_low,ex_high)/(ex_high-ex_low);


    graphh->SetBinContent(k+1,y/evalfit);
    graphh->SetBinError(k+1,error/evalfit);

  }

  return graphh;

}


void ClearBinsAbove(TH1* h, Double_t maxpt) {
  Int_t maxptbin = h->GetXaxis()->FindBin(maxpt);
  for (Int_t iBin = maxptbin; iBin < h->GetNbinsX()+1; iBin++) {
    h->SetBinContent(iBin,0.);
  }
}

TH1F* AddHistosDiffRange(TH1F* h1, TH1F* h2) {
  for (Int_t i = 0; i <= h1->GetNbinsX(); i++) {
    Double_t pt = h1->GetXaxis()->GetBinCenter(i);
    Int_t ptbinh2 = h2->GetXaxis()->FindBin(pt);
    // cout << "bin h1 = " << i <<", pt h1 = " << pt << ", bin h2 = " << ptbinh2 << ", pt h2 = " <<  h2->GetBinContent(ptbinh2) << endl;
    if(ptbinh2 < 1) continue;
    double sum = h2->GetBinContent(ptbinh2) + h1->GetBinContent(i);
    h1->SetBinContent(i, sum);
  }
  return h1;
}

//__________________________________________________


void anaPlots(TString inputFile)
{
  TFile *fIn  = TFile::Open(inputFile.Data());
  for (Int_t ith_PIDscenario = 1; ith_PIDscenario <= nPIDscenarios; ith_PIDscenario++) {


  // select which PID cenario (written in the input file) should be analyzed
  TString pathRecPIDScenario = Form("reconstructed/PIDscenario_%i", ith_PIDscenario);
  TString pathGenPIDScenario = Form("generated/PIDscenario_%i", ith_PIDscenario);

  double BField;
  TString collSystem;
  Bool_t findBField;
  Bool_t findCollSystem;
  findBField = inputFile.Contains("B=2kG");
  if (findBField) BField = 0.2;
  findBField = inputFile.Contains("B=5kG");
  if (findBField) BField = 0.5;
  findBField = inputFile.Contains("B=20kG");
  if (findBField) BField = 2.0;
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


    // only used if additional files need to be read in
    TFile *fReadChPi;
    TH1F* hPaper_chPi_Pt_firstCentralSection;
    TH1F* hPaper_chPi_Pt_secondCentralSection;
    TH1F* hPaper_chPi_Pt_fullCentralSection;
    TH1F* hPaper_chPi_Pt_firstCentralSection_stat_Err;
    TH1F* hPaper_chPi_Pt_secondCentralSection_stat_Err;
    TH1F* hPaper_chPi_Pt_fullCentralSection_stat_Err;
    TH1F* hPaper_chPi_Pt_firstCentralSection_sys_Err;
    TH1F* hPaper_chPi_Pt_secondCentralSection_sys_Err;
    TH1F* hPaper_chPi_Pt_fullCentralSection_sys_Err;
    TH1F* hPaper_chPi_Pt_firstCentralSection_sys_unc;
    TH1F* hPaper_chPi_Pt_secondCentralSection_sys_unc;
    TH1F* hPaper_chPi_Pt_fullCentralSection_sys_unc;
    TH1F* hCombStatErr;
    TH1F* hCombSysErr;
    TH1F* hCombSysUnc;
    TH1F* hCombStatErr_rebin;
    TH1F* hCombSysErr_rebin;
    TH1F* hCombSysUnc_rebin;

    TFile *fReadSingleESpec;
    TH1F* hSingleElefromHF_Cocktail;

    TFile *fReadHEP_paper_HFtoE;
    TH1F* hPaper_HFtoe_Pt;
    TH1F* hPaper_HFtoe_Pt_stat;
    TH1F* hPaper_HFtoe_Pt_syst;


    if (bUseAdditionalFiles) {
      if(bReadPionPt){
        fReadChPi = TFile::Open("/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/efficiency/data/PtePi_PbPb_HEPData-ins1759506-v1-Table_1.root");
        Int_t centralSection = 1; // 0-5% central   x+1 = 5-10%
        if (bSemiCentral) centralSection = 5;  // 30-40% central   x+1 = 40-50%
        hPaper_chPi_Pt_firstCentralSection = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i",centralSection));
        hPaper_chPi_Pt_firstCentralSection_stat_Err = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i_e1",centralSection));
        hPaper_chPi_Pt_firstCentralSection_sys_Err = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i_e2",centralSection));
        hPaper_chPi_Pt_firstCentralSection_sys_unc = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i_e3",centralSection));
        hPaper_chPi_Pt_secondCentralSection = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i",centralSection+1));
        hPaper_chPi_Pt_secondCentralSection_stat_Err = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i_e1",centralSection+1));
        hPaper_chPi_Pt_secondCentralSection_sys_Err = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i_e2",centralSection+1));
        hPaper_chPi_Pt_secondCentralSection_sys_unc = (TH1F*) fReadChPi->Get(Form("Table 1/Hist1D_y%i_e3",centralSection+1));

        hPaper_chPi_Pt_fullCentralSection = (TH1F*) hPaper_chPi_Pt_firstCentralSection->Clone();
        hCombStatErr = (TH1F*) hPaper_chPi_Pt_firstCentralSection_stat_Err->Clone();
        hCombSysErr = (TH1F*) hPaper_chPi_Pt_firstCentralSection_sys_Err->Clone();
        hCombSysUnc = (TH1F*) hPaper_chPi_Pt_firstCentralSection_sys_unc->Clone();
        hPaper_chPi_Pt_fullCentralSection->Add(hPaper_chPi_Pt_secondCentralSection,1);
        hCombStatErr->Add(hPaper_chPi_Pt_secondCentralSection_stat_Err,1);
        hCombSysErr->Add(hPaper_chPi_Pt_secondCentralSection_sys_Err,1);
        hCombSysUnc->Add(hPaper_chPi_Pt_secondCentralSection_sys_unc,1);

        hPaper_chPi_Pt_fullCentralSection->Scale(1./2.);
        hCombStatErr->Scale(1./2.);
        hCombSysErr->Scale(1./2.);
        hCombSysUnc->Scale(1./2.);

        ScaleBinWidth(hCombStatErr, kFALSE);
        ScaleBinWidth(hCombSysErr, kFALSE);
        ScaleBinWidth(hCombSysUnc, kFALSE);
        ScaleBinWidth(hPaper_chPi_Pt_fullCentralSection, kFALSE);

        hCombStatErr_rebin = (TH1F*) hCombStatErr->Rebin(nbinspt_proj_2D_01_10,"hCombStatErr_rebin",&pte_bin_2D_01_10[0]);
        hCombSysErr_rebin = (TH1F*) hCombSysErr->Rebin(nbinspt_proj_2D_01_10,"hCombSysErr_rebin",&pte_bin_2D_01_10[0]);
        hCombSysUnc_rebin = (TH1F*) hCombSysUnc->Rebin(nbinspt_proj_2D_01_10,"hCombSysUnc_rebin",&pte_bin_2D_01_10[0]);
        // hPaper_chPi_Pt_fullCentralSection_stat_Err = (TH1F*) hPaper_chPi_Pt_fullCentralSection->Clone();
        // hPaper_chPi_Pt_fullCentralSection_sys_Err = (TH1F*) hPaper_chPi_Pt_fullCentralSection->Clone();
        // hPaper_chPi_Pt_fullCentralSection_sys_unc = (TH1F*) hPaper_chPi_Pt_fullCentralSection->Clone();
        hPaper_chPi_Pt_fullCentralSection_stat_Err = (TH1F*) hPaper_chPi_Pt_fullCentralSection->Rebin(nbinspt_proj_2D_01_10,"hPaper_chPi_Pt_fullCentralSection_stat_Err",&pte_bin_2D_01_10[0]);
        hPaper_chPi_Pt_fullCentralSection_sys_Err = (TH1F*) hPaper_chPi_Pt_fullCentralSection->Rebin(nbinspt_proj_2D_01_10,"hPaper_chPi_Pt_fullCentralSection_sys_Err",&pte_bin_2D_01_10[0]);
        hPaper_chPi_Pt_fullCentralSection_sys_unc = (TH1F*) hPaper_chPi_Pt_fullCentralSection->Rebin(nbinspt_proj_2D_01_10,"hPaper_chPi_Pt_fullCentralSection_sys_unc",&pte_bin_2D_01_10[0]);

        ScaleBinWidth(hCombStatErr_rebin, kTRUE);
        ScaleBinWidth(hCombSysErr_rebin, kTRUE);
        ScaleBinWidth(hCombSysUnc_rebin, kTRUE);
        ScaleBinWidth(hPaper_chPi_Pt_fullCentralSection_stat_Err, kTRUE);
        ScaleBinWidth(hPaper_chPi_Pt_fullCentralSection_sys_Err, kTRUE);
        ScaleBinWidth(hPaper_chPi_Pt_fullCentralSection_sys_unc, kTRUE);

        for (size_t iBin = 1; iBin < hPaper_chPi_Pt_fullCentralSection_stat_Err->GetNbinsX(); iBin++) {
          double iBinStatErr = hCombStatErr_rebin->GetBinContent(iBin);
          double iBinSysErr = hCombSysErr_rebin->GetBinContent(iBin);
          double iBinSysUnc = hCombSysUnc_rebin->GetBinContent(iBin);
          hPaper_chPi_Pt_fullCentralSection_stat_Err->SetBinError(iBin,iBinStatErr);
          hPaper_chPi_Pt_fullCentralSection_sys_Err->SetBinError(iBin,iBinSysErr);
          hPaper_chPi_Pt_fullCentralSection_sys_unc->SetBinError(iBin,iBinSysUnc);
        }

        makeHistNice(hPaper_chPi_Pt_fullCentralSection_stat_Err,kBlack);
        makeHistNice(hPaper_chPi_Pt_fullCentralSection_sys_Err,kBlack);
        makeHistNice(hPaper_chPi_Pt_fullCentralSection_sys_unc,kBlack);
        hPaper_chPi_Pt_fullCentralSection_sys_unc->SetFillColor(kGray+1);
      }


      if(bReadSingleElectronSepc){
        fReadSingleESpec = TFile::Open("/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/efficiency/data/Single_electron_spectrum.root");
         TH1F* hbeauty_bTOe = (TH1F*) fReadSingleESpec->Get("eformbeauty_be");
         TH1F* hbeauty_bTOcTOe = (TH1F*) fReadSingleESpec->Get("eformbeauty_bce");
         TH1F* hbeauty_cTOe = (TH1F*) fReadSingleESpec->Get("eformCharm_EPS09");

        hSingleElefromHF_Cocktail = (TH1F*) hbeauty_bTOe->Clone();
        hSingleElefromHF_Cocktail->Add(hbeauty_bTOcTOe);
//        hSingleElefromHF_Cocktail->Reset();
        hSingleElefromHF_Cocktail->Add(hbeauty_cTOe);

        // contribution from cocktail already should be e+ + e-
        // multiply the cocktil also by a factor of 2 (like the HFtoe data below) to have a smooth transition.
        // use it as an extrapolation, and not as a given expectation
         hSingleElefromHF_Cocktail->Scale(2.);

        makeHistNice(hSingleElefromHF_Cocktail,kRed);
        hSingleElefromHF_Cocktail->Scale(1./1.6); //normalize for eta
        if(bSemiCentral) hSingleElefromHF_Cocktail->Scale(250.8/1444.); // scale for 30-50 central events Ncoll_semicentral/Ncoll_central
//        ScaleBinWidth(hSingleElefromHF_Cocktail,kTRUE);
        invariantYield(hSingleElefromHF_Cocktail);


        // paper: yield of electrons from HF (c,b -> e)
        TString sPtHFtoeData = "/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/efficiency/data/PtHFtoe_PbPb_HEPData-ins1759860-v1-Table_2.1.root";
        TString fTable = "Table 2.1";
        if (bSemiCentral) {
          sPtHFtoeData = "/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/efficiency/data/PtHFtoe_PbPb_HEPData-ins1759860-v1-Table_2.2_30-50central.root";
          fTable = "Table 2.2";
        }
        fReadHEP_paper_HFtoE = TFile::Open(sPtHFtoeData.Data());
        hPaper_HFtoe_Pt = (TH1F*) fReadHEP_paper_HFtoE->Get(Form("%s/Hist1D_y1",fTable.Data()));
        TH1F* hPaper_stat_Err = (TH1F*) fReadHEP_paper_HFtoE->Get(Form("%s/Hist1D_y1_e1",fTable.Data()));
        TH1F* hPaper_sys_Err = (TH1F*) fReadHEP_paper_HFtoE->Get(Form("%s/Hist1D_y1_e2",fTable.Data()));
        hPaper_stat_Err->Scale(2.);             // The data is (e+ + e-)/2 to only have electron contribution. Multiply by two to have e+ + e-
        hPaper_sys_Err->Scale(2.);              // The data is (e+ + e-)/2 to only have electron contribution. Multiply by two to have e+ + e-
        hPaper_HFtoe_Pt->Scale(2.);             // The data is (e+ + e-)/2 to only have electron contribution. Multiply by two to have e+ + e-

        hPaper_HFtoe_Pt_stat = (TH1F*) hPaper_HFtoe_Pt->Clone();
        hPaper_HFtoe_Pt_syst = (TH1F*) hPaper_HFtoe_Pt->Clone();
        for (size_t iBin = 1; iBin < hPaper_HFtoe_Pt->GetNbinsX(); iBin++) {
          double iBinStatErr = hPaper_stat_Err->GetBinContent(iBin);
          double iBinSysErr = hPaper_sys_Err->GetBinContent(iBin);
          hPaper_HFtoe_Pt_stat->SetBinError(iBin,iBinStatErr);
          hPaper_HFtoe_Pt_syst->SetBinError(iBin,iBinSysErr);

          makeHistNice(hPaper_HFtoe_Pt_stat,kBlack);
          makeHistNice(hPaper_HFtoe_Pt_syst,kBlack);
          hPaper_HFtoe_Pt_syst->SetFillColor(kGray+1);
        }

      }


    }


  // Track histograms

  TH3F* hRec_Track_Pt_Eta_Phi_BeforeSmearing = (TH3F*) fIn->Get("hBeforeSmearing_Pt_Eta_Phi_rec");
  TH3F* hRec_Track_Pt_Eta_Phi_AfterSmearing = (TH3F*) fIn->Get("hAfterSmearing_Pt_Eta_Phi_rec");

  TH3F* hRec_Track_Pt_Eta_Phi_AfterKineCuts = (TH3F*) fIn->Get(Form("%s/hAfterKineCuts_Pt_Eta_Phi_rec_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Rec_Pt_Eta_Phi_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_ElePos_Pi0_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get("hTrack_ElePos_Pi0_Rec_Pt_Eta_Phi_beforePID");
  // TH3F* hRec_Track_ElePos_LF_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get("hTrack_ElePos_LF_Rec_Pt_Eta_Phi_beforePID");
  // TH3F* hRec_Track_ElePos_HF_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get("hTrack_ElePos_HF_Rec_Pt_Eta_Phi_beforePID");

  // TH3F* hRec_Track_ElePos_Pi0_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Pi0_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_ElePos_LF_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_LF_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_ElePos_HF_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_HF_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_ElePos_Pi0_Pt_Eta_Phi_afterPID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Pi0_Rec_Pt_Eta_Phi_afterPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_ElePos_LF_Pt_Eta_Phi_afterPID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_LF_Rec_Pt_Eta_Phi_afterPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_ElePos_HF_Pt_Eta_Phi_afterPID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_HF_Rec_Pt_Eta_Phi_afterPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH3F* hRec_Track_ElePos_Pi0_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_LF_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_HF_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_Pi0_Pt_Eta_Phi_afterPID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_MotherPi0_Rec_Pt_Eta_Phi_afterPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_LF_Pt_Eta_Phi_afterPID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_MotherLF_Rec_Pt_Eta_Phi_afterPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_HF_Pt_Eta_Phi_afterPID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_MotherHF_Rec_Pt_Eta_Phi_afterPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_ElePos_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_Muon_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_Muon_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hRec_Track_Pion_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_Pion_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_Kaon_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_Kaon_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_Proton_Pt_Eta_Phi_beforePID = (TH3F*) fIn->Get(Form("%s/hTrack_Proton_Rec_Pt_Eta_Phi_beforePID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_Pion_Pt_Eta_Phi_afterPSPID = (TH3F*) fIn->Get(Form("%s/hTrack_Pion_Rec_Pt_Eta_Phi_afterPSPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hRec_Track_ElePos_Pt_Eta_Phi_afterPSPID = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Rec_Pt_Eta_Phi_afterPSPID_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
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
  TH3F* hGen_Track_Pion_Pt_Rap_Phi = (TH3F*) fIn->Get("hTrack_Pion_Gen_Pt_Rap_Phi");
  // TH3F* hGen_Track_Pion_Pt_Rap_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Pion_Gen_Pt_Rap_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hGen_Track_Kaon_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Kaon_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hGen_Track_Proton_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_Proton_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_ElePos_LF_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_LF_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  TH3F* hGen_Track_ElePos_HF_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_HF_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hGen_Track_ElePos_MotherPi0_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_MotherPi0_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH3F* hGen_Track_ElePos_MotherPi0_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_Pi0_Gen_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));

  // TH3F* hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts = (TH3F*) fIn->Get("hTrack_ElePos_GenSmeared_Pt_Eta_Phi_beforeKineCuts");
  TH3F* hGenSmeared_Track_ElePos_Pt_Eta_Phi = (TH3F*) fIn->Get(Form("%s/hTrack_ElePos_GenSmeared_Pt_Eta_Phi_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));

  // TH1F* hTracks_ElePos_Gen_inPFvector_Pt = (TH1F*) fIn->Get(Form("%s/hTracks_ElePos_Gen_inPFvector_Pt_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));
  // TH1F* hTracks_ElePos_Rec_inPFvector_Pt = (TH1F*) fIn->Get(Form("%s/hTracks_ElePos_Rec_inPFvector_Pt_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH2F* hMeeOpAngle_RecULS_bPF = (TH2F*) fIn->Get(Form("%s/hMeeOpAngle_RecULS_bPF_pi0Dalitz_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH2F* hMeeOpAngle_RecULS_aPF = (TH2F*) fIn->Get(Form("%s/hMeeOpAngle_RecULS_aPF_pi0Dalitz_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH2F* hMeeOpAngle_RecULS_bPF = (TH2F*) fIn->Get(Form("%s/hMeeOpAngle_RecULS_bPF_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  // TH2F* hMeeOpAngle_RecULS_aPF = (TH2F*) fIn->Get(Form("%s/hMeeOpAngle_RecULS_aPF_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TF1* funcPSparam = new TF1("funcPSparam","0.8*(1.-exp(-1.6*(x-0.05)))",0,4);


  std::vector<TH2F*> vecTOF_PIDplots;
  std::vector<TH2F*> vecRICH_PIDplots;
  TH2F* hRec_TOF_BetaP_beforeSmearing;
  TH2F* hRec_TOF_BetaP_afterSmearing;
  TH2F* hRec_RICH_CherenkovAngleP_beforeSmearing;
  TH2F* hRec_RICH_CherenkovAngleP_afterSmearing;

  TH2 *hNsigmaP_TOF_trueElec[5];
  TH2 *hNsigmaP_TOF_trueMuon[5];
  TH2 *hNsigmaP_TOF_truePion[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueElec[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_trueMuon[5];
  TH2 *hNsigmaP_afterPIDcuts_TOF_truePion[5];

  TH2 *hNsigmaP_RICH[5];
  TH2 *hNsigmaP_RICH_trueElec[5];
  TH2 *hNsigmaP_RICH_trueMuon[5];
  TH2 *hNsigmaP_RICH_truePion[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueElec[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_trueMuon[5];
  TH2 *hNsigmaP_afterPIDcuts_RICH_truePion[5];


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
    // hNsigmaP_afterPIDcuts_TOF[i] = (TH2F*) fIn->Get(Form("hNsigmaP_%s_afterPIDcuts_TOF", pname[i]));
    hNsigmaP_afterPIDcuts_TOF_trueElec[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_trueElec_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
    hNsigmaP_afterPIDcuts_TOF_trueMuon[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_trueMuon_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
    hNsigmaP_afterPIDcuts_TOF_truePion[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_TOF_truePion_sce%i", pathRecPIDScenario.Data(),pname[i],ith_PIDscenario));
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
    hNsigmaP_afterPIDcuts_RICH[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_trueElec[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_trueElec_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_trueMuon[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_trueMuon_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
    hNsigmaP_afterPIDcuts_RICH_truePion[i] = (TH2F*) fIn->Get(Form("%s/hNsigmaP_%s_afterPIDcuts_RICH_truePion_sce%i",pathRecPIDScenario.Data(), pname[i], ith_PIDscenario));
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

  //Pair histograms M,pT
  TH3F* hMPtDCA_ULS_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH3F* hMPtDCA_ULS_rec_MCpidEle = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_ULS_rec_MCpidEle_charmTOe = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_charmTOe_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_ULS_rec_MCpidEle_beautyTOe = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_beautyTOe_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_ULS_rec_MCpidEle_hfTOe = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_hfTOe_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_ULS_rec_MCpidEle_lfTOee = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_lfTOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_ULS_rec_MCpidEle_ccTOee = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_ccTOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_ULS_rec_MCpidEle_bbTOee = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_bbTOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_ULS_rec_MCpidEle_hfTOee = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_rec_MCpidEle_hfTOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH3F* hMPtDCA_ULS_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_ULS_gen_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));


  TH3F* hMPtDCA_LS_rec = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_MCpidEle = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_MCpidEle_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_lfTOe    = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_lfTOe_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_charmTOe  = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_charmTOe_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_beautyTOe = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_beautyTOe_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_hfTOe     = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_hfTOe_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_Pi0TOee   = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_Pi0TOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_lfTOee_selectPDG = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_lfTOee_selectPDG_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_hfTOee = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_hfTOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_ccTOee = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_ccTOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));
  TH3F* hMPtDCA_LS_rec_bbTOee = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_rec_bbTOee_sce%i",pathRecPIDScenario.Data(),ith_PIDscenario));

  TH3F* hMPtDCA_LS_gen = (TH3F*) fIn->Get(Form("%s/hMPtDCA_LS_gen_sce%i",pathGenPIDScenario.Data(),ith_PIDscenario));


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

  TH1F* ptRecTrackElePos_Pi0_beforePID = (TH1F*) hRec_Track_ElePos_Pi0_Pt_Eta_Phi_beforePID->ProjectionX("Rec ElePos Pi0 beforePIDCuts ptTrack");
  TH1F* ptRecTrackElePos_LF_beforePID  = (TH1F*) hRec_Track_ElePos_LF_Pt_Eta_Phi_beforePID->ProjectionX("Rec ElePos LF beforePIDCuts ptTrack");
  TH1F* ptRecTrackElePos_HF_beforePID  = (TH1F*) hRec_Track_ElePos_HF_Pt_Eta_Phi_beforePID->ProjectionX("Rec ElePos HF beforePIDCuts ptTrack");

  TH1F* ptRecTrackElePos_Pi0_afterPID = (TH1F*) hRec_Track_ElePos_Pi0_Pt_Eta_Phi_afterPID->ProjectionX("Rec ElePos Pi0 afterPIDCuts ptTrack");
  TH1F* ptRecTrackElePos_LF_afterPID  = (TH1F*) hRec_Track_ElePos_LF_Pt_Eta_Phi_afterPID->ProjectionX("Rec ElePos LF afterPIDCuts ptTrack");
  TH1F* ptRecTrackElePos_HF_afterPID  = (TH1F*) hRec_Track_ElePos_HF_Pt_Eta_Phi_afterPID->ProjectionX("Rec ElePos HF afterPIDCuts ptTrack");

  TH2F* ptEtaRecTrackElePos = (TH2F*) hRec_Track_ElePos_Pt_Eta_Phi->Project3D("yx o");

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
  TH1F* ptRecProtonTrack  = (TH1F*) hRec_Track_Proton_Pt_Eta_Phi->ProjectionX("Rec Proton ptTrack");

  TH1F* ptRecTrackElePos_beforePID  = (TH1F*) hRec_Track_ElePos_Pt_Eta_Phi_beforePID->ProjectionX("Rec ElePos afterPIDCuts ptTrack before PID");
  // TH1F* ptRecTrackElePos_afterPSPID  = (TH1F*) hRec_Track_ElePos_Pt_Eta_Phi_afterPSPID->ProjectionX("Rec ElePos afterPIDCuts ptTrack after PS PID");
  // TH1F* ptRecMuonTrack_beforePID    = (TH1F*) hRec_Track_Muon_Pt_Eta_Phi_beforePID->ProjectionX("Rec Muon ptTrack before PID");
  TH1F* ptRecPionTrack_beforePID    = (TH1F*) hRec_Track_Pion_Pt_Eta_Phi_beforePID->ProjectionX("Rec Pion ptTrack before PID");
  // TH1F* ptRecKaonTrack_beforePID    = (TH1F*) hRec_Track_Kaon_Pt_Eta_Phi_beforePID->ProjectionX("Rec Kaon ptTrack before PID");
  // TH1F* ptRecProtonTrack_beforePID  = (TH1F*) hRec_Track_Proton_Pt_Eta_Phi_beforePID->ProjectionX("Rec Proton ptTrack before PID");
  // TH1F* ptRecPionTrack_afterPSPID    = (TH1F*) hRec_Track_Pion_Pt_Eta_Phi_afterPSPID->ProjectionX("Rec Pion ptTrack after PS PID");


  TH1F* ptGenTrackElePos_beforeKineCuts  = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionX("Gen ElePos ptTrack beforeKineCuts");
  TH1F* etaGenTrackElePos_beforeKineCuts = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionY("Gen ElePos etaTrack beforeKineCuts");
  TH1F* phiGenTrackElePos_beforeKineCuts = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionZ("Gen ElePos phiTrack beforeKineCuts");

  TH1F* ptGenTrackElePos  = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi->ProjectionX("Gen ElePos ptTrack");
  TH1F* etaGenTrackElePos = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi->ProjectionY("Gen ElePos etaTrack");
  TH1F* phiGenTrackElePos = (TH1F*) hGen_Track_ElePos_Pt_Eta_Phi->ProjectionZ("Gen ElePos phiTrack");

  TH2F* ptEtaGenTrackElePos = (TH2F*) hGen_Track_ElePos_Pt_Eta_Phi->Project3D("yx o");

  TH1F* ptGenTrackElePos_LF  = (TH1F*) hGen_Track_ElePos_LF_Pt_Eta_Phi->ProjectionX("Gen ElePos from LF ptTrack");
  TH1F* ptGenTrackElePos_HF  = (TH1F*) hGen_Track_ElePos_HF_Pt_Eta_Phi->ProjectionX("Gen ElePos from HF ptTrack");
  // TH1F* ptGenTrackElePos_Pi0  = (TH1F*) hGen_Track_ElePos_MotherPi0_Pt_Eta_Phi->ProjectionX("Gen ElePos from Pi0 ptTrack");


  TH1F* ptGenTrackAll   = (TH1F*) hGen_Track_All_Pt_Eta_Phi->ProjectionX("Gen All ptTrack");
  TH1F* ptGenTrackMuon  = (TH1F*) hGen_Track_Muon_Pt_Eta_Phi->ProjectionX("Gen Muon ptTrack");
  TH1F* ptGenTrackPion  = (TH1F*) hGen_Track_Pion_Pt_Eta_Phi->ProjectionX("Gen Pion ptTrack");
  // TH1F* ptGenTrackKaon  = (TH1F*) hGen_Track_Kaon_Pt_Eta_Phi->ProjectionX("Gen Kaon ptTrack");
  // TH1F* ptGenTrackProton  = (TH1F*) hGen_Track_Proton_Pt_Eta_Phi->ProjectionX("Gen Proton ptTrack");
  TH1F* ptGenTrackPionRapSel  = (TH1F*) hGen_Track_Pion_Pt_Rap_Phi->ProjectionX("Gen Pion ptTrack RapSel");



  // TH1F* ptGenSmearedTrackElePos_beforeKineCuts  = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionX("GenSmeared ElePos ptTrack beforeKineCuts");
  // TH1F* etaGenSmearedTrackElePos_beforeKineCuts = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionY("GenSmeared ElePos etaTrack beforeKineCuts");
  // TH1F* phiGenSmearedTrackElePos_beforeKineCuts = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi_beforeKineCuts->ProjectionZ("GenSmeared ElePos phiTrack beforeKineCuts");

  TH1F* ptGenSmearedTrackElePos  = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi->ProjectionX("GenSmeared ElePos ptTrack");
  TH1F* etaGenSmearedTrackElePos = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi->ProjectionY("GenSmeared ElePos etaTrack");
  TH1F* phiGenSmearedTrackElePos = (TH1F*) hGenSmeared_Track_ElePos_Pt_Eta_Phi->ProjectionZ("GenSmeared ElePos phiTrack");


    TH1F* hParticlesFIT = (TH1F*) fIn->Get("nParticlesFIT");
    Int_t nentries=0;
    Int_t nbins = hParticlesFIT->GetXaxis()->GetNbins();
    Int_t centralNinetyPercent = hParticlesFIT->GetEntries() * 0.9;
    Int_t centralsixtyPercent = hParticlesFIT->GetEntries() * 0.6;
    Int_t centralfortyPercent = hParticlesFIT->GetEntries() * 0.4;
    cout  << nbins << endl;
    for (size_t ibin = 1; ibin < nbins ; ibin++) {
      if(nentries >= centralNinetyPercent) continue; // look for the lowest 90% to find the lower egde of the 10% cental events
      nentries = nentries + hParticlesFIT->GetBinContent(ibin);
      if(nentries >= centralNinetyPercent) cout << "nentries has reached 10% at " << hParticlesFIT->GetXaxis()->GetBinCenter(ibin) << endl;
    }
    nentries=0;
    for (size_t ibin = 1; ibin < nbins ; ibin++) {
      if(nentries >= centralsixtyPercent) continue; // look for the lowest 60% to find the lower egde of the 40% cental events
      nentries = nentries + hParticlesFIT->GetBinContent(ibin);
      if(nentries >= centralsixtyPercent) cout << "nentries has reached 40% at " << hParticlesFIT->GetXaxis()->GetBinCenter(ibin) << endl;
    }
    nentries=0;
    for (size_t ibin = 1; ibin < nbins ; ibin++) {
      if(nentries >= centralfortyPercent) continue; // look for the lowest 40% to find the lower egde of the 60% cental events
      nentries = nentries + hParticlesFIT->GetBinContent(ibin);
      if(nentries >= centralfortyPercent) cout << "nentries has reached 60% at " << hParticlesFIT->GetXaxis()->GetBinCenter(ibin) << endl;
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
    hdNdeta_midrap->GetYaxis()->SetTitle("1/#it{N}_{ev}d#it{N}_{ch}/d#eta");
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


    // TH1F* proj_recULS_MCpidEle_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle->ProjectionX("proj_recULS_MCpidEle_Mee");
    // TH1F* proj_recULS_MCpidEle_Ptee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle->ProjectionY("proj_recULS_MCpidEle_Ptee");

    TH1F* proj_recULS_MCpidEle_charmTOe_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle_charmTOe->ProjectionX("proj_recULS_Mee_charmTOe");
    TH1F* proj_recULS_MCpidEle_beautyTOe_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle_beautyTOe->ProjectionX("proj_recULS_Mee_beautyTOe");
    TH1F* proj_recULS_MCpidEle_hfTOe_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle_hfTOe->ProjectionX("proj_recULS_Mee_hfTOe");
    TH1F* proj_recULS_MCpidEle_lfTOee_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle_lfTOee->ProjectionX("proj_recULS_Mee_lfTOee");

    TH1F* proj_recULS_MCpidEle_ccTOee_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle_ccTOee->ProjectionX("proj_recULS_Mee_ccTOee");
    TH1F* proj_recULS_MCpidEle_bbTOee_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle_bbTOee->ProjectionX("proj_recULS_Mee_bbTOee");
    TH1F* proj_recULS_MCpidEle_hfTOee_Mee = (TH1F*) hMPtDCA_ULS_rec_MCpidEle_hfTOee->ProjectionX("proj_recULS_Mee_hfTOee");

    TH1F* proj_recLS_Mee = (TH1F*) hMPtDCA_LS_rec->ProjectionX("proj_recLS_Mee");
    TH1F* proj_recLS_Ptee = (TH1F*) hMPtDCA_LS_rec->ProjectionY("proj_recLS_Ptee");

    printf(" dcacut %i\n", DCAcut);


    // DCA cut on pairs
    if(DCAcut) hMPtDCA_LS_rec_MCpidEle->GetZaxis()->SetRangeUser(0.,1.2);
    if(DCAcut) hMPtDCA_LS_rec_lfTOee_selectPDG->GetZaxis()->SetRangeUser(0.,1.2);
    if(DCAcut) hMPtDCA_LS_rec_hfTOe->GetZaxis()->SetRangeUser(0.,1.2);
    if(DCAcut) hMPtDCA_LS_rec_hfTOee->GetZaxis()->SetRangeUser(0.,1.2);

    TH1F* proj_recLS_MCpidEle_Mee = (TH1F*) hMPtDCA_LS_rec_MCpidEle->Project3D("x");
    TH1F* proj_recLS_MCpidEle_Ptee = (TH1F*) hMPtDCA_LS_rec_MCpidEle->Project3D("y");
    TH1F* proj_recLS_MCpidEle_DCA = (TH1F*) hMPtDCA_LS_rec_MCpidEle->Project3D("z");
    proj_recLS_MCpidEle_Mee->SetTitle("proj_recLS_MCpidEle_Mee");
    proj_recLS_MCpidEle_Ptee->SetTitle("proj_recLS_MCpidEle_Ptee");
    proj_recLS_MCpidEle_DCA->SetTitle("proj_recLS_MCpidEle_DCA");

    TH1F* proj_recLS_Ctoe_Mee  = (TH1F*) hMPtDCA_LS_rec_charmTOe->ProjectionX("proj_recLS_Ctoe_Mee");
    TH1F* proj_recLS_Ctoe_Ptee = (TH1F*) hMPtDCA_LS_rec_charmTOe->ProjectionY("proj_recLS_Ctoe_Ptee");
    TH1F* proj_recLS_Ctoe_DCA  = (TH1F*) hMPtDCA_LS_rec_charmTOe->ProjectionZ("proj_recLS_Ctoe_DCA");

    TH1F* proj_recLS_Btoe_Mee  = (TH1F*) hMPtDCA_LS_rec_beautyTOe->ProjectionX("proj_recLS_Btoe_Mee");
    TH1F* proj_recLS_Btoe_Ptee = (TH1F*) hMPtDCA_LS_rec_beautyTOe->ProjectionY("proj_recLS_Btoe_Ptee");
    TH1F* proj_recLS_Btoe_DCA  = (TH1F*) hMPtDCA_LS_rec_beautyTOe->ProjectionZ("proj_recLS_Btoe_DCA");

    TH1F* proj_recLS_LFtoe_Mee  = (TH1F*) hMPtDCA_LS_rec_lfTOe->Project3D("x");
    TH1F* proj_recLS_LFtoe_Ptee = (TH1F*) hMPtDCA_LS_rec_lfTOe->Project3D("y");
    TH1F* proj_recLS_LFtoe_DCA  = (TH1F*) hMPtDCA_LS_rec_lfTOe->Project3D("z");
    proj_recLS_LFtoe_Mee->SetTitle("proj_recLS_LFtoe_Mee");
    proj_recLS_LFtoe_Ptee->SetTitle("proj_recLS_LFtoe_Ptee");
    proj_recLS_LFtoe_DCA->SetTitle("proj_recLS_LFtoe_DCA");

    TH1F* proj_recLS_HFtoe_Mee  = (TH1F*) hMPtDCA_LS_rec_hfTOe->Project3D("x");
    TH1F* proj_recLS_HFtoe_Ptee = (TH1F*) hMPtDCA_LS_rec_hfTOe->Project3D("y");
    TH1F* proj_recLS_HFtoe_DCA  = (TH1F*) hMPtDCA_LS_rec_hfTOe->Project3D("z");
    proj_recLS_HFtoe_Mee->SetTitle("proj_recLS_HFtoe_Mee");
    proj_recLS_HFtoe_Ptee->SetTitle("proj_recLS_HFtoe_Ptee");
    proj_recLS_HFtoe_DCA->SetTitle("proj_recLS_HFtoe_DCA");

    TH1F* proj_recLS_Pi0toee_Mee  = (TH1F*) hMPtDCA_LS_rec_Pi0TOee->Project3D("x");
    TH1F* proj_recLS_Pi0toee_Ptee = (TH1F*) hMPtDCA_LS_rec_Pi0TOee->Project3D("y");
    TH1F* proj_recLS_Pi0toee_DCA  = (TH1F*) hMPtDCA_LS_rec_Pi0TOee->Project3D("z");
    proj_recLS_Pi0toee_Mee->SetTitle("proj_recLS_Pi0toee_Mee");
    proj_recLS_Pi0toee_Ptee->SetTitle("proj_recLS_Pi0toee_Ptee");
    proj_recLS_Pi0toee_DCA->SetTitle("proj_recLS_Pi0toee_DCA");

    TH1F* proj_recLS_LFtoee_Mee  = (TH1F*) hMPtDCA_LS_rec_lfTOee_selectPDG->Project3D("x");
    TH1F* proj_recLS_LFtoee_Ptee = (TH1F*) hMPtDCA_LS_rec_lfTOee_selectPDG->Project3D("y");
    TH1F* proj_recLS_LFtoee_DCA  = (TH1F*) hMPtDCA_LS_rec_lfTOee_selectPDG->Project3D("z");
    proj_recLS_LFtoee_Mee->SetTitle("proj_recLS_LFtoee_Mee");
    proj_recLS_LFtoee_Ptee->SetTitle("proj_recLS_LFtoee_Ptee");
    proj_recLS_LFtoee_DCA->SetTitle("proj_recLS_LFtoee_DCA");

    TH1F* proj_recLS_HFtoee_Mee  = (TH1F*) hMPtDCA_LS_rec_hfTOee->Project3D("x");
    TH1F* proj_recLS_HFtoee_Ptee = (TH1F*) hMPtDCA_LS_rec_hfTOee->Project3D("y");
    TH1F* proj_recLS_HFtoee_DCA  = (TH1F*) hMPtDCA_LS_rec_hfTOee->Project3D("z");
    proj_recLS_HFtoee_Mee->SetTitle("proj_recLS_HFtoee_Mee");
    proj_recLS_HFtoee_Ptee->SetTitle("proj_recLS_HFtoee_Ptee");
    proj_recLS_HFtoee_DCA->SetTitle("proj_recLS_HFtoee_DCA");

    TH1F* proj_recLS_ccTOee_Mee  = (TH1F*) hMPtDCA_LS_rec_ccTOee->ProjectionX("proj_recLS_cctoee_Mee");
    TH1F* proj_recLS_bbTOee_Mee  = (TH1F*) hMPtDCA_LS_rec_bbTOee->ProjectionX("proj_recLS_bbtoee_Mee");


    TH1F* proj_genULS_Mee = (TH1F*) hMPtDCA_ULS_gen->ProjectionX("proj_genULS_Mee");
    TH1F* proj_genULS_Ptee = (TH1F*) hMPtDCA_ULS_gen->ProjectionY("proj_genULS_Ptee");


    TH1F* proj_genLS_Mee = (TH1F*) hMPtDCA_LS_gen->ProjectionX("proj_genLS_Mee");
    TH1F* proj_genLS_Ptee = (TH1F*) hMPtDCA_LS_gen->ProjectionY("proj_genLS_Ptee");


    proj_recULS_Mee->Scale(1./nEventsCent);
    proj_recULS_Ptee->Scale(1./nEventsCent);

    proj_recLS_Mee->Scale(1./nEventsCent);
    proj_recLS_Ptee->Scale(1./nEventsCent);

    proj_genULS_Mee->Scale(1./nEventsCent);
    proj_genULS_Ptee->Scale(1./nEventsCent);
    proj_genLS_Mee->Scale(1./nEventsCent);
    proj_genLS_Ptee->Scale(1./nEventsCent);

    double etacoverage;
    if(!bSemiCentral) etacoverage = 2*etaCutsVec08[ith_PIDscenario-1];
    else if(bSemiCentral) etacoverage = 2*etaCutsVec175[ith_PIDscenario-1];

    // normlise by Nevents, dy = 1.6 (for |eta|<0.8)
    proj_recULS_MCpidEle_charmTOe_Mee->Scale(1/(etacoverage*nEventsCent));
    proj_recULS_MCpidEle_beautyTOe_Mee->Scale(1/(etacoverage*nEventsCent));
    proj_recULS_MCpidEle_hfTOe_Mee->Scale(1/(etacoverage*nEventsCent));
    proj_recULS_MCpidEle_lfTOee_Mee->Scale(1/(etacoverage*nEventsCent));
    proj_recULS_MCpidEle_ccTOee_Mee->Scale(1/(etacoverage*nEventsCent));
    proj_recULS_MCpidEle_bbTOee_Mee->Scale(1/(etacoverage*nEventsCent));
    proj_recULS_MCpidEle_hfTOee_Mee->Scale(1/(etacoverage*nEventsCent));

    proj_recLS_Ctoe_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_Ctoe_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_Ctoe_DCA->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_Btoe_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_Btoe_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_Btoe_DCA->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_LFtoe_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_LFtoe_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_LFtoe_DCA->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_HFtoe_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_HFtoe_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_HFtoe_DCA->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_MCpidEle_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_MCpidEle_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_MCpidEle_DCA->Scale(1./(etacoverage*nEventsCent));

    proj_recLS_Pi0toee_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_Pi0toee_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_Pi0toee_DCA->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_LFtoee_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_LFtoee_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_LFtoee_DCA->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_HFtoee_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_HFtoee_Ptee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_HFtoee_DCA->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_ccTOee_Mee->Scale(1./(etacoverage*nEventsCent));
    proj_recLS_bbTOee_Mee->Scale(1./(etacoverage*nEventsCent));


    normalizeToBinWidth(proj_recULS_Mee);
    normalizeToBinWidth(proj_recULS_Ptee);

    normalizeToBinWidth(proj_recLS_Mee);
    normalizeToBinWidth(proj_recLS_Ptee);

    normalizeToBinWidth(proj_genULS_Mee);
    normalizeToBinWidth(proj_genULS_Ptee);
    normalizeToBinWidth(proj_genLS_Mee);
    normalizeToBinWidth(proj_genLS_Ptee);

    // normalizeToBinWidth(proj_recLS_MCpidEle_Ptee);
    // normalizeToBinWidth(proj_recLS_LFtoee_Ptee);
    // normalizeToBinWidth(proj_recLS_HFtoe_Ptee);


    make3HistNice(proj_recLS_MCpidEle_Mee,kBlack);
    make3HistNice(proj_recLS_MCpidEle_Ptee,kBlack);
    make3HistNice(proj_recLS_MCpidEle_DCA,kBlack);
    make3HistNice(proj_recLS_Ctoe_Mee,kOrange+2);
    make3HistNice(proj_recLS_Ctoe_Ptee,kOrange+2);
    make3HistNice(proj_recLS_Ctoe_DCA,kOrange+2);
    make3HistNice(proj_recLS_Btoe_Mee,kMagenta+1);
    make3HistNice(proj_recLS_Btoe_Ptee,kMagenta+1);
    make3HistNice(proj_recLS_Btoe_DCA,kMagenta+1);
    make3HistNice(proj_recLS_LFtoe_Mee,kGreen+1);
    make3HistNice(proj_recLS_LFtoe_Ptee,kGreen+1);
    make3HistNice(proj_recLS_LFtoe_DCA,kGreen+1);
    make3HistNice(proj_recLS_HFtoe_Mee,kRed+2);
    make3HistNice(proj_recLS_HFtoe_Ptee,kRed+2);
    make3HistNice(proj_recLS_HFtoe_DCA,kRed+2);
    make3HistNice(proj_recLS_Pi0toee_Mee,kRed);
    make3HistNice(proj_recLS_Pi0toee_Ptee,kRed);
    make3HistNice(proj_recLS_Pi0toee_DCA,kRed);
    make3HistNice(proj_recLS_LFtoee_Mee,kCyan+1);
    make3HistNice(proj_recLS_LFtoee_Ptee,kCyan+1);
    make3HistNice(proj_recLS_LFtoee_DCA,kCyan+1);
    make3HistNice(proj_recLS_HFtoee_Mee,kRed+2);
    make3HistNice(proj_recLS_HFtoee_Ptee,kRed+2);
    make3HistNice(proj_recLS_HFtoee_DCA,kRed+2);
    makeHistNice(proj_recLS_ccTOee_Mee,kOrange+2);
    makeHistNice(proj_recLS_bbTOee_Mee,kMagenta+1);

    makeHistNice(proj_recULS_MCpidEle_charmTOe_Mee,kOrange+2);
    makeHistNice(proj_recULS_MCpidEle_beautyTOe_Mee,kMagenta+1);
    makeHistNice(proj_recULS_MCpidEle_hfTOe_Mee,kRed+2);
    makeHistNice(proj_recULS_MCpidEle_lfTOee_Mee,kBlue+1);
    makeHistNice(proj_recULS_MCpidEle_ccTOee_Mee,kOrange+2);
    makeHistNice(proj_recULS_MCpidEle_bbTOee_Mee,kMagenta+1);
    makeHistNice(proj_recULS_MCpidEle_hfTOee_Mee,kRed+2);




 TH1F* ptRecTrackBeforeSmearing_rebin = (TH1F*) ptRecTrackBeforeSmearing->Rebin(nbinspt_proj_2D,"ptRecTrackBeforeSmearing_rebin",&pte_bin_2D[0]);
 TH1F* ptRecTrackAfterSmearing_rebin = (TH1F*) ptRecTrackAfterSmearing->Rebin(nbinspt_proj_2D,"ptRecTrackAfterSmearing_rebin",&pte_bin_2D[0]);
 TH1F* ptRecTrackAfterKineCuts_rebin = (TH1F*) ptRecTrackAfterKineCuts->Rebin(nbinspt_proj_2D,"ptRecTrackAfterKineCuts_rebin",&pte_bin_2D[0]);

 TH1F* ptRecTrackElePos_rebin = (TH1F*) ptRecTrackElePos->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_rebin",&pte_bin_2D[0]);
 TH1F* ptGenTrackElePos_rebin = (TH1F*) ptGenTrackElePos->Rebin(nbinspt_proj_2D,"ptGenTrackElePos_rebin",&pte_bin_2D[0]);
 TH1F* ptGenTrackAll_rebin = (TH1F*) ptGenTrackAll->Rebin(nbinspt_proj_2D,"ptGenTrackAll_rebin",&pte_bin_2D[0]);
 TH1F* ptGenTrackMuon_rebin = (TH1F*) ptGenTrackMuon->Rebin(nbinspt_proj_2D,"ptGenTrackMuon_rebin",&pte_bin_2D[0]);
 TH1F* ptGenTrackPion_rebin = (TH1F*) ptGenTrackPion->Rebin(nbinspt_proj_2D,"ptGenTrackPion_rebin",&pte_bin_2D[0]);
 // TH1F* ptGenTrackKaon_rebin = (TH1F*) ptGenTrackKaon->Rebin(nbinspt_proj_2D,"ptGenTrackKaon_rebin",&pte_bin_2D[0]);
 // TH1F* ptGenTrackProton_rebin = (TH1F*) ptGenTrackProton->Rebin(nbinspt_proj_2D,"ptGenTrackProton_rebin",&pte_bin_2D[0]);
 TH1F* ptGenSmearedTrackElePos_rebin = (TH1F*) ptGenSmearedTrackElePos->Rebin(nbinspt_proj_2D,"ptGenSmearedTrackElePos_rebin",&pte_bin_2D[0]);
 TH1F* ptGenTrackElePos_beforeKineCuts_rebin = (TH1F*) ptGenTrackElePos_beforeKineCuts->Rebin(nbinspt_proj_2D,"ptGenTrackElePos_beforeKineCuts_rebin",&pte_bin_2D[0]);
 // TH1F* ptGenSmearedTrackElePos_beforeKineCuts_rebin = (TH1F*) ptGenSmearedTrackElePos_beforeKineCuts->Rebin(nbinspt_proj_2D,"ptGenSmearedTrackElePos_beforeKineCuts_rebin",&pte_bin_2D[0]);
 TH1F* ptGenTrackElePos_LF_rebin = (TH1F*) ptGenTrackElePos_LF->Rebin(nbinspt_proj_2D,"ptGenTrackElePos_LF_rebin",&pte_bin_2D[0]);
 TH1F* ptGenTrackElePos_HF_rebin = (TH1F*) ptGenTrackElePos_HF->Rebin(nbinspt_proj_2D,"ptGenTrackElePos_HF_rebin",&pte_bin_2D[0]);
 // TH1F* ptGenTrackElePos_Pi0_rebin = (TH1F*) ptGenTrackElePos_Pi0->Rebin(nbinspt_proj_2D,"ptGenTrackElePos_Pi0_rebin",&pte_bin_2D[0]);



 TH1F* ptRecTrackElePos_Pi0_beforePID_rebin = (TH1F*) ptRecTrackElePos_Pi0_beforePID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_Pi0_beforePID_rebin",&pte_bin_2D[0]);
 TH1F* ptRecTrackElePos_LF_beforePID_rebin = (TH1F*) ptRecTrackElePos_LF_beforePID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_LF_beforePID_rebin",&pte_bin_2D[0]);
 TH1F* ptRecTrackElePos_HF_beforePID_rebin = (TH1F*) ptRecTrackElePos_HF_beforePID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_HF_beforePID_rebin",&pte_bin_2D[0]);

 TH1F* ptRecTrackElePos_Pi0_afterPID_rebin = (TH1F*) ptRecTrackElePos_Pi0_afterPID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_Pi0_afterPID_rebin",&pte_bin_2D[0]);
 TH1F* ptRecTrackElePos_LF_afterPID_rebin = (TH1F*) ptRecTrackElePos_LF_afterPID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_LF_afterPID_rebin",&pte_bin_2D[0]);
 TH1F* ptRecTrackElePos_HF_afterPID_rebin = (TH1F*) ptRecTrackElePos_HF_afterPID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_HF_afterPID_rebin",&pte_bin_2D[0]);

 TH1F* ptAllRecTrack_rebin = (TH1F*) ptAllRecTrack->Rebin(nbinspt_proj_2D,"ptAllRecTrack_rebin",&pte_bin_2D[0]);
 TH1F* ptRecMuonTrack_rebin = (TH1F*) ptRecMuonTrack->Rebin(nbinspt_proj_2D,"ptRecMuonTrack_rebin",&pte_bin_2D[0]);
 TH1F* ptRecPionTrack_rebin = (TH1F*) ptRecPionTrack->Rebin(nbinspt_proj_2D,"ptRecPionTrack_rebin",&pte_bin_2D[0]);
 TH1F* ptRecKaonTrack_rebin = (TH1F*) ptRecKaonTrack->Rebin(nbinspt_proj_2D,"ptRecKaonTrack_rebin",&pte_bin_2D[0]);
 TH1F* ptRecProtonTrack_rebin = (TH1F*) ptRecProtonTrack->Rebin(nbinspt_proj_2D,"ptRecProtonTrack_rebin",&pte_bin_2D[0]);

 TH1F* ptRecTrackElePos_beforePID_rebin = (TH1F*) ptRecTrackElePos_beforePID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_beforePID_rebin",&pte_bin_2D[0]);
 // TH1F* ptRecTrackElePos_afterPSPID_rebin = (TH1F*) ptRecTrackElePos_afterPSPID->Rebin(nbinspt_proj_2D,"ptRecTrackElePos_afterPSPID_rebin",&pte_bin_2D[0]);
 // TH1F* ptRecMuonTrack_beforePID_rebin = (TH1F*) ptRecMuonTrack_beforePID->Rebin(nbinspt_proj_2D,"ptRecMuonTrack_beforePID_rebin",&pte_bin_2D[0]);
 TH1F* ptRecPionTrack_beforePID_rebin = (TH1F*) ptRecPionTrack_beforePID->Rebin(nbinspt_proj_2D,"ptRecPionTrack_beforePID_rebin",&pte_bin_2D[0]);
 // TH1F* ptRecPionTrack_afterPSPID_rebin = (TH1F*) ptRecPionTrack_afterPSPID->Rebin(nbinspt_proj_2D,"ptRecPionTrack_afterPSPID_rebin",&pte_bin_2D[0]);


 TH1F* ptRecNegTrack_rebin = (TH1F*) ptRecNegTrack->Rebin(nbinspt_proj_2D,"ptRecNegTrack_rebin",&pte_bin_2D[0]);
 TH1F* ptRecPosTrack_rebin = (TH1F*) ptRecPosTrack->Rebin(nbinspt_proj_2D,"ptRecPosTrack_rebin",&pte_bin_2D[0]);


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
  TH1F* ptEffElePos_PreShower;
  TH1F* ptEffElePos_LF_afterPF;
  TH1F* ptEffElePos_HF_afterPF;

  TH2F* ptEtaEffElePos;
  TH2F* mPtPairEffElePos;

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
    ptEffElePosGen->Divide(ptEffElePosGen,ptGenTrackElePos_rebin,1,1,"B");
    ptEffElePosGenSmeared = (TH1F*) ptRecTrackElePos_rebin->Clone("eff_pT_singleElePos_GenSmeared");
    ptEffElePosGenSmeared->Divide(ptEffElePosGenSmeared,ptGenSmearedTrackElePos_rebin,1,1,"B");

    ptEffElePosGen_woPID = (TH1F*) ptRecTrackAfterKineCuts_rebin->Clone("eff_pT_singleElePos_Gen_woPID");
    ptEffElePosGen_woPID->Divide(ptEffElePosGen_woPID,ptGenTrackElePos_rebin,1,1,"B");

    // single track efficiency TH2 as function of pt and eta
    ptEtaEffElePos = (TH2F*) ptEtaRecTrackElePos->Clone("eff_ptEta_singleElePos");
    ptEtaEffElePos->Divide(ptEtaEffElePos,ptEtaGenTrackElePos,1,1,"B");

    ptEffElePos_LF_afterPF = (TH1F*) ptRecTrackElePos_LF_afterPID_rebin->Clone("eff_pT_singleElePos_LF_afterPF");
    ptEffElePos_LF_afterPF->Divide(ptEffElePos_LF_afterPF,ptGenTrackElePos_LF_rebin,1,1,"B");

    ptEffElePos_HF_afterPF = (TH1F*) ptRecTrackElePos_HF_afterPID_rebin->Clone("eff_pT_singleElePos_HF_afterPF");
    ptEffElePos_HF_afterPF->Divide(ptEffElePos_HF_afterPF,ptGenTrackElePos_HF_rebin,1,1,"B");


   // pair efficiencies
   // // ptPairEffULS = (TH1F*) proj_recULS_Ptee->Clone();
   // ptPairEffULS = (TH1F*) proj_recULS_MCpidEle_Ptee->Clone();
   // // massPairEffULS = (TH1F*) proj_recULS_Mee->Clone();
   // massPairEffULS = (TH1F*) proj_recULS_MCpidEle_Mee->Clone();
   // ptPairEffULS->Sumw2();
   // massPairEffULS->Sumw2();
   // ptPairEffULS->Divide(ptPairEffULS,proj_genULS_Ptee,1,1,"B");
   // massPairEffULS->Divide(massPairEffULS,proj_genULS_Mee,1,1,"B");

   mPtPairEffElePos = (TH2F*) hMPtDCA_ULS_rec_MCpidEle->Clone("pairEff_mPt_ElePos");
   // mPtPairEffElePos->Sumw2();
   mPtPairEffElePos->Divide(mPtPairEffElePos,hMPtDCA_ULS_gen,1,1,"B");

   ptEffElePos_PreShower = (TH1F*) ptRecTrackElePos_rebin->Clone("eff_pT_ElePos_PreShower");
   ptEffElePos_PreShower->Divide(ptEffElePos_PreShower,ptGenTrackElePos_rebin,1,1,"B");

  }


  TH1F* hMuonRejectionFactorPt;
  TH1F* hPionRejectionFactorPt;
  // TH1F* hKaonRejectionFactorPt;
  // TH1F* hProtonRejectionFactorPt;
  TH1F* hTotalRejectionFactorPt;
  TH1F* ptRecContamination;
  ptRecContamination = (TH1F*)  ptAllRecTrack_rebin->Clone("pt_Total_contamination");
  ptRecContamination->Add(ptRecTrackElePos_rebin,-1);

  hMuonRejectionFactorPt = (TH1F*) ptGenTrackMuon_rebin->Clone("pT_Muon_Rejectionfactor");
  hPionRejectionFactorPt = (TH1F*) ptGenTrackPion_rebin->Clone("pT_Pion_Rejectionfactor");
  // hKaonRejectionFactorPt = (TH1F*) ptGenTrackKaon_rebin->Clone("pT_Kaon_Rejectionfactor");
  // hProtonRejectionFactorPt = (TH1F*) ptGenTrackProton_rebin->Clone("pT_Proton_Rejectionfactor");
  hTotalRejectionFactorPt = (TH1F*) ptGenTrackAll_rebin->Clone("pT_Total_Rejectionfactor");
  hTotalRejectionFactorPt->Add(ptGenTrackElePos_rebin,-1);


  hMuonRejectionFactorPt->Divide(hMuonRejectionFactorPt, ptRecMuonTrack_rebin,1,1,"B");
  hPionRejectionFactorPt->Divide(hPionRejectionFactorPt, ptRecPionTrack_rebin,1,1,"B");
  // hKaonRejectionFactorPt->Divide(hKaonRejectionFactorPt, ptRecKaonTrack_rebin,1,1,"B");
  // hProtonRejectionFactorPt->Divide(hProtonRejectionFactorPt, ptRecProtonTrack_rebin,1,1,"B");
  hTotalRejectionFactorPt->Divide(hTotalRejectionFactorPt, ptRecContamination,1,1,"B");


  make3HistNice(hMuonRejectionFactorPt, kBlue+1);
  make3HistNice(hPionRejectionFactorPt, kRed+2);
  // make3HistNice(hKaonRejectionFactorPt, kGreen+2);
  // make3HistNice(hProtonRejectionFactorPt, kOrange+2);
  make3HistNice(hTotalRejectionFactorPt, kBlack);
  hTotalRejectionFactorPt->SetMarkerStyle(24);
  hTotalRejectionFactorPt->SetMarkerSize(0.9);



  TH1F* hTotalPureContaminationRecPt;
  TH1F* hTotalPureContaminationRecEta;
  TH1F* hTotalPureContaminationRecPhi;
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
  TH1F* hProtonContaminationRecPt;



if (bPlotTrackContamination) {
  hTotalRecTrackPt = (TH1F*)  ptAllRecTrack_rebin->Clone();
  hTotalRecTrackEta = (TH1F*) etaAllRecTrack->Clone();
  hTotalRecTrackPhi = (TH1F*) phiAllRecTrack->Clone();

  hTotalPureContaminationRecPt = (TH1F*)  ptAllRecTrack_rebin->Clone("TotalContaminationPt");
  hTotalPureContaminationRecEta = (TH1F*) etaAllRecTrack->Clone("TotalContaminationEta");
  hTotalPureContaminationRecPhi = (TH1F*) phiAllRecTrack->Clone("TotalContaminationPhi");
  hMuonContaminationRecPt = (TH1F*)    ptRecMuonTrack_rebin->Clone();
  hMuonContaminationRecEta = (TH1F*)   etaRecMuonTrack->Clone();
  hMuonContaminationRecPhi = (TH1F*)   phiRecMuonTrack->Clone();
  hPionContaminationRecPt = (TH1F*)    ptRecPionTrack_rebin->Clone();
  hPionContaminationRecEta = (TH1F*)   etaRecPionTrack->Clone();
  hPionContaminationRecPhi = (TH1F*)   phiRecPionTrack->Clone();
  hKaonContaminationRecPt = (TH1F*)    ptRecKaonTrack_rebin->Clone();
  hProtonContaminationRecPt = (TH1F*)  ptRecProtonTrack_rebin->Clone();



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


  hMuonContaminationRecPt->Divide(hMuonContaminationRecPt, hTotalRecTrackPt,1,1,"B");
  hMuonContaminationRecEta->Divide(hMuonContaminationRecEta, hTotalRecTrackEta,1,1,"B");
  hMuonContaminationRecPhi->Divide(hMuonContaminationRecPhi, hTotalRecTrackPhi,1,1,"B");
  hPionContaminationRecPt->Divide(hPionContaminationRecPt, hTotalRecTrackPt,1,1,"B");
  hPionContaminationRecEta->Divide(hPionContaminationRecEta, hTotalRecTrackEta,1,1,"B");
  hPionContaminationRecPhi->Divide(hPionContaminationRecPhi, hTotalRecTrackPhi,1,1,"B");
  hKaonContaminationRecPt->Divide(hKaonContaminationRecPt, hTotalRecTrackPt,1,1,"B");
  hProtonContaminationRecPt->Divide(hProtonContaminationRecPt, hTotalRecTrackPt,1,1,"B");


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


  make3HistNice(hMuonContaminationRecPt, kBlue+1);
  make3HistNice(hMuonContaminationRecEta, kBlue+1);
  make3HistNice(hMuonContaminationRecPhi, kBlue+1);
  make3HistNice(hPionContaminationRecPt, kRed+2);
  make3HistNice(hPionContaminationRecEta, kRed+2);
  make3HistNice(hPionContaminationRecPhi, kRed+2);
  make3HistNice(hKaonContaminationRecPt, kGreen+2);
  make3HistNice(hProtonContaminationRecPt, kOrange+2);
  make3HistNice(hTotalPureContaminationRecPt, kBlack);
  make3HistNice(hTotalPureContaminationRecEta, kBlack);
  make3HistNice(hTotalPureContaminationRecPhi, kBlack);


  hTotalPureContaminationRecPt->SetMarkerStyle(24);
  hTotalPureContaminationRecPt->SetMarkerSize(0.9);
  hTotalPureContaminationRecEta->SetMarkerStyle(24);
  hTotalPureContaminationRecEta->SetMarkerSize(0.9);
  hTotalPureContaminationRecPhi->SetMarkerStyle(24);
  hTotalPureContaminationRecPhi->SetMarkerSize(0.9);

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

  makeHistNice(ptRecTrackElePos_Pi0_beforePID_rebin,kCyan+2);
  makeHistNice(ptRecTrackElePos_LF_beforePID_rebin,kBlue+3);
  makeHistNice(ptRecTrackElePos_HF_beforePID_rebin,kRed+3);
  ptRecTrackElePos_Pi0_beforePID_rebin->SetMarkerStyle(24);
  ptRecTrackElePos_LF_beforePID_rebin->SetMarkerStyle(24);
  ptRecTrackElePos_HF_beforePID_rebin->SetMarkerStyle(24);
  makeHistNice(ptRecTrackElePos_Pi0_afterPID_rebin,kCyan);
  makeHistNice(ptRecTrackElePos_LF_afterPID_rebin,kBlue+1);
  makeHistNice(ptRecTrackElePos_HF_afterPID_rebin,kRed+1);

  makeHistNice(ptGenTrackElePos_HF_rebin,kBlue);
  ptGenTrackElePos_HF_rebin->SetMarkerStyle(24);

  // makeHistNice(ptGenTrackElePos_Pi0_rebin, kOrange+1);

  normalizeToBinWidth(ptGenTrackElePos_rebin);
  normalizeToBinWidth(ptGenSmearedTrackElePos_rebin);
  normalizeToBinWidth(ptRecTrackElePos_rebin);
  //
  normalizeToBinWidth(ptGenTrackElePos_beforeKineCuts_rebin);
  // normalizeToBinWidth(ptGenSmearedTrackElePos_beforeKineCuts_rebin);

  normalizeToBinWidth(ptRecTrackElePos_Pi0_beforePID_rebin);
  normalizeToBinWidth(ptRecTrackElePos_LF_beforePID_rebin);
  normalizeToBinWidth(ptRecTrackElePos_HF_beforePID_rebin);
  normalizeToBinWidth(ptRecTrackElePos_Pi0_afterPID_rebin);
  normalizeToBinWidth(ptRecTrackElePos_LF_afterPID_rebin);
  normalizeToBinWidth(ptRecTrackElePos_HF_afterPID_rebin);


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

  normalizeToBinWidth(ptRecTrackBeforeSmearing_rebin);
  normalizeToBinWidth(ptRecTrackAfterSmearing_rebin);
  normalizeToBinWidth(ptRecTrackAfterKineCuts_rebin);




  if (bPlotEfficiency) {
    makeHistNice(ptEffElePosGen,kBlue+1);
    makeHistNice(ptEffElePosGenSmeared,kRed+2);
    makeHistNice(ptEffElePosGen_woPID,kViolet+2);

    makeHistNice(ptEffElePos_PreShower,kBlue+1);
    makeHistNice(ptEffElePos_HF_afterPF,kBlue+1);
    makeHistNice(ptEffElePos_LF_afterPF,kBlue+1);

  }

  makeHistNice(proj_recULS_Mee,kBlack);
  makeHistNice(proj_recULS_Ptee,kBlack);


  makeHistNice(proj_recLS_Mee,kRed+2);
  makeHistNice(proj_recLS_Ptee,kRed+2);



  // calculate ULS-LS (Signal) for charm,beauty,hf->e
  TH1F* hSignal_charmTOe_Mee;
  TH1F* hSignal_beautyTOe_Mee;
  TH1F* hSignal_HFTOe_Mee;
  TH1F* hSignal_LFTOee_Mee;
  TH1F* hSignal_HFTOee_Mee;
  TH1F* hSignal_ccTOee_Mee;
  TH1F* hSignal_bbTOee_Mee;
  TH1F* hSignal_ccTOee_Mee_rebin;
  TH1F* hSignal_bbTOee_Mee_rebin;
  TH1F* hSignal_HFTOee_Mee_rebin;
  TH1F* hSignal_LFTOee_Mee_rebin;
  TH1F* hCocktailCCtoee;
  TH1F* hCocktailBBtoee;
  TH1F* hCocktailPi0toee;
  TH1F* hCocktailEtatoee;
  TH1F* hCocktailEtaPrimetoee;
  TH1F* hCocktailOmegatoee;
  TH1F* hCocktailPhitoee;
  TH1F* hTotalLFCocktailtoee;

  if(bPlotLFHFcontributions){
    hSignal_charmTOe_Mee = (TH1F*) proj_recULS_MCpidEle_charmTOe_Mee->Clone();
    hSignal_beautyTOe_Mee = (TH1F*) proj_recULS_MCpidEle_beautyTOe_Mee->Clone();
    hSignal_HFTOe_Mee = (TH1F*) proj_recULS_MCpidEle_hfTOe_Mee->Clone();
    hSignal_charmTOe_Mee->Add(proj_recLS_Ctoe_Mee,-1.);
    hSignal_beautyTOe_Mee->Add(proj_recLS_Btoe_Mee,-1.);
    hSignal_HFTOe_Mee->Add(proj_recLS_HFtoe_Mee,-1.);
    ScaleBinWidth(hSignal_charmTOe_Mee,kTRUE);
    ScaleBinWidth(hSignal_beautyTOe_Mee,kTRUE);
    ScaleBinWidth(hSignal_HFTOe_Mee,kTRUE);

    hSignal_LFTOee_Mee = (TH1F*) proj_recULS_MCpidEle_lfTOee_Mee->Clone();
    hSignal_LFTOee_Mee->Add(proj_recLS_LFtoee_Mee,-1.);
    makeHistNice(hSignal_LFTOee_Mee,kBlue+1);


    hSignal_ccTOee_Mee  = (TH1F*) proj_recULS_MCpidEle_ccTOee_Mee->Clone();
    hSignal_bbTOee_Mee  = (TH1F*) proj_recULS_MCpidEle_bbTOee_Mee->Clone();
    hSignal_HFTOee_Mee  = (TH1F*) proj_recULS_MCpidEle_hfTOee_Mee->Clone();
    hSignal_ccTOee_Mee->Add(proj_recLS_ccTOee_Mee,-1.);
    hSignal_bbTOee_Mee->Add(proj_recLS_bbTOee_Mee,-1.);
    hSignal_HFTOee_Mee->Add(proj_recLS_HFtoee_Mee,-1.);

    Double_t pt_bin_HFTOee[]  = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.5,3.0,3.5,4.0};
    Int_t nbinspt_HFTOee  = sizeof(pt_bin_HFTOee)/sizeof(*pt_bin_HFTOee) -1;
    hSignal_ccTOee_Mee_rebin = (TH1F*) hSignal_ccTOee_Mee->Rebin(nbinspt_HFTOee,"ptRecULS_ccTOee_rebin",&pt_bin_HFTOee[0]);
    hSignal_bbTOee_Mee_rebin = (TH1F*) hSignal_bbTOee_Mee->Rebin(nbinspt_HFTOee,"ptRecULS_bbTOee_rebin",&pt_bin_HFTOee[0]);
    hSignal_HFTOee_Mee_rebin = (TH1F*) hSignal_HFTOee_Mee->Rebin(nbinspt_HFTOee,"ptRecULS_HFTOee_rebin",&pt_bin_HFTOee[0]);
    hSignal_LFTOee_Mee_rebin = (TH1F*) hSignal_LFTOee_Mee->Rebin(nbinspt_HFTOee,"ptRecULS_LFTOee_rebin",&pt_bin_HFTOee[0]);
    // ScaleBinWidth(hSignal_ccTOee_Mee_rebin,kTRUE);
    // ScaleBinWidth(hSignal_bbTOee_Mee_rebin,kTRUE);
    // ScaleBinWidth(hSignal_HFTOee_Mee_rebin,kTRUE);

    // load cocktail files
    if(bUseAdditionalFiles && bCocktailFile){
      TFile* fCocktailLFHF = TFile::Open("/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/efficiency/data/cocktails_check.root");
      hCocktailCCtoee = (TH1F*) fCocktailLFHF->Get("ccbar_witheff");
      hCocktailBBtoee = (TH1F*) fCocktailLFHF->Get("bbbar_witheff");
      hCocktailPi0toee = (TH1F*) fCocktailLFHF->Get("pi0_witheff");
      hCocktailEtatoee = (TH1F*) fCocktailLFHF->Get("eta_witheff");
      hCocktailEtaPrimetoee = (TH1F*) fCocktailLFHF->Get("etaprim_witheff");
      hCocktailOmegatoee = (TH1F*) fCocktailLFHF->Get("omega_witheff");
      hCocktailPhitoee = (TH1F*) fCocktailLFHF->Get("phi_witheff");
      hTotalLFCocktailtoee = (TH1F*) hCocktailPi0toee->Clone();
      hTotalLFCocktailtoee->Add(hCocktailEtatoee);
      hTotalLFCocktailtoee->Add(hCocktailEtaPrimetoee);
      hTotalLFCocktailtoee->Add(hCocktailOmegatoee);
      hTotalLFCocktailtoee->Add(hCocktailPhitoee);
      makeHistNice(hCocktailCCtoee,kMagenta+3);
      makeHistNice(hCocktailBBtoee,kOrange+4);
      makeHistNice(hCocktailPi0toee,kRed-2);
      makeHistNice(hCocktailEtatoee,kGreen);
      makeHistNice(hCocktailEtaPrimetoee,kGreen+3);
      makeHistNice(hCocktailOmegatoee,kBlue-3);
      makeHistNice(hCocktailPhitoee,kYellow+2);
      makeHistNice(hTotalLFCocktailtoee,kBlack);
    }
  }

  // Plot the ratio of electrons/pions
  // TH1F* hElePionRatio_beforePID;
  // TH1F* hMuonEleRatio_beforePID;
  // TH1F* hElePionRatio_afterPID;
  // TH1F* hMuonEleRatio_afterPID;
  // if(bElePionRatios){
  //   normalizeToBinWidth(ptRecMuonTrack_rebin);
  //   normalizeToBinWidth(ptRecPionTrack_rebin);
  //   normalizeToBinWidth(ptRecTrackElePos_beforePID_rebin);
  //   normalizeToBinWidth(ptRecMuonTrack_beforePID_rebin);
  //   normalizeToBinWidth(ptRecPionTrack_beforePID_rebin);
  //
  //   // before PID
  //   hElePionRatio_beforePID = (TH1F*) ptRecTrackElePos_beforePID_rebin->Clone();
  //   hMuonEleRatio_beforePID = (TH1F*) ptRecMuonTrack_beforePID_rebin->Clone();
  //   // after PID
  //   hElePionRatio_afterPID = (TH1F*) ptRecTrackElePos_rebin->Clone();
  //   hMuonEleRatio_afterPID = (TH1F*) ptRecMuonTrack_rebin->Clone();
  //
  //   hElePionRatio_beforePID->Divide(ptRecPionTrack_beforePID_rebin);
  //   hElePionRatio_afterPID ->Divide(ptRecPionTrack_rebin);
  //   hMuonEleRatio_beforePID->Divide(ptRecTrackElePos_beforePID_rebin);
  //   hMuonEleRatio_afterPID ->Divide(ptRecTrackElePos_rebin);
  //
  //   makeHistNice(hElePionRatio_beforePID, kBlack);
  //   makeHistNice(hElePionRatio_afterPID, kBlue+1);
  //   makeHistNice(hMuonEleRatio_beforePID, kRed+2);
  //   makeHistNice(hMuonEleRatio_afterPID, kOrange+1);
  // }



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

  auto legEff_PreShower = new TLegend(0.2,0.8,0.5,0.9);
  legEff_PreShower->SetBorderSize(0);
  legEff_PreShower->SetFillStyle(0);
  legEff_PreShower->SetTextSize(0.04);
  legEff_PreShower->AddEntry(ptEffElePos_PreShower,"rec/gen","p");
  legEff_PreShower->AddEntry(funcPSparam,"PS eff fit","l");

  auto legEff_LF_afterPF = new TLegend(0.6,0.8,0.9,0.9);
  legEff_LF_afterPF->SetBorderSize(0);
  legEff_LF_afterPF->SetFillStyle(0);
  legEff_LF_afterPF->SetTextSize(0.04);
  legEff_LF_afterPF->AddEntry(ptEffElePos_LF_afterPF,"LF efficiency ","p");

  auto legEff_HF_afterPF = new TLegend(0.6,0.8,0.9,0.9);
  legEff_HF_afterPF->SetBorderSize(0);
  legEff_HF_afterPF->SetFillStyle(0);
  legEff_HF_afterPF->SetTextSize(0.04);
  legEff_HF_afterPF->AddEntry(ptEffElePos_HF_afterPF,"HF efficiency ","p");


  auto legChPionPt = new TLegend(legPosTrack[0]-0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legChPionPt->SetBorderSize(0);
  legChPionPt->SetFillStyle(0);
  legChPionPt->SetTextFont(43);
  legChPionPt->SetTextSize(24);
  legChPionPt->AddEntry(ptGenTrackPionRapSel,"fast sim. #pi^{+} + #pi^{-}","p");
  legChPionPt->AddEntry(hPaper_chPi_Pt_fullCentralSection_stat_Err,"paper #pi^{+} + #pi^{-}","p");


  auto legSigCBHFtoe = new TLegend(legPosTrack[0]-0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legSigCBHFtoe->SetBorderSize(0);
  legSigCBHFtoe->SetFillStyle(0);
  legSigCBHFtoe->SetTextFont(43);
  legSigCBHFtoe->SetTextSize(24);
  legSigCBHFtoe->AddEntry(hSignal_HFTOe_Mee,"ULS-LS HF #rightarrow e","p");
  legSigCBHFtoe->AddEntry(hSignal_charmTOe_Mee,"ULS-LS c #rightarrow e","p");
  legSigCBHFtoe->AddEntry(hSignal_beautyTOe_Mee,"ULS-LS b #rightarrow e","p");

  auto legHFtoe = new TLegend(legPosTrack[0]-0.15,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legHFtoe->SetBorderSize(0);
  legHFtoe->SetFillStyle(0);
  legHFtoe->SetTextFont(43);
  legHFtoe->SetTextSize(24);
  legHFtoe->SetHeader("HF #rightarrow e");
  legHFtoe->AddEntry(ptGenTrackElePos_HF_rebin,"fast sim.","l");
  legHFtoe->AddEntry(hSingleElefromHF_Cocktail,"extrapolation","l");
  legHFtoe->AddEntry(hPaper_HFtoe_Pt_stat,"PLB 804 (2020) 135377","pe1");





auto legSigLFtoee = new TLegend(legPosTrack[0],legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
legSigLFtoee->SetBorderSize(0);
legSigLFtoee->SetFillStyle(0);
legSigLFtoee->SetTextFont(43);
legSigLFtoee->SetTextSize(24);
legSigLFtoee->AddEntry(hSignal_LFTOee_Mee,"ULS-LS LF #rightarrow ee","p");
legSigLFtoee->AddEntry(hTotalLFCocktailtoee,"cocktail total LF","l");
legSigLFtoee->AddEntry(hCocktailPi0toee,"cocktail #pi^{0} #rightarrow ee","l");
legSigLFtoee->AddEntry(hCocktailEtatoee,"cocktail #eta #rightarrow ee","l");
legSigLFtoee->AddEntry(hCocktailEtaPrimetoee,"cocktail #eta' #rightarrow ee","l");
legSigLFtoee->AddEntry(hCocktailOmegatoee,"cocktail #omega #rightarrow ee","l");
legSigLFtoee->AddEntry(hCocktailPhitoee,"cocktail #phi #rightarrow ee","l");

  auto legSigCBHFtoee = new TLegend(legPosTrack[0],legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legSigCBHFtoee->SetBorderSize(0);
  legSigCBHFtoee->SetFillStyle(0);
  legSigCBHFtoee->SetTextFont(43);
  legSigCBHFtoee->SetTextSize(24);
  legSigCBHFtoee->AddEntry(hSignal_HFTOee_Mee,"ULS-LS HF #rightarrow ee","p");
  legSigCBHFtoee->AddEntry(hSignal_ccTOee_Mee,"ULS-LS cc #rightarrow ee","p");
  legSigCBHFtoee->AddEntry(hSignal_bbTOee_Mee,"ULS-LS bb #rightarrow ee","p");
  legSigCBHFtoee->AddEntry(hCocktailCCtoee,"cocktail cc #rightarrow ee","l");
  legSigCBHFtoee->AddEntry(hCocktailBBtoee,"cocktail bb #rightarrow ee","l");


  auto legTrackPi0LFHF = new TLegend(legPosTrack[0],legPosTrack[1]-0.1,legPosTrack[2]+0.05,legPosTrack[3]);
  legTrackPi0LFHF->SetBorderSize(0);
  legTrackPi0LFHF->SetFillStyle(0);
  legTrackPi0LFHF->SetTextFont(43);
  legTrackPi0LFHF->SetTextSize(25);
  legTrackPi0LFHF->AddEntry((TObject*)0,"open: before PID","");
  legTrackPi0LFHF->AddEntry((TObject*)0,"solid: after PID","");
  legTrackPi0LFHF->AddEntry(ptRecTrackElePos_rebin,"Sum, after PID","p");
  legTrackPi0LFHF->AddEntry(ptRecTrackElePos_LF_afterPID_rebin,"LF->e","p");
  legTrackPi0LFHF->AddEntry(ptRecTrackElePos_Pi0_afterPID_rebin,"#pi^{0}->e","p");
  legTrackPi0LFHF->AddEntry(ptRecTrackElePos_HF_afterPID_rebin,"HF->e","p");

  auto legTrackLFHFtoe = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legTrackLFHFtoe->SetBorderSize(0);
  legTrackLFHFtoe->SetFillStyle(0);
  legTrackLFHFtoe->SetTextFont(43);
  legTrackLFHFtoe->SetTextSize(25);
  legTrackLFHFtoe->AddEntry(proj_recLS_Ctoe_Mee,"c->e","p");
  legTrackLFHFtoe->AddEntry(proj_recLS_Btoe_Mee,"b->e","p");
  legTrackLFHFtoe->AddEntry(proj_recLS_HFtoe_Mee,"HF->e","p");

  auto legTrackLFHFtoee = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legTrackLFHFtoee->SetBorderSize(0);
  legTrackLFHFtoee->SetFillStyle(0);
  legTrackLFHFtoee->SetTextFont(43);
  legTrackLFHFtoee->SetTextSize(25);
  legTrackLFHFtoee->AddEntry(proj_recLS_LFtoee_Mee,"LF->ee","p");
  legTrackLFHFtoee->AddEntry(proj_recLS_HFtoee_Mee,"HF->ee","p");


  auto legTrackLFtoeeHFtoe = new TLegend(legPosTrack[0]-0.15,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legTrackLFtoeeHFtoe->SetBorderSize(0);
  legTrackLFtoeeHFtoe->SetFillStyle(0);
  legTrackLFtoeeHFtoe->SetTextFont(43);
  legTrackLFtoeeHFtoe->SetTextSize(17);
  legTrackLFtoeeHFtoe->AddEntry(proj_recLS_MCpidEle_Mee,"LS, Sum, two electrons","p");
  legTrackLFtoeeHFtoe->AddEntry(proj_recLS_LFtoe_Mee,"LS, at least one LF leg","p");
  legTrackLFtoeeHFtoe->AddEntry(proj_recLS_Pi0toee_Mee,"LS, both legs from #pi^{0}","p");
  legTrackLFtoeeHFtoe->AddEntry(proj_recLS_LFtoee_Mee,"LS, both legs from LF","p");
  legTrackLFtoeeHFtoe->AddEntry(proj_recLS_HFtoe_Mee,"LS, at least one HF leg","p");

  auto legGenGenSmearRecElePos_top = new TLegend(legPosTrack[0]+0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  legGenGenSmearRecElePos_top->SetBorderSize(0);
  legGenGenSmearRecElePos_top->SetFillStyle(0);
  legGenGenSmearRecElePos_top->AddEntry(ptGenTrackElePos,"gen tracks","p");
  legGenGenSmearRecElePos_top->AddEntry(ptGenSmearedTrackElePos,"gen smeared tracks","p");
  legGenGenSmearRecElePos_top->AddEntry(ptRecTrackElePos,"rec tracks","p");

  auto legGenGenSmearRecElePos_Plus = new TLegend(legPosBottomLeft[0]+0.07,legPosBottomLeft[1],legPosBottomLeft[2],legPosBottomLeft[3]);
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
  legPIDContamination->AddEntry(hKaonContaminationRecPt,"Kaon","p");
  legPIDContamination->AddEntry(hProtonContaminationRecPt,"Proton","p");


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

  // auto legRatiobeforePID = new TLegend(legPosTrack[0]-0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  // legRatiobeforePID->SetBorderSize(0);
  // legRatiobeforePID->SetFillStyle(0);
  // legRatiobeforePID->SetTextSize(25);
  // legRatiobeforePID->SetTextFont(43);
  // legRatiobeforePID->AddEntry(hElePionRatio_beforePID,"before PID e/#pi","p");
  // legRatiobeforePID->AddEntry(hMuonEleRatio_beforePID,"before PID #mu/e","p");
  //
  // auto legRatioafterPID = new TLegend(legPosTrack[0]-0.05,legPosTrack[1],legPosTrack[2]+0.05,legPosTrack[3]);
  // legRatioafterPID->SetBorderSize(0);
  // legRatioafterPID->SetFillStyle(0);
  // legRatioafterPID->SetTextSize(25);
  // legRatioafterPID->SetTextFont(43);
  // legRatioafterPID->AddEntry(hElePionRatio_afterPID,"after PID e/#pi","p");
  // legRatioafterPID->AddEntry(hMuonEleRatio_afterPID,"after PID #mu/e","p");



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
  TLegendEntry *entryInfo5=legendInfo->AddEntry("Layout",Form("Layout v1, |#it{#eta}_{e}| < 0.8, #it{B} = %g T",BField),"");
  // TLegendEntry *entryInfo5=legendInfo->AddEntry("Layout",Form("iTOF 19cm, |#it{#eta}_{e}| < 1.1, #it{B} = %g T",BField),"");

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

  auto legendInfoTrackEff = new TLegend(0.35, 0.77, 0.87, 0.97);
  legendInfoTrackEff->AddEntry(entryInfo0,entryInfo0->GetLabel(),"");
  legendInfoTrackEff->AddEntry(entryInfo1,entryInfo1->GetLabel(),"");
  legendInfoTrackEff->AddEntry(entryInfo5,entryInfo5->GetLabel(),"");
  legendInfoTrackEff->SetBorderSize(0);
  legendInfoTrackEff->SetFillColorAlpha(0, 0.0);
  // legendInfoTrackEff->SetTextSize(0.025);
  legendInfoTrackEff->SetTextFont(43);
  legendInfoTrackEff->SetTextSize(24);

  auto legendInfoPairEff = (TLegend*)legendInfoTrackEff->Clone();
  legendInfoPairEff->AddEntry(entryInfo3,entryInfo3->GetLabel(),"");





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

  // auto cTrackPFVectorPt = new TCanvas("cTrackPFVectorPt","cTrackPFVectorPt",800,800);
  // cTrackPFVectorPt->SetLogy();
  // cTrackPFVectorPt->SetTopMargin(0.03);
  // cTrackPFVectorPt->SetRightMargin(0.03);
  // cTrackPFVectorPt->SetLeftMargin(0.13);
  // makeHistNice(hTracks_ElePos_Gen_inPFvector_Pt,kBlack);
  // makeHistNice(hTracks_ElePos_Rec_inPFvector_Pt,kRed);
  //
  // hTracks_ElePos_Rec_inPFvector_Pt->GetYaxis()->SetTitle("N Tracks");
  // hTracks_ElePos_Rec_inPFvector_Pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  // hTracks_ElePos_Rec_inPFvector_Pt->GetXaxis()->SetRangeUser(0.0,0.5);
  // hTracks_ElePos_Gen_inPFvector_Pt->Scale(1.,"width");
  // hTracks_ElePos_Rec_inPFvector_Pt->Scale(1.,"width");
  // hTracks_ElePos_Rec_inPFvector_Pt->Draw("hist p e1 ");
  // hTracks_ElePos_Gen_inPFvector_Pt->Draw("hist p e1 same");
  //
  // cTrackPFVectorPt->SaveAs("./plots/PtTrack_prefilterVec_GenRec.png");

  TLine* line1 = new TLine(0.1,1.,10.,1.);
  line1->SetLineStyle(2);


  // plot Pt of single Track charged Pions
  ptGenTrackPionRapSel->Scale(1./(nEventsCent));
  normalizeToBinWidth(ptGenTrackPionRapSel);
  make3HistNice(ptGenTrackPionRapSel,kBlue+1);
  TH1F* hRatioPlot = (TH1F*) RatioToTheoryy(ptGenTrackPionRapSel, hPaper_chPi_Pt_fullCentralSection_stat_Err);

  auto cChPiPt = new TCanvas("cChPiPt","cChPiPt",800,800);
  // cChPiPt->SetLogy();
  cChPiPt->SetTopMargin(0.03);
  cChPiPt->SetRightMargin(0.03);
  cChPiPt->SetLeftMargin(0.13);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetLogx();
    pad1->SetLogy();
    pad1->Draw();
  cChPiPt->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0);
    // pad2->SetLogy();
    pad2->SetLogx();
    pad2->Draw();
    pad2->cd();
    pad1->SetTopMargin(pad1->GetTopMargin()*0.5);
    pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
    pad2->SetBottomMargin(pad2->GetBottomMargin()*2.0);
  pad1->cd();
  ptGenTrackPionRapSel->GetYaxis()->SetTitle("1/(#it{N}_{ev})d^{2}N/(d#it{p}_{T}dy)");
  ptGenTrackPionRapSel->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  ptGenTrackPionRapSel->GetXaxis()->SetRangeUser(0.1,10.0);
  ptGenTrackPionRapSel->GetYaxis()->SetRangeUser(0.01,1000000.0);
  ptGenTrackPionRapSel->Draw("hist p e1");
  // hPaper_chPi_Pt_firstCentralSection->Draw("hist p E1 same");
  // hPaper_chPi_Pt_secondCentralSection->Draw("hist p E1 same");
  hPaper_chPi_Pt_fullCentralSection_sys_unc->Draw("hist p E3 same");
  hPaper_chPi_Pt_fullCentralSection_sys_Err->Draw("hist p E2 same");
  hPaper_chPi_Pt_fullCentralSection_stat_Err->Draw("hist p E1 same");
  legChPionPt->Draw("same");
  pad2->cd();
  hRatioPlot->GetYaxis()->SetTitle("Ratio paper/gen");
  hRatioPlot->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hRatioPlot->GetXaxis()->SetRangeUser(0.1,10.0);
  hRatioPlot->GetXaxis()->SetLabelSize(0.1);
  hRatioPlot->GetYaxis()->SetLabelSize(0.1);
  hRatioPlot->GetXaxis()->SetTitleSize(0.1);
  hRatioPlot->GetYaxis()->SetTitleSize(0.1);
  hRatioPlot->GetXaxis()->SetTitleOffset(0.8);
  hRatioPlot->GetYaxis()->SetTitleOffset(0.4);
  hRatioPlot->Draw("c hist");
  line1->Draw("same");

  cChPiPt->SaveAs("./plots/ChPionPt.png");

  if(bwriteFileLFHFweights){
    TString centralityRange = "0-10central";
    if(bSemiCentral) centralityRange = "30-50central";
    TFile *fOutputWeightsLF = TFile::Open(Form("./corrWeights/lfe_weights_%s.root",centralityRange.Data()),"RECREATE");
    hRatioPlot->SetName("lftoe_weights_paper_data_ratio");
    hRatioPlot->Write();
    fOutputWeightsLF->Close();
  }






  if(bPlotLFHFcontributions) {
    // Plot single electrons from HF as function of pt, compaere cocktail, my study and paper
    // HF -> e
    auto cPt_HFtoe = new TCanvas("cPt_HFtoe","cPt_HFtoe",800,800);
    // cPt_HFtoe->SetLogy();
    cPt_HFtoe->SetTopMargin(0.03);
    cPt_HFtoe->SetRightMargin(0.03);
    cPt_HFtoe->SetLeftMargin(0.13);
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
      pad1->SetBottomMargin(0);
      pad1->SetLogx();
      pad1->SetLogy();
      pad1->Draw();
    cPt_HFtoe->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
      pad2->SetTopMargin(0);
      // pad2->SetLogy();
      pad2->SetLogx();
      pad2->Draw();
      pad2->cd();
      pad1->SetTopMargin(pad1->GetTopMargin()*0.5);
      pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
      pad2->SetBottomMargin(pad2->GetBottomMargin()*2.0);
    pad1->cd();
    ptGenTrackElePos_HF_rebin->Scale(1./(etacoverage * nEventsCent));
    invariantYield(ptGenTrackElePos_HF_rebin);
    normalizeToBinWidth(ptGenTrackElePos_HF_rebin);
    ptGenTrackElePos_HF_rebin->GetYaxis()->SetTitle("1/(2#pi#it{p}_{T}#it{N}_{ev})d^{2}N/(d#it{p}_{T}dy)");
    ptGenTrackElePos_HF_rebin->GetXaxis()->SetTitle("#it{p}_{T,ee} (GeV/#it{c}^{2})");
    ptGenTrackElePos_HF_rebin->GetXaxis()->SetRangeUser(0.1,10.);
    ptGenTrackElePos_HF_rebin->GetYaxis()->SetRangeUser(0.0000001,100);
    ptGenTrackElePos_HF_rebin->Draw("p e1 hist");
    hSingleElefromHF_Cocktail->Draw("hist p e1 same");
    hPaper_HFtoe_Pt_syst->Draw("p E2 same");
    hPaper_HFtoe_Pt_stat->Draw("p E1 same");
    // hPaper_HFtoe_Pt->Draw("same p");
    legHFtoe->Draw("same");
    pad2->cd();
    TH1F* hRatioCocktailPaper = (TH1F*) RatioToTheoryy(hSingleElefromHF_Cocktail,hPaper_HFtoe_Pt_stat);
    TH1F* hRatioRecHFPaper = (TH1F*) RatioToTheoryy(ptGenTrackElePos_HF_rebin,hPaper_HFtoe_Pt_stat);
    TH1F* hdummy = new TH1F("hdummy","", 1, 0.1, 10.); makeHistNice(hdummy,kBlack);
    TH1F* hRatioCocktailFastSim = (TH1F*) RatioToTheoryy(ptGenTrackElePos_HF_rebin,hSingleElefromHF_Cocktail);
    hRatioCocktailFastSim->GetXaxis()->SetRangeUser(0.,0.6);
    makeHistNice(hRatioCocktailFastSim,kBlue);
    hRatioCocktailFastSim->SetTitle("hftoe_weights_cocktail_data_ratio");
    hRatioCocktailFastSim->SetName("hftoe_weights_cocktail_data_ratio");
    hRatioRecHFPaper->SetTitle("hftoe_weights_paper_data_ratio");
    hRatioRecHFPaper->SetName("hftoe_weights_paper_data_ratio");
    // hRatioRecHFPaper->Add(hRatioCocktailFastSim);
    makeHistNice(hRatioCocktailPaper,kRed);
    makeHistNice(hRatioRecHFPaper,kBlue);
    hdummy->GetYaxis()->SetTitle("Ratio PLB 804/X");
    hdummy->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hdummy->GetXaxis()->SetRangeUser(0.1,10.0);
    hdummy->GetYaxis()->SetRangeUser(0.,2.);
    hdummy->GetXaxis()->SetLabelSize(0.1);
    hdummy->GetYaxis()->SetLabelSize(0.1);
    hdummy->GetXaxis()->SetTitleSize(0.1);
    hdummy->GetYaxis()->SetTitleSize(0.1);
    hdummy->GetXaxis()->SetTitleOffset(0.8);
    hdummy->GetYaxis()->SetTitleOffset(0.4);
    hdummy->Draw("");
    line1->Draw("same");
    hRatioCocktailFastSim->Draw("p e1 same");
    hRatioRecHFPaper->Draw("p e1 same");
    hRatioCocktailPaper->Draw("p e1 same");
    cPt_HFtoe->SaveAs("./plots/Pte_HFtoe+Cocktail+Paper.png");

    if(bwriteFileLFHFweights){
      TString centralityRange = "0-10central";
      if(bSemiCentral) centralityRange = "30-50central";
      TFile *fOutputWeights = TFile::Open(Form("./corrWeights/hfe_weights_%s.root",centralityRange.Data()),"RECREATE");
      hRatioCocktailFastSim->Write(); // below 0.6 GeV/c
      hRatioRecHFPaper->Write();  // above 0.6 GeV/c
      fOutputWeights->Close();
    }



    // // plot the ratio (π0->e  +  c,b->e)/(π+ + π-) , using generated π0->e and data of c,b->e and π+ + π-
    // auto cRatioPi0cbToeOverChPion = new TCanvas("cRatioPi0cbToeOverChPion","cRatioPi0cbToeOverChPion",800,800);
    // cRatioPi0cbToeOverChPion->SetTopMargin(0.03);
    // cRatioPi0cbToeOverChPion->SetRightMargin(0.03);
    //
    // TH1F* ptGenTrackElePos_HF_lowPt;
    //
    // pad1->SetBottomMargin(0.);
    // pad1->SetLeftMargin(0.13);
    // pad1->SetRightMargin(0.02);
    // pad1->SetLogx(false);
    // pad1->SetLogy();
    // pad1->SetTicks();
    // pad1->Draw();
    // cRatioPi0cbToeOverChPion->cd();
    //
    // pad2->SetTopMargin(0);
    // pad2->SetLeftMargin(0.13);
    // pad2->SetRightMargin(0.02);
    // pad2->SetLogx(false);
    // pad2->SetTicks();
    // pad2->Draw();
    // pad2->cd();
    // pad1->SetTopMargin(0.02);
    // pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
    // pad2->SetBottomMargin(0.25);
    // pad1->cd();
    //
    //
    //
    // // for c,b->e data only down to pt > 0.6 GeV/c, below use corrected fast simulation
    // ptGenTrackElePos_Pi0_rebin->Scale(1./(etacoverage * nEventsCent));
    // normalizeToBinWidth(ptGenTrackElePos_Pi0_rebin);
    // // invariantYield(ptGenTrackElePos_Pi0_rebin);
    // // hPaper_HFtoe_Pt_stat->GetXaxis()->SetRangeUser(0.1,10.);
    // hPaper_HFtoe_Pt_stat->GetXaxis()->SetRangeUser(0.1,4.);
    // ptGenTrackElePos_HF_lowPt = (TH1F*) ptGenTrackElePos_HF_rebin->Clone();
    // ptGenTrackElePos_HF_lowPt->GetXaxis()->SetRangeUser(0.,0.5);
    // ClearBinsAbove(ptGenTrackElePos_HF_lowPt,0.5);
    // reverseInvariantYield(ptGenTrackElePos_HF_lowPt);
    // reverseInvariantYield(hPaper_HFtoe_Pt_stat);
    //
    // TH1F* hNominator = (TH1F*) ptGenTrackElePos_Pi0_rebin->Clone();
    // hNominator->Add(ptGenTrackElePos_HF_lowPt,1.);
    // AddHistosDiffRange(hNominator,hPaper_HFtoe_Pt_stat);
    // hNominator->SetMarkerColor(kRed+1);
    // hNominator->SetLineColor(kRed+1);
    // // invariantYield(hPaper_chPi_Pt_fullCentralSection_stat_Err);
    //
    // hPaper_chPi_Pt_fullCentralSection_stat_Err->SetMarkerColor(kGreen);
    // hPaper_chPi_Pt_fullCentralSection_stat_Err->SetLineColor(kGreen);
    // TH1F* hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej1000 = (TH1F*) hPaper_chPi_Pt_fullCentralSection_stat_Err->Clone();
    // TH1F* hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000 = (TH1F*) hPaper_chPi_Pt_fullCentralSection_stat_Err->Clone();
    // makeHistNice(hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej1000, kGreen+2);
    // makeHistNice(hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000, kGreen+4);
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej1000->Scale(1./1000.);
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->Scale(1./5000.);
    // ptGenTrackElePos_HF_lowPt->SetMarkerStyle(20);
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->SetMinimum(0.0000001);
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->SetMaximum(1000000);
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->GetYaxis()->SetTitleOffset(1.15);
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->GetYaxis()->SetTitle("1/#it{N}_{ev} d^{2}N/(d#it{p}_{T}dy)");
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->GetXaxis()->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->GetXaxis()->SetRangeUser(0.1,4.);
    //
    // TF1 *fPSeff = new TF1("effPS_inverse","0.8*(1.-exp(-1.6*(x-0.05)))",0,10);
    // TF1 *fPSeff_inverse = new TF1("effPS_inverse","1./( 0.8*(1.-exp(-1.6*(x-0.05))) )",0,10);
    // TH1F* hNominator_scaled = (TH1F*) hNominator->Clone();
    // hNominator_scaled->Divide(fPSeff_inverse);
    // hNominator_scaled->SetMarkerColor(kRed+2);
    // hNominator_scaled->SetLineColor(kRed+2);
    //
    // auto legPi0cbToeOverChPion = new TLegend(0.5,0.75,0.8,0.96);
    // legPi0cbToeOverChPion->SetBorderSize(0);
    // legPi0cbToeOverChPion->SetFillStyle(0);
    // legPi0cbToeOverChPion->SetTextFont(43);
    // legPi0cbToeOverChPion->SetTextSize(20);
    // // legPi0cbToeOverChPion->AddEntry(ptGenTrackElePos_Pi0_rebin,"fast sim, #pi^{0} #rightarrow e","p");
    // // legPi0cbToeOverChPion->AddEntry(ptGenTrackElePos_HF_lowPt,"fast sim, c,b #rightarrow e, #it{p}_{T}<0.5 GeV/#it{c}","p");
    // // legPi0cbToeOverChPion->AddEntry(hPaper_HFtoe_Pt_stat,"paper, c,b #rightarrow e, #it{p}_{T}>0.5 GeV/#it{c}","p");
    // legPi0cbToeOverChPion->AddEntry(hPaper_chPi_Pt_fullCentralSection_stat_Err,"paper, #pi^{+} + #pi^{-}","p");
    // legPi0cbToeOverChPion->AddEntry(hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej1000,"paper, #pi^{+} + #pi^{-} #times (1/1000)","p");
    // legPi0cbToeOverChPion->AddEntry(hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000,"paper, #pi^{+} + #pi^{-} #times (1/5000)","p");
    // legPi0cbToeOverChPion->AddEntry(hNominator,"#pi^{0} #rightarrow e + c,b #rightarrow e","p");
    // legPi0cbToeOverChPion->AddEntry(hNominator_scaled,"eff_{e^{#pm}} #times (#pi^{0} #rightarrow e + c,b #rightarrow e)","p");
    //
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000->Draw("hist p e1");
    // hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej1000->Draw("hist p e1 same");
    // hPaper_chPi_Pt_fullCentralSection_stat_Err->Draw("hist p e1 same");
    // // ptGenTrackElePos_Pi0_rebin->Draw("hist p e1 same");
    // // hPaper_HFtoe_Pt_stat->Draw("hist p e1 same");
    // // ptGenTrackElePos_HF_lowPt->Draw("hist p e1 same");
    // hNominator->Draw("same hist p e1");
    // hNominator_scaled->Draw("same hist p e1");
    // fPSeff->Draw("same hist c");
    // legPi0cbToeOverChPion->Draw("same");
    //
    // pad2->cd();
    // // pad2->SetLogy();
    //
    // // TH1F* hRatioPi0cbToeOverChPion = (TH1F*) RatioToTheoryy(hPaper_chPi_Pt_fullCentralSection_stat_Err,hNominator);
    // TH1F* hRatioPi0cbToeOverChPion = (TH1F*) RatioToTheoryy(hNominator_scaled,hPaper_chPi_Pt_fullCentralSection_stat_Err);
    // TH1F* hRatioPi0cbToeOverChPion_PionRej1000 = (TH1F*) RatioToTheoryy(hNominator_scaled,hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej1000);
    // TH1F* hRatioPi0cbToeOverChPion_PionRej5000 = (TH1F*) RatioToTheoryy(hNominator_scaled,hPaper_chPi_Pt_fullCentralSection_stat_Err_PionRej5000);
    // // hRatioPi0cbToeOverChPion->GetXaxis()->SetRangeUser(0.1,10.);
    // hRatioPi0cbToeOverChPion->GetXaxis()->SetRangeUser(0.1,4.);
    // // hRatioPi0cbToeOverChPion->GetYaxis()->SetRangeUser(0.,100.);
    // // TString sYRatioLabel = "#frac{#pi^{0} #rightarrow e + c,b #rightarrow e}{#pi^{+} + #pi^{-}}";
    // TString sYRatioLabel = "#frac{#pi^{+} + #pi^{-}}{#pi^{0} #rightarrow e + c,b #rightarrow e}";
    // TString sXRatioLabel = "#it{p}_{T,e} (GeV/#it{c})";
    // // TLine* line10 = new TLine(0.,0.1,10.,0.1);
    // TLine* line10 = new TLine(0.1,0.1,4.,0.1);
    // line10->SetLineStyle(2);
    // makeRatioNice(hRatioPi0cbToeOverChPion,sXRatioLabel,sYRatioLabel);
    // hRatioPi0cbToeOverChPion->GetXaxis()->SetTitleOffset(1.08);
    // hRatioPi0cbToeOverChPion->GetYaxis()->SetTitleOffset(0.63);
    // // hRatioPi0cbToeOverChPion->GetYaxis()->SetRangeUser(0.001,3000.);
    // hRatioPi0cbToeOverChPion->GetYaxis()->SetRangeUser(0.,0.8);
    // hRatioPi0cbToeOverChPion->Draw();
    // hRatioPi0cbToeOverChPion_PionRej1000->Draw("same");
    // hRatioPi0cbToeOverChPion_PionRej5000->Draw("same");
    // line10->Draw("same");
    // cRatioPi0cbToeOverChPion->SaveAs("./plots/RatioPi0cbToeOverChPion.png");




    auto cPtElePosPi0LFHF = new TCanvas("cPtElePosPi0LFHF","cPtElePosPi0LFHF",800,800);
    cPtElePosPi0LFHF->SetTopMargin(0.03);
    cPtElePosPi0LFHF->SetRightMargin(0.03);
    cPtElePosPi0LFHF->SetLeftMargin(0.13);
    // TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
      pad1->SetBottomMargin(0);
      pad1->SetLogy();
      pad1->Draw();
    cPtElePosPi0LFHF->cd();
    // TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
      pad2->SetTopMargin(0);
      // pad2->SetLogy();
      pad2->Draw();
      pad2->cd();
      // pad1->SetTopMargin(pad1->GetTopMargin()*0.5);
      // pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
      // pad2->SetBottomMargin(pad2->GetBottomMargin()*2.0);
    pad1->cd();
    ptRecTrackElePos_rebin->GetYaxis()->SetTitle("dN/d#it{p}_{T,e}");
    ptRecTrackElePos_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptRecTrackElePos_rebin->GetXaxis()->SetRangeUser(0.06,4.0);
    ptRecTrackElePos_rebin->GetYaxis()->SetRangeUser(1.0,100000000);
    ptRecTrackElePos_rebin->SetMarkerColor(kBlack);
    ptRecTrackElePos_rebin->SetLineColor(kBlack);
    ptRecTrackElePos_rebin->Draw("hist p e1");
    ptRecTrackElePos_LF_afterPID_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_LF_beforePID_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_Pi0_afterPID_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_Pi0_beforePID_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_HF_afterPID_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_HF_beforePID_rebin->Draw("hist p e1 same");
    legendInfo->Draw("same");
    legTrackPi0LFHF->Draw("same");
    pad2->cd();
    TH1F* ratioTotalHF = (TH1F*) ptRecTrackElePos_LF_afterPID_rebin->Clone();
    ratioTotalHF->Add(ptRecTrackElePos_HF_afterPID_rebin);
    ratioTotalHF->Divide(ptRecTrackElePos_HF_afterPID_rebin);
    // ratioTotalHF->GetYaxis()->SetRangeUser(0.5,300);
    ratioTotalHF->GetXaxis()->SetRangeUser(0.06,4.0);
    ratioTotalHF->GetYaxis()->SetTitle("Ratio (LF+HF)/HF");
    ratioTotalHF->GetXaxis()->SetLabelSize(0.1);
    ratioTotalHF->GetYaxis()->SetLabelSize(0.1);
    ratioTotalHF->GetXaxis()->SetTitleSize(0.1);
    ratioTotalHF->GetYaxis()->SetTitleSize(0.1);
    ratioTotalHF->GetXaxis()->SetTitleOffset(0.8);
    ratioTotalHF->GetYaxis()->SetTitleOffset(0.4);
    ratioTotalHF->Draw();
    cPtElePosPi0LFHF->SaveAs("./plots/PtElePosPi0LFHF.png");

    ptRecTrackElePos_rebin->SetLineColor(kBlue+1);
    ptRecTrackElePos_rebin->SetMarkerColor(kBlue+1);

    auto cMeePteeDCACharmBeautyHFtoe = new TCanvas("cMeePteeDCACharmBeautyHFtoe","cMeePteeDCACharmBeautyHFtoe",1350,600);
    cMeePteeDCACharmBeautyHFtoe->Divide(3,1);
    for (size_t i = 1; i < 4; i++) {
      cMeePteeDCACharmBeautyHFtoe->cd(i)->SetTopMargin(0.03);
      cMeePteeDCACharmBeautyHFtoe->cd(i)->SetRightMargin(0.03);
    }
    cMeePteeDCACharmBeautyHFtoe->cd(1);
    cMeePteeDCACharmBeautyHFtoe->cd(1)->SetLogy();
    proj_recLS_HFtoe_Mee->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_HFtoe_Mee->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->e}/d#it{m}_{ee}dy");
    proj_recLS_HFtoe_Mee->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    proj_recLS_HFtoe_Mee->Draw("hist p e1");
    proj_recLS_Ctoe_Mee->Draw("hist p e1 same");
    proj_recLS_Btoe_Mee->Draw("hist p e1 same");
    legTrackLFHFtoe->Draw("same");

    cMeePteeDCACharmBeautyHFtoe->cd(2);
    cMeePteeDCACharmBeautyHFtoe->cd(2)->SetLogy();
    proj_recLS_HFtoe_Ptee->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_HFtoe_Ptee->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->e}/d#it{p}_{T,ee}dy");
    proj_recLS_HFtoe_Ptee->GetXaxis()->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    proj_recLS_HFtoe_Ptee->GetXaxis()->SetRangeUser(0.0,4.0);
    proj_recLS_HFtoe_Ptee->Draw("hist p e1");
    proj_recLS_Ctoe_Ptee->Draw("hist p e1 same");
    proj_recLS_Btoe_Ptee->Draw("hist p e1 same");

    cMeePteeDCACharmBeautyHFtoe->cd(3);
    cMeePteeDCACharmBeautyHFtoe->cd(3)->SetLogy();
    proj_recLS_HFtoe_DCA->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_HFtoe_DCA->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->e}/dDCAdy");
    proj_recLS_HFtoe_DCA->GetXaxis()->SetTitle("DCA_{ee} (cm)");
    proj_recLS_HFtoe_DCA->GetXaxis()->SetRangeUser(0.0,10.0);
    proj_recLS_HFtoe_DCA->Draw("hist p e1");
    proj_recLS_Ctoe_DCA->Draw("hist p e1 same");
    proj_recLS_Btoe_DCA->Draw("hist p e1 same");

    cMeePteeDCACharmBeautyHFtoe->SaveAs("./plots/MeePteeDCA_CharmBeautyHFtoe.png");



    auto cMeePteeDCACharmBeautyHFtoee = new TCanvas("cMeePteeDCACharmBeautyHFtoee","cMeePteeDCACharmBeautyHFtoee",1350,600);
    cMeePteeDCACharmBeautyHFtoee->SetLogy();
    cMeePteeDCACharmBeautyHFtoee->Divide(3,1);
    for (size_t i = 1; i < 4; i++) {
      cMeePteeDCACharmBeautyHFtoee->cd(i)->SetTopMargin(0.03);
      cMeePteeDCACharmBeautyHFtoee->cd(i)->SetRightMargin(0.03);
    }
    cMeePteeDCACharmBeautyHFtoee->cd(1);
    cMeePteeDCACharmBeautyHFtoee->cd(1)->SetLogy();
    proj_recLS_LFtoee_Mee->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_LFtoee_Mee->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->ee}/d#it{m}_{ee}dy");
    proj_recLS_LFtoee_Mee->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    proj_recLS_LFtoee_Mee->Draw("hist p e1");
    proj_recLS_HFtoee_Mee->Draw("hist p e1 same");
    legTrackLFHFtoee->Draw("same");

    cMeePteeDCACharmBeautyHFtoee->cd(2);
    cMeePteeDCACharmBeautyHFtoee->cd(2)->SetLogy();
    proj_recLS_LFtoee_Ptee->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_LFtoee_Ptee->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->ee}/d#it{p}_{T,ee}dy");
    proj_recLS_LFtoee_Ptee->GetXaxis()->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    proj_recLS_LFtoee_Ptee->GetXaxis()->SetRangeUser(0.0,4.0);
    proj_recLS_LFtoee_Ptee->Draw("hist p e1");
    proj_recLS_HFtoee_Ptee->Draw("hist p e1 same");

    cMeePteeDCACharmBeautyHFtoee->cd(3);
    cMeePteeDCACharmBeautyHFtoee->cd(3)->SetLogy();
    proj_recLS_LFtoee_DCA->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_LFtoee_DCA->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->ee}/dDCAdy");
    proj_recLS_LFtoee_DCA->GetXaxis()->SetTitle("DCA_{ee} (cm)");
    proj_recLS_LFtoee_DCA->GetXaxis()->SetRangeUser(0.0,10.0);
    proj_recLS_LFtoee_DCA->Draw("hist p e1");
    proj_recLS_HFtoee_DCA->Draw("hist p e1 same");

    cMeePteeDCACharmBeautyHFtoee->SaveAs("./plots/MeePteeDCA_LFHFtoee.png");



    auto cMeePteeDCALFtoeeHFtoe = new TCanvas("cMeePteeDCALFtoeeHFtoe","cMeePteeDCALFtoeeHFtoe",1350,600);
    cMeePteeDCALFtoeeHFtoe->SetLogy();
    cMeePteeDCALFtoeeHFtoe->Divide(3,1);
    for (size_t i = 1; i < 4; i++) {
      cMeePteeDCALFtoeeHFtoe->cd(i)->SetTopMargin(0.03);
      cMeePteeDCALFtoeeHFtoe->cd(i)->SetRightMargin(0.03);
    }
    cMeePteeDCALFtoeeHFtoe->cd(1);
    cMeePteeDCALFtoeeHFtoe->cd(1)->SetLogy();
    proj_recLS_MCpidEle_Mee->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_MCpidEle_Mee->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->ee}/d#it{m}_{ee}dy");
    proj_recLS_MCpidEle_Mee->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    proj_recLS_MCpidEle_Mee->Draw("hist p e1");
    proj_recLS_LFtoe_Mee->Draw("hist p e1 same");
    proj_recLS_Pi0toee_Mee->Draw("hist p e1 same");
    proj_recLS_LFtoee_Mee->Draw("hist p e1 same");
    proj_recLS_HFtoe_Mee->Draw("hist p e1 same");
    legTrackLFtoeeHFtoe->Draw("same");

    cMeePteeDCALFtoeeHFtoe->cd(2);
    cMeePteeDCALFtoeeHFtoe->cd(2)->SetLogy();
    proj_recLS_MCpidEle_Ptee->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_MCpidEle_Ptee->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->ee}/d#it{p}_{T,ee}dy");
    proj_recLS_MCpidEle_Ptee->GetXaxis()->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    proj_recLS_MCpidEle_Ptee->GetXaxis()->SetRangeUser(0.0,4.0);
    proj_recLS_MCpidEle_Ptee->Draw("hist p e1");
    proj_recLS_LFtoe_Ptee->Draw("hist p e1 same");
    proj_recLS_LFtoee_Ptee->Draw("hist p e1 same");
    proj_recLS_HFtoe_Ptee->Draw("hist p e1 same");

    cMeePteeDCALFtoeeHFtoe->cd(3);
    cMeePteeDCALFtoeeHFtoe->cd(3)->SetLogy();
    proj_recLS_MCpidEle_DCA->GetYaxis()->SetTitleOffset(1.);
    proj_recLS_MCpidEle_DCA->GetYaxis()->SetTitle("1/#it{N}_{ev}d^{2}#it{N}_{X->ee}/dDCAdy");
    proj_recLS_MCpidEle_DCA->GetXaxis()->SetTitle("DCA_{ee} (cm)");
    proj_recLS_MCpidEle_DCA->GetXaxis()->SetRangeUser(0.0,10.0);
    proj_recLS_MCpidEle_DCA->Draw("hist p e1");
    proj_recLS_LFtoe_DCA->Draw("hist p e1 same");
    proj_recLS_LFtoee_DCA->Draw("hist p e1 same");
    proj_recLS_HFtoe_DCA->Draw("hist p e1 same");

    cMeePteeDCALFtoeeHFtoe->SaveAs("./plots/MeePteeDCA_LFtoeeHFtoe.png");

    auto cPF_opAngle_mass = new TCanvas("cTrackPFVectorPt","cTrackPFVectorPt",800,800);
    // auto cPF_opAngle_mass = new TCanvas("cTrackPFVectorPt","cTrackPFVectorPt",1900,800);
    // cPF_opAngle_mass->Divide(2,1);
    // for (size_t i = 1; i < 4; i++) {
    //   cPF_opAngle_mass->cd(i)->SetTopMargin(0.03);
    //   cPF_opAngle_mass->cd(i)->SetRightMargin(0.03);
    // }
    // cPF_opAngle_mass->cd(1);
    // cPF_opAngle_mass->cd(1)->SetLogz();
    cPF_opAngle_mass->cd(1)->SetLogz();
    makeHistNiceTH2(hMeeOpAngle_RecULS_bPF,kBlack);
    makeHistNiceTH2(hMeeOpAngle_RecULS_aPF,kBlack);
    // hMeeOpAngle_RecULS_bPF->GetYaxis()->SetTitleOffset(1.);
    hMeeOpAngle_RecULS_bPF->GetYaxis()->SetTitle("#it{#omega}_{ee} (rad)");
    hMeeOpAngle_RecULS_bPF->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    hMeeOpAngle_RecULS_bPF->GetXaxis()->SetRangeUser(0.0,0.15);
    hMeeOpAngle_RecULS_bPF->Draw("colz");
    // cPF_opAngle_mass->cd(2);
    // cPF_opAngle_mass->cd(2)->SetLogz();
    // hMeeOpAngle_RecULS_aPF->GetYaxis()->SetTitleOffset(1.);
    // hMeeOpAngle_RecULS_aPF->GetYaxis()->SetTitle("opening angle (rad)");
    // hMeeOpAngle_RecULS_aPF->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    // hMeeOpAngle_RecULS_aPF->GetXaxis()->SetRangeUser(0.0,0.15);
    // hMeeOpAngle_RecULS_aPF->Draw("colz");
    cPF_opAngle_mass->SaveAs(Form("./plots/PF_BeforeAfter_opAngle_Mass_sce%i.png",ith_PIDscenario));


    // generate root file with LS background mee ptee and  LF->ee , HF->e and Sum
    if (bGenerateBackground) {
      TString sDCAcut;
      if(DCAcut) sDCAcut = "wDCAcut";
      else if(!DCAcut) sDCAcut = "woDCAcut";
      // TFile *fBackground_weights = TFile::Open("./plots/MeePtee_LFtoeeHFtoe_noweight.root","RECREATE");
      TFile *fBackground_weights = TFile::Open(Form("./plots/plottingMacros/MeePtee_recLS_wPreFilter_%s_sce%i.root",sDCAcut.Data(),ith_PIDscenario),"RECREATE");
      proj_recLS_MCpidEle_Mee->SetName("proj_recLS_MCpidEle_Mee");
      proj_recLS_LFtoee_Mee->SetName("proj_recLS_LFtoee_Mee");
      proj_recLS_HFtoe_Mee->SetName("proj_recLS_HFtoe_Mee");
      proj_recLS_HFtoee_Mee->SetName("proj_recLS_HFtoee_Mee");
      proj_recLS_MCpidEle_Ptee->SetName("proj_recLS_MCpidEle_Ptee");
      proj_recLS_LFtoee_Ptee->SetName("proj_recLS_LFtoee_Ptee");
      proj_recLS_HFtoe_Ptee->SetName("proj_recLS_HFtoe_Ptee");
      proj_recLS_HFtoee_Ptee->SetName("proj_recLS_HFtoee_Ptee");
      proj_recLS_MCpidEle_Mee->SetMarkerStyle(24);
      proj_recLS_LFtoee_Mee->SetMarkerStyle(24);
      proj_recLS_HFtoe_Mee->SetMarkerStyle(24);
      proj_recLS_HFtoee_Mee->SetMarkerStyle(24);
      proj_recLS_MCpidEle_Ptee->SetMarkerStyle(24);
      proj_recLS_LFtoee_Ptee->SetMarkerStyle(24);
      proj_recLS_HFtoe_Ptee->SetMarkerStyle(24);
      proj_recLS_HFtoee_Ptee->SetMarkerStyle(24);
      proj_recLS_MCpidEle_Mee->Write();
      proj_recLS_LFtoee_Mee->Write();
      proj_recLS_HFtoe_Mee->Write();
      proj_recLS_HFtoee_Mee->Write();
      proj_recLS_MCpidEle_Ptee->Write();
      proj_recLS_LFtoee_Ptee->Write();
      proj_recLS_HFtoe_Ptee->Write();
      proj_recLS_HFtoee_Ptee->Write();
      hMeeOpAngle_RecULS_bPF->Write("beforePF_opAngle_mass");
      hMeeOpAngle_RecULS_aPF->Write("afterPF_opAngle_mass");
      fBackground_weights->Close();
      printf(" MeePtee_recLS_wPreFilter_%s_sce%i.root has been written. \n",sDCAcut.Data(),ith_PIDscenario);
    }


    // Plot comparison between cocktail and ULS-LS (charm, beauty, hf)
    // LF, HF -> e
    auto cSignalCBHFtoe = new TCanvas("cSignalCBHFtoe","cSignalCBHFtoe",800,800);
    // cChPiPt->SetLogy();
    cSignalCBHFtoe->SetTopMargin(0.03);
    cSignalCBHFtoe->SetRightMargin(0.03);
    cSignalCBHFtoe->SetLeftMargin(0.13);
    // cSignalCBHFtoe->SetLogx();
    cSignalCBHFtoe->SetLogy();
    hSignal_HFTOe_Mee->GetYaxis()->SetTitle("1/(#it{N}_{ev})d^{2}N/(d#it{p}_{T}dy)");
    hSignal_HFTOe_Mee->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    hSignal_HFTOe_Mee->Draw("hist p E1 ");
    hSignal_charmTOe_Mee->Draw("hist p E1 same");
    hSignal_beautyTOe_Mee->Draw("hist p E1 same");
    legSigCBHFtoe->Draw("same");
    cSignalCBHFtoe->SaveAs("./plots/Mee_CBHFtoe+Cocktail.png");

    // LF -> ee
    auto cSignalLFtoee= new TCanvas("cSignalLFtoee","cSignalLFtoee",800,800);
    // cChPiPt->SetLogy();
    cSignalLFtoee->SetTopMargin(0.03);
    cSignalLFtoee->SetRightMargin(0.03);
    cSignalLFtoee->SetLeftMargin(0.13);
    // cSignalLFtoee->SetLogx();
    cSignalLFtoee->SetLogy();
    hSignal_LFTOee_Mee_rebin->GetYaxis()->SetTitle("1/(#it{N}_{ev})d^{2}N/(d#it{p}_{T}dy)");
    hSignal_LFTOee_Mee_rebin->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    // hSignal_LFTOee_Mee->GetYaxis()->SetRangeUser(0.00001,1.);
    // proj_recULS_MCpidEle_lfTOee_Mee->Draw("c");
    // proj_recLS_LFtoee_Mee->Draw("c same");
    ScaleBinWidth(hSignal_LFTOee_Mee_rebin,kTRUE);
    hSignal_LFTOee_Mee_rebin->Draw("c");
    hCocktailPi0toee->Draw("c same");
    hCocktailEtatoee->Draw("c same");
    hCocktailEtaPrimetoee->Draw("c same");
    hCocktailOmegatoee->Draw("c same");
    hCocktailPhitoee->Draw("c same");
    hTotalLFCocktailtoee->Draw("c same");
    legSigLFtoee->Draw("same");
    cSignalLFtoee->SaveAs("./plots/Mee_LFtoee+Cocktail.png");

    // HF -> ee
    auto cSignalCBHFtoee= new TCanvas("cSignalCBHFtoee","cSignalCBHFtoee",800,800);
    // cChPiPt->SetLogy();
    cSignalCBHFtoee->SetTopMargin(0.03);
    cSignalCBHFtoee->SetRightMargin(0.03);
    cSignalCBHFtoee->SetLeftMargin(0.13);
    // cSignalCBHFtoee->SetLogx();
    // cSignalCBHFtoee->SetLogy();
    hSignal_HFTOee_Mee->GetYaxis()->SetTitle("1/(#it{N}_{ev})d^{2}N/(d#it{p}_{T}dy)");
    hSignal_HFTOee_Mee->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    // hSignal_HFTOee_Mee->GetYaxis()->SetRangeUser(0.00001,1.);
    hSignal_HFTOee_Mee->Draw("c ");
    hSignal_ccTOee_Mee->Draw("c same");
    hSignal_bbTOee_Mee->Draw("c same");
    hCocktailCCtoee->Draw("c same");
    hCocktailBBtoee->Draw("c same");
    legSigCBHFtoee->Draw("same");
    cSignalCBHFtoee->SaveAs("./plots/Mee_CCBBHFtoee+Cocktail.png");


  }





  if (bPlotEfficiency) {

    auto cEffAddedElePosPt = new TCanvas("cEffAddedElePosPt","cEffAddedElePosPt",800,800);
    // cEffAddedElePosPt->SetLogy();
    cEffAddedElePosPt->SetTopMargin(0.03);
    cEffAddedElePosPt->SetRightMargin(0.03);
    cEffAddedElePosPt->SetLeftMargin(0.13);
    ptEffElePosGen->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffElePosGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // ptEffElePosGen->GetXaxis()->SetRangeUser(0.0,0.4);
    ptEffElePosGen->GetXaxis()->SetRangeUser(0.0,4.0);
    // ptEffElePosGen->SetMaximum(1.1*std::max({ptEffElePosGen->GetMaximum(),ptEffElePosGenSmeared->GetMaximum()}));
    ptEffElePosGen->SetMaximum(1.2);
    // ptEffElePosGen->SetMaximum(1.1*std::max({ptEffElePosGen->GetMaximum()}));
    // ptEffElePosGenSmeared->SetMaximum(1.1*std::max({ptEffElePosGenSmeared->GetMaximum()}));
    // ptEffElePosGen->SetMinimum(0.);
    // ptEffElePosGen->GetYaxis()->SetRangeUser(0.,1.);
    ptEffElePosGen->Draw("hist p e1");
    ptEffElePosGenSmeared->Draw("hist p e1 same");
    legEff_GenGenSmear->Draw("same");
    textBField->Draw("same");
    // textPIDScenario->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_GenGenSmear.png");
    ptEffElePosGen->GetYaxis()->SetRangeUser(0.,1.);
    ptEffElePosGen->Draw("hist p e1 ");
    legEff_Gen->Draw("same");
    textBField->Draw("same");
    // textPIDScenario->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_Gen.png");
    ptEffElePosGen_woPID->Draw("hist p e1 ");
    legEff_Gen_woPID->Draw("same");
    textBField->Draw("same");
    // textPIDScenario->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_Gen_woPID.png");
    ptEffElePos_PreShower->GetYaxis()->SetRangeUser(0.,1.);
    ptEffElePos_PreShower->Draw("hist p e1");
    ptEffElePos_PreShower->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffElePos_PreShower->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    funcPSparam->Draw("same");
    legEff_PreShower->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_PreShower.png");
    ptEffElePos_LF_afterPF->GetXaxis()->SetRangeUser(0.,4.);
    ptEffElePos_LF_afterPF->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffElePos_LF_afterPF->Draw("hist p e1");
    legEff_LF_afterPF->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_LF_afterPF.png");
    ptEffElePos_HF_afterPF->GetXaxis()->SetRangeUser(0.,4.);
    ptEffElePos_HF_afterPF->GetYaxis()->SetTitle("#it{p}_{T}^{rec}/#it{p}_{T}^{gen}");
    ptEffElePos_HF_afterPF->Draw("hist p e1");
    legEff_HF_afterPF->Draw("same");
    cEffAddedElePosPt->SaveAs("./plots/Eff_ElePos_Pt_HF_afterPF.png");



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
    // ptGenTrackElePos_rebin->GetXaxis()->SetRangeUser(0.,0.5);
    ptGenTrackElePos_rebin->Draw("hist p e1");
    ptGenSmearedTrackElePos_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_rebin->Draw("hist p e1 same");
    textBField->Draw("same");
    // textPIDScenario->Draw("same");
    legGenGenSmearRecElePos_Plus->Draw("same");
    ptGenTrackElePos_beforeKineCuts_rebin->Draw("hist p e1 same");
    // ptGenSmearedTrackElePos_beforeKineCuts_rebin->Draw("hist p e1 same");
    ptRecTrackBeforeSmearing_rebin->Draw("hist p e1 same");
    ptRecTrackAfterSmearing_rebin->Draw("hist p e1 same");
    ptRecTrackAfterKineCuts_rebin->Draw("hist p e1 same");
    ptRecTrackElePos_rebin->Draw("hist p e1 same");
    cGenGenSmearRecElePosPt->SaveAs("./plots/GenGenSmearRecElePos_Pt_0-0.5+woPIDcuts.png");



    auto cEffElePosPtEta = new TCanvas("cEffElePosPtEta","cEffElePosPtEta",800,800);
    gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
    // cEffElePosPtEta->SetLogz();
    cEffElePosPtEta->SetLogx();
    cEffElePosPtEta->SetTopMargin(0.03);
    cEffElePosPtEta->SetRightMargin(0.13);
    cEffElePosPtEta->SetLeftMargin(0.13);
    ptEtaEffElePos->SetTitle("");
    ptEtaEffElePos->GetYaxis()->SetTitle("#eta");
    ptEtaEffElePos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    ptEtaEffElePos->GetZaxis()->SetTitle("track efficiency");
    ptEtaEffElePos->GetZaxis()->SetTitleOffset(1.2);
    ptEtaEffElePos->GetZaxis()->RotateTitle(kTRUE);
    // ptEtaEffElePos->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEtaEffElePos->GetXaxis()->SetRangeUser(0.0,0.6);
    ptEtaEffElePos->GetYaxis()->SetRangeUser(-1.2,1.8);
    ptEtaEffElePos->GetZaxis()->SetRangeUser(0,1.0);
    ptEtaEffElePos->Draw("Colz");
    legendInfoTrackEff->Draw("same");
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


    auto cPairEffElePosMPt = new TCanvas("cPairEffElePosMPt","cPairEffElePosMPt",800,800);
    gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
    // cPairEffElePosMPt->SetLogz();
    // cPairEffElePosMPt->SetLogy();
    cPairEffElePosMPt->SetTopMargin(0.03);
    cPairEffElePosMPt->SetRightMargin(0.13);
    cPairEffElePosMPt->SetLeftMargin(0.13);
    mPtPairEffElePos->SetTitle("");
    mPtPairEffElePos->GetYaxis()->SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    mPtPairEffElePos->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    mPtPairEffElePos->GetZaxis()->SetTitle("pair efficiency");
    mPtPairEffElePos->GetZaxis()->RotateTitle(kTRUE);
    // mPtPairEffElePos->GetYaxis()->SetRangeUser(0.,4.);
    mPtPairEffElePos->GetXaxis()->SetRangeUser(0.0,0.7);
    mPtPairEffElePos->GetYaxis()->SetRangeUser(0.0,0.7);
    // mPtPairEffElePos->GetZaxis()->SetRangeUser(0,1.0);
    mPtPairEffElePos->Draw("Colz");
    legendInfoPairEff->Draw("same");
    cPairEffElePosMPt->SaveAs("./plots/PairEffElePosMPt.png");

  }

  auto cBeforeAfterSmearing = new TCanvas("cBeforeAfterSmearing","cBeforeAfterSmearing",1800,800);

  cBeforeAfterSmearing->Divide(3,1);
  for (size_t i = 1; i < 4; i++) {
    cBeforeAfterSmearing->cd(i)->SetTopMargin(0.03);
    cBeforeAfterSmearing->cd(i)->SetRightMargin(0.03);
  }
  cBeforeAfterSmearing->cd(1)->SetLogy();
  ptRecTrackBeforeSmearing->GetXaxis()->SetRangeUser(0.,4.);
  // ptRecTrackBeforeSmearing_rebin->GetXaxis()->SetRangeUser(0.,0.4);
  // ptRecTrackBeforeSmearing->SetMinimum(0.01);
  ptRecTrackBeforeSmearing_rebin->SetMaximum(1.1*std::max({ptRecTrackBeforeSmearing_rebin->GetMaximum(),ptRecTrackAfterSmearing_rebin->GetMaximum(),ptRecTrackElePos_rebin->GetMaximum()}));
  ptRecTrackBeforeSmearing_rebin->Draw("pe1");
  ptRecTrackAfterSmearing_rebin->Draw("pe1 same");
  ptRecTrackAfterKineCuts_rebin->Draw("pe1 same");

  ptRecTrackElePos_rebin->Draw("pe1 same");
  legSmearLabel->Draw("same");
  cBeforeAfterSmearing->cd(2);
  etaRecTrackBeforeSmearing->GetXaxis()->SetRangeUser(-4.5,4.5);
  etaRecTrackBeforeSmearing->SetMinimum(0.);
  etaRecTrackBeforeSmearing->SetMaximum(1.1*std::max({etaRecTrackBeforeSmearing->GetMaximum(),etaRecTrackAfterSmearing->GetMaximum(),etaRecTrackElePos->GetMaximum()}));
  etaRecTrackBeforeSmearing->Draw("axis");
  etaRecTrackBeforeSmearing->Draw("pe1 same");
  etaRecTrackAfterSmearing->Draw("pe1 same");
  etaRecTrackAfterKineCuts->Draw("pe1 same");

  etaRecTrackElePos->Draw("pe1 same");
  cBeforeAfterSmearing->cd(3);
  phiRecTrackBeforeSmearing->SetMinimum(0.);
  phiRecTrackBeforeSmearing->SetMaximum(1.1*std::max({phiRecTrackBeforeSmearing->GetMaximum(),phiRecTrackAfterSmearing->GetMaximum(),phiRecTrackElePos->GetMaximum()}));
  phiRecTrackBeforeSmearing->Draw("axis");
  phiRecTrackBeforeSmearing->Draw("pe1 same");
  phiRecTrackAfterSmearing->Draw("pe1 same");
  phiRecTrackAfterKineCuts->Draw("pe1 same");

  phiRecTrackElePos->Draw("pe1 same");
  cBeforeAfterSmearing->SaveAs("./plots/RecPtEtaPhiBeforeAfterSmearing.png");



  // if(bElePionRatios){
  //   auto cRatiobeforePID = new TCanvas("cRatiobeforePID","cRatiobeforePID",800,800);
  //   cRatiobeforePID->SetTopMargin(0.03);
  //   cRatiobeforePID->SetRightMargin(0.03);
  //   cRatiobeforePID->SetLeftMargin(0.13);
  //   cRatiobeforePID->SetLogx();
  //   cRatiobeforePID->SetLogy();
  //   // hbeforePIDRatio_beforePID->GetYaxis()->SetRangeUser(0.01,10.);
  //   hElePionRatio_beforePID->SetMaximum(10*std::max({hElePionRatio_beforePID->GetMaximum(),hMuonEleRatio_beforePID->GetMaximum()}));
  //   hElePionRatio_beforePID->SetMinimum(0.001+0.1*std::min({hElePionRatio_beforePID->GetMinimum(),hMuonEleRatio_beforePID->GetMinimum()}));
  //   hElePionRatio_beforePID->Draw("");
  //   hMuonEleRatio_beforePID->Draw("same");
  //   legRatiobeforePID->Draw("same");
  //   cRatiobeforePID->SetTicks();
  //   cRatiobeforePID->SaveAs("./plots/ratiobeforePID.png");
  //
  //   auto cRatioafterPID = new TCanvas("cRatioafterPID","cRatioafterPID",800,800);
  //   cRatioafterPID->SetTopMargin(0.03);
  //   cRatioafterPID->SetRightMargin(0.03);
  //   cRatioafterPID->SetLeftMargin(0.13);
  //   cRatioafterPID->SetLogx();
  //   cRatioafterPID->SetLogy();
  //   // hafterPIDRatio_beforePID->GetYaxis()->SetRangeUser(0.01,10.);
  //   hElePionRatio_afterPID->SetMaximum(10*std::max({hElePionRatio_afterPID->GetMaximum(),hMuonEleRatio_afterPID->GetMaximum()}));
  //   hElePionRatio_afterPID->SetMinimum(0.001+0.1*std::min({hElePionRatio_afterPID->GetMinimum(),hMuonEleRatio_afterPID->GetMinimum()}));
  //   hElePionRatio_afterPID->Draw("");
  //   hMuonEleRatio_afterPID->Draw("same");
  //   legRatioafterPID->Draw("same");
  //   cRatioafterPID->SetTicks();
  //   cRatioafterPID->SaveAs("./plots/ratioafterPID.png");
  // }



  // TCanvas *cChPions_beforeafterPSPID = new  TCanvas("cChPions_beforeafterPSPID","cChPions_beforeafterPSPID",670*1.3,700*1.3);
  //
  // // TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  //   pad1->SetBottomMargin(0);
  //   pad1->SetTopMargin(0.03);
  //   pad1->SetRightMargin(0.03);
  //   pad1->SetLogx(false);
  //   pad1->SetLogy();
  //   pad1->SetTicks();
  //   pad1->Draw();
  // cChPions_beforeafterPSPID->cd();
  // cChPions_beforeafterPSPID->cd(1);
  // // TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  //   pad2->SetTopMargin(0);
  //   // pad2->SetLogy();
  //   pad2->SetLogx(false);
  //   pad2->Draw();
  //   pad2->cd();
  //   pad2->SetTicks();
  //   // pad1->SetTopMargin(pad1->GetTopMargin()*0.5);
  //   // pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
  //   pad2->SetRightMargin(0.03);
  //   // pad2->SetBottomMargin(pad2->GetBottomMargin()*2.5);
  // pad1->cd();
  //
  // normalizeToBinWidth(ptRecPionTrack_beforePID_rebin);
  // normalizeToBinWidth(ptRecPionTrack_afterPSPID_rebin);
  // makeHistNice(ptRecPionTrack_beforePID_rebin, kBlack);
  // makeHistNice(ptRecPionTrack_afterPSPID_rebin, kRed+1);
  // ptRecPionTrack_beforePID_rebin->GetYaxis()->SetTitle("Counts");
  //
  // ptRecPionTrack_beforePID_rebin->Draw("hist p e");
  // ptRecPionTrack_afterPSPID_rebin->Draw("hist p e same");
  //
  // pad2->cd();
  // TH1F* hratioPlotSum;
  // hratioPlotSum = (TH1F*) ptRecPionTrack_afterPSPID_rebin->Clone();
  // hratioPlotSum->Divide(ptRecPionTrack_beforePID_rebin);
  //
  // TString ptXLabel = "#it{p}_{T} (GeV/#it{c})";
  // TString yRatioLabel = "after PID/before PID";
  // makeRatioNice(hratioPlotSum,ptXLabel,yRatioLabel);
  // hratioPlotSum->GetYaxis()->SetTitleOffset(0.5);
  // hratioPlotSum->Draw("hist p e");
  //
  // cChPions_beforeafterPSPID->SaveAs("./plots/ChPions_beforeafterPSPID.png");


  //
  // TCanvas *PtElePos_beforeafterPSPID = new  TCanvas("PtElePos_beforeafterPSPID","PtElePos_beforeafterPSPID",670*1.3,700*1.3);
  //
  // // TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  //   pad1->SetBottomMargin(0);
  //   pad1->SetTopMargin(0.03);
  //   pad1->SetRightMargin(0.03);
  //   pad1->SetLogx(false);
  //   pad1->SetLogy();
  //   pad1->SetTicks();
  //   pad1->Draw();
  // PtElePos_beforeafterPSPID->cd();
  // PtElePos_beforeafterPSPID->cd(1);
  // // TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  //   pad2->SetTopMargin(0);
  //   // pad2->SetLogy();
  //   pad2->SetLogx(false);
  //   pad2->Draw();
  //   pad2->cd();
  //   pad2->SetTicks();
  //   // pad1->SetTopMargin(pad1->GetTopMargin()*0.5);
  //   // pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
  //   pad2->SetRightMargin(0.03);
  //   // pad2->SetBottomMargin(pad2->GetBottomMargin()*2.5);
  // pad1->cd();
  //
  // normalizeToBinWidth(ptRecTrackElePos_beforePID_rebin);
  // // normalizeToBinWidth(ptRecTrackElePos_afterPSPID_rebin);
  // makeHistNice(ptRecTrackElePos_beforePID_rebin, kBlack);
  // // makeHistNice(ptRecTrackElePos_afterPSPID_rebin, kRed+1);
  // ptRecTrackElePos_beforePID_rebin->GetYaxis()->SetTitle("Counts");
  //
  // ptRecTrackElePos_beforePID_rebin->Draw("hist p e");
  // // ptRecTrackElePos_afterPSPID_rebin->Draw("hist p e same");
  //
  // pad2->cd();
  // hratioPlotSum;
  // // hratioPlotSum = (TH1F*) ptRecTrackElePos_afterPSPID_rebin->Clone();
  // hratioPlotSum->Divide(ptRecTrackElePos_beforePID_rebin);
  // TF1 *fPSeff = new TF1("effPS_inverse","0.8*(1.-exp(-1.6*(x-0.05)))",0,10);
  //
  // makeRatioNice(hratioPlotSum,ptXLabel,yRatioLabel);
  // hratioPlotSum->GetYaxis()->SetTitleOffset(0.5);
  // hratioPlotSum->Draw("hist p e");
  // fPSeff->Draw("same");
  //
  //
  // PtElePos_beforeafterPSPID->SaveAs("./plots/PtElePos_beforeafterPSPID.png");






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
    hTotalPureContaminationRecEta->Draw("same hist p e1");
    cContaminationPtEtaPhi->cd(3);
    hMuonContaminationRecPhi->GetYaxis()->SetTitle("contamination");
    hMuonContaminationRecPhi->GetXaxis()->SetTitle("#varphi (rad)");
    hMuonContaminationRecPhi->SetMaximum(0.1);
    hMuonContaminationRecPhi->SetMinimum(0.);
    hMuonContaminationRecPhi->GetXaxis()->SetRangeUser(-7.0,7.0);
    hMuonContaminationRecPhi->Draw("hist p e1");
    hPionContaminationRecPhi->Draw("same hist p e1");
    hTotalPureContaminationRecPhi->Draw("same hist p e1");
    cContaminationPtEtaPhi->SaveAs("./plots/PID_Contamiantion.png");



    auto cContaminationPt = new TCanvas("cContaminationPt","cContaminationPt",800,800);
    cContaminationPt->SetTopMargin(0.03);
    cContaminationPt->SetRightMargin(0.03);
    cContaminationPt->SetLeftMargin(0.13);
    hMuonContaminationRecPt->GetYaxis()->SetTitle("contamination");
    hMuonContaminationRecPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    // hMuonContaminationRecPt->GetXaxis()->SetRangeUser(0.0,0.4);
    hMuonContaminationRecPt->GetXaxis()->SetRangeUser(0.0,4.0);
    hMuonContaminationRecPt->SetMaximum(1.3);
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
    cContaminationPt->SaveAs("./plots/PID_Pt_Contamiantion.png");
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
  // hPionRejectionFactorPt->SetMaximum(10*std::max({hMuonRejectionFactorPt->GetMaximum(),hPionRejectionFactorPt->GetMaximum()/*,hTotalRejectionFactorPt->GetMaximum()*/}));
  // hPionRejectionFactorPt->SetMaximum(1.1);
  // hPionRejectionFactorPt->SetMinimum(0.1*std::min({hMuonRejectionFactorPt->GetMinimum(),hPionRejectionFactorPt->GetMinimum()/*,hTotalRejectionFactorPt->GetMaximum()*/}));
  // hPionRejectionFactorPt->SetMaximum(100000000000);
  hPionRejectionFactorPt->SetMaximum(100000000);
  hPionRejectionFactorPt->SetMinimum(0.1);
  hPionRejectionFactorPt->GetXaxis()->SetRangeUser(0.0,4.0);
  // hPionRejectionFactorPt->GetXaxis()->SetRangeUser(0.0,0.4);
  hMuonRejectionFactorPt->SetMinimum(0.1);
  hPionRejectionFactorPt->Draw("hist p e1");
  hMuonRejectionFactorPt->Draw("same hist p e1");
  // hKaonRejectionFactorPt->Draw("same hist p e1");
  // hProtonRejectionFactorPt->Draw("same hist p e1");
  // hTotalRejectionFactorPt->Draw("same hist p e1");
  legRejFactor->Draw("same");
  legendInfo->Draw("same");
  // textBField_conta->Draw("same");
  // textPIDScenario_conta->Draw("same");
  cRejectionFactorPt->SetTicks();
  cRejectionFactorPt->SaveAs("./plots/PID_Pt_RejectionFactor.png");
  // hPionRejectionFactorPt->SetMinimum(1000);
  // hPionRejectionFactorPt->SetMaximum(100000000);
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
      cNSigmaSeparatePID_TOF->SaveAs(Form("./plots/PID_histograms/NSigma%s_TOF_SeparatePiEleMuPID_noCuts.png",namesEleMuPi[iNSig]));

      auto cNSigmaSeparatePID_TOF_afterCuts = new TCanvas("cNSigmaPionSeparatePID_TOF_afterCuts","cNSigmaPionSeparatePID_TOF_afterCuts",800,800);
      cNSigmaSeparatePID_TOF_afterCuts->SetTopMargin(0.03);
      cNSigmaSeparatePID_TOF_afterCuts->SetRightMargin(0.03);
      cNSigmaSeparatePID_TOF_afterCuts->SetLeftMargin(0.13);
      cNSigmaSeparatePID_TOF_afterCuts->SetLogx();
      hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->GetXaxis()->SetRangeUser(0.0,4.0);
      // hNsigmaP_afterPIDcuts_TOF_trueElec[iNSig]->GetXaxis()->SetRangeUser(0.0,0.4);
      hNsigmaP_afterPIDcuts_TOF_trueElec[iNSig]->GetYaxis()->SetRangeUser(-15.,25.);
      // hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->SetMaximum(1.1);
      // hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->SetMinimum(0.);
      hNsigmaP_afterPIDcuts_TOF_trueElec[iNSig]->Draw("");
      hNsigmaP_afterPIDcuts_TOF_truePion[iNSig]->Draw("same");
      hNsigmaP_afterPIDcuts_TOF_trueMuon[iNSig]->Draw("same");
      legPIDSeparateColorTOF->Draw("same");
      legendInfoTOF->Draw("same");
      cNSigmaSeparatePID_TOF_afterCuts->SetTicks();
      cNSigmaSeparatePID_TOF_afterCuts->SaveAs(Form("./plots/PID_histograms/NSigma%s_TOF_SeparatePiEleMuPID_afterCuts.png",namesEleMuPi[iNSig]));



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
      cNSigmaSeparatePID_RICH->SaveAs(Form("./plots/PID_histograms/NSigma%s_RICH_SeparatePiEleMuPID_noCuts.png",namesEleMuPi[iNSig]));


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
      cNSigmaSeparatePID_RICH_afterCuts->SaveAs(Form("./plots/PID_histograms/NSigma%s_RICH_SeparatePiEleMuPID_afterCuts.png",namesEleMuPi[iNSig]));



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

      auto cprojMeePtee = new TCanvas("cprojMeePtee","cprojMeePtee",1800,800);
      cprojMeePtee->Divide(2,1);
      for (size_t i = 1; i < 4; i++) {
        cprojMeePtee->cd(i)->SetTopMargin(0.03);
        cprojMeePtee->cd(i)->SetRightMargin(0.03);
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
      cprojMeePtee->cd(1);
      cprojMeePtee->cd(1)->SetLogy();
      proj_recULS_Mee->GetYaxis()->SetTitle("Yield");
      proj_recULS_Mee->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
      // proj_recULS_Mee->RebinX(2);
      // proj_recLS_Mee->RebinX(2);
      proj_recULS_Mee->GetXaxis()->SetRangeUser(0.0,3.0);
      proj_recULS_Mee->Draw("hist p e1");
      proj_recLS_Mee->Draw("same hist p e1");
      legendInfo->Draw("same");
      cprojMeePtee->cd(2);
      cprojMeePtee->cd(2)->SetLogy();
      proj_recULS_Ptee->GetYaxis()->SetTitle("Yield");
      proj_recULS_Ptee->GetXaxis()->SetTitle("p_{T,ee} (GeV/#it{c})");
      // proj_recULS_Ptee->RebinX(2);
      proj_recULS_Ptee->GetXaxis()->SetRangeUser(0.0,4.0);
      proj_recULS_Ptee->Draw("hist p e1");
      proj_recLS_Ptee->Draw("same hist p e1");
      legendPtee->Draw("same");

      cprojMeePtee->SaveAs(Form("./plots/projMeePteeee_sce%i.png",ith_PIDscenario));

      // auto cMeePteeTH2 = new TCanvas("cMeePteeTH2","cMeePteeTH2",800,800);
      // cMeePteeTH2->SetTopMargin(0.03);
      // cMeePteeTH2->SetRightMargin(0.03);
      // cMeePteeTH2->SetLeftMargin(0.13);
      // hMPtDCA_ULS_rec->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      // hMPtDCA_ULS_rec->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
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
      // proj_recULS_Mee->GetXaxis()->SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
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


  // Generate further root files with new plots

  TString nameEffRootFile = inputFile.Data();
  // TString endDirName = "12.9.21_ptcut0.2_0.8eta_6PFstudies_woPF_oTOFRICH_iTOFpfacc_diffOpAngleMassCuts";
  TString endDirName = "18.9.21_ptcut0.2_0.8eta_prefilter_update_0-10cent";
  // TString endDirName = "19.9.21_ptcut0.2_1.75eta_prefilter_update_30-50cent";
  if (inputFile.Contains("./data")) nameEffRootFile.ReplaceAll("./data/prod/anaEEstudy.", "");
  if (inputFile.Contains("../grid/second_job")) nameEffRootFile.ReplaceAll("../grid/second_job/data/prod/anaEEstudy.", "");
  // if (inputFile.Contains("B2_100k")) nameEffRootFile.ReplaceAll(Form("../grid/output/B2_100k_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
  // if (inputFile.Contains("B5_100k")) nameEffRootFile.ReplaceAll(Form("../grid/output/B5_100k_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
  // if (inputFile.Contains("B2_2M_502TeV")) nameEffRootFile.ReplaceAll(Form("../grid/output/B2_2M_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
  // if (inputFile.Contains("B5_2M_502TeV")) nameEffRootFile.ReplaceAll(Form("../grid/output/B5_2M_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
  // if (inputFile.Contains("B=0.5_2M_502TeV")) nameEffRootFile.ReplaceAll("/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/efficiency/data/prod/anaEEstudy.", "");
  if (inputFile.Contains("B5_100k_502TeV")) nameEffRootFile.ReplaceAll(Form("/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/grid/output/B5_100k_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
  if (inputFile.Contains("B5_200k_502TeV")) nameEffRootFile.ReplaceAll(Form("/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/grid/output/B5_200k_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
  // if (inputFile.Contains("B20_100k_502TeV")) nameEffRootFile.ReplaceAll(Form("/data/feisenhut/DelphesO2/ALICE3-LoI-LMee/grid/output/B20_100k_502TeV_%s/anaEEstudy.",endDirName.Data()), "");
  nameEffRootFile.ReplaceAll(".root", "");
  TFile *fOut = TFile::Open(Form("./data/efficiencies/TrackEff_%s_PIDscenario%i.root",nameEffRootFile.Data(),ith_PIDscenario),"RECREATE");
  //
  if (bPlotEfficiency) {
    ptEffElePosGen->SetTitle("eff_Track_Rec/Generated");
    ptEffElePosGenSmeared->SetTitle("eff_Track_Rec/GeneratedSmeared");
    ptEffElePosGen->GetXaxis()->SetRangeUser(0.0,10.0);
    ptEffElePosGenSmeared->GetXaxis()->SetRangeUser(0.0,10.0);
    ptEffElePosGen->Write();
    ptEffElePosGenSmeared->Write();
    ptEffElePosGenSmeared->SetTitle("eff_HF_Track_Rec/Generated");
    ptEffElePos_LF_afterPF->Write();
    ptEffElePos_HF_afterPF->Write();

    ptEtaEffElePos->SetTitle("eff_pteta_2D");
    ptEtaEffElePos->GetXaxis()->SetRangeUser(0.0,4.0);
    ptEtaEffElePos->GetYaxis()->SetRangeUser(-1.2,1.2);
    ptEtaEffElePos->Write();

  }

  if (bPlotTrackContamination) {
    hTotalPureContaminationRecPt->SetTitle("Total Contamination");
    hTotalPureContaminationRecPt->GetXaxis()->SetRangeUser(0.0,4.0);
    hTotalPureContaminationRecPt->Write();
  }
  fOut->Close();

}
}
