TH1D* convertCocktailHistogram(TH1D *cocktailHisto, Bool_t doPtee);
TH1D* convertRappHistogram(TH1D *rappHisto, Bool_t doPtee);
TH1D* convertCharmHistogram(TH1D *charmHisto, Bool_t doPtee);
void plotFinalPlotsLowB_Systems (
				 Int_t iSystem = 0,
				 Bool_t doPlot = kTRUE,
				 Bool_t Run5   = kFALSE,
				 Bool_t RICH   = kFALSE,
				 Bool_t doPtee = kTRUE,
				 Int_t meebin = 0,
				 Double_t nEventsLHC = 2.500000000000E9,
				 Bool_t doIPcut = kFALSE,
				 Bool_t writeFile = kTRUE,
				 Bool_t doPlain  = 0, // no Poisson Sample & Bkg Error
				 Bool_t doNoSyst = 1, // no Systematic errors
				 Bool_t doEff = 1,    // efficiency implementation
				 Double_t conversionBGScaling = 1.0 // scaling from BG w/o conversions 
				 ){


  // ----------------------------------------------------------------------
  // STYLE
  // ----------------------------------------------------------------------
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);// no box around histogram title
  gStyle->SetOptStat(0);// statistics yes, less, no (1, 11, 0)
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);




  
  // ----------------------------------------------------------------------
  // CONFIGURATION
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // Some numbers 
  Double_t dNchdy_55    = 1930.; // to scale Rapp histos to 5.5TeV
  Double_t dNchdeta_55  = 1750.; // to scale Averbeck histos to 5.5TeV
  Double_t dNchdeta_276 = 1207.; // to downscale Averbeck histos (0-20%) first (estimate from PR, ALICE paper: cent. dep. of mult.)
  Double_t Ncoll        = 1625.; // to scale charm (estimate from Harry: Ncoll = sigma_pp * T_AA = 69.2 mb * 23.48 mb^-1 = 1625
  Double_t NunitsEta    = 1.7; // to scale cocktail histos
  Double_t NunitsEtaExp = 1.6; // to scale ALICE histos 
  Double_t MatchRhos    = 1.35; // to scale Rapp histos (first version: 2.8)
  // ----------------------------------------------------------------------
  
  // ----------------------------------------------------------------------
  // Collision system set up
  const Int_t nSystems = 4;
  TString sSystem[nSystems] = {"Pb","Xe","Kr","Ar"};
  TString sSystem2[nSystems] = {"Pb-Pb","Xe-Xe","Kr-Kr","Ar-Ar"};
  Double_t A[nSystems] = {208.,129.,78.,40.};
  Double_t lumiFactor[nSystems]     = {1.,7.8,44.,646.};
  Double_t sigmaHadFactor[nSystems] = {7.8,5.67,4.06,2.6};

   if(iSystem >= nSystems){
    Printf("This system is not known...");
    return;
  }
  else{
    Printf("Process system %d = %s (%s)",iSystem,sSystem[iSystem].Data(),sSystem2[iSystem].Data());
  }
  
  if(doIPcut){
    Printf("Not implemented yet");
    return;
  }

  Printf("Scale number of events for Pb-Pb with %f",lumiFactor[iSystem]*sigmaHadFactor[iSystem]/sigmaHadFactor[0]);
  nEventsLHC *= lumiFactor[iSystem]*sigmaHadFactor[iSystem]/sigmaHadFactor[0];
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // prompt and non-prompt efficiencies
  Double_t PrimEff  = 1.;
  Double_t CharmEff = 1.;
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // Titles and names
  char MeeOrPtee[10]  = "Mee";
  char titleXaxis[30] = "#it{M}_{ee} (GeV/#it{c }^{2})";
  char titleYaxis[80] = "1/#it{N}_{evt} d#it{N }^{2}/d#it{M}_{ee}d#it{y} (GeV/#it{c }^{2})^{-1}";
  if (doPtee) {
    sprintf(MeeOrPtee,  "Ptee");
    sprintf(titleXaxis, "#it{p}_{T,ee} (GeV/#it{c})");
    sprintf(titleYaxis, "d#it{N}^{2}/d#it{p}_{T,ee}d#it{y} (GeV/#it{c})^{-1}");
  }
  Int_t    nLogEventsLHC = (Int_t) TMath::Log10(nEventsLHC); // (Int_t) not needed but ok here
  Int_t    nIntEventsLHC = nEventsLHC / TMath::Power(10, nLogEventsLHC); // do not use (Int_t) here..!
  Double_t nDouEventsLHC = nEventsLHC / TMath::Power(10, nLogEventsLHC);

  const Int_t nMeeRanges         = 6;
  Double_t MeeRanges[nMeeRanges] = {0.15, 0.3, 0.7, 1.1, 1.5, 3.0};
  Double_t MeeRangeMin           = 0.15;
  Double_t MeeRangeMax           = 3.0;
  
  char projXorY[10] = "px";
  if (doPtee) {
    if(meebin >= nMeeRanges){
      Printf("You chose the wrong pt range");
      return;
    }
    if(meebin == 0){
      MeeRangeMin = MeeRanges[0];
      MeeRangeMax = MeeRanges[nMeeRanges-1];
    }
    else{
      MeeRangeMin = MeeRanges[meebin-1];
      MeeRangeMax = MeeRanges[meebin];
    }
    sprintf(projXorY, "py");
  }
  else { meebin = 0; } // Mee is currently done with just one integrated Pt bin.

  if(doPtee)
    Printf("We are doing projXorY=%s with meebin=%d (%.2f<M<%.2f) ",projXorY,meebin,MeeRangeMin,MeeRangeMax);
  else
    Printf("We are doing projXorY=%s with meebin=%d",projXorY,meebin);
  
  char ITSup_opt[10];
  char canvname_opt[64];
  char reichelt_opt[64];
  char filename_reichelt[128];
  Double_t ScaleReicheltHistos_becauseIPcut = 1.; // only for Current ITS, see below!
  Double_t xmaxMee  = 2.5;
  Double_t ymindNdx = 1E-4 *PrimEff;
  Double_t ymaxdNdx = 1000. *PrimEff;// from Npp360
  if (doPtee) {
    xmaxMee         = 1.5; // gives 3...
    ymindNdx        = 1E-4 *PrimEff;
    ymaxdNdx        = 50. *PrimEff;// from Npp360
  }
  cout<<"We plot from "<<ymindNdx<<" "<<ymaxdNdx<<endl;

  Double_t yminSoverB = 1E-3;	// from Npp360
  Double_t ymaxSoverB = 1E+2;  // from Npp1
  Double_t yminSignif = 1E-4;   // from Npp1
  Double_t ymaxSignif = 1E+1;	// from Npp360
  // ----------------------------------------------------------------------




  
  // ----------------------------------------------------------------------
  // GET INPUTS AND PRODUCE RESULTS
  // ---------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // Input files with all needed hiostograms (from plotAnapartTreeWithCocktail)
  TCanvas *cInputQA = new TCanvas("cInputQA","cInputQA",1200,500);
  cInputQA->Divide(3,1);

  TString inputfilename = "";
  
  if(doIPcut){ 
    Printf("Not implemented yet");
    return;
  }
  else{//no additional DCA cut
    if(Run5){
      if(!doPtee){
	if(RICH)
	  inputfilename = Form("inputFiles/AnaPartTreeOut_RecRun5Pt200_RICH_Results.root");
	else
	  inputfilename = Form("inputFiles/AnaPartTreeOut_RecRun5Pt200_Results.root");
      }
      else{
	if(RICH)
	  inputfilename = Form("inputFiles/AnaPartTreeOut_RecRun5Pt200_RICH_pT_Results.root");
	else
	  inputfilename = Form("inputFiles/AnaPartTreeOut_RecRun5Pt200_pT_Results.root");
      }
    }
    else{
      if(!doPtee)
	inputfilename = Form("inputFiles/AnaPartTreeOut_RecRun3Pt200_Results.root");
      else{
	inputfilename = Form("inputFiles/AnaPartTreeOut_RecRun3Pt200_pT_Results.root");
      }
    }
  }
  TFile *fIn = TFile::Open(inputfilename,"READ");
  TH1D* hSignalwoRho   = NULL;
  TH1D* hSignalAll     = NULL;
  TH1D* hSignal        = NULL;
  TH1D* hSBQA          = NULL;
  TH1D* hBG            = NULL;

  // just get one of the output histograms for the proper settings
  if(!doPtee){
    hSignalwoRho   = (TH1D*)fIn->Get("fHistBGnC")->Clone("hSignalwoRho");
    hSignalAll     = (TH1D*)fIn->Get("fHistInvMassPt_BG_px")->Clone("hSignalAll");
    hSignal        = (TH1D*)fIn->Get("fHistBGnC")->Clone("hSignal");
    hSBQA          = (TH1D*)fIn->Get("fHistBGnC")->Clone("hSBQA");
  }
  else{
    hSignalwoRho   = (TH1D*)fIn->Get("fHistInvMassPt_BG_py")->Clone("hSignalwoRho");
    hSignalAll     = (TH1D*)fIn->Get("fHistInvMassPt_BG_py")->Clone("hSignalAll");
    hSignal        = (TH1D*)fIn->Get("fHistInvMassPt_BG_py")->Clone("hSignal");
    hSBQA          = (TH1D*)fIn->Get("fHistInvMassPt_BG_py")->Clone("hSBQA");
  }
  
  // we use ideal BG (no conversions) here, since
  // a) in "Eff" generated files there are no conversions and
  // b) we do not know the conversion quite well for a Run5 setup
  hBG = (TH1D*)fIn->Get("fHistBGnC")->Clone("hBG");// ideal BG (no conversions)
  hBG->Scale(1./(hBG->GetBinWidth(1)));// not done in plotAnaPartTree.C

  // scaling for conversions in BG (to be determined externally so far), could be collision system etc. dependent
  Printf("Scale BG input (without conversions) with %f to match expectation with conversions",conversionBGScaling);
  hBG->Scale(conversionBGScaling);
	     
  // scaling for collision system
  Double_t scaleBG = TMath::Power(A[iSystem]/A[0],2.); 
  Printf("Scale BG input for Pb-Pb with %f",scaleBG);
  hBG->Scale(scaleBG);
  
  // just get one of clone the BG histogram for the proper histogram settings (content reset later)
  hSignalwoRho   = (TH1D*)hBG->Clone("hSignalwoRho");
  hSignalAll     = (TH1D*)hBG->Clone("hSignalAll");
  hSignal        = (TH1D*)hBG->Clone("hSignal");
  hSBQA          = (TH1D*)hBG->Clone("hSBQA");
  
  // cocktail input
  TString filenameCocktail = "~/MAC_201909/ALICE//dileptons/Upgrade/ITS/analysis_Patrick/DileptonUpgradeStudiesPatrick/Input_Averbeck/cocktail_10M_pbpb_0020_upgrade_2760GeV_v1.root";
  TFile *fCocktail = TFile::Open(filenameCocktail,"READ");
  
  Double_t scaleCocktail = TMath::Power(A[iSystem]/A[0],1.); 
  Printf("Scale cocktail input for Pb-Pb with %f",scaleCocktail);

  Int_t binMaxPt  = 0;
  Int_t binMinMee = 0;
  Int_t binMaxMee = 0;
  Double_t binWidthC = 0.;

  TH2D* hPion2D     = (TH2D*)fCocktail->Get("pteevsmeePion");
  TH2D* hEta2D      = (TH2D*)fCocktail->Get("pteevsmeeEta");
  TH2D* hEtaPrime2D = (TH2D*)fCocktail->Get("pteevsmeeEtaprime");
  TH2D* hOmega2D    = (TH2D*)fCocktail->Get("pteevsmeeOmega");
  TH2D* hRho2D      = (TH2D*)fCocktail->Get("pteevsmeeRho");
  TH2D* hPhi2D      = (TH2D*)fCocktail->Get("pteevsmeePhi");

  TCanvas *cQA = new TCanvas("cQA","cQA",1200,900);
  cQA->Divide(3,2);
  cQA->cd(1)->SetLogz();
  hPion2D->Draw("colz");
  cQA->cd(2)->SetLogz();
  hEta2D->Draw("colz");
  cQA->cd(3)->SetLogz();
  hEtaPrime2D->Draw("colz");
  cQA->cd(4)->SetLogz();
  hOmega2D->Draw("colz");
  cQA->cd(5)->SetLogz();
  hRho2D->Draw("colz");
  cQA->cd(6)->SetLogz();
  hPhi2D->Draw("colz");
    
  TH1D* hPion       = NULL; 
  TH1D* hEta        = NULL; 
  TH1D* hEtaPrime   = NULL;
  TH1D* hOmega      = NULL; 
  TH1D* hRho        = NULL; 
  TH1D* hPhi        = NULL; 

  if(!doPtee){
    binWidthC   = hPion2D->GetYaxis()->GetBinWidth(1);
    binMaxPt    = hPion2D->GetYaxis()->FindBin(3.0-0.00001);  
    hPion       = (TH1D*)convertCocktailHistogram((TH1D*)hPion2D->ProjectionX("hPion",1,binMaxPt),doPtee);
    hEta        = (TH1D*)convertCocktailHistogram((TH1D*)hEta2D->ProjectionX("hEta",1,binMaxPt),doPtee);
    hEtaPrime   = (TH1D*)convertCocktailHistogram((TH1D*)hEtaPrime2D->ProjectionX("hEtaPrime",1,binMaxPt),doPtee);
    hOmega      = (TH1D*)convertCocktailHistogram((TH1D*)hOmega2D->ProjectionX("hOmega",1,binMaxPt),doPtee);
    hRho        = (TH1D*)convertCocktailHistogram((TH1D*)hRho2D->ProjectionX("hRho",1,binMaxPt),doPtee);
    hPhi        = (TH1D*)convertCocktailHistogram((TH1D*)hPhi2D->ProjectionX("hPhi",1,binMaxPt),doPtee);
  }
  else{
    binWidthC    = hPion2D->GetXaxis()->GetBinWidth(1);
    binMinMee    = hPion2D->GetXaxis()->FindBin(MeeRangeMin+0.00001);
    binMaxMee    = hPion2D->GetXaxis()->FindBin(MeeRangeMax-0.00001);
    cout<<binMinMee<<" "<<binMaxMee<<endl;
    hPion       = (TH1D*)convertCocktailHistogram((TH1D*)hPion2D->ProjectionY("hPion",binMinMee,binMaxMee),doPtee);
    hEta        = (TH1D*)convertCocktailHistogram((TH1D*)hEta2D->ProjectionY("hEta",binMinMee,binMaxMee),doPtee);
    hEtaPrime   = (TH1D*) convertCocktailHistogram((TH1D*)hEtaPrime2D->ProjectionY("hEtaPrime",binMinMee,binMaxMee),doPtee);
    hOmega      = (TH1D*) convertCocktailHistogram((TH1D*)hOmega2D->ProjectionY("hOmega",binMinMee,binMaxMee),doPtee);
    hRho        = (TH1D*) convertCocktailHistogram((TH1D*)hRho2D->ProjectionY("hRho",binMinMee,binMaxMee),doPtee);
    hPhi        = (TH1D*) convertCocktailHistogram((TH1D*)hPhi2D->ProjectionY("hPhi",binMinMee,binMaxMee),doPtee);
  }

  Printf("Also scale with bin width of 2D cocktail histo (done previously before) = %f",binWidthC);
  hPion->Scale(scaleCocktail*binWidthC);
  hEta->Scale(scaleCocktail*binWidthC);
  hEtaPrime->Scale(scaleCocktail*binWidthC);
  hOmega->Scale(scaleCocktail*binWidthC);
  hRho->Scale(scaleCocktail*binWidthC);
  hPhi->Scale(scaleCocktail*binWidthC);

  // Rapp input
  TString filenameRapp = "~/MAC_201909/ALICE//dileptons/Upgrade/ITS/theory/Rapp/hist2D_Rapp_v4_1930_update_multqt.root";
  TFile *fRapp = TFile::Open(filenameRapp,"READ");

  TH2D* fHistRhoVac2D = (TH2D*)fRapp->Get("fH2D_vacSF");
  TH2D* fHistRhoMed2D = (TH2D*)fRapp->Get("fH2D_medSF");
  TH2D* fHistRhoDrop2D = (TH2D*)fRapp->Get("fH2D_dropSF");
  TH2D* fHistQGP2D = (TH2D*)fRapp->Get("fH2D_QGP");

  TH1D* hRhoVacRR  = NULL; 
  TH1D* hRhoDropRR = NULL; 
  TH1D* hRhoRR     = NULL; 
  TH1D* hQGP       = NULL; 

  if(!doPtee){
    binMaxPt = fHistRhoMed2D->GetYaxis()->FindBin(3.0-0.00001);
    hRhoVacRR  = convertRappHistogram((TH1D*)fHistRhoVac2D->ProjectionX("fHistRhoVacRR",1,binMaxPt),doPtee);
    hRhoDropRR = convertRappHistogram((TH1D*)fHistRhoDrop2D->ProjectionX("fHistRhoDropRR",1,binMaxPt),doPtee);
    hRhoRR     = convertRappHistogram((TH1D*)fHistRhoMed2D->ProjectionX("fHistRhoRR",1,binMaxPt),doPtee);
    hQGP       = convertRappHistogram((TH1D*)fHistQGP2D->ProjectionX("fHistQGP",1,binMaxPt),doPtee);
  }
  else{
    binMinMee    = hPion2D->GetYaxis()->FindBin(MeeRangeMin+0.00001);
    binMaxMee    = hPion2D->GetYaxis()->FindBin(MeeRangeMax-0.00001);
    hRhoVacRR  = convertRappHistogram((TH1D*)fHistRhoVac2D->ProjectionY("fHistRhoVacRR",binMinMee,binMaxMee),doPtee);
    hRhoDropRR = convertRappHistogram((TH1D*)fHistRhoDrop2D->ProjectionY("fHistRhoDropRR",binMinMee,binMaxMee),doPtee);
    hRhoRR     = convertRappHistogram((TH1D*)fHistRhoMed2D->ProjectionY("fHistRhoRR",binMinMee,binMaxMee),doPtee);
    hQGP       = convertRappHistogram((TH1D*)fHistQGP2D->ProjectionY("fHistQGP",binMinMee,binMaxMee),doPtee);
  }
  
  hRhoVacRR->Scale(fHistRhoVac2D->GetYaxis()->GetBinWidth(1));//this is also done in Patrick's projection macro
  hRhoDropRR->Scale(fHistRhoDrop2D->GetYaxis()->GetBinWidth(1));//this is also done in Patrick's projection macro
  hRhoRR->Scale(fHistRhoMed2D->GetYaxis()->GetBinWidth(1));//this is also done in Patrick's projection macro
  hQGP->Scale(fHistQGP2D->GetYaxis()->GetBinWidth(1));//this is also done in Patrick's projection macro

  Double_t scaleRapp = TMath::Power(A[iSystem]/A[0],1.4); 
  Printf("Scale Rapp input for Pb-Pb with %f",scaleRapp);

  hRhoVacRR->Scale(scaleRapp);
  hRhoDropRR->Scale(scaleRapp);
  hRhoRR->Scale(scaleRapp);
  hQGP->Scale(scaleRapp);

  // Albericas input:
  //	hpt_mass2;1     mass vs pt,|y|<0.84,charm
  //	hpt_mass3;1     mass vs pt,|y|<0.84,charm, effi weight
  //	hpt_mass4;1     mass vs pt,|y|<0.84,charm, effi and pt weight
  //	hpt_mass5;1     mass vs pt,|y|<0.84,charm, pt > 0.4 GeV/c
  //	hpt_mass6;1     mass vs pt,|y|<0.84,charm, pt > 0.2 GeV/c
  TString filenameCharm = "inputFiles/charm_diele_v3_hpt_mass6_proj.root";
  TFile* fCharm = new TFile(filenameCharm, "READ");
  TH1D* hCharm    = convertCharmHistogram((TH1D*) fCharm->Get(Form("meeCharm_%s%d", projXorY, meebin)),doPtee);
  hCharm->Smooth(1000); // remove stat fluctuations

  Double_t scaleHF = TMath::Power(A[iSystem]/A[0],4./3.); 
  Printf("Scale HF input for Pb-Pb with %f",scaleHF);

  hCharm->Scale(scaleHF);

  // pair efficiency
  TString filenameEff = "";
  if(!Run5){
    if(doPtee){
      if(meebin==0)
	filenameEff = "inputFiles/pairEfficiency_pT_Study7.root";
      else
	filenameEff = Form("inputFiles/pairEfficiency_pT_meebin%d_Study7.root",meebin);
    }
    else
      filenameEff = "inputFiles/pairEfficiency_Study7.root";
  }
  else{
    if(RICH){
      if(doPtee){
	if(meebin==0)
	  filenameEff = "inputFiles/pairEfficiency_pT_Study9.root";
	else
	  filenameEff = Form("inputFiles/pairEfficiency_pT_meebin%d_Study9.root",meebin);
      }
      else
	filenameEff = "inputFiles/pairEfficiency_Study9.root";
    }
    else{
      if(doPtee){
	if(meebin==0)
	  filenameEff = "inputFiles/pairEfficiency_pT_Study8.root";
	else
	  filenameEff = Form("inputFiles/pairEfficiency_pT_meebin%d_Study8.root",meebin);
      }
      else
	filenameEff = "inputFiles/pairEfficiency_Study8.root";
    }
  }

  Printf("We take this efficiency file: %s",filenameEff.Data());
  
  TFile *fEff = TFile::Open(filenameEff,"READ");
  TH1D* meePairEff  = NULL;
  //meePairEff  = (TH1D*) fEff->Get(Form("fHistEffFit"));//gives more stable results, but problems if values go to ZERO (Run 5 w/o RICH)
  meePairEff  = (TH1D*) fEff->Get(Form("fHistEffPrimNum"));
 
  if(meePairEff){
    TCanvas* cEff = new TCanvas("cEff","cEff",600,500);
    cEff->cd();
    meePairEff->Draw();
  }
  else{
    Printf("Efficiency file not found...");
    return;
  }
 
  // scalings
  hPion     ->Scale(dNchdeta_55 / dNchdeta_276 / NunitsEta * PrimEff);
  hEta      ->Scale(dNchdeta_55 / dNchdeta_276 / NunitsEta * PrimEff);
  hEtaPrime ->Scale(dNchdeta_55 / dNchdeta_276 / NunitsEta * PrimEff);
  hRho      ->Scale(dNchdeta_55 / dNchdeta_276 / NunitsEta * PrimEff);
  hOmega    ->Scale(dNchdeta_55 / dNchdeta_276 / NunitsEta * PrimEff);
  hPhi      ->Scale(dNchdeta_55 / dNchdeta_276 / NunitsEta * PrimEff);
  hQGP      ->Scale(dNchdy_55 * MatchRhos * PrimEff);
  hRhoVacRR ->Scale(dNchdy_55 * MatchRhos * PrimEff);
  hRhoDropRR ->Scale(dNchdy_55 * MatchRhos * PrimEff);
  hRhoRR    ->Scale(dNchdy_55 * MatchRhos * PrimEff);
  hCharm    ->Scale(Ncoll / NunitsEta * CharmEff);

  hBG->Scale(1. / NunitsEtaExp);
  //hCharm    ->Scale();
  //hBottom   ->Scale();

  if(doEff){
    
    Double_t bincont  = 1.;
    Double_t paireffi = 1.;
    Int_t binPairEff  = 1;
    for (Int_t ix=1; ix<=hRhoRR->GetNbinsX(); ix++) {
      
      binPairEff = meePairEff->FindBin(hRhoRR->GetBinCenter(ix));
      paireffi   = meePairEff->GetBinContent(binPairEff);
      
      bincont  = hPion->GetBinContent(ix) * paireffi;
      hPion->SetBinContent(ix, bincont);
      
      bincont  = hEta->GetBinContent(ix) * paireffi;
      hEta->SetBinContent(ix, bincont);
      
      bincont  = hEtaPrime->GetBinContent(ix) * paireffi;
      hEtaPrime->SetBinContent(ix, bincont);
      
      bincont  = hRho->GetBinContent(ix) * paireffi;
      hRho->SetBinContent(ix, bincont);
      
      bincont  = hOmega->GetBinContent(ix) * paireffi;
      hOmega->SetBinContent(ix, bincont);
      
      bincont  = hPhi->GetBinContent(ix) * paireffi;
      hPhi->SetBinContent(ix, bincont);
      
      bincont  = hRhoVacRR->GetBinContent(ix) * paireffi;
      hRhoVacRR->SetBinContent(ix, bincont);
    
      bincont  = hRhoDropRR->GetBinContent(ix) * paireffi;
      hRhoDropRR->SetBinContent(ix, bincont);
    
      bincont  = hRhoRR->GetBinContent(ix) * paireffi;
      hRhoRR->SetBinContent(ix, bincont);
      
      bincont  = hQGP->GetBinContent(ix) * paireffi;
      hQGP->SetBinContent(ix, bincont);
      
      bincont  = hCharm->GetBinContent(ix) * paireffi;
      hCharm->SetBinContent(ix, bincont);
      
      // bincont  = hBottom->GetBinContent(ix) * paireffi;
      // hBottom->SetBinContent(ix, bincont);
    }
  }


  // build total signal
  hSignalwoRho->Reset();
  hSignalwoRho->Add(hPion);
  hSignalwoRho->Add(hEta);
  hSignalwoRho->Add(hEtaPrime);
  hSignalwoRho->Add(hOmega);
  hSignalwoRho->Add(hPhi);

  hSignalAll->Reset();
  hSignalAll->Add(hSignalwoRho);
  hSignalAll->Add(hRhoRR);
  hSignalAll->Add(hQGP);
  hSignalAll->Add(hCharm);
  //hSignalAll->Add(hBottom);

  hSignal->Reset();
  hSignal->Add(hSignalwoRho);
  hSignal->Add(hRhoRR);
  hSignal->Add(hQGP);
  hSignal->Add(hCharm);
  //hSignal->Add(hBottom);
 
  cInputQA->cd(1)->SetLogy();
  hSignalAll->GetXaxis()->SetRangeUser(0.,xmaxMee);
  hSignalAll->SetFillColor(0);
  hSignalAll->SetLineWidth(2);
  hSignalAll->SetLineColor(1);
  hSignalAll->Draw("P");
  hBG->SetLineColor(1);
  hBG->SetMarkerColor(1);
  hBG->Draw("same");


  cInputQA->cd(2)->SetLogy();
  hRhoVacRR->GetXaxis()->SetRangeUser(0,xmaxMee);
  hRhoVacRR->GetYaxis()->SetRangeUser(1e-6,1);
  hRhoVacRR->SetFillColor(0);
  hRhoVacRR->SetLineWidth(2);
  hRhoVacRR->SetLineColor(1);
  hRhoVacRR->Draw("L");
  hRho->SetFillColor(0);
  hRho->SetLineWidth(2);
  hRho->SetLineColor(4);
  hRho->Draw("Lsame");
 

  // S/B and significance
  TH1D* hSB           = (TH1D*)fIn->Get("fHistSBT");
  hSBQA->Reset();
  hSBQA->Add(hSignalwoRho);
  hSBQA->Add(hRhoRR);
  hSBQA->Add(hQGP);
  hSBQA->Add(hCharm);
  //hSBQA->Add(hBottom);
  hSBQA->Divide(hBG);

  TH1D* hSignificance = (TH1D*)hSignalAll->Clone("hSignificance");
  TH1D* hSigDen = (TH1D*)hSignalAll->Clone("hSigDen");
  hSigDen->Add(hBG,2.); // significance = S/sqr(S+B) or S+2B ? -->  use LS for EoI and YR
  for(Int_t i = 0; i < hSigDen->GetNbinsX(); i++){
    hSigDen->SetBinContent(i+1,TMath::Sqrt(hSigDen->GetBinContent(i+1)));
    hSigDen->SetBinError(i+1,TMath::Sqrt(hSigDen->GetBinError(i+1))); // FIXME: that's wrong!!!
  }
  hSignificance->Divide(hSigDen);
  for(Int_t i = 0; i < hSigDen->GetNbinsX(); i++){
    hSB->SetBinError(i+1,0.00000001); // FIXME: that's wrong!!!
    hSBQA->SetBinError(i+1,0.000000001); // FIXME: that's wrong!!!
    hSignificance->SetBinError(i+1,0.00000001); // FIXME: that's wrong!!!
  }

  cInputQA->cd(3)->SetLogy();
  hSBQA->GetXaxis()->SetRangeUser(0.,xmaxMee);
  hSBQA->SetMarkerColor(4);
  hSBQA->SetLineColor(4);
  hSBQA->Draw();
  hSB->Draw("same");
  hSignificance->SetLineColor(8);
  hSignificance->SetMarkerColor(8);
  hSignificance->Draw("same");

  cInputQA->cd(1)->SetLogy();
  hSignalAll->GetYaxis()->SetRangeUser(1e-4,10);
  hSignalAll->SetTitle("Signal and Background");
  hSignalAll->Draw("P");
  hBG->Draw("same");
  cInputQA->cd(2)->SetLogy();
  hSBQA->GetYaxis()->SetRangeUser(1e-3,100);
  hSBQA->SetTitle("S/B");
  hSBQA->Draw("");
  cInputQA->cd(3)->SetLogy();
  hSignificance->GetYaxis()->SetRangeUser(1e-4,10);
  hSignificance->SetTitle("Significance");
  hSignificance->Draw("");
  // ----------------------------------------------------------------------
  
  // ----------------------------------------------------------------------
  // Produce Sample Spectrum with Uncertainty using Poisson Distribution according to Significance
  cout << "Now Sample Spectrum" << endl;
  
  TRandom* rand3 = new TRandom3();
  rand3->SetSeed(0);						// initialize random generator
  
  TH1D* hSignalStatError = (TH1D*) hSignal->Clone("hSignalStatError");
	
  Double_t BinWidthMee = hSignal->GetBinWidth(1);
  Double_t PoissonMeas = 0;
  Double_t PoissonErr  = 0;
  Double_t x_val = 0;
  Double_t y_val = 0;
  Double_t x_err = 0;
  Double_t y_err = 0;
  Double_t signif_m = 0;
  for (Int_t ix=1; ix<=hSignal->GetNbinsX(); ix++){
    
    x_val = hSignal->GetBinCenter(ix);
    y_val = hSignal->GetBinContent(ix);
    
    signif_m = hSignificance->GetBinContent(ix);
    signif_m = signif_m * TMath::Sqrt(nEventsLHC) * TMath::Sqrt(BinWidthMee);
    if (signif_m <= 0) {
      Printf("Significance smaller than 0 --> set to small");
      signif_m = 1E-9;
    }
		
    PoissonMeas = rand3->PoissonD(signif_m * signif_m);
    //PoissonErr  = TMath::Sqrt(PoissonMeas) / (signif_m * signif_m); // original version from patrick
    PoissonErr  = TMath::Sqrt(PoissonMeas) / PoissonMeas;  // update from Patrick (11.01.2017)
    //cout<<TMath::Sqrt(PoissonMeas) / (signif_m * signif_m)<<" -> "<<PoissonErr<<endl;

    hSignalStatError->SetBinContent(ix, y_val * PoissonMeas / (signif_m * signif_m));
    hSignalStatError->SetBinError(ix, y_val * PoissonErr);
  
  }
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // Calculate Excess Spectrum
  // Errors are statistical of Poisson sample
	
  Double_t zeroErr[1000] = {0.};
  hCharm->SetError(zeroErr);
  //hBottom->SetError(zeroErr);
  hSignalwoRho->SetError(zeroErr);
  hSignalAll->SetError(zeroErr);
  hSignal->SetError(zeroErr);

  TH1D* hExcess = (TH1D*)hSignalStatError->Clone("hExcess");
  hExcess->Add(hCharm, -1);
  //hExcess->Add(hBottom, -1);
  hExcess->Add(hSignalwoRho, -1);

  
  // ----------------------------------------------------------------------

  // ----------------------------------------------------------------------
  // Calculate Systematic Errors due to Subtraction
  // - of Combinatorial Background (band) in Signal Spectrum
  // - of Combinatorial Background (band) in Excess Spectrum
  // - of Charm + Cocktail (boxes) in Excess Spectrum
  cout << "Now Systematic Errors" << endl;
  
  Double_t GlobalErr_Bottom   = 0.2; // NOT USED HERE

  Double_t GlobalErr_Charm    = 0.02;// (assume measurement (DCA) with 2% uncertainty in ALICE 3)
  if(!Run5)GlobalErr_Charm    = 0.02;// (assume measurement (DCA) with 2% uncertainty in ITS 3, 15% in YR (external))

  Double_t GlobalErr_Cocktail = 0.02;// (assume measurement with 2% uncertainty in ALICE 3)
  if(!Run5)GlobalErr_Cocktail = 0.02;// (assume measurement with 2% uncertainty in ITS 3, 10% in YR)

  Double_t GlobalErr_Signal   = 0.02;// (tracking + PID, use 1% per leg)
  if(!Run5)GlobalErr_Signal   = 0.02;// (2% leg uncertainty in ITS 3, 5% in YR)

  Double_t LocalErr_CombBkg_R = 0.0002;// (as in YR)
  // use a constant syst. error (fom Run 2 Pb-Pb, Carsten's studies on R factor)
  //LocalErr_CombBkg_R = 0.;

  Double_t LocalErr_CombBkg_LS     = 0.;// (YR version)
  Double_t LocalErr_CombBkg_Sig    = 0.;// (YR version)
  Double_t LocalErr_CombBkg_Excess = 0.;// (post-YR version)

  
  if(doNoSyst){
    GlobalErr_Bottom   = 0.;
    GlobalErr_Charm    = 0.;
    GlobalErr_Cocktail = 0.;
    GlobalErr_Signal   = 0.;
    LocalErr_CombBkg_R    = 0.;
    LocalErr_CombBkg_LS   = 0.;
    LocalErr_CombBkg_Sig  = 0.;
  }
  
  // to check what we use for signal unceratinties Run1/2
  Double_t GlobalErr_Signal_mass1 = 0.07;
  Double_t GlobalErr_Signal_mass2 = 0.35;
  Double_t GlobalErr_Signal_mass3 = 0.15;

  Double_t GlobalErr_CocktailRun1_mass1 = 0.31;
  Double_t GlobalErr_CocktailRun1_mass2 = 0.48;
  Double_t GlobalErr_CocktailRun1_mass3 = 0.48;
  
  Int_t skipFirstPoint = 0;
  
  TGraphErrors* fGrSystErr_Signal_CombBkg = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErr_Signal_CombBkg->SetName("fGrSystErr_Signal_CombBkg");
  TGraphErrors* fGrSystErr_Excess_CombBkg = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErr_Excess_CombBkg->SetName("fGrSystErr_Excess_CombBkg");
  TGraphErrors* fGrSystErr_Charm_Cocktail = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErr_Charm_Cocktail->SetName("fGrSystErr_Charm_Cocktail");
  TGraphErrors* fGrSystErr_Signal = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErr_Signal->SetName("fGrSystErr_Signal");

  TGraphErrors* fGrSystErrOnly_Signal_CombBkg = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Signal_CombBkg->SetName("fGrSystErrOnly_Signal_CombBkg");
  TGraphErrors* fGrSystErrOnly_Signal_S = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Signal_S->SetName("fGrSystErrOnly_Signal_S");
  TGraphErrors* fGrSystErrOnly_Signal_LS = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Signal_LS->SetName("fGrSystErrOnly_Signal_LS");
  TGraphErrors* fGrSystErrOnly_Signal_R = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Signal_R->SetName("fGrSystErrOnly_Signal_R");
  TGraphErrors* fGrSystErrOnly_Signal = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Signal->SetName("fGrSystErrOnly_Signal");
  TGraphErrors* fGrSystErrOnly_Cocktail = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Cocktail->SetName("fGrSystErrOnly_Cocktail");
  TGraphErrors* fGrSystErrOnly_CocktailRun1 = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_CocktailRun1->SetName("fGrSystErrOnly_CocktailRun1");

  TGraphErrors* fGrSystErrOnly_Excess_CombBkg = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Excess_CombBkg->SetName("fGrSystErrOnly_Excess_CombBkg");
  TGraphErrors* fGrSystErrOnly_Charm_Cocktail = new TGraphErrors( (Int_t) hSignalStatError->GetNbinsX() -skipFirstPoint );
  fGrSystErrOnly_Charm_Cocktail->SetName("fGrSystErrOnly_Charm_Cocktail");
  
  TGraphErrors* fGrFlowPtSignifStatBkgCharm = NULL;
  TGraphErrors* fGrFlowPtSignifStatBkg = NULL;
  TGraphErrors* fGrFlowPtSignifStat = NULL;
  
  Double_t y_errExcess   = 0;
  Double_t y_errStat     = 0;
  Double_t y_errCharm    = 0;
  Double_t y_errBottom   = 0;
  Double_t y_errCocktail = 0;
  Double_t y_errSignal   = 0;
  Double_t soverb_m = 0;
 
  TF1 *fDeltaRfactor  = new TF1("fDeltaRfactor","[0]+[1]/x+[2]/(x*x)",0,1.5);
  TF1 *fDeltaRfactor2 = new TF1("fDeltaRfactor2","[0]+[1]/x+[2]/(x*x)",0,1.5);
  fDeltaRfactor->SetParameters(-0.00103727,-7.79264e-05,0);
  fDeltaRfactor2->SetParameters(-0.00101783,-6.44731e-05,0);

  cout<<"Delta R factor at 0.3 GeV/c2 = "<<(fDeltaRfactor->Eval(0.3))*100<<"% and "<<(fDeltaRfactor2->Eval(0.3))*100<<"% ==> "<<TMath::Sqrt(fDeltaRfactor->Eval(0.3)*fDeltaRfactor->Eval(0.3)+fDeltaRfactor2->Eval(0.3)*fDeltaRfactor2->Eval(0.3))*100<<"%"<<endl;
  cout<<"Delta R factor at 0.5 GeV/c2 = "<<(fDeltaRfactor->Eval(0.5))*100<<"% and "<<(fDeltaRfactor2->Eval(0.5))*100<<"% ==> "<<TMath::Sqrt(fDeltaRfactor->Eval(0.5)*fDeltaRfactor->Eval(0.5)+fDeltaRfactor2->Eval(0.5)*fDeltaRfactor2->Eval(0.5))*100<<"%"<<endl;
  cout<<"Delta R factor at 1.0 GeV/c2 = "<<(fDeltaRfactor->Eval(1.0))*100<<"% and "<<(fDeltaRfactor2->Eval(1.0))*100<<"% ==> "<<TMath::Sqrt(fDeltaRfactor->Eval(1.0)*fDeltaRfactor->Eval(1.0)+fDeltaRfactor2->Eval(1.0)*fDeltaRfactor2->Eval(1.0))*100<<"%"<<endl;
  
  for (Int_t ix=1+skipFirstPoint; ix<=hSignalStatError->GetNbinsX(); ix++){

      // for Signal
      x_val  = hSignalStatError->GetBinCenter(ix);
      y_val  = hSignalStatError->GetBinContent(ix);
      y_errStat  = hSignalStatError->GetBinError(ix);
      
      // Comb Bkg
      // abs. Error in Signal: S * dS/S ~= S * dB/S = S * dB/B * B/S
      //soverb_m = hSB->GetBinContent(ix);//original
      soverb_m = hSBQA->GetBinContent(ix);//modified background
      if (soverb_m <= 0) {
	soverb_m = 0.000001;
      }
      
      x_err = 0.0125;
      LocalErr_CombBkg_Sig = 0.;
      //LocalErr_CombBkg_R   = 0.;
      LocalErr_CombBkg_LS  = 0.;
      
      LocalErr_CombBkg_Sig = TMath::Sqrt(LocalErr_CombBkg_R*LocalErr_CombBkg_R + LocalErr_CombBkg_LS*LocalErr_CombBkg_LS)/soverb_m;// BG uncertainties
      LocalErr_CombBkg_Excess = TMath::Sqrt(LocalErr_CombBkg_R*LocalErr_CombBkg_R + LocalErr_CombBkg_LS*LocalErr_CombBkg_LS)/soverb_m;// BG uncertainties
      LocalErr_CombBkg_Sig = TMath::Sqrt(LocalErr_CombBkg_Sig*LocalErr_CombBkg_Sig + GlobalErr_Signal*GlobalErr_Signal);//signal uncertainties



      y_err = y_val * LocalErr_CombBkg_Sig;
      y_errExcess = y_val * LocalErr_CombBkg_Excess;
	    
      fGrSystErr_Signal_CombBkg->SetPoint(ix-1-skipFirstPoint, x_val, y_val);
      fGrSystErr_Signal_CombBkg->SetPointError(ix-1-skipFirstPoint, x_err, y_err);
      fGrSystErrOnly_Signal_LS->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      fGrSystErrOnly_Signal_LS->SetPointError(ix-1-skipFirstPoint, x_err, LocalErr_CombBkg_LS/soverb_m);
      fGrSystErrOnly_Signal_R->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      fGrSystErrOnly_Signal_R->SetPointError(ix-1-skipFirstPoint, x_err, LocalErr_CombBkg_R/soverb_m);
      fGrSystErrOnly_Signal_S->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      fGrSystErrOnly_Signal_S->SetPointError(ix-1-skipFirstPoint, x_err, GlobalErr_Signal);
      fGrSystErrOnly_Signal_CombBkg->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      fGrSystErrOnly_Signal_CombBkg->SetPointError(ix-1-skipFirstPoint, x_err, y_err/y_val);
      
      // Signal (to check what we apply for Run1/2)
      x_err = 0.006; // for visualization
      if (doPtee) { 
	x_err *= 2; 
      } // because of the larger binwidth in Ptee spectra.
      if(x_val<0.5)
	y_errSignal     = y_val * GlobalErr_Signal_mass1;
      else if(x_val<2.8)
	y_errSignal     = y_val * GlobalErr_Signal_mass2;
      else
	y_errSignal     = y_val * GlobalErr_Signal_mass3;
      fGrSystErr_Signal->SetPoint(ix-1, x_val, y_val);
      fGrSystErr_Signal->SetPointError(ix-1, x_err,y_errSignal);
      fGrSystErrOnly_Signal->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      fGrSystErrOnly_Signal->SetPointError(ix-1-skipFirstPoint, x_err, y_errSignal/ y_val);
      
      // for Excess
      x_val  = hExcess->GetBinCenter(ix);
      y_val  = hExcess->GetBinContent(ix);
      
      // Comb Bkg
      // the same absolute Errors as above are applied
      x_err = 0.0125; // for visualization
      if (doPtee) { 
	x_err *= 2; 
      } // because of the larger binwidth in Ptee spectra.

      LocalErr_CombBkg_Excess = y_errExcess/y_val;
      LocalErr_CombBkg_Excess = TMath::Sqrt(LocalErr_CombBkg_Excess*LocalErr_CombBkg_Excess + GlobalErr_Signal*GlobalErr_Signal);//signal uncertainties
      y_errExcess = LocalErr_CombBkg_Excess*y_val;
      
      fGrSystErr_Excess_CombBkg->SetPoint(ix-1-skipFirstPoint, x_val, y_val);
      fGrSystErr_Excess_CombBkg->SetPointError(ix-1-skipFirstPoint, x_err, y_errExcess);
      fGrSystErrOnly_Excess_CombBkg->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      fGrSystErrOnly_Excess_CombBkg->SetPointError(ix-1-skipFirstPoint, x_err, y_errExcess/y_val);
      
      // Charm + Cocktail
      x_err = 0.006; // for visualization
      if (doPtee) { 
	x_err *= 2; 
      } // because of the larger binwidth in Ptee spectra.
      y_errCharm    = hCharm->GetBinContent(ix) * GlobalErr_Charm;
      //y_errBottom   = hBottom->GetBinContent(ix) * GlobalErr_Bottom;
      y_errCocktail = hSignalwoRho->GetBinContent(ix) * GlobalErr_Cocktail;
      y_err = y_errCharm + y_errBottom + y_errCocktail;
      fGrSystErr_Charm_Cocktail->SetPoint(ix-1, x_val, y_val);
      fGrSystErr_Charm_Cocktail->SetPointError(ix-1, x_err, y_err);
      fGrSystErrOnly_Charm_Cocktail->SetPoint(ix-1, x_val, 0);
      fGrSystErrOnly_Charm_Cocktail->SetPointError(ix-1, x_err, y_err/y_val);
      fGrSystErrOnly_Cocktail->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      fGrSystErrOnly_Cocktail->SetPointError(ix-1-skipFirstPoint, x_err, GlobalErr_Cocktail);

      fGrSystErrOnly_CocktailRun1->SetPoint(ix-1-skipFirstPoint, x_val, 0);
      if(x_val<0.5)
	fGrSystErrOnly_CocktailRun1->SetPointError(ix-1-skipFirstPoint, x_err, GlobalErr_CocktailRun1_mass1);
      else if(x_val<2.8)
	fGrSystErrOnly_CocktailRun1->SetPointError(ix-1-skipFirstPoint, x_err, GlobalErr_CocktailRun1_mass2);
      else
	fGrSystErrOnly_CocktailRun1->SetPointError(ix-1-skipFirstPoint, x_err, GlobalErr_CocktailRun1_mass3);

  }

    // Systematic error comparison
  TCanvas *cSystComp = new TCanvas("cSystComp","cSystComp",1200,500);
  cSystComp->Divide(2,1);
  
  cSystComp->cd(1);
  fGrSystErrOnly_Signal_CombBkg->SetTitle("Systematic uncertainty (Background/Signal)");
  fGrSystErrOnly_Signal_CombBkg->GetXaxis()->SetTitle(titleXaxis);
  fGrSystErrOnly_Signal_CombBkg->GetXaxis()->SetRangeUser(0,2.5);
  fGrSystErrOnly_Signal_CombBkg->GetYaxis()->SetRangeUser(0,1);  
  fGrSystErrOnly_Signal_CombBkg->SetFillColor(kGray+1);
  fGrSystErrOnly_Signal_CombBkg->SetLineColor(kGray+1);
  fGrSystErrOnly_Signal_CombBkg->SetFillStyle(3244);
  fGrSystErrOnly_Signal_CombBkg->Draw("A2");
  fGrSystErrOnly_Signal_S->SetFillColor(kGreen+1);
  fGrSystErrOnly_Signal_S->SetLineColor(kGreen+1);
  fGrSystErrOnly_Signal_S->SetFillStyle(3244);
  fGrSystErrOnly_Signal_S->Draw("2");
   // TO BE REMOVED
  // fGrSystErrOnly_Signal_LS->SetFillColor(kBlue+1);
  // fGrSystErrOnly_Signal_LS->SetLineColor(kBlue+1);
  // fGrSystErrOnly_Signal_LS->SetFillStyle(3245);
  // //fGrSystErrOnly_Signal_LS->Draw("2");
  fGrSystErrOnly_Signal_R->SetFillColor(kOrange);
  fGrSystErrOnly_Signal_R->SetLineColor(kOrange);
  fGrSystErrOnly_Signal_R->SetFillStyle(3246);
  fGrSystErrOnly_Signal_R->Draw("2");

  fGrSystErrOnly_Signal->SetFillColor(kRed);
  fGrSystErrOnly_Signal->SetLineColor(kRed);
  fGrSystErrOnly_Signal->SetFillStyle(3244);
  //fGrSystErrOnly_Signal->Draw("Csame");

  // Legend
  Double_t leg_textsize = 0.029;
  TLegend *legSyst=new TLegend(0.45,0.66,0.88,0.89);
  legSyst->SetFillColor(0);
  legSyst->SetBorderSize(0);
  legSyst->SetTextSize(leg_textsize); 

  legSyst->AddEntry(fGrSystErrOnly_Signal_CombBkg,"from Sig. + BG (quadr. sum)","f");
  
   // TO BE REMOVED
  //legSyst->AddEntry(fGrSystErrOnly_Signal_LS,"from Sig. + BG (LS)","f");
  legSyst->AddEntry(fGrSystErrOnly_Signal_R,"from BG (R factor)","f");
  //legSyst->AddEntry(fGrSystErrOnly_Signal,"Run1, Pb-Pb 0-10%","l");
  
  legSyst->Draw();


  cSystComp->cd(2);
  fGrSystErrOnly_Excess_CombBkg->SetTitle("Systematic uncertainty (Excess)");
  fGrSystErrOnly_Excess_CombBkg->GetXaxis()->SetTitle(titleXaxis);
  fGrSystErrOnly_Excess_CombBkg->GetXaxis()->SetRangeUser(0,2.5);
  fGrSystErrOnly_Excess_CombBkg->GetYaxis()->SetRangeUser(0,1);  
  fGrSystErrOnly_Excess_CombBkg->SetFillColor(kGray+1);
  fGrSystErrOnly_Excess_CombBkg->SetLineColor(kGray+1);
  fGrSystErrOnly_Excess_CombBkg->SetFillStyle(3244);
  fGrSystErrOnly_Excess_CombBkg->Draw("A2");
  fGrSystErrOnly_Charm_Cocktail->SetFillColor(kOrange);
  fGrSystErrOnly_Charm_Cocktail->SetLineColor(kOrange);
  fGrSystErrOnly_Charm_Cocktail->SetFillStyle(3246);
  fGrSystErrOnly_Charm_Cocktail->Draw("2");

   // Legend
  TLegend *legSyst2=new TLegend(0.45,0.66,0.88,0.89);
  legSyst2->SetFillColor(0);
  legSyst2->SetBorderSize(0);
  legSyst2->SetTextSize(leg_textsize); 

  legSyst2->AddEntry(fGrSystErrOnly_Excess_CombBkg,"from Signal + BG","f");
  legSyst2->AddEntry(fGrSystErrOnly_Charm_Cocktail,"from cocktail and charm","f");
  legSyst2->Draw();

  if(meebin==0)
    cSystComp->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_Syst_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
  else
    cSystComp->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_Syst_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
   // ----------------------------------------------------------------------


  // ----------------------------------------------------------------------
  // Plotting
  if(!doPlot)
    return;

  cout << "Now Plotting" << endl;

  // ----------------------------------------------------------------------
  // Signals
  
  char s_cent[10]  = "5TeV_cent";
  Int_t canvWidth  = 670;	// 670x700 ok to get vertical pdf when importing to ppt.
  Int_t canvHeight = 700;
  Int_t canvXshift = 50;
  TCanvas* fCanvSignalsVsMee = new TCanvas(Form("cSignalsVs%s_%s_%dE%d_bin%d%s", MeeOrPtee, s_cent, nIntEventsLHC, nLogEventsLHC, meebin, canvname_opt), Form("Signals vs %s", MeeOrPtee), canvXshift, 0, canvWidth, canvHeight);
  fCanvSignalsVsMee->SetLeftMargin(0.13);
  fCanvSignalsVsMee->cd(1)->SetLogy();
  
  TH1D* fH1D_axis_mee = (TH1D*) hRho->Clone("fH1D_axis_mee");
  fH1D_axis_mee->SetAxisRange(0., xmaxMee, "X");
  fH1D_axis_mee->SetLineColor(0);
  fH1D_axis_mee->SetMarkerColor(0);
  fH1D_axis_mee->SetTitle(Form(";%s;%s", titleXaxis, titleYaxis));
  fH1D_axis_mee->SetTitleOffset(1.5,"Y");
  fH1D_axis_mee->SetMinimum(ymindNdx);
  fH1D_axis_mee->SetMaximum(ymaxdNdx);
  fH1D_axis_mee->DrawCopy("axis");
  
  // Systematic Uncertainties
  fGrSystErr_Signal_CombBkg->SetFillColor(kGreen+1);
  fGrSystErr_Signal_CombBkg->SetFillStyle(3244);
  if (!doPlain) fGrSystErr_Signal_CombBkg->Draw("2same");

  // Sum
  hSignalAll->SetLineColor(kBlack);
  hSignalAll->DrawCopy("csame");
  hSignalStatError->SetMarkerStyle(20); //(20 or 24);
  hSignalStatError->SetMarkerColor(kBlack); 
  hSignalStatError->SetLineColor(kBlack);
  hSignalStatError->GetXaxis()->SetRangeUser(0,xmaxMee);
  if (!doPlain) hSignalStatError->DrawCopy("p same");
  
  // Rapp
  hRhoRR->SetFillColor(0);
  hRhoRR->SetLineColor(kRed);
  hQGP->SetFillColor(0);
  hQGP->SetLineColor(kOrange+1);
  hRhoRR->DrawCopy("csame");
  hQGP->DrawCopy("csame");
  
  TH1D* fH1D_Rapp_Sum = (TH1D*) hRhoRR->Clone("fH1D_Rapp_Sum");
  fH1D_Rapp_Sum->SetFillColor(0);
  fH1D_Rapp_Sum->SetLineColor(kRed+1);
  fH1D_Rapp_Sum->Add(hQGP);

  TH1D* fH1D_Rapp_Sum_drop = (TH1D*) hRhoDropRR->Clone("fH1D_Rapp_Sum_drop");
  fH1D_Rapp_Sum_drop->SetFillColor(0);
  fH1D_Rapp_Sum_drop->SetLineColor(kRed);
  fH1D_Rapp_Sum_drop->SetLineStyle(2);
  fH1D_Rapp_Sum_drop->SetLineWidth(2);
  fH1D_Rapp_Sum_drop->Add(hQGP);

  TH1D* fH1D_Rapp_Sum_vac = (TH1D*) hRhoVacRR->Clone("fH1D_Rapp_Sum_vac");
  fH1D_Rapp_Sum_vac->SetFillColor(0);
  fH1D_Rapp_Sum_vac->SetLineColor(kBlue+3);
  fH1D_Rapp_Sum_vac->SetLineStyle(2);
  fH1D_Rapp_Sum_vac->SetLineWidth(2);
  fH1D_Rapp_Sum_vac->Add(hQGP);
  

  // Cocktail sum
  hSignalwoRho->SetFillColor(0);
  hSignalwoRho->SetLineColor(kCyan+1);
  hSignalwoRho->SetMarkerColor(kCyan+1);
  hSignalwoRho->Draw("c same"); //("csamehist");

  // Charm and bottom
  hCharm->SetFillColor(0);
  hCharm->SetLineColor(kMagenta);
  hCharm->SetMarkerColor(kMagenta);
  hCharm->Draw("c same"); //("csamehist");
  // hBottom->SetFillColor(0);
  // hBottom->SetLineColor(kMagenta);
  // hBottom->SetMarkerColor(kMagenta);
  // hBottom->Draw("c same"); //("csamehist");

  // Legend
  leg_textsize = 0.029;
  TLegend *leg1=new TLegend(0.46,0.59,0.88,0.82);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(leg_textsize); 

  leg1->AddEntry(hSignalAll,"Sum","l");
  leg1->AddEntry(hRhoRR,"Rapp in-medium SF","l");
  leg1->AddEntry(hQGP,"Rapp QGP","l");
  
  if (!doPlain) leg1->AddEntry(hSignalwoRho, Form("cocktail w/o #rho (#pm %1.0f%%)", 100.*GlobalErr_Cocktail), "l");
  else          leg1->AddEntry(hSignalwoRho, Form("cocktail w/o #rho"), "l");
  if (!doPlain) leg1->AddEntry(hCharm,Form("c#bar{c} #rightarrow ee (#pm %1.0f%%)", 100.*GlobalErr_Charm), "l");
  else          leg1->AddEntry(hCharm,Form("c#bar{c} #rightarrow ee"), "l");
  
  if (!doPlain) leg1->AddEntry(hSignalStatError, Form("L_{int}=%.0f nb^{-1} 'measured'",nEventsLHC*1.2e-9), "ple");
  if (!doPlain && !doNoSyst) leg1->AddEntry(fGrSystErr_Signal_CombBkg,Form("Syst. err. sig.+ bkg."), "f");
  leg1->Draw();

  TLatex *   logo = new TLatex(0.16,0.84,"ALICE Upgrade Simulation");
  logo->SetNDC();
  logo->SetTextSize(0.04);
  logo->SetTextFont(42);
  logo->DrawClone();
  
  TPaveText* fPave = new TPaveText(0.16, 0.52, .475, .81, "NDC");
  fPave->AddText(Form("%s #sqrt{#it{s}_{NN}} = 5.5 TeV",sSystem2[iSystem].Data()));
  if (!doPlain) fPave->AddText(Form("0 - 10%%, %3.1fE%d events", nDouEventsLHC, nLogEventsLHC));
  else          fPave->AddText(Form("0 - 10%%"));
    if(!Run5)
    fPave->AddText("Run3");
  else
    fPave->AddText("Run5+");
  fPave->AddText("|#eta_{e}| < 0.8");
  fPave->AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}"); //to be done, change all pT cuts in analysis chain to lowest possible value (needs to be done for all predictions consistently)
  if(doPtee)
    fPave->AddText(Form("%.2f < #it{M}_{ee} (GeV/#it{c }^{2}) < %.2f", MeeRangeMin, MeeRangeMax));
  fPave->SetTextAlign(11);
  fPave->SetTextSize(leg_textsize);
  fPave->SetBorderSize(0);
  fPave->SetFillStyle(0);
  fPave->DrawClone();

  if(meebin==0)
    fCanvSignalsVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_Mee_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
  else
    fCanvSignalsVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_Mee_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));

  cInputQA->Update();

  if(meebin==0){
    cInputQA->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_inputQA_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
  }
  else{
    cInputQA->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_inputQA_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
  }

  // ----------------------------------------------------------------------
  // Excess
  
  TCanvas* fCanvExcessVsMee = new TCanvas(Form("cExcessVs%s_%s_%dE%d_bin%d%s", MeeOrPtee, s_cent, nIntEventsLHC, nLogEventsLHC, meebin, canvname_opt), Form("Excess Signal vs %s", MeeOrPtee), canvXshift+700, 0, canvWidth, canvHeight);
  fCanvExcessVsMee->SetLeftMargin(0.13);
  fCanvExcessVsMee->cd(1)->SetLogy();
  
  Int_t reduceMaxFac = 1;
  fH1D_axis_mee->SetMinimum(ymindNdx);
  fH1D_axis_mee->SetMaximum(ymaxdNdx/reduceMaxFac);
  fH1D_axis_mee->DrawCopy("axis");
  
  // Rapp
  hRhoRR->DrawCopy("chistsame");
  hQGP->DrawCopy("csame");
  fH1D_Rapp_Sum->SetFillColor(0);
  fH1D_Rapp_Sum->SetLineColor(kBlack);
  fH1D_Rapp_Sum->DrawCopy("csame");
  
  // Systematic Uncertainties
  fGrSystErr_Excess_CombBkg->SetFillColor(kGreen+1);
  fGrSystErr_Excess_CombBkg->SetFillStyle(3244);
  fGrSystErr_Excess_CombBkg->Draw("2 same");
  fGrSystErr_Charm_Cocktail->SetFillColor(kMagenta);
  fGrSystErr_Charm_Cocktail->SetFillStyle(3105);  //(3001); not viewable in Illustrator
  fGrSystErr_Charm_Cocktail->Draw("2 same");
  
  // Excess Spectrum
  hExcess->SetLineColor(kBlack);
  hExcess->SetMarkerColor(kBlack);
  hExcess->SetMarkerStyle(25);
  hExcess->GetXaxis()->SetRangeUser(0,xmaxMee);
  hExcess->DrawCopy("same");
  
  // Legend
  TLegend *leg2=new TLegend(0.46,0.59,0.88,0.82);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(leg_textsize); 
  leg2->AddEntry(fH1D_Rapp_Sum,"Rapp Sum","l");
  leg2->AddEntry(hRhoRR,"Rapp in-medium SF","l");
  leg2->AddEntry(hQGP,"Rapp QGP","l");
  leg2->AddEntry(hExcess,Form("L_{int}=%.0f nb^{-1} 'measured'", nEventsLHC*1.2e-9), "ple");
  if(!doNoSyst) leg2->AddEntry(fGrSystErr_Excess_CombBkg,Form("Syst. err. sig. + bkg."), "f");
  if(!doNoSyst) leg2->AddEntry(fGrSystErr_Charm_Cocktail,Form("Syst. err. c#bar{c} + cocktail"), "f");
  leg2->Draw();
  logo->DrawClone();
  fPave->DrawClone();
  fPave->DrawClone();

  if(meebin==0){
    fCanvExcessVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_Excess_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
  }
  else{
    fCanvExcessVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_Excess_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));

  }

 
  // ----------------------------------------------------------------------

  //------------------------------
  // Excess - 2
  //------------------------------
  
  TCanvas* fCanvExcessVsMee2 = new TCanvas(Form("cExcessVs%s_%s_%dE%d_bin%d%s_2", MeeOrPtee, s_cent, nIntEventsLHC, nLogEventsLHC, meebin, canvname_opt), Form("Excess Signal vs %s", MeeOrPtee), canvXshift+700, 0, canvWidth, canvHeight);
  fCanvExcessVsMee2->SetLeftMargin(0.13);
  fCanvExcessVsMee2->cd(1)->SetLogy();
  fH1D_axis_mee->DrawCopy("axis");

  // Rapp
  fH1D_Rapp_Sum->SetFillColor(0);
  fH1D_Rapp_Sum->SetLineColor(kBlack);
  fH1D_Rapp_Sum->DrawCopy("csame");
  fH1D_Rapp_Sum_vac->DrawCopy("csame");
  fH1D_Rapp_Sum_drop->DrawCopy("csame");
  
  // Systematic Uncertainties
  fGrSystErr_Excess_CombBkg->SetFillColor(kGreen+1);
  fGrSystErr_Excess_CombBkg->SetFillStyle(3244);
  fGrSystErr_Excess_CombBkg->Draw("2 same");
  fGrSystErr_Charm_Cocktail->SetFillColor(kMagenta);
  fGrSystErr_Charm_Cocktail->SetFillStyle(3105);  //(3001); not viewable in Illustrator
  fGrSystErr_Charm_Cocktail->Draw("2 same");
  
  // Excess Spectrum
  hExcess->SetLineColor(kBlack);
  hExcess->SetMarkerColor(kBlack);
  hExcess->SetMarkerStyle(25);
  hExcess->GetXaxis()->SetRangeUser(0,xmaxMee);
  hExcess->DrawCopy("same");

  // for(Int_t iBin = 0; iBin < hExcess->GetNbinsX(); iBin++){
  //   cout<<iBin<<" "<<hExcess->GetBinCenter(iBin+1)<<" "<<hExcess->GetBinContent(iBin+1)<<"    "<<hExcess->GetBinError(iBin+1)/hExcess->GetBinContent(iBin+1)<<endl;
  // }
  
  // Legend
  TLegend *leg3=new TLegend(0.46,0.59,0.88,0.82);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(leg_textsize); 
  leg3->AddEntry(fH1D_Rapp_Sum,"Rapp Sum (broad #rho)","l");
  leg3->AddEntry(fH1D_Rapp_Sum_drop,"Rapp Sum (drop #rho)","l");
  leg3->AddEntry(fH1D_Rapp_Sum_vac,"Rapp Sum (vac #rho)","l");
  leg3->AddEntry(hExcess,Form("L_{int}=%.0f nb^{-1} 'measured'", nEventsLHC*1.2e-9), "ple");
  if(!doNoSyst) leg3->AddEntry(fGrSystErr_Excess_CombBkg,Form("Syst. err. sig. + bkg."), "f");
  if(!doNoSyst) leg3->AddEntry(fGrSystErr_Charm_Cocktail,Form("Syst. err. c#bar{c} + cocktail"), "f");
  leg3->Draw();
  
  fPave->DrawClone();
  fPave->DrawClone();
  logo->DrawClone();

  if(meebin==0)
    fCanvExcessVsMee2->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_Excess2_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
  else
    fCanvExcessVsMee2->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_Excess2_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));


  if(!doPtee){
	
	//------------------------------
	// Polarization
	//------------------------------
        Double_t NBinsAngularDistribution = 80.; // 8 x 10 in cos(theta) and phi (CS frame, see NA 60 paper: PRL 102, 222301 (2009))
	cout << "Now Polarization (with "<<NBinsAngularDistribution<<" bins in angular distribution)" << endl; // (if " << doPtee << " == 1)" << endl;
	Bool_t checkPolarization  = 0;
	
	// for Polarization vs Ptee:
	const Int_t npolarizationPtBins = 3;
	Double_t polarizationPtBins[npolarizationPtBins+1] = {0.1, 0.8, 1.5, 3.};
		
	// for Polarization vs Mee:
	const Int_t npolarizationMeeBins = 5;
	Double_t polarizationMeeBins[npolarizationMeeBins+1] = {0.15, 0.4, 0.6, 0.9, 1.5,2.5};
	
	Int_t nbins = 0;
	if (doPtee) nbins = npolarizationPtBins;
	else        nbins = npolarizationMeeBins;
	
	const Int_t npolarizationBins = nbins;
	Double_t polarizationBins[npolarizationBins+1];
	for (Int_t ipolarization=0; ipolarization<npolarizationBins+1; ipolarization++)
	{
		if (doPtee) polarizationBins[ipolarization] = polarizationPtBins[ipolarization];
		else        polarizationBins[ipolarization] = polarizationMeeBins[ipolarization];
	}
	
	Double_t polarizationX[npolarizationBins];
	Double_t polarizationY[npolarizationBins];
	Double_t polarizationXerr[npolarizationBins];
	Double_t polarizationYerrStat[npolarizationBins];
	Double_t polarizationYerrBkg[npolarizationBins];
	Double_t polarizationYerrCharm[npolarizationBins];
	
	Int_t nPointsInBin = 0;
	Int_t nLostPoints  = 0;
	Double_t ySumPtBin = 0;
	Double_t ySignifStat2PtBin  = 0;
	Double_t ySignifBkg2PtBin   = 0;
	Double_t ySignifCharm2PtBin = 0;
	Double_t y_err_stat  = 0;
	Double_t y_err_bkg   = 0;
	Double_t y_err_charm = 0;
	
//		Int_t idebug = 0;
//		cout << "  debug " << idebug << endl; idebug++;
	
	for (Int_t ipolarization=0; ipolarization<npolarizationBins; ipolarization++)
	{
		nPointsInBin=0;
		nLostPoints=0;
		ySumPtBin=0;
		ySignifStat2PtBin=0;
		ySignifBkg2PtBin=0;
		ySignifCharm2PtBin=0;
		
		skipFirstPoint=0;
		for (Int_t ix=1+skipFirstPoint; ix<=hExcess->GetNbinsX(); ix++) 
		{
			x_val  = hExcess->GetBinCenter(ix);
			y_val  = hExcess->GetBinContent(ix);
			y_err_stat  = hExcess->GetBinError(ix);
			y_err_bkg   = fGrSystErr_Excess_CombBkg->GetEY()[ix-1];
			y_err_charm = fGrSystErr_Charm_Cocktail->GetEY()[ix-1];
			
			if (   (fGrSystErr_Excess_CombBkg->GetX()[ix-1] != x_val)
					|| (fGrSystErr_Charm_Cocktail->GetX()[ix-1] != x_val) )
			{
				cout << " ERROR in Polarization: Graph from Sampling has other binning as histogram!" << endl;
				break;
			}
			
			if (x_val < polarizationBins[ipolarization]) continue;
			if (x_val > polarizationBins[ipolarization+1]) break;
//			cout << " x_val, y_val, y_err_stat = " << x_val << " , " << y_val << " , " << y_err_stat << endl;
//			cout << "   y_err_bkg, y_err_charm = " << y_err_bkg << " , " << y_err_charm << endl;
			
			nPointsInBin++;
			if (y_val <= 0) {
				nLostPoints++;
				continue;
			}
			
			ySumPtBin  += y_val;
			ySignifStat2PtBin  += TMath::Power(y_val/y_err_stat, 2);
			ySignifBkg2PtBin   += TMath::Power(y_val/y_err_bkg, 2);
			ySignifCharm2PtBin += TMath::Power(y_val/y_err_charm, 2);
		}
		
//		cout << " nPointsInBin, nLostPoints = " << nPointsInBin << " , " << nLostPoints << endl;
		if (nPointsInBin-nLostPoints) { // != 0
			ySumPtBin  *= (Double_t) nPointsInBin / (nPointsInBin-nLostPoints);
		}
		// this correction is an overestimate:
		Double_t err_corr_fac = 1.; // = (Double_t) nPointsInBin / (nPointsInBin-nLostPoints);
		ySignifStat2PtBin  *= err_corr_fac;
		ySignifBkg2PtBin   *= err_corr_fac;
		ySignifCharm2PtBin *= err_corr_fac;
		
		polarizationX[ipolarization] = (polarizationBins[ipolarization+1] - polarizationBins[ipolarization]) / 2 + polarizationBins[ipolarization];
		//			cout << "  before div" << endl;
		polarizationY[ipolarization] = ySumPtBin / nPointsInBin;
		polarizationXerr[ipolarization]	= (polarizationBins[ipolarization+1] - polarizationBins[ipolarization]) / 2;
		if (nPointsInBin-nLostPoints) { // != 0
			polarizationYerrStat[ipolarization]	 = polarizationY[ipolarization] / TMath::Sqrt(ySignifStat2PtBin);
			polarizationYerrBkg[ipolarization]	 = polarizationY[ipolarization] / TMath::Sqrt(ySignifBkg2PtBin); //  + polarizationYerrStat[ipolarization];
			polarizationYerrCharm[ipolarization] = polarizationY[ipolarization] / TMath::Sqrt(ySignifCharm2PtBin); // + polarizationYerrBkg[ipolarization];
		}
		else {
			polarizationYerrStat[ipolarization]	 = 0.;
			polarizationYerrBkg[ipolarization]	 = 0.;
			polarizationYerrCharm[ipolarization] = 0.;
		}
//			cout << polarizationX[ipolarization] << endl;
//  		        cout << polarizationY[ipolarization] << endl;
//			cout << polarizationXerr[ipolarization] << endl;
//            		cout << polarizationYerrStat[ipolarization] << endl;
//			cout << polarizationYerrBkg[ipolarization] << endl;
//			cout << polarizationYerrCharm[ipolarization] << endl;
//       		cout<<endl;	
	}
	
	
	
	TGraphErrors* fGrPolarizationPtErrStat = new TGraphErrors(npolarizationBins, polarizationX, polarizationY, polarizationXerr, polarizationYerrStat);
	TGraphErrors* fGrPolarizationPtErrStatBkg = new TGraphErrors(npolarizationBins, polarizationX, polarizationY, polarizationXerr, polarizationYerrBkg);
	TGraphErrors* fGrPolarizationPtErrStatBkgCharm = new TGraphErrors(npolarizationBins, polarizationX, polarizationY, polarizationXerr, polarizationYerrCharm);
	
	fGrPolarizationPtErrStatBkgCharm->SetMarkerStyle(21);
	fGrPolarizationPtErrStatBkgCharm->SetMarkerColor(kMagenta);
	fGrPolarizationPtErrStatBkgCharm->SetLineColor(kMagenta);
	fGrPolarizationPtErrStatBkgCharm->SetLineWidth(4);
	fGrPolarizationPtErrStatBkg->SetMarkerStyle(21);
	fGrPolarizationPtErrStatBkg->SetMarkerColor(kGreen+1);
	fGrPolarizationPtErrStatBkg->SetLineColor(kGreen+1);
	fGrPolarizationPtErrStatBkg->SetLineWidth(3);
	fGrPolarizationPtErrStat->SetMarkerStyle(21);
	fGrPolarizationPtErrStat->SetMarkerColor(kBlack);
	fGrPolarizationPtErrStat->SetLineColor(kBlack);
	fGrPolarizationPtErrStat->SetLineWidth(2);
	if (checkPolarization) {
		fCanvExcessVsMee->cd();
		//fGrPolarizationPtErrStatBkgCharm->DrawClone("p same");
		//fGrPolarizationPtErrStatBkg->DrawClone("p same");
		fGrPolarizationPtErrStat->DrawClone("p same");
	}
	
	TGraphErrors* fGrPolarizationPtSignifStatBkgCharm = (TGraphErrors*) fGrPolarizationPtErrStatBkgCharm->Clone("fGrPolarizationPtSignifStatBkgCharm");
	TGraphErrors* fGrPolarizationPtSignifStatBkg = (TGraphErrors*) fGrPolarizationPtErrStatBkg->Clone("fGrPolarizationPtSignifStatBkg");
	TGraphErrors* fGrPolarizationPtSignifStat = (TGraphErrors*) fGrPolarizationPtErrStat->Clone("fGrPolarizationPtSignifStat");
	
	// ----------------------------------------------------------------------
	// Plot Polarization
	cout << " plotting Polarization" << endl;

	Double_t canvYstretchPolarization = 0.8;
	
	TCanvas* fCanvPolarizationVsMee = new TCanvas(Form("cPolarizationVs%s_%s_%dE%d_bin%d%s", MeeOrPtee, s_cent, nIntEventsLHC, nLogEventsLHC, meebin, canvname_opt), Form("Polarization Background vs %s", MeeOrPtee), canvXshift+300, 0, canvWidth*canvYstretchPolarization, canvHeight*canvYstretchPolarization);
	fCanvPolarizationVsMee->SetLeftMargin(0.13);
	fCanvPolarizationVsMee->cd(1)->SetLogy();
	
	for (Int_t ipolarization=0; ipolarization<npolarizationBins; ipolarization++)
	{
		x_val = fGrPolarizationPtSignifStatBkgCharm->GetX()[ipolarization];
		x_err = fGrPolarizationPtSignifStatBkgCharm->GetEX()[ipolarization];
		y_val = fGrPolarizationPtSignifStatBkgCharm->GetY()[ipolarization];
		y_err = fGrPolarizationPtSignifStatBkgCharm->GetEY()[ipolarization];
		
		if (!y_val) {
			cout << "y_val not defined! continue." << endl;
			continue;
		}
		
		if (y_val <= 0) {cout << "y_val <=0" << endl; y_val = 0.; }
		if (y_err <= 0) {cout << "y_err <=0" << endl; y_err = 1.; }
		y_val = TMath::Sqrt(NBinsAngularDistribution) / ( y_val / y_err ); // divide signal in NBinsAngularDistribution bins for polarization measurement, since this is relative error it is in the numerator (if I am not wrong)
		fGrPolarizationPtSignifStatBkgCharm->SetPoint(ipolarization, x_val, y_val);
		fGrPolarizationPtSignifStatBkgCharm->SetPointError(ipolarization, x_err, 0);
		
		y_val = fGrPolarizationPtSignifStatBkg->GetY()[ipolarization];
		y_err = fGrPolarizationPtSignifStatBkg->GetEY()[ipolarization];
		cout << " x_val, y_val, y_err = " << x_val << " , " << y_val << " , " << y_err << endl;
		cout << " sqrt(y_val*N_{evt}),N_{evt}*y_err_stat = "<< TMath::Sqrt(y_val*nEventsLHC) << " , " << nEventsLHC*y_err << endl;
		if (y_val <= 0) {cout << "y_val <=0" << endl; y_val = 0.; }
		if (y_err <= 0) {cout << "y_err <=0" << endl; y_err = 1.; }
		y_val = TMath::Sqrt(NBinsAngularDistribution) / ( y_val / y_err ); // divide signal in NBinsAngularDistribution bins for polarization measurement, since this is relative error it is in the numerator (if I am not wrong)
		fGrPolarizationPtSignifStatBkg->SetPoint(ipolarization, x_val, y_val);
		fGrPolarizationPtSignifStatBkg->SetPointError(ipolarization, x_err, 0);
		
		y_val = fGrPolarizationPtSignifStat->GetY()[ipolarization];
		y_err = fGrPolarizationPtSignifStat->GetEY()[ipolarization];
		if (y_val <= 0) {cout << "y_val <=0" << endl; y_val = 0.; }
		if (y_err <= 0) {cout << "y_err <=0" << endl; y_err = 1.; }
		y_val = TMath::Sqrt(NBinsAngularDistribution) / ( y_val / y_err ); // divide signal in NBinsAngularDistribution bins for polarization measurement, since this is relative error it is in the numerator (if I am not wrong)
		fGrPolarizationPtSignifStat->SetPoint(ipolarization, x_val, y_val);
		fGrPolarizationPtSignifStat->SetPointError(ipolarization, x_err, 0);
	}
	
	fH1D_axis_mee->SetMinimum(1e-3);
	fH1D_axis_mee->SetMaximum(10);
	fH1D_axis_mee->SetTitle(Form(";%s;%s", titleXaxis, "Avg. rel. error per bin in polarization variables"));
	fH1D_axis_mee->DrawCopy("axis");
	
	//fGrPolarizationPtSignifStatBkgCharm->DrawClone("p same");
	//fGrPolarizationPtSignifStatBkg->DrawClone("p same");
	fGrPolarizationPtSignifStat->DrawClone("p same");
	
	TLegend* fLegendPolarizationSignif = new TLegend(.55, .75, .88, .88, "");
	fLegendPolarizationSignif->SetTextSize(leg_textsize * 1.25);
	fLegendPolarizationSignif->SetBorderSize(0);
	fLegendPolarizationSignif->SetFillColor(0);
	fLegendPolarizationSignif->AddEntry(fGrPolarizationPtSignifStat, "Statistical", "lp");
	if(!doNoSyst) fLegendPolarizationSignif->AddEntry(fGrPolarizationPtSignifStatBkg, "Comb. Bkg.", "lp");
	if(!doNoSyst) fLegendPolarizationSignif->AddEntry(fGrPolarizationPtSignifStatBkgCharm, "Charm+Cocktail", "lp");
	fLegendPolarizationSignif->Draw();

	if(meebin==0)
	  fCanvPolarizationVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_Polarization_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
	else
	  fCanvPolarizationVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_Polarization_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
	
	cout << "Polarization done" << endl;



	// ----------------------------------------------------------------------
	// Flow
	cout << "Now Flow" << endl; // (if " << doPtee << " == 1)" << endl;
	Bool_t checkFlow = 0;
	
	// for Flow vs Ptee:
	const Int_t nflowPtBins = 3;
	Double_t flowPtBins[nflowPtBins+1] = {0.1, 0.8, 1.5, 3.};
		
	// for Flow vs Mee:
	const Int_t nflowMeeBins = 5;
	Double_t flowMeeBins[nflowMeeBins+1] = {0.15, 0.3, 0.7, 0.9, 1.5, 2.5};
	
	nbins = 0;
	if (doPtee) nbins = nflowPtBins;
	else        nbins = nflowMeeBins;
	
	const Int_t nflowBins = nbins;
	Double_t flowBins[nflowBins+1];
	for (Int_t iflow=0; iflow<nflowBins+1; iflow++)
	{
		if (doPtee) flowBins[iflow] = flowPtBins[iflow];
		else        flowBins[iflow] = flowMeeBins[iflow];
	}
	
	Double_t flowX[nflowBins];
	Double_t flowY[nflowBins];
	Double_t flowXerr[nflowBins];
	Double_t flowYerrStat[nflowBins];
	Double_t flowYerrBkg[nflowBins];
	Double_t flowYerrCharm[nflowBins];
	
	nPointsInBin = 0;
	nLostPoints  = 0;
	ySumPtBin = 0;
	ySignifStat2PtBin  = 0;
	ySignifBkg2PtBin   = 0;
	ySignifCharm2PtBin = 0;
	y_err_stat  = 0;
	y_err_bkg   = 0;
	y_err_charm = 0;
	
//		Int_t idebug = 0;
//		cout << "  debug " << idebug << endl; idebug++;
	
	for (Int_t iflow=0; iflow<nflowBins; iflow++)
	{
		nPointsInBin=0;
		nLostPoints=0;
		ySumPtBin=0;
		ySignifStat2PtBin=0;
		ySignifBkg2PtBin=0;
		ySignifCharm2PtBin=0;
		
		skipFirstPoint=0;
		for (Int_t ix=1+skipFirstPoint; ix<=hExcess->GetNbinsX(); ix++) 
		{
			x_val  = hExcess->GetBinCenter(ix);
			y_val  = hExcess->GetBinContent(ix);
			y_err_stat  = hExcess->GetBinError(ix);
			y_err_bkg   = fGrSystErr_Excess_CombBkg->GetEY()[ix-1];
			y_err_charm = fGrSystErr_Charm_Cocktail->GetEY()[ix-1];
			
			if (   (fGrSystErr_Excess_CombBkg->GetX()[ix-1] != x_val)
					|| (fGrSystErr_Charm_Cocktail->GetX()[ix-1] != x_val) )
			{
				cout << " ERROR in Flow: Graph from Sampling has other binning as histogram!" << endl;
				break;
			}
			
			if (x_val < flowBins[iflow]) continue;
			if (x_val > flowBins[iflow+1]) break;
//			cout << " x_val, y_val, y_err_stat = " << x_val << " , " << y_val << " , " << y_err_stat << endl;
//			cout << "   y_err_bkg, y_err_charm = " << y_err_bkg << " , " << y_err_charm << endl;
			
			nPointsInBin++;
			if (y_val <= 0) {
				nLostPoints++;
				continue;
			}
			
			ySumPtBin  += y_val;
			ySignifStat2PtBin  += TMath::Power(y_val/y_err_stat, 2);
			ySignifBkg2PtBin   += TMath::Power(y_val/y_err_bkg, 2);
			ySignifCharm2PtBin += TMath::Power(y_val/y_err_charm, 2);
		}
		
//		cout << " nPointsInBin, nLostPoints = " << nPointsInBin << " , " << nLostPoints << endl;
		if (nPointsInBin-nLostPoints) { // != 0
			ySumPtBin  *= (Double_t) nPointsInBin / (nPointsInBin-nLostPoints);
		}
		// this correction is an overestimate:
		Double_t err_corr_fac = 1.; // = (Double_t) nPointsInBin / (nPointsInBin-nLostPoints);
		ySignifStat2PtBin  *= err_corr_fac;
		ySignifBkg2PtBin   *= err_corr_fac;
		ySignifCharm2PtBin *= err_corr_fac;
		
		flowX[iflow] = (flowBins[iflow+1] - flowBins[iflow]) / 2 + flowBins[iflow];
		//			cout << "  before div" << endl;
		flowY[iflow] = ySumPtBin / nPointsInBin;
		flowXerr[iflow]	= (flowBins[iflow+1] - flowBins[iflow]) / 2;
		if (nPointsInBin-nLostPoints) { // != 0
			flowYerrStat[iflow]	 = flowY[iflow] / TMath::Sqrt(ySignifStat2PtBin);
			flowYerrBkg[iflow]	 = flowY[iflow] / TMath::Sqrt(ySignifBkg2PtBin); //  + flowYerrStat[iflow];
			flowYerrCharm[iflow] = flowY[iflow] / TMath::Sqrt(ySignifCharm2PtBin); // + flowYerrBkg[iflow];
		}
		else {
			flowYerrStat[iflow]	 = 0.;
			flowYerrBkg[iflow]	 = 0.;
			flowYerrCharm[iflow] = 0.;
		}
//			cout << flowX[iflow] << endl;
//			cout << flowY[iflow] << endl;
//			cout << flowXerr[iflow] << endl;
//			cout << flowYerrStat[iflow] << endl;
//			cout << flowYerrBkg[iflow] << endl;
//			cout << flowYerrCharm[iflow] << endl;
		
	}
	
	
	
	TGraphErrors* fGrFlowPtErrStat = new TGraphErrors(nflowBins, flowX, flowY, flowXerr, flowYerrStat);
	TGraphErrors* fGrFlowPtErrStatBkg = new TGraphErrors(nflowBins, flowX, flowY, flowXerr, flowYerrBkg);
	TGraphErrors* fGrFlowPtErrStatBkgCharm = new TGraphErrors(nflowBins, flowX, flowY, flowXerr, flowYerrCharm);
	
	fGrFlowPtErrStatBkgCharm->SetMarkerStyle(21);
	fGrFlowPtErrStatBkgCharm->SetMarkerColor(kMagenta);
	fGrFlowPtErrStatBkgCharm->SetLineColor(kMagenta);
	fGrFlowPtErrStatBkgCharm->SetLineWidth(4);
	fGrFlowPtErrStatBkg->SetMarkerStyle(21);
	fGrFlowPtErrStatBkg->SetMarkerColor(kGreen+1);
	fGrFlowPtErrStatBkg->SetLineColor(kGreen+1);
	fGrFlowPtErrStatBkg->SetLineWidth(3);
	fGrFlowPtErrStat->SetMarkerStyle(21);
	fGrFlowPtErrStat->SetMarkerColor(kBlack);
	fGrFlowPtErrStat->SetLineColor(kBlack);
	fGrFlowPtErrStat->SetLineWidth(2);
	if (checkFlow) {
		fCanvExcessVsMee->cd();
		//fGrFlowPtErrStatBkgCharm->DrawClone("p same");
		//fGrFlowPtErrStatBkg->DrawClone("p same");
		fGrFlowPtErrStat->DrawClone("p same");
	}
	
	fGrFlowPtSignifStatBkgCharm = (TGraphErrors*) fGrFlowPtErrStatBkgCharm->Clone("fGrFlowPtSignifStatBkgCharm");
	fGrFlowPtSignifStatBkg = (TGraphErrors*) fGrFlowPtErrStatBkg->Clone("fGrFlowPtSignifStatBkg");
	fGrFlowPtSignifStat = (TGraphErrors*) fGrFlowPtErrStat->Clone("fGrFlowPtSignifStat");
  
	// ----------------------------------------------------------------------
	// Plot Flow
	cout << " plotting Flow" << endl;

	Double_t canvYstretchFlow = 0.8;
	
	TCanvas* fCanvFlowVsMee = new TCanvas(Form("cFlowVs%s_%s_%dE%d_bin%d%s", MeeOrPtee, s_cent, nIntEventsLHC, nLogEventsLHC, meebin, canvname_opt), Form("Flow Background vs %s", MeeOrPtee), canvXshift+300, 0, canvWidth*canvYstretchFlow, canvHeight*canvYstretchFlow);
	fCanvFlowVsMee->SetLeftMargin(0.13);
	fCanvFlowVsMee->cd(1)->SetLogy();
	
	for (Int_t iflow=0; iflow<nflowBins; iflow++)
	{
		x_val = fGrFlowPtSignifStatBkgCharm->GetX()[iflow];
		x_err = fGrFlowPtSignifStatBkgCharm->GetEX()[iflow];
		y_val = fGrFlowPtSignifStatBkgCharm->GetY()[iflow];
		y_err = fGrFlowPtSignifStatBkgCharm->GetEY()[iflow];
//			cout << " x_val, y_val, y_err = " << x_val << " , " << y_val << " , " << y_err << endl;
		
		if (!y_val) {
			cout << "y_val not defined! continue." << endl;
			continue;
		}
		
		if (y_val <= 0) {cout << "y_val <=0" << endl; y_val = 0.; }
		if (y_err <= 0) {cout << "y_err <=0" << endl; y_err = 1.; }
		y_val = 0.7 / ( y_val / y_err );//0.7 = 1/sqrt(2) [dividing signal in two half = in- and out-of-plane]
		fGrFlowPtSignifStatBkgCharm->SetPoint(iflow, x_val, y_val);
		fGrFlowPtSignifStatBkgCharm->SetPointError(iflow, x_err, 0);
		
		y_val = fGrFlowPtSignifStatBkg->GetY()[iflow];
		y_err = fGrFlowPtSignifStatBkg->GetEY()[iflow];
		if (y_val <= 0) {cout << "y_val <=0" << endl; y_val = 0.; }
		if (y_err <= 0) {cout << "y_err <=0" << endl; y_err = 1.; }
		y_val = 0.7 / ( y_val / y_err );//0.7 = 1/sqrt(2) [dividing signal in two half = in- and out-of-plane]
		fGrFlowPtSignifStatBkg->SetPoint(iflow, x_val, y_val);
		fGrFlowPtSignifStatBkg->SetPointError(iflow, x_err, 0);
		
		y_val = fGrFlowPtSignifStat->GetY()[iflow];
		y_err = fGrFlowPtSignifStat->GetEY()[iflow];
		if (y_val <= 0) {cout << "y_val <=0" << endl; y_val = 0.; }
		if (y_err <= 0) {cout << "y_err <=0" << endl; y_err = 1.; }
		y_val = 0.7 / ( y_val / y_err );//0.7 = 1/sqrt(2) [dividing signal in two half = in- and out-of-plane]
		fGrFlowPtSignifStat->SetPoint(iflow, x_val, y_val);
		fGrFlowPtSignifStat->SetPointError(iflow, x_err, 0);
	}
	
	fH1D_axis_mee->SetMinimum(1e-3);
	fH1D_axis_mee->SetMaximum(10);
//	fH1D_axis_mee->SetTitle(Form(";%s;%s", titleXaxis, "dN / #sigma_{dN}"));
//	fH1D_axis_mee->SetTitle(Form(";%s;%s", titleXaxis, "#sigma_{v2} = 0.7 / (dN/#sigma_{dN})"));
	fH1D_axis_mee->SetTitle(Form(";%s;%s", titleXaxis, "#sigma_{v2}"));
	fH1D_axis_mee->DrawCopy("axis");
	
	//fGrFlowPtSignifStatBkgCharm->DrawClone("p same");
	//fGrFlowPtSignifStatBkg->DrawClone("p same");
	fGrFlowPtSignifStat->DrawClone("p same");
	
	TLegend* fLegendFlowSignif = new TLegend(.55, .75, .88, .88, "");
	fLegendFlowSignif->SetTextSize(leg_textsize * 1.25);
	fLegendFlowSignif->SetBorderSize(0);
	fLegendFlowSignif->SetFillColor(0);
	fLegendFlowSignif->AddEntry(fGrFlowPtSignifStat, "Statistical", "lp");
	if(!doNoSyst) fLegendFlowSignif->AddEntry(fGrFlowPtSignifStatBkg, "Comb. Bkg.", "lp");
	if(!doNoSyst) fLegendFlowSignif->AddEntry(fGrFlowPtSignifStatBkgCharm, "Charm+Cocktail", "lp");
	fLegendFlowSignif->Draw();

	if(meebin==0)
	  fCanvFlowVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_Flow_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
	else
	  fCanvFlowVsMee->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_Flow_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
	
	cout << "Flow done" << endl;

  }
	
	
  // ----------------------------------------------------------------------
  // QGP Temperature Fit
  Double_t canvYstretch = 0.8;
  TString str1 = "";

  x_val = 0;
  y_val = 0;
  Double_t y_err_stat = 0;
  Double_t y_err_bkg = 0;
  Double_t y_err_charm = 0;

  Int_t fi = 0;
  Double_t T_fit[7];
  Double_t T_norm[7]; 

  cout << "Now QGP T Fit"<< endl;
  if (!doPtee){

      // new histograms contain statistical errors from Poisson sample:
      TH1D* fH1DExcess_stat           = (TH1D*) hExcess->Clone("fH1DExcess_stat");
      TH1D* fH1DExcess_systCharmMin   = (TH1D*) hExcess->Clone("fH1DExcess_systCharmMin");
      TH1D* fH1DExcess_systCharmMax   = (TH1D*) hExcess->Clone("fH1DExcess_systCharmMax");
      TH1D* fH1DExcess_systCombBkgMin = (TH1D*) hExcess->Clone("fH1DExcess_systCombBkgMin");
      TH1D* fH1DExcess_systCombBkgMax = (TH1D*) hExcess->Clone("fH1DExcess_systCombBkgMax");      
      
      skipFirstPoint=0;
      for (Int_t ix=1+skipFirstPoint; ix<=hExcess->GetNbinsX(); ix++) 
	{
	  x_val  = hExcess->GetBinCenter(ix);
	  y_val  = hExcess->GetBinContent(ix);
	  y_err_stat  = hExcess->GetBinError(ix);
	  y_err_bkg   = fGrSystErr_Excess_CombBkg->GetEY()[ix-1];
	  y_err_charm = fGrSystErr_Charm_Cocktail->GetEY()[ix-1];
	  
	  if (   (fGrSystErr_Excess_CombBkg->GetX()[ix-1] != x_val)
		 || (fGrSystErr_Charm_Cocktail->GetX()[ix-1] != x_val) )
	    {
	      cout << " ERROR in T Fit: Graph from Sampling has other binning as histogram!" << endl;
	      break;
	    }
	  
	  fH1DExcess_systCharmMin->SetBinContent(ix, y_val - y_err_charm);
	  fH1DExcess_systCharmMax->SetBinContent(ix, y_val + y_err_charm);
	  fH1DExcess_systCombBkgMin->SetBinContent(ix, y_val - y_err_bkg);
	  fH1DExcess_systCombBkgMax->SetBinContent(ix, y_val + y_err_bkg);
	  
	  
	}
      
      TCanvas* fCanvQGPTVsMee = new TCanvas(Form("cQGPTVs%s_%s_%dE%d_bin%d%s", MeeOrPtee, s_cent, nIntEventsLHC, nLogEventsLHC, meebin, canvname_opt), Form("QGP T Fit vs %s", MeeOrPtee), canvXshift+100, 0, canvWidth*canvYstretch, canvHeight*canvYstretch);
      fCanvQGPTVsMee->Divide(3,2);
      
      Double_t fitmin = 1.1;//ITS TDR
      Double_t fitmax = 2.0;//ITS TDR
      
      TF1* fcnFitT = new TF1("fcnFitTRappQGP", "[1] * TMath::Power(x,3/2) * TMath::Exp(-x/[0])", fitmin, fitmax);
      fcnFitT->SetParameter(0, 0.500); // start fitting with 350 MeV
      fcnFitT->SetParLimits(0, 0.1, 2.0);
      fcnFitT->SetParameter(1, 1e-3);
     
      // real QGP slope:
      Float_t qgpplotmin=0.000001;
      Float_t qgpplotmax=100;
      Float_t qgpplotminX=0.0;
      Float_t qgpplotmaxX=2.5;
      
      fCanvQGPTVsMee->cd(1)->SetLogy(1);
      fH1D_Rapp_Sum->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fH1D_Rapp_Sum->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1D_Rapp_Sum->SetDrawOption("c"); // changes drawn style from hist to curved line.
      fH1D_Rapp_Sum->Fit(fcnFitT, "IQL","",fitmin,fitmax);
      fH1D_Rapp_Sum->Draw();
      hRhoRR->DrawCopy("chistsame");
      hQGP->DrawCopy("csame");
      
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      cout << "INPUT: T = " << T_fit[fi-1] << " GeV/c2 (Rapp QGP)" << endl;

      fCanvQGPTVsMee->cd(4)->SetLogy(1);
      fH1DExcess_stat->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fH1DExcess_stat->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_stat->Fit(fcnFitT, "RIQ");
      fH1DExcess_stat->Draw();
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      T_fit[fi] = fcnFitT->GetParError(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      cout << "FITRESULT (stat): T = " << T_fit[fi-2] << "+- " << T_fit[fi-1] <<" GeV/c2  --> " << T_fit[fi-1]/T_fit[fi-2]*100 << " % " << endl;

      fCanvQGPTVsMee->cd(2)->SetLogy(1);
      fGrSystErr_Charm_Cocktail->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fGrSystErr_Charm_Cocktail->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCharmMin->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fH1DExcess_systCharmMin->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCharmMin->Fit(fcnFitT, "RIQ");
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      fGrSystErr_Charm_Cocktail->Draw("2");
      fH1DExcess_systCharmMin->DrawCopy("same");

      fCanvQGPTVsMee->cd(5)->SetLogy(1);
      fGrSystErr_Charm_Cocktail->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fGrSystErr_Charm_Cocktail->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCharmMax->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fH1DExcess_systCharmMax->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCharmMax->Fit(fcnFitT, "RIQ");
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      cout << "FITRESULT (syst, cc): T = " << T_fit[fi-2] << " - " << T_fit[fi-1] <<" GeV/c2 " << endl;
      fGrSystErr_Charm_Cocktail->Draw("2");
      fH1DExcess_systCharmMax->DrawCopy("same");

      fCanvQGPTVsMee->cd(3)->SetLogy(1);
      fGrSystErr_Excess_CombBkg->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fGrSystErr_Excess_CombBkg->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCombBkgMin->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fH1DExcess_systCombBkgMin->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCombBkgMin->Fit(fcnFitT, "RIQ");
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      fGrSystErr_Excess_CombBkg->Draw("2");
      fH1DExcess_systCombBkgMin->DrawCopy("same");

      fCanvQGPTVsMee->cd(6)->SetLogy(1);
      fGrSystErr_Excess_CombBkg->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fGrSystErr_Excess_CombBkg->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCombBkgMax->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fH1DExcess_systCombBkgMax->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_systCombBkgMax->Fit(fcnFitT, "RIQ");
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      cout << "FITRESULT (syst, BG): T = " << T_fit[fi-2] << " - " << T_fit[fi-1] <<" GeV/c2 " << endl;
      fGrSystErr_Excess_CombBkg->Draw("2");
      fH1DExcess_systCombBkgMax->DrawCopy("same");
  }

   else{

      // new histograms contain statistical errors from Poisson sample:
      TH1D* fH1DExcess_stat           = (TH1D*) hExcess->Clone("fH1DExcess_stat");  
      
      TCanvas* fCanvQGPTVsPt = new TCanvas(Form("cQGPTVs%s_%s_%dE%d_bin%d%s", MeeOrPtee, s_cent, nIntEventsLHC, nLogEventsLHC, meebin, canvname_opt), Form("QGP T Fit vs %s", MeeOrPtee), canvXshift+100, 0, canvWidth*canvYstretch, canvHeight*canvYstretch);
      Double_t fitmin = 0.9;
      Double_t fitmax = 1.5;
      
      // simple fit
      TF1* fcnFitT = new TF1("fcnFitT", "exp(-x/[0])*[1]", fitmin, fitmax);
      fcnFitT->SetParameter(0, 0.500); 
      fcnFitT->SetParLimits(0, 0.1, 2.0);
      fcnFitT->SetParameter(1, 1e-3); 
        
      // real QGP slope:
      Float_t qgpplotmin=0.0000001;
      Float_t qgpplotmax=1.0;
      Float_t qgpplotminX=0.0;
      Float_t qgpplotmaxX=1.5;
      
      fCanvQGPTVsPt->cd()->SetLogy(1);
      fH1DExcess_stat->GetXaxis()->SetRangeUser(qgpplotminX,qgpplotmaxX);
      fH1DExcess_stat->GetYaxis()->SetRangeUser(qgpplotmin,qgpplotmax);
      fH1DExcess_stat->Fit(fcnFitT, "RIQ");
      fH1DExcess_stat->Draw();
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      T_fit[fi] = fcnFitT->GetParameter(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      T_fit[fi] = fcnFitT->GetParError(0); T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      cout << "FITRESULT (stat): T = " << T_fit[fi-2] << "+- " << T_fit[fi-1] <<" GeV/c2  --> " << T_fit[fi-1]/T_fit[fi-2]*100 << " % " << endl;

      T_fit[fi] = 0.; T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      T_fit[fi] = 0.; T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      T_fit[fi] = 0.; T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;
      T_fit[fi] = 0.; T_norm[fi] = T_fit[fi]/T_fit[0]; fi++;

   }

  TCanvas *cTemp = new TCanvas("cTemp","cTemp",500,500);
    
  TGraphErrors* fGrErr_T_Stat        = new TGraphErrors( 1 );
  TGraphErrors* fGrSystErr_T_CombBkg = new TGraphErrors( 1 );
  TGraphErrors* fGrSystErr_T_Charm   = new TGraphErrors( 1 );

  fGrErr_T_Stat->SetPoint(0, 1, T_norm[1]);
  fGrErr_T_Stat->SetPointError(0, 0, T_norm[2]);
  fGrSystErr_T_Charm->SetPoint(0, 1, (T_norm[3]+T_norm[4])/2.);
  fGrSystErr_T_Charm->SetPointError(0, 0.1,TMath::Abs(T_norm[3]-T_norm[4])/2.);
  fGrSystErr_T_CombBkg->SetPoint(0, 1, (T_norm[5]+T_norm[6])/2.);
  fGrSystErr_T_CombBkg->SetPointError(0, 0.1,TMath::Abs(T_norm[5]-T_norm[6])/2.);

  fGrErr_T_Stat->GetYaxis()->SetTitle("T_{fit}/T_{real}");
  fGrErr_T_Stat->SetTitle("Temperature fit");
  fGrErr_T_Stat->GetYaxis()->SetRangeUser(0.0,2.0);
  fGrErr_T_Stat->SetLineColor(kBlack);
  fGrErr_T_Stat->SetMarkerColor(kBlack);
  fGrErr_T_Stat->SetMarkerStyle(25);
  fGrErr_T_Stat->Draw("AP");
  fGrSystErr_T_CombBkg->SetFillColor(kGreen+1);
  fGrSystErr_T_CombBkg->SetFillStyle(3244);
  fGrSystErr_T_CombBkg->Draw("2same");
  fGrSystErr_T_Charm->SetFillColor(kMagenta);
  fGrSystErr_T_Charm->SetFillStyle(3105);  //(3001); not viewable in Illustrator
  fGrSystErr_T_Charm->Draw("2same");
  fGrErr_T_Stat->Draw("Psame");


  if(meebin==0)
    cTemp->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_Temperature_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));
  else
    cTemp->SaveAs(Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_Temperature_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.png",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling));

  
  // ----------------------------------------------------------------------
  // Write output file
  if(writeFile){

    TString fileOutName = Form("./finalPlotsLowB_Systems_preliminary_pT%d_%s_histos_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.root",doPtee,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling);
    if(meebin!=0)
      fileOutName = Form("./finalPlotsLowB_Systems_preliminary_pT%d_meebin%d_%s_histos_Run5_%d_RICH_%d_IPcut%d_events%.0f_Conv%.1f.root",doPtee,meebin,sSystem[iSystem].Data(),Run5,RICH,doIPcut,nEventsLHC,conversionBGScaling);
    Printf("Write to file %s",fileOutName.Data()); 
    
    TFile *fOut = TFile::Open(fileOutName,"RECREATE");
    fOut->cd();

     
    hSignalAll->Write();
    hSignalStatError->Write();
    hExcess->Write();
    hBG->Write();
    hSBQA->Write();
    hSignificance->Write();
    hSignalwoRho->Write();
    hRhoRR->Write();
    hQGP->Write();
    hCharm->Write();

    fGrSystErr_Signal_CombBkg->SetName("fGrSystErr_Signal_CombBkg");
    fGrSystErr_Excess_CombBkg->SetName("fGrSystErr_Excess_CombBkg");
    fGrSystErr_Charm_Cocktail->SetName("fGrSystErr_Charm_Cocktail");

    fGrSystErr_Signal_CombBkg->Write();
    fGrSystErr_Excess_CombBkg->Write();
    fGrSystErr_Charm_Cocktail->Write();
    
    fGrErr_T_Stat->SetName("fGrErr_T_Stat");
    fGrSystErr_T_CombBkg->SetName("fGrSystErr_T_CombBkg");
    fGrSystErr_T_Charm->SetName("fGrSystErr_T_Charm");
  
    fGrErr_T_Stat->Write();
    fGrSystErr_T_CombBkg->Write();
    fGrSystErr_T_Charm->Write();

    fH1D_Rapp_Sum->Write();
    fH1D_Rapp_Sum_vac->Write();
    fH1D_Rapp_Sum_drop->Write();

    if(!doPtee){
      fGrFlowPtSignifStat->Write();
    }
    
    fOut->Close();
  }
}



//____________________________________________________________//
TH1D* convertCocktailHistogram(TH1D *cocktailHisto, Bool_t doPtee) {

  cout<<"Cocktail histo "<<cocktailHisto<<" with "<<cocktailHisto->GetNbinsX()<<"bins from "<<cocktailHisto->GetXaxis()->GetBinCenter(1)<<" to "<<cocktailHisto->GetXaxis()->GetBinCenter(cocktailHisto->GetNbinsX())<<" with bin width "<<cocktailHisto->GetXaxis()->GetBinWidth(1)<<endl;


  TH1D* newCocktailHisto = NULL;
  if(!doPtee)
    newCocktailHisto = new TH1D(Form("%sNew",cocktailHisto->GetName()),Form("%s",cocktailHisto->GetTitle()),160,0,4);
  else
    newCocktailHisto = new TH1D(Form("%sNew",cocktailHisto->GetName()),Form("%s",cocktailHisto->GetTitle()),160,0,8);
  
  for(Int_t i = 0; i < 160; i++){
    if(i<120){
      if(doPtee && i>=60){
	newCocktailHisto->SetBinContent(1+i,cocktailHisto->GetBinContent(60));
	newCocktailHisto->SetBinError(1+i,0.);
      }
      else{
	newCocktailHisto->SetBinContent(1+i,cocktailHisto->GetBinContent(1+i));
	newCocktailHisto->SetBinError(1+i,cocktailHisto->GetBinError(1+i));
      }
    }
    else{
      newCocktailHisto->SetBinContent(1+i,0.);
      newCocktailHisto->SetBinError(1+i,0.);
    }
  }
  
  return newCocktailHisto;

}

//____________________________________________________________//
TH1D* convertRappHistogram(TH1D *rappHisto, Bool_t doPtee) {

  TString histName = rappHisto->GetName();
    
  cout<<"Rapp histo "<<histName<<" with "<<rappHisto->GetNbinsX()<<" bins from "<<rappHisto->GetXaxis()->GetBinCenter(1)<<" to "<<rappHisto->GetXaxis()->GetBinCenter(rappHisto->GetNbinsX())<<" with bin width "<<rappHisto->GetXaxis()->GetBinWidth(1)<<endl;

  // extrapolation to 2.5
  Double_t fitmin = 0.9;
  Double_t fitmax = 1.5;
  TF1* fcnFitT = new TF1("fcnFitT", "exp(-x/[0])*[1]", fitmin, fitmax);
  fcnFitT->SetParameter(0, 0.500); // start fitting with 500 MeV
  fcnFitT->SetParLimits(0, 0.1, 2.0);
  fcnFitT->SetParameter(1, 1e-3);

  TF1* fcnFitTRappQGP = new TF1("fcnFitTRappQGP", "[0] * TMath::Power(x,3/2) * TMath::Exp(-x/[1])", fitmin, fitmax);
  fcnFitTRappQGP->SetParameter(1, 0.500); // start fitting with 500 MeV
  fcnFitTRappQGP->SetParLimits(1, 0.1, 2.0);
  fcnFitTRappQGP->SetParameter(0, 1e-3);

  TF1* fcnFitTRappHad = new TF1("fcnFitTRappHad", "[0] * TMath::Power(x,3/4) * TMath::Exp(-x/[1])", fitmin, fitmax);
  fcnFitTRappHad->SetParameter(1, 0.500); // start fitting with 500 MeV
  fcnFitTRappHad->SetParLimits(1, 0.1, 2.0);
  fcnFitTRappHad->SetParameter(0, 1e-3);
  
  rappHisto->Fit(fcnFitT, "IL","0",fitmin,fitmax);
  rappHisto->Fit(fcnFitTRappQGP, "IL","0",fitmin,fitmax);
  rappHisto->Fit(fcnFitTRappHad, "IL","0",fitmin,fitmax);

  TH1D* newRappHisto = NULL;
  if(!doPtee)
    newRappHisto = new TH1D(Form("%sNew",rappHisto->GetName()),Form("%s",rappHisto->GetTitle()),160,0,4);
  else
    newRappHisto = new TH1D(Form("%sNew",rappHisto->GetName()),Form("%s",rappHisto->GetTitle()),160,0,8);

  for(Int_t i = 0; i < 160; i++){
    
    if(i<60){
      newRappHisto->SetBinContent(1+i,rappHisto->GetBinContent(1+i));
      newRappHisto->SetBinError(1+i,0.);
     }
    else{
      if(!doPtee){
	// as presented by Patrick in LMee PAG (23.03.2017)
	if(histName.Contains("fHistRho"))
	  newRappHisto->SetBinContent(1+i,fcnFitTRappHad->Eval(newRappHisto->GetBinCenter(1+i)));
	else
	  newRappHisto->SetBinContent(1+i,fcnFitTRappQGP->Eval(newRappHisto->GetBinCenter(1+i)));
	// simplified method as presented by MW in ITS plenary session
	//newRappHisto->SetBinContent(1+i,fcnFitT->Eval(newRappHisto->GetBinCenter(1+i))); 
	newRappHisto->SetBinError(1+i,0);
      }
      else{
	newRappHisto->SetBinContent(1+i,rappHisto->GetBinContent(60));
	newRappHisto->SetBinError(1+i,0.);
      }
    }
  }
  
  return newRappHisto;
  
}


//____________________________________________________________//
TH1D* convertCharmHistogram(TH1D *charmHisto, Bool_t doPtee) {

  cout<<"Charm histo "<<charmHisto<<" with "<<charmHisto->GetNbinsX()<<"bins from "<<charmHisto->GetXaxis()->GetBinCenter(1)<<" to "<<charmHisto->GetXaxis()->GetBinCenter(charmHisto->GetNbinsX())<<" with bin width "<<charmHisto->GetXaxis()->GetBinWidth(1)<<endl;

  TH1D* newCharmHisto = NULL;
  if(!doPtee)
    newCharmHisto = new TH1D(Form("%sNew",charmHisto->GetName()),Form("%s",charmHisto->GetTitle()),160,0,4);
  else
    newCharmHisto = new TH1D(Form("%sNew",charmHisto->GetName()),Form("%s",charmHisto->GetTitle()),160,0,8);

  
  for(Int_t i = 0; i < 160; i++){
    newCharmHisto->SetBinContent(1+i,charmHisto->GetBinContent(1+i));
    newCharmHisto->SetBinError(1+i,charmHisto->GetBinError(1+i));
  }

  return newCharmHisto;

}
