void plotFinalPlotsTemperatureComparison(Bool_t bAsymmetric = kFALSE, Int_t iStudy = 0){

  gROOT->SetStyle("Plain");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  Int_t col[]      = {kBlack,kRed,kOrange-3,kGreen+3,kOrange-1,kOrange+7,kGreen-3,kMagenta+3,kAzure+7,kCyan+2,kPink-8,kYellow+2,kBlue+2,kMagenta-6};

  const Int_t nFileMax = 5;
  Int_t nFile = 0;
  TString sFile[nFileMax];
  TString sFile2[nFileMax];
  TString sFileSyst[nFileMax];
  TString sFileSyst2[nFileMax];
  for(Int_t iFile = 0; iFile < nFileMax; iFile++){
    sFile[iFile] = "";
    sFile2[iFile] = "";
    sFileSyst[iFile] = "";
    sFileSyst2[iFile] = "";
  }

  TFile *fFile[nFileMax];
  TFile *fFileSyst[nFileMax];
  TFile *fFileSyst2[nFileMax];
  TGraphErrors* fHistStat;
  TGraphErrors* fHistSystCharm;
  TGraphErrors* fHistSystBG;
  TGraphErrors* fHistSystStat;
  TGraphErrors* fHistNormStat;
  TGraphErrors* fHistNormSystCharm;
  TGraphErrors* fHistNormSystBG;
  TGraphErrors* fHistNormSystStat;

  TGraphErrors* fGrErr_T_Stat;
  TGraphErrors* fGrSystErr_T_CombBkg;
  TGraphErrors* fGrSystErr_T_Charm;
  TGraphErrors* fGrErr_T_SystStat;

  Double_t x[nFileMax];
  Double_t xE[nFileMax];

  Double_t yStat[nFileMax];
  Double_t ySystCharm[nFileMax];
  Double_t ySystBG[nFileMax];
  Double_t ySystStat[nFileMax];

  Double_t yEStat[nFileMax];
  Double_t yESystCharm[nFileMax];
  Double_t yESystBG[nFileMax];
  Double_t yESystStat[nFileMax];

  Double_t yNormStat[nFileMax];
  Double_t yNormSystCharm[nFileMax];
  Double_t yNormSystBG[nFileMax];
  Double_t yNormSystStat[nFileMax];

  Double_t yNormEStat[nFileMax];
  Double_t yNormESystCharm[nFileMax];
  Double_t yNormESystBG[nFileMax];
  Double_t yNormESystStat[nFileMax];

  Double_t     nEventsLHC = 2500000000.;
  

  if(iStudy ==  0){

    nFile = 4;
    
    sFile[0]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.0.root";
    sFile[1]     = "./finalPlotsLowB_Systems_preliminary_pT0_Xe_histos_Run5_1_RICH_1_IPcut0_events14175000000_Conv1.0.root";
    sFile[2]     = "./finalPlotsLowB_Systems_preliminary_pT0_Kr_histos_Run5_1_RICH_1_IPcut0_events57256410256_Conv1.0.root";
    sFile[3]     = "./finalPlotsLowB_Systems_preliminary_pT0_Ar_histos_Run5_1_RICH_1_IPcut0_events538333333333_Conv1.0.root";
            
    sFileSyst[0]  = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.0.root";
    sFileSyst[1]  = "./finalPlotsLowB_Systems_preliminary_pT0_Xe_histos_Run5_1_RICH_1_IPcut0_events14175000000_Conv1.0.root";
    sFileSyst[2]  = "./finalPlotsLowB_Systems_preliminary_pT0_Kr_histos_Run5_1_RICH_1_IPcut0_events57256410256_Conv1.0.root";
    sFileSyst[3]  = "./finalPlotsLowB_Systems_preliminary_pT0_Ar_histos_Run5_1_RICH_1_IPcut0_events538333333333_Conv1.0.root";

    sFile2[0] = "Pb-Pb";
    sFile2[1] = "Xe-Xe";
    sFile2[2] = "Kr-Kr";
    sFile2[3] = "Ar-Ar";
    
  }


  else if(iStudy ==  1){

    nFile = 5;

    sFile[0]     = "~/MAC_201909/ALICE/dileptons/Upgrade/ITS/analysis/out/lowB_HF/finalPlots/finalPlotsLowB_ITS3_EoI_updateRFactorSyst_histos_ITSCyl0_IPcut1_events2500000000.root";
    sFile[1]     = "~/MAC_201909/ALICE/dileptons/Upgrade/ITS/analysis/out/lowB_HF/finalPlots/finalPlotsLowB_ITS3_EoI_updateRFactorSyst_histos_ITSCyl1_IPcut1_events2500000000.root";
    sFile[2]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_0_RICH_0_IPcut0_events2500000000_Conv1.0.root";
    sFile[3]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_0_IPcut0_events2500000000_Conv1.0.root";
    sFile[4]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.0.root";

            
    sFileSyst[0]     = "~/MAC_201909/ALICE/dileptons/Upgrade/ITS/analysis/out/lowB_HF/finalPlots/finalPlotsLowB_ITS3_EoI_updateRFactorSyst_histos_ITSCyl0_IPcut1_events2500000000000.root";
    sFileSyst[1]     = "~/MAC_201909/ALICE/dileptons/Upgrade/ITS/analysis/out/lowB_HF/finalPlots/finalPlotsLowB_ITS3_EoI_updateRFactorSyst_histos_ITSCyl1_IPcut1_events2500000000000.root";
    sFileSyst[2]  = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_0_RICH_0_IPcut0_events2500000000_Conv1.0.root";
    sFileSyst[3]  = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_0_IPcut0_events2500000000_Conv1.0.root";
    sFileSyst[4]  = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.0.root";

    sFile2[0] = "ITS3 (LoI)";
    sFile2[1] = "ITS3 (LoI)";
    sFile2[2] = "ITS3 (improved)";
    sFile2[3] = "ALICE 3";
    sFile2[4] = "ALICE 3 (+RICH)";

  }
  else if(iStudy ==  2){

    nFile = 5;
  
    sFile[0]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.0.root";
    sFile[1]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.2.root";
    sFile[2]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.5.root";
    sFile[3]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.8.root";
    sFile[4]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv2.0.root";

    sFileSyst[0]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.0.root";
    sFileSyst[1]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.2.root";
    sFileSyst[2]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.5.root";
    sFileSyst[3]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.8.root";
    sFileSyst[4]     = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv2.0.root";
   

    sFile2[0] = "0% Conversions";
    sFile2[1] = "25% ";
    sFile2[2] = "50%";
    sFile2[3] = "75%";
    sFile2[4] = "100%";

  }
  

  else{
    Printf("Study = %d not known",iStudy);
    return;
  }


  Int_t    nLogEventsLHC = (Int_t) TMath::Log10(nEventsLHC); // (Int_t) not needed but ok here
  Int_t    nIntEventsLHC = nEventsLHC / TMath::Power(10, nLogEventsLHC); // do not use (Int_t) here..!
  Double_t nDouEventsLHC = nEventsLHC / TMath::Power(10, nLogEventsLHC);

    
  for(Int_t iFile = 0; iFile < nFile; iFile++){

    if(iFile==0){
      fHistStat      = new TGraphErrors(nFile);
      fHistSystCharm = new TGraphErrors(nFile);
      fHistSystBG    = new TGraphErrors(nFile);
      fHistSystStat  = new TGraphErrors(nFile);
      fHistNormStat      = new TGraphErrors(nFile);
      fHistNormSystCharm = new TGraphErrors(nFile);
      fHistNormSystBG    = new TGraphErrors(nFile);
      fHistNormSystStat  = new TGraphErrors(nFile);
    }
    
    fFile[iFile] = TFile::Open(sFile[iFile].Data(),"READ");
    if(!fFile[iFile]){
      cerr<<"File "<<sFile[iFile].Data()<<" not found"<<endl;
      return;
    }
    //fFile[iFile]->ls();

    fFileSyst[iFile] = TFile::Open(sFileSyst[iFile].Data(),"READ");
    if(!fFileSyst[iFile]){
      cerr<<"File (syst) "<<sFileSyst[iFile].Data()<<" not found"<<endl;
      return;
    }
    //fFileSyst[iFile]->ls();

    if(sFileSyst2[iFile]!=""){
      fFileSyst2[iFile] = TFile::Open(sFileSyst2[iFile].Data(),"READ");
      if(!fFileSyst2[iFile]){
	cerr<<"File (syst) "<<sFileSyst2[iFile].Data()<<" not found"<<endl;
	return;
      }
      //fFileSyst2[iFile]->ls();
    }
    
    fGrErr_T_Stat        = (TGraphErrors*)fFile[iFile]->Get("fGrErr_T_Stat");
    if(sFileSyst2[iFile]!=""){
      fGrSystErr_T_CombBkg = (TGraphErrors*)fFileSyst2[iFile]->Get("fGrSystErr_T_CombBkg");
    }else{
      fGrSystErr_T_CombBkg = (TGraphErrors*)fFileSyst[iFile]->Get("fGrSystErr_T_CombBkg");
    }
    fGrSystErr_T_Charm   = (TGraphErrors*)fFileSyst[iFile]->Get("fGrSystErr_T_Charm");
    fGrErr_T_SystStat    = (TGraphErrors*)fFileSyst[iFile]->Get("fGrErr_T_Stat");

    fGrErr_T_Stat->GetPoint(0,x[iFile],yStat[iFile]);
    fGrSystErr_T_Charm->GetPoint(0,x[iFile],ySystCharm[iFile]);
    fGrSystErr_T_CombBkg->GetPoint(0,x[iFile],ySystBG[iFile]);
    fGrErr_T_SystStat->GetPoint(0,x[iFile],ySystStat[iFile]);
    
    yEStat[iFile]      = fGrErr_T_Stat->GetErrorY(0);
    yESystCharm[iFile] = fGrSystErr_T_Charm->GetErrorY(0);
    yESystBG[iFile]    = fGrSystErr_T_CombBkg->GetErrorY(0);
    yESystStat[iFile]  = fGrErr_T_SystStat->GetErrorY(0);

    if(bAsymmetric){
      yNormStat[iFile]       = yStat[iFile]/yStat[iFile];
      yNormEStat[iFile]      = yEStat[iFile]/yStat[iFile];
      yNormSystCharm[iFile]  = ySystCharm[iFile]/ySystStat[iFile];
      yNormESystCharm[iFile] = yESystCharm[iFile]/ySystStat[iFile];
      yNormSystBG[iFile]     = ySystBG[iFile]/ySystStat[iFile];
      yNormESystBG[iFile]    = yESystBG[iFile]/ySystStat[iFile];
      yNormSystStat[iFile]   = ySystStat[iFile]/ySystStat[iFile];
      yNormESystStat[iFile]  = yESystStat[iFile]/ySystStat[iFile];
    }
    else{
      yNormStat[iFile]       = yStat[iFile]/yStat[iFile];
      yNormEStat[iFile]      = yEStat[iFile]/yStat[iFile];
      yNormSystCharm[iFile]  = ySystCharm[iFile]/ySystCharm[iFile];
      yNormESystCharm[iFile] = yESystCharm[iFile]/ySystCharm[iFile];
      yNormSystBG[iFile]     = ySystBG[iFile]/ySystBG[iFile];
      yNormESystBG[iFile]    = yESystBG[iFile]/ySystBG[iFile];
      yNormSystStat[iFile]   = ySystStat[iFile]/ySystStat[iFile];
      yNormESystStat[iFile]  = yESystStat[iFile]/ySystStat[iFile];
    }
    

    x[iFile] = iFile + 0.5;
    
    // cout<<x[iFile]<<" "<<yStat[iFile]<<" "<<ySystCharm[iFile]<<" "<<ySystBG[iFile]<<" "<<ySystStat[iFile]<<endl;
    // cout<<" --> "<<yEStat[iFile]<<" "<<yESystCharm[iFile]<<" "<<yESystBG[iFile]<<" "<<yESystStat[iFile]<<endl<<endl;
    // cout<<x[iFile]<<" "<<yNormStat[iFile]<<" "<<yNormSystCharm[iFile]<<" "<<yNormSystBG[iFile]<<" "<<yNormSystStat[iFile]<<endl;
    // cout<<" --> "<<yNormEStat[iFile]<<" "<<yNormESystCharm[iFile]<<" "<<yNormESystBG[iFile]<<" "<<yNormESystStat[iFile]<<endl<<endl;
    cout<<" Stat. "<<yNormEStat[iFile]*100<<"% --- ratio =  "<<yNormEStat[iFile]/yNormEStat[0]*100<<"% ==> "<<yNormEStat[0]/yNormEStat[iFile]<<endl;
    // cout<<" Syst. (Charm) "<<(1-(yNormSystCharm[iFile]-yNormESystCharm[iFile]))*100<<" - "<<(yNormSystCharm[iFile]+yNormESystCharm[iFile]-1)*100<<"% --- ratio =  "<<(1-(yNormSystCharm[iFile]-yNormESystCharm[iFile]))/(1-(yNormSystCharm[0]-yNormESystCharm[0]))*100<<" - "<<(yNormSystCharm[iFile]+yNormESystCharm[iFile]-1)/(yNormSystCharm[0]+yNormESystCharm[0]-1)*100<<"% ==> "<<(yNormSystCharm[0]+yNormESystCharm[0]-1)/(yNormSystCharm[iFile]+yNormESystCharm[iFile]-1)<<endl;
    // cout<<" Syst. (BG)    "<<(1-(yNormSystBG[iFile]-yNormESystBG[iFile]))*100<<" - "<<(yNormSystBG[iFile]+yNormESystBG[iFile]-1)*100<<"% --- ratio =  "<<(1-(yNormSystBG[iFile]-yNormESystBG[iFile]))/(1-(yNormSystBG[0]-yNormESystBG[0]))*100<<" - "<<(yNormSystBG[iFile]+yNormESystBG[iFile]-1)/(yNormSystBG[0]+yNormESystBG[0]-1)*100<<"% ==> "<<(yNormSystBG[0]+yNormESystBG[0]-1)/(yNormSystBG[iFile]+yNormESystBG[iFile]-1)<<endl<<endl<<endl;

    fHistStat->SetPoint(iFile,x[iFile],yStat[iFile]);
    fHistStat->SetPointError(iFile,0,yEStat[iFile]);

    fHistSystBG->SetPoint(iFile,x[iFile],ySystBG[iFile]);
    fHistSystBG->SetPointError(iFile,0.1,yESystBG[iFile]);

    fHistSystCharm->SetPoint(iFile,x[iFile],ySystCharm[iFile]);
    fHistSystCharm->SetPointError(iFile,0.07,yESystCharm[iFile]);
    
    fHistSystStat->SetPoint(iFile,x[iFile]+0.05,ySystStat[iFile]);
    fHistSystStat->SetPointError(iFile,0,yESystStat[iFile]);

    fHistNormStat->SetPoint(iFile,x[iFile],yNormStat[iFile]);
    fHistNormStat->SetPointError(iFile,0,yNormEStat[iFile]);

    fHistNormSystBG->SetPoint(iFile,x[iFile],yNormSystBG[iFile]);
    fHistNormSystBG->SetPointError(iFile,0.1,yNormESystBG[iFile]);

    fHistNormSystCharm->SetPoint(iFile,x[iFile],yNormSystCharm[iFile]);
    fHistNormSystCharm->SetPointError(iFile,0.07,yNormESystCharm[iFile]);
    
    fHistNormSystStat->SetPoint(iFile,x[iFile]+0.05,yNormSystStat[iFile]);
    fHistNormSystStat->SetPointError(iFile,0,yNormESystStat[iFile]);
    
  }


  TCanvas *c = new TCanvas("c","c",900,900);
  TCanvas *cNorm = new TCanvas("cNorm","cNorm",900,900);
  
  c->cd();
  
  TH1F* hPad = new TH1F("hPad","",nFile,0,nFile);
  hPad->GetYaxis()->SetRangeUser(0.9,1.2);
  for(Int_t iFile = 0; iFile < nFile; iFile++){
    hPad->GetXaxis()->SetBinLabel(iFile+1,sFile2[iFile]);
  }
  hPad->GetXaxis()->SetLabelSize(0.05);
  hPad->GetYaxis()->SetTitle("#it{T}_{fit}/#it{T}_{real}");
  hPad->GetYaxis()->SetTitleOffset(1.2);
  hPad->SetTitle("Temperature fit");
  hPad->Draw();

  fHistSystBG->SetFillColor(kGreen+1);
  fHistSystBG->SetFillStyle(3244);
  fHistSystBG->Draw("2same");
  fHistSystCharm->SetFillColor(kMagenta);
  fHistSystCharm->SetFillStyle(3105);  //(3001); not viewable in Illustrator
  fHistSystCharm->Draw("2same");
  fHistStat->SetLineColor(kBlack);
  fHistStat->SetMarkerColor(kBlack);
  fHistStat->SetMarkerStyle(21);
  fHistStat->Draw("Psame");
  // fHistSystStat->SetLineColor(kGray+2);
  // fHistSystStat->SetMarkerColor(kGray+2);
  // fHistSystStat->SetMarkerStyle(20);
  // fHistSystStat->Draw("Psame");

  Double_t leg_textsize = 0.029;

  TLatex *   logo = new TLatex(0.4,0.84,"ALICE Upgrade Simulation");
  logo->SetNDC();
  logo->SetTextSize(0.04);
  //logo->SetTextFont(42);
  logo->Draw();


  // Legend
  TLegend *leg2=new TLegend(0.46,0.59,0.88,0.82,Form(""),"NDC");
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(leg_textsize); 
  leg2->AddEntry(fHistStat,Form("#it{B} = 0.2 T"), "ple");
  leg2->AddEntry(fHistStat,Form("1.1 < #it{M}_{ee,fit} (GeV/#it{c}^{2}) < 2.0"), "");
  leg2->AddEntry(fHistSystBG,Form("Syst. err. sig. + bkg."), "f");
  leg2->AddEntry(fHistSystCharm,Form("Syst. err. c#bar{c} + cocktail"), "f");
  leg2->Draw();

  gPad->SetTicks(0,1);

  TLine *lineOne = new TLine(0.,1.,2,1.);
  lineOne->SetLineWidth(2);
  lineOne->SetLineStyle(2);
  lineOne->SetLineColor(1);
  lineOne->Draw();
  
  cNorm->cd();
  
  TH1F* hPadNorm = new TH1F("hPadNorm","",nFile,0,nFile);
  hPadNorm->GetYaxis()->SetRangeUser(0.9,1.2);
  for(Int_t iFile = 0; iFile < nFile; iFile++){
    hPadNorm->GetXaxis()->SetBinLabel(iFile+1,sFile2[iFile]);
  }
  hPadNorm->GetXaxis()->SetLabelSize(0.05);
  hPadNorm->GetYaxis()->SetTitleOffset(1.2);
  hPadNorm->GetYaxis()->SetTitle("#it{T}_{fit}/#it{T}_{fit}");
  hPadNorm->SetTitle("Temperature fit");
  hPadNorm->Draw();
  
  fHistNormSystBG->SetFillColor(kGreen+1);
  fHistNormSystBG->SetFillStyle(3244);
  fHistNormSystBG->Draw("2same");
  fHistNormSystCharm->SetFillColor(kMagenta);
  fHistNormSystCharm->SetFillStyle(3105);  //(3001); not viewable in Illustrator
  fHistNormSystCharm->Draw("2same");
  fHistNormStat->SetLineColor(kBlack);
  fHistNormStat->SetMarkerColor(kBlack);
  fHistNormStat->SetMarkerStyle(21);
  fHistNormStat->Draw("Psame");
  // fHistNormSystStat->SetLineColor(kGray+2);
  // fHistNormSystStat->SetMarkerColor(kGray+2);
  // fHistNormSystStat->SetMarkerStyle(20);
  // fHistNormSystStat->Draw("Psame");

  logo->Draw();
  leg2->Draw();
  
  c->SaveAs(Form("./finalPlotsTemperatureComparison_study%d.png",iStudy));
  cNorm->SaveAs(Form("./finalPlotsTemperatureComparison_norm_study%d.png",iStudy));
}
