void setupHisto(TH1 *currentGraph=0,
		int myMarkerStyle=8,
		int myMarkerColor=1,
		float myMarkerSize=1,
		int myLineStyle=1,
		int myLineColor=1,
		float myLineWidth=1,
		int myFillStyle =1001,
		int myFillColor =1);

void setupGraph(TGraph *currentGraph=0,
		int myMarkerStyle=8,
		int myMarkerColor=1,
		float myMarkerSize=1,
		int myLineStyle=1,
		int myLineColor=1,
		float myLineWidth=1,
		int myFillStyle =1001,
		int myFillColor =1);
  
void plotFinalPlotsComparison(Int_t iStudy = 0){

  gStyle->SetOptStat(0);

  Int_t col[]      = {kBlack,kRed,kOrange-3,kGreen+3,kOrange-1,kOrange+7,kGreen-3,kMagenta+3,kAzure+7,kCyan+2,kPink-8,kYellow+2,kBlue+2,kMagenta-6};

  const Int_t nFileMax = 5;
  Int_t nFile = 0;
  TString sFile[nFileMax];
  TString sFile2[nFileMax];
  for(Int_t iFile = 0; iFile < nFileMax; iFile++){
    sFile[iFile] = "";
    sFile2[iFile] = "";
  }

  
  if(iStudy ==  0){

    nFile = 4;
    sFile[0] = "./finalPlotsLowB_Systems_preliminary_pT0_Pb_histos_Run5_1_RICH_1_IPcut0_events2500000000_Conv1.0.root";
    sFile[1] = "./finalPlotsLowB_Systems_preliminary_pT0_Xe_histos_Run5_1_RICH_1_IPcut0_events14175000000_Conv1.0.root";
    sFile[2] = "./finalPlotsLowB_Systems_preliminary_pT0_Kr_histos_Run5_1_RICH_1_IPcut0_events57256410256_Conv1.0.root";
    sFile[3] = "./finalPlotsLowB_Systems_preliminary_pT0_Ar_histos_Run5_1_RICH_1_IPcut0_events538333333333_Conv1.0.root";
            
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

    sFile2[0] = "0% Conversions";
    sFile2[1] = "25% ";
    sFile2[2] = "50%";
    sFile2[3] = "75%";
    sFile2[4] = "100%";

  }
  else
    return;

 
  const Int_t nType  = 6;
  const Int_t nType2 = 3;
  TString sType[nType] = {
    "hSBQA", 
    "hSignificance",
    "hSignalAll",
    "hBG",
    "hSignalwoRho",
    "hpt_mass_charm_px0New"
  };
  TString sType2[nType] = {
    "S/B",
    "significance",
    "signal",
    "BG",
    "cocktail",
    "c#bar{c}"
  };
  TString sType3[nType2] = {
    "fGrSystErr_Signal_CombBkg", 
    "fGrSystErr_Excess_CombBkg",
    "fGrSystErr_Charm_Cocktail"
  };
  TString sType4[nType2] = {
    "Syst. uncert. BG (cocktail+charm+signal)",
    "Syst. uncert. BG (signal)",
    "Syst. uncert. Charm + Cocktail (signal)"
  };

  TFile *fFile[nFileMax];
  TH1D* fHist[nType][nFileMax];
  TH1D* fHistRatio[nType][nFileMax];
  TGraphErrors* fGSyst[nType2][nFileMax];
  TGraphErrors* fGSystRelative[nType2][nFileMax];
  
  TCanvas *c = new TCanvas("c","c",1200,900);
  c->Divide(3,2);
  TCanvas *cRatio = new TCanvas("cRatio","cRatio",1200,900);
  cRatio->Divide(3,2);
  TCanvas *d = new TCanvas("d","d",1200,500);
  d->Divide(3,1);
  
  TLegend *l = new TLegend(0.7,0.6,0.85,0.85,"","brNDC");
  l->SetLineColor(0);
  TLegend *l2 = new TLegend(0.7,0.6,0.85,0.85,"","brNDC");
  l2->SetLineColor(0);
  
  for(Int_t iFile = 0; iFile < nFile; iFile++){
    
    fFile[iFile] = TFile::Open(sFile[iFile].Data(),"READ");
    if(!fFile[iFile]){
      cerr<<"File "<<sFile[iFile].Data()<<" not found"<<endl;
      return;
    }
    //fFile[iFile]->ls();
    for(Int_t iType = 0; iType < nType; iType++){
      fHist[iType][iFile] = (TH1D*)fFile[iFile]->Get(Form("%s",sType[iType].Data()));
      if(!fHist[iType][iFile]){
       	cerr<<"Histogram "<<sType[iType].Data()<<" not found in file "<<sFile[iFile].Data()<<endl;
       	return;
      }
      fHist[iType][iFile]->GetXaxis()->SetRangeUser(0,1.5);
      setupHisto(fHist[iType][iFile],20,col[iFile],1.0,1,col[iFile]);
      fHistRatio[iType][iFile] = (TH1D*)fHist[iType][iFile]->Clone(Form("%s_ratio",fHist[iType][iFile]->GetName()));
      fHistRatio[iType][iFile]->Divide(fHist[iType][0]);
      fHist[iType][iFile]->SetTitle(Form("%s",sType2[iType].Data()));
      fHistRatio[iType][iFile]->SetTitle(Form("%s ratio",sType2[iType].Data()));
    }
    for(Int_t iType = 0; iType < nType2; iType++){

      fGSyst[iType][iFile] = (TGraphErrors*)fFile[iFile]->Get(Form("%s",sType3[iType].Data()));
      if(!fGSyst[iType][iFile]){
    	cerr<<"TGraphErrors "<<sType3[iType].Data()<<" not found in file "<<sFile[iFile].Data()<<endl;
    	return;
      }
      fGSyst[iType][iFile]->GetXaxis()->SetRangeUser(0,1.5);
      setupGraph(fGSyst[iType][iFile],20,col[iFile],1.0,1,col[iFile]);
      fGSyst[iType][iFile]->SetTitle(Form("%s",sType4[iType].Data()));
      fGSystRelative[iType][iFile] = (TGraphErrors*)fGSyst[iType][iFile]->Clone(Form("%s_relative",fGSyst[iType][iFile]->GetName()));
      for(Int_t iBin = 0; iBin < fGSystRelative[iType][iFile]->GetN(); iBin++){
	Double_t x, y, xE, yE;
	fGSyst[iType][iFile]->GetPoint(iBin,x,y);
	xE = fGSyst[iType][iFile]->GetErrorX(iBin);
	yE = fGSyst[iType][iFile]->GetErrorY(iBin);
	fGSystRelative[iType][iFile]->SetPoint(iBin,x,yE/y);
	fGSystRelative[iType][iFile]->SetPointError(iBin,xE,0);
	
      }
    }
  }


  for(Int_t iFile = 0; iFile < nFile; iFile++){
    for(Int_t iType = 0; iType < nType; iType++){

      
      if(iType==0)
  	l->AddEntry(fHist[iType][iFile],sFile2[iFile].Data(),"lp");

      
      c->cd(iType+1)->SetLogy();
      if(iFile==0){
  	fHist[iType][iFile]->Draw("P");
      }
      else{
  	fHist[iType][iFile]->Draw("Psame");
      }
      if(iFile==nFile-1){
  	l->Draw();
      }

      cRatio->cd(iType+1);//->SetLogy();
      if(iFile==0){
  	fHistRatio[iType][iFile]->GetYaxis()->SetRangeUser(0.0,2.0);
  	fHistRatio[iType][iFile]->Draw("P");
      }
      else{
  	fHistRatio[iType][iFile]->Draw("Psame");
      }
      if(iFile==nFile-1){
  	l->Draw();
      }
    }
  }


  for(Int_t iFile = 0; iFile < nFile; iFile++){
    for(Int_t iType = 0; iType < nType2; iType++){

      
      if(iType==0)
  	l2->AddEntry(fGSyst[iType][iFile],sFile2[iFile].Data(),"lp");
      
      d->cd(iType+1)->SetLogy();
      if(iFile==0){
	fGSystRelative[iType][iFile]->GetYaxis()->SetRangeUser(1e-2,1e2);
  	fGSystRelative[iType][iFile]->Draw("AP");
      }
      else{
  	fGSystRelative[iType][iFile]->Draw("P");
      }
      if(iFile==nFile-1){
  	l2->Draw();
      }
    }
  }

  c->SaveAs(Form("./finalPlotsComparison_study%d.png",iStudy));
  cRatio->SaveAs(Form("./finalPlotsComparison_ratios_study%d.png",iStudy));
  d->SaveAs(Form("./finalPlotsComparison_syst_study%d.png",iStudy));

}

//____________________________________________________________//
void setupHisto(TH1 *currentGraph=0,
		int myMarkerStyle=8,
		int myMarkerColor=1,
		float myMarkerSize=1,
		int myLineStyle=1,
		int myLineColor=1,
		float myLineWidth=1,
		int myFillStyle =1001,
		int myFillColor =1) {
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
//   currentGraph->Set();
}
//____________________________________________________________//
void setupGraph(TGraph *currentGraph=0,
		int myMarkerStyle=8,
		int myMarkerColor=1,
		float myMarkerSize=1,
		int myLineStyle=1,
		int myLineColor=1,
		float myLineWidth=1,
		int myFillStyle =1001,
		int myFillColor =1) {
  currentGraph->SetMarkerStyle(myMarkerStyle);
  currentGraph->SetMarkerColor(myMarkerColor);
  currentGraph->SetMarkerSize(myMarkerSize);
  currentGraph->SetLineColor(myLineColor);
  currentGraph->SetLineStyle(myLineStyle);
  currentGraph->SetLineWidth(myLineWidth);
  currentGraph->SetFillStyle(myFillStyle);
  currentGraph->SetFillColor(myFillColor);
//   currentGraph->Set();
}


