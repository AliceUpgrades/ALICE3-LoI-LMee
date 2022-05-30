void DrawEfficiency01(){ 
  TFile *file = new TFile("findable.root", "READ");
  TH1D *fGenerated = (TH1D*) file->Get("hGenerated");
    
  TTree *fTree = (TTree*) file->Get("fTree");
  
  TH1D *fAnalysis = new TH1D("fAnalysis", "", 200,0,20);
  TH1D *fReco = new TH1D("fReco", "", 200,0,20);
  TH1D *fReco1 = new TH1D("fReco1", "", 200,0,20);
  TH1D *fReco6 = new TH1D("fReco6", "", 200,0,20);
  TH1D *fReco6cont = new TH1D("fReco6cont", "", 200,0,20);
  TH1D *fReco6tracked = new TH1D("fReco6tracked", "", 200,0,20);
  TH1D *fReco6analysis = new TH1D("fReco6analysis", "", 200,0,20);
  
  Float_t fPtMC, fEtaMC;
  Int_t fBachHits, fNegHits, fPosHits;
  
  fTree->SetBranchAddress("fPtMC", &fPtMC);
  fTree->SetBranchAddress("fEtaMC", &fEtaMC);
  
  fTree->SetBranchAddress("fBachHits", &fBachHits);
  fTree->SetBranchAddress("fPosHits", &fPosHits);
  fTree->SetBranchAddress("fNegHits", &fNegHits);
  
  Bool_t fBachTracked, fPosTracked, fNegTracked;
  
  fTree->SetBranchAddress("fBachTracked", &fBachTracked);
  fTree->SetBranchAddress("fPosTracked", &fPosTracked);
  fTree->SetBranchAddress("fNegTracked", &fNegTracked);
  
  Bool_t fCascadeReco;
  
  fTree->SetBranchAddress("fCascadeReco", &fCascadeReco);
  
  Bool_t fBachHitLayer[12];
  Bool_t fPosHitLayer[12];
  Bool_t fNegHitLayer[12];
  
  for(Int_t ih=0; ih<12; ih++) fTree->SetBranchAddress( Form("fBachHits%i", ih), &(fBachHitLayer[ih]));
  for(Int_t ih=0; ih<12; ih++) fTree->SetBranchAddress( Form("fPosHits%i", ih), &(fPosHitLayer[ih]));
  for(Int_t ih=0; ih<12; ih++) fTree->SetBranchAddress( Form("fNegHits%i", ih), &(fNegHitLayer[ih]));
  
  
  for(Long_t ii=0; ii<fTree->GetEntries(); ii++){
    fTree->GetEntry(ii);
    fReco->Fill(fPtMC);
    if(fBachHits>=1&&fPosHits>=1&&fNegHits>=1) fReco1->Fill(fPtMC);
    if(fBachHits>=6&&fPosHits>=6&&fNegHits>=6) fReco6->Fill(fPtMC);
    
    if(fBachTracked&&fPosTracked&&fNegTracked) fReco6tracked->Fill(fPtMC);
    if(fCascadeReco) fReco6analysis->Fill(fPtMC);
    
    //check contiguous
    Bool_t lBach6=kFALSE;
    Bool_t lPos6=kFALSE;
    Bool_t lNeg6=kFALSE;
    for(Int_t ih=0; ih<7; ih++){
      Int_t lCondBach=1;
      Int_t lCondPos=1;
      Int_t lCondNeg=1;
      for(Int_t ic=0; ic<6; ic++){
        lCondBach *= fBachHitLayer[ih+ic];
        lCondPos *= fPosHitLayer[ih+ic];
        lCondNeg *= fNegHitLayer[ih+ic];
      }
      if(lCondBach>0) lBach6=kTRUE;
      if(lCondPos>0) lPos6=kTRUE;
      if(lCondNeg>0) lNeg6=kTRUE;
    }
    
    if( lBach6 && lPos6 && lNeg6 ) fReco6cont->Fill(fPtMC);
  }
  
  
  fReco->Sumw2();
  fReco1->Sumw2();
  fReco6->Sumw2();
  fReco6cont->Sumw2();
  fReco6tracked->Sumw2();
  fReco6analysis->Sumw2();
  fGenerated->Sumw2();
  
  fReco->Divide(fGenerated);
  fReco1->Divide(fGenerated);
  fReco6->Divide(fGenerated);
  fReco6cont->Divide(fGenerated);
  fReco6tracked->Divide(fGenerated);
  fReco6analysis->Divide(fGenerated);
  
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->SetTicks(1,1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.14);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.04);
  
  fReco->GetYaxis()->SetTitle("Acceptance #times efficiency");
  fReco->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fReco->GetXaxis()->SetTitleSize(0.055);
  fReco->GetYaxis()->SetTitleSize(0.055);
  fReco->GetXaxis()->SetLabelSize(0.042);
  fReco->GetYaxis()->SetLabelSize(0.042);
  fReco->GetXaxis()->SetRangeUser(0,4);
  fReco->GetYaxis()->SetRangeUser(-0.25,0.8);
  
  fReco->SetMarkerStyle(24);
  fReco->SetMarkerColor(kBlack);
  fReco->SetLineColor(kBlack);
  
  fReco1->SetMarkerStyle(20);
  fReco1->SetMarkerSize(0.6);
  fReco1->SetMarkerColor(kBlack);
  fReco1->SetLineColor(kBlack);
  
  fReco->Draw();
  fReco6->SetMarkerStyle(20);
  fReco6->SetMarkerColor(kRed);
  fReco6->SetLineColor(kRed);
  
  fReco6cont->SetMarkerStyle(20);
  fReco6cont->SetMarkerColor(kBlue);
  fReco6cont->SetLineColor(kBlue);
  
  fReco6tracked->SetMarkerStyle(47);
  fReco6tracked->SetMarkerSize(1.3);
  fReco6tracked->SetMarkerColor(kPink+3);
  fReco6tracked->SetLineColor(kPink+3);
  
  fReco6analysis->SetMarkerStyle(20);
  fReco6analysis->SetMarkerSize(0.5);
  fReco6analysis->SetMarkerColor(kGreen+1);
  fReco6analysis->SetLineColor(kGreen+1);
  
  fReco1->Draw("same");
  fReco6->Draw("same");
  fReco6cont->Draw("same");
  fReco6tracked->Draw("same");
  fReco6analysis->Draw("same");
  
  TLegend *leg = new TLegend (0.184, 0.154783, 0.534, 0.307391);
  leg->SetBorderSize(0);
  leg->AddEntry(fReco, "All in decay channel", "lp");
  leg->AddEntry(fReco1, "1 hit minimum", "lp");
  leg->AddEntry(fReco6, "6 hits minimum", "lp");
  leg->AddEntry(fReco6cont, "6 hits no holes", "lp");
  leg->AddEntry(fReco6tracked, "6-hit tracker", "lp");
  leg->AddEntry(fReco6analysis, "6-hit candidate", "pl");
  leg->Draw();
  
} 
