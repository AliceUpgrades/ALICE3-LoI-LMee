void preparePreshowerEff(){

  // Create histograms that are then used by Delphes analysis:
  // These should be replaced later with the input
  // from the dedicated Preshower studies
  auto preshowerEffElectrons = new TH1F("preshowerEffElectrons","Preshower efficiency for electrons; p (GeV/c); efficiency",200,0,10);
  auto preshowerEffPions = new TH1F("preshowerEffPions","Preshower efficiency for pions; p (GeV/c); efficiency",200,0,10);

  // Simple efficiency onset parametrizations
  auto f_cauchy1 = new TF1("f_cauchy", "[2]*(TMath::ATan((x - [0]) / [1]) / TMath::Pi() + 0.5)+[3]", 0., 10.);
  auto f_cauchy2 = new TF1("f_cauchy", "[2]*(0.5-TMath::ATan((x - [0]) / [1]) / TMath::Pi())+[3]", 0., 10.);

  // draw the parametrizations
  TCanvas* c = new TCanvas("c","c",900,500);
  c->Divide(2,1);
  
  // Simple electron parametrization
  f_cauchy1->SetParameter(0,0.5);
  f_cauchy1->SetParameter(1,0.1);
  f_cauchy1->SetParameter(2,0.9);
  f_cauchy1->SetParameter(3,0.);
  f_cauchy1->SetTitle("Simple electron parametrization");
  c->cd(1);
  f_cauchy1->Draw("C");

  // Simple pion parametrization
  f_cauchy2->SetParameter(0,0.5);
  f_cauchy2->SetParameter(1,0.1);
  f_cauchy2->SetParameter(2,0.005);
  f_cauchy2->SetParameter(3,0.0005);
  f_cauchy2->SetTitle("Simple pion parametrization");
  c->cd(2);
  f_cauchy2->Draw("C");

  // Fill histograms that are then used by Delphes analysis
  for(Int_t i = 0; i < 200; i++){
    preshowerEffElectrons->Fill(i*0.05,f_cauchy1->Eval((i*0.05+0.025)));
    preshowerEffPions->Fill(i*0.05,f_cauchy2->Eval((i*0.05+0.025)));
  }

  // check histograms
  c->cd(1);
  preshowerEffElectrons->Draw("same hist");
  c->cd(2);
  preshowerEffPions->Draw("same hist");

  auto fOut = TFile::Open("preshowerEff.root","RECREATE");
  preshowerEffElectrons->Write();
  preshowerEffPions->Write();
  fOut->Close();
  
}
