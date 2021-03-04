#include "genericTOF.h"

//===========================================================================
// Simple macro to extract electron efficiency
// (from Roberto Preghnella, modified by Michael Weber)
//
// in ROOT:
// > .L pidEfficiency.C
// > drawElectron()
// > drawElectronOptimized()
//===========================================================================


genericTOF TOF;
TF1 fres("fres", "TMath::Gaus(x, [0], [1], true)");

void setup()
{
  TOF.setName("TOF");
  TOF.setType(genericTOF::kBarrel);
  TOF.setSigma(0.020);       // [ns]
  TOF.setLength(200.);       // [cm]
  TOF.setRadius(100.);       // [cm]
  TOF.setMagneticField(0.2); // [T]  
}

double
pidEfficiency(double p, double eta, double imass, double jmass, double nsigma) 
{
  if (!TOF.isHit(eta, p)) return 0.;

  auto length = TOF.trackLength(eta);
  auto iflight = TOF.timeOfFlight(length, imass, p);
  auto jflight = TOF.timeOfFlight(length, jmass, p);
  auto sigma = TOF.getSigma();
  
  fres.SetParameter(0, jflight);
  fres.SetParameter(1, sigma);
  
  return fres.Integral(iflight - nsigma * sigma, iflight + nsigma * sigma);

}

void
drawElectron()
{
  setup();
  
  auto c = new TCanvas("c", "c", 800, 800);
  c->DrawFrame(0., 0., 10., 1.1, "electron ID (3#sigma cut);#it{p} (GeV/#it{c});efficiency");
  auto fel = new TF1("fel", "pidEfficiency(x, 0., 0.00051099891, 0.00051099891, 3.)", 0., 10.);
  auto fmu = new TF1("fmu", "pidEfficiency(x, 0., 0.00051099891, 0.10565800,    3.)", 0., 10.);
  auto fpi = new TF1("fpi", "pidEfficiency(x, 0., 0.00051099891, 0.13957000,    3.)", 0., 10.);
  auto fka = new TF1("fka", "pidEfficiency(x, 0., 0.00051099891, 0.49367700,    3.)", 0., 10.);
  auto fpr = new TF1("fpr", "pidEfficiency(x, 0., 0.00051099891, 0.93827200,    3.)", 0., 10.);
  fel->SetLineColor(1);
  fmu->SetLineColor(2);
  fpi->SetLineColor(4);
  fka->SetLineColor(6);
  fpr->SetLineColor(8);
  fel->Draw("same");
  fmu->Draw("same");
  fpi->Draw("same");
  fka->Draw("same");
  fpr->Draw("same");

  TLegend *l1 = new TLegend(0.65,0.15,0.85,0.45,"","brNDC");
  l1->SetLineColor(0);
  l1->SetFillColor(0);
  l1->AddEntry(fel,"Electrons","lp");
  l1->AddEntry(fmu,"Muons","lp");
  l1->AddEntry(fpi,"Pions","lp");
  l1->AddEntry(fka,"Kaons","lp");
  l1->AddEntry(fpr,"Protons","lp");
  l1->Draw();
  
}

void
drawPion()
{
  setup();
  
  auto c = new TCanvas("c", "c", 800, 800);
  c->DrawFrame(0., 0., 10., 1.1, "electron ID (3#sigma cut);#it{p} (GeV/#it{c});efficiency");
  auto fel = new TF1("fel", "pidEfficiency(x, 0.,  0.13957000, 0.00051099891, 3.)", 0., 10.);
  auto fmu = new TF1("fmu", "pidEfficiency(x, 0.,  0.13957000, 0.10565800,    3.)", 0., 10.);
  auto fpi = new TF1("fpi", "pidEfficiency(x, 0.,  0.13957000, 0.13957000,    3.)", 0., 10.);
  auto fka = new TF1("fka", "pidEfficiency(x, 0.,  0.13957000, 0.49367700,    3.)", 0., 10.);
  auto fpr = new TF1("fpr", "pidEfficiency(x, 0.,  0.13957000, 0.93827200,    3.)", 0., 10.);
  fel->SetLineColor(1);
  fmu->SetLineColor(2);
  fpi->SetLineColor(4);
  fka->SetLineColor(6);
  fpr->SetLineColor(8);
  fel->Draw("same");
  fmu->Draw("same");
  fpi->Draw("same");
  fka->Draw("same");
  fpr->Draw("same");

  TLegend *l1 = new TLegend(0.65,0.15,0.85,0.45,"","brNDC");
  l1->SetLineColor(0);
  l1->SetFillColor(0);
  l1->AddEntry(fel,"Electrons","lp");
  l1->AddEntry(fmu,"Muons","lp");
  l1->AddEntry(fpi,"Pions","lp");
  l1->AddEntry(fka,"Kaons","lp");
  l1->AddEntry(fpr,"Protons","lp");
  l1->Draw();
  
  
}

void
drawElectronOptimized(double pionEff = 0.001)
{
  setup();
  
  auto c = new TCanvas("c", "c", 800, 800);
  c->DrawFrame(0., 0., 10., 1.1, Form("electron ID (%.1f%% pion efficiency);#it{p} (GeV/#it{c});efficiency",100.*pionEff));

  auto hNSigmaThr = new TH1F("hNSigmaThr","",1000,0.,10.);
  
  TF1* fpi[5000];
  TF1* fel[5000];

  double nsigma = 0.;
  double mom    = 0.;
  double eff    = 0.;
  bool filled   = false;
  
  for(Int_t i = 0; i < 5000; i++){

    nsigma = (1.+(double)i)*0.001;
    fpi[i] = new TF1("fpi", "pidEfficiency(x, 0., 0.00051099891, 0.13957000,[0])", 0., 10.);
    fel[i] = new TF1("fel", "pidEfficiency(x, 0., 0.00051099891, 0.00051099891,[0])", 0., 10.);
    fpi[i]->SetParameter(0,nsigma);
    fel[i]->SetParameter(0,nsigma);

  }

  for(Int_t j = 0; j < 1000; j++){

    mom = (double)j*0.010;
    filled = false;
    for(Int_t i = 0; i < 5000; i++){
      eff = fpi[i]->Eval(mom);
      if (eff>=pionEff){
	hNSigmaThr->SetBinContent(j+1,fel[i-1]->Eval(mom));
	filled = true;
	break;
      }
    }
    if(filled==false)
      hNSigmaThr->SetBinContent(j+1,1.);
  }
  hNSigmaThr->SetLineWidth(2);
  hNSigmaThr->SetLineColor(2);
  hNSigmaThr->Draw("same");

  TLegend *l1 = new TLegend(0.65,0.75,0.85,0.85,"","brNDC");
  l1->SetLineColor(0);
  l1->SetFillColor(0);
  l1->AddEntry(hNSigmaThr,"Electrons","lp");
  l1->Draw();
  

  auto fOut = TFile::Open("electronEff.root","RECREATE");
  hNSigmaThr->Write();
  fOut->Close();

  
}
