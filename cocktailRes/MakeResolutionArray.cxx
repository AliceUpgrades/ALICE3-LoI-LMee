#include "TFile.h"
#include "TH2.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"

#include <iostream>

TFile *fIn(0x0), *fOut(0x0);
TString fNameIn, fNameOut;

std::vector<Double_t> fP_bins_vec;
std::vector<Double_t> fOpAngle_bins_vec;

//_______________________________________________________________________________________________
void Init_ppAt5TeV() {
 
  fNameIn   = "/mnt/Work/ALICE3/Input/cocktailInput.rmin100.2kG.root";
  fNameOut  = "resolution_test_2kG.root";
 
  std::cout << fNameIn << " " << fNameOut << std::endl;



  // Preliminary binning for the maps
  Double_t p_bin[] = {0.000000,
    0.03,
    0.04,
    0.05,
    0.06,
    0.07,
    0.08,
    0.09,
    0.1,
    0.11,
    0.12,
    0.13,
    0.14,
    0.15,
    0.16,
    0.17,
    0.18,
    0.19,
    0.2,
    0.21,
    0.22,
    0.23,
    0.24,
    0.25,
    0.26,
    0.27,
    0.28,
    0.29,
    0.3,
    0.320000,
    0.340000,
    0.360000,
    0.380000,
    0.400000,
    0.420000,
    0.440000,
    0.460000,
    0.500000,
    0.540000,
    0.580000,
    0.620000,
    0.660000,
    0.700000,
    0.760000,
    0.820000,
    0.880000,
    0.940000,
    1.000000,
    1.100000,
    1.200000,
    1.300000,
    1.400000,
    1.600000,
    1.800000,
    2.,
    2.5,
    3.,
    4.};
  
  
  Int_t n_p_bin = sizeof(p_bin) / sizeof(*p_bin) - 1;
  for (int i = 0; i < n_p_bin + 1; ++i)
    fP_bins_vec.push_back(p_bin[i]);

  Double_t opAng_bin[] = {0.,   0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
                          0.08, 0.10, 0.12, 0.16, 0.20, 0.30, 0.40,
                          0.6,  0.8,  1.0,  1.2,  1.4,  1.6,  1.8,
                          2.0,  2.2,  2.4,  2.6,  2.8,  3.0,  3.2};
  Int_t n_opAng_bin = sizeof(opAng_bin) / sizeof(*opAng_bin) - 1;
  for (int i = 0; i < n_opAng_bin + 1; ++i)
    fOpAngle_bins_vec.push_back(opAng_bin[i]);
}
// taken from plotting framework
//_______________________________________________________________________________________________
TH2D *Rebin2DHistogram(TH2F &hIn, Int_t n_bins_x, Double_t *bins_x, Int_t n_bins_y, Double_t *bins_y) {

  TString sname(hIn.GetName()); // Otherwise potential memory leak
  TH2D *hOut = new TH2D("hOut", hIn.GetTitle(), n_bins_x, bins_x, n_bins_y, bins_y);
  hOut->GetXaxis()->SetTitle(hIn.GetXaxis()->GetTitle());
  hOut->GetYaxis()->SetTitle(hIn.GetYaxis()->GetTitle());
  hOut->GetZaxis()->SetTitle(hIn.GetZaxis()->GetTitle());
  Double_t bin_content_temp = 0;
  Double_t bin_center_mee = -1;
  Double_t bin_center_ptee = -1;
  for (Int_t i_mee = 0; i_mee <= hIn.GetNbinsX() + 1; ++i_mee) { // Overflow?
    for (Int_t j_ptee = 0; j_ptee <= hIn.GetNbinsY() + 1; ++j_ptee) {
      bin_content_temp = hIn.GetBinContent(i_mee, j_ptee);
      bin_center_mee = hIn.GetXaxis()->GetBinCenter(i_mee);
      bin_center_ptee = hIn.GetYaxis()->GetBinCenter(j_ptee);
      hOut->Fill(bin_center_mee, bin_center_ptee, bin_content_temp);
    }
  }
  // Setting the bin error manually to avoid conflicts between root5 and root6
  for (Int_t i_mee = 0; i_mee <= hOut->GetNbinsX() + 1; ++i_mee) { // Overflow?
    for (Int_t j_ptee = 0; j_ptee <= hOut->GetNbinsY() + 1; ++j_ptee) {
      Double_t bin_content = hOut->GetBinContent(i_mee, j_ptee);
      hOut->SetBinError(i_mee, j_ptee, TMath::Sqrt(bin_content));
    }
  }

 
  hOut->SetName(sname.Data());
  return hOut;
}
//_______________________________________________________________________________________________
void MakeArray2Sum(TString hName, TString arrName, TString arrNamee, std::vector<Double_t> x_bins_vec) {

  TH2F *hRes2D_or = static_cast<TH2F *>(fIn->Get(Form("%s_pos",hName.Data())));
  hRes2D_or->SetName(arrNamee.Data());
  TH2F *hRes2D_bis = static_cast<TH2F *>(fIn->Get(Form("%s_ele",hName.Data())));
  hRes2D_or->Add(hRes2D_bis);

  //
  // Rebinning:
  //
  // for x, take the binning from the vector
  Int_t n_x_bins = x_bins_vec.size() - 1;
  // Double_t* x_bins = x_bins_vec.data(); // c++11
  Double_t *x_bins = &(x_bins_vec[0]);

  // for y, keep the original binning
  Int_t n_y_bins = hRes2D_or->GetYaxis()->GetNbins();
  Double_t y_bins[10000];
  hRes2D_or->GetYaxis()->GetLowEdge(y_bins); // Fills the binning of the original histogram into the array.
  //    cout << " y_bins[1] = " << y_bins[1] << " y_bins[99] = " << y_bins[99]
  //    << " y_bins[100] = " << y_bins[100] << endl;
  y_bins[n_y_bins] = hRes2D_or->GetYaxis()->GetBinLowEdge(n_y_bins + 1); // fill last bin edge manually, because sometimes this is
  // wrong with GetLowEdge()...

  TH2D *hRes2D = Rebin2DHistogram(*hRes2D_or, n_x_bins, x_bins, n_y_bins, y_bins);

  //
  // Slices:
  //
  TObjArray *resArr = new TObjArray();
  resArr->SetName(arrName.Data());
  resArr->SetOwner(kTRUE);
  resArr->Add(hRes2D);

  std::cout << "fill resolution array: " << resArr->GetName() << std::flush;
  Int_t nbinsx = hRes2D->GetXaxis()->GetNbins();
  for (Int_t bini = 1; bini <= nbinsx; bini++) {
    TString newname = Form("bin%d", bini);
    TH1D *hSlice = new TH1D(
        *((TH1D *)hRes2D->ProjectionY(newname.Data(), bini, bini, "e")));
    hSlice->SetMarkerStyle(20);
    resArr->Add(hSlice);
  }
  std::cout << "   ... done. Number of entries: " << resArr->GetLast()
            << std::endl;

  fOut->cd();
  resArr->Write(resArr->GetName(), TObject::kSingleKey);
}

//_______________________________________________________________________________________________
void MakeArray2(TString hName, TString arrName, TString arrNamee, std::vector<Double_t> x_bins_vec) {

  TH2F *hRes2D_or = static_cast<TH2F *>(fIn->Get(hName.Data()));
  hRes2D_or->SetName(arrNamee.Data());

  //
  // Rebinning:
  //
  // for x, take the binning from the vector
  Int_t n_x_bins = x_bins_vec.size() - 1;
  // Double_t* x_bins = x_bins_vec.data(); // c++11
  Double_t *x_bins = &(x_bins_vec[0]);

  // for y, keep the original binning
  Int_t n_y_bins = hRes2D_or->GetYaxis()->GetNbins();
  Double_t y_bins[10000];
  hRes2D_or->GetYaxis()->GetLowEdge(y_bins); // Fills the binning of the original histogram into the array.
  //    cout << " y_bins[1] = " << y_bins[1] << " y_bins[99] = " << y_bins[99]
  //    << " y_bins[100] = " << y_bins[100] << endl;
  y_bins[n_y_bins] = hRes2D_or->GetYaxis()->GetBinLowEdge(n_y_bins + 1); // fill last bin edge manually, because sometimes this is
  // wrong with GetLowEdge()...

  TH2D *hRes2D = Rebin2DHistogram(*hRes2D_or, n_x_bins, x_bins, n_y_bins, y_bins);

  //
  // Slices:
  //
  TObjArray *resArr = new TObjArray();
  resArr->SetName(arrName.Data());
  resArr->SetOwner(kTRUE);
  resArr->Add(hRes2D);

  std::cout << "fill resolution array: " << resArr->GetName() << std::flush;
  Int_t nbinsx = hRes2D->GetXaxis()->GetNbins();
  for (Int_t bini = 1; bini <= nbinsx; bini++) {
    TString newname = Form("bin%d", bini);
    TH1D *hSlice = new TH1D(
        *((TH1D *)hRes2D->ProjectionY(newname.Data(), bini, bini, "e")));
    hSlice->SetMarkerStyle(20);
    resArr->Add(hSlice);
  }
  std::cout << "   ... done. Number of entries: " << resArr->GetLast()
            << std::endl;

  fOut->cd();
  resArr->Write(resArr->GetName(), TObject::kSingleKey);
}

int MakeResolutionArray() {
  
  Init_ppAt5TeV();
  
  //...
  fIn = new TFile(fNameIn.Data(), "READ");
  fOut = new TFile(fNameOut.Data(), "RECREATE");

  MakeArray2Sum("hist_deltaPtRel", "RelPtResArrCocktail", "PtGen_DeltaPtOverPtGen", fP_bins_vec);
  MakeArray2Sum("hist_deltaEta", "EtaResArrVsPt", "PtGen_DeltaEta", fP_bins_vec);
  MakeArray2("hist_deltaPhi_ele", "PhiEleResArrVsPt", "PtGen_DeltaPhi_Ele", fP_bins_vec);
  MakeArray2("hist_deltaPhi_pos", "PhiPosResArrVsPt","PtGen_DeltaPhi_Pos", fP_bins_vec);
 
 
  fOut->Close();

  return 0;
}
