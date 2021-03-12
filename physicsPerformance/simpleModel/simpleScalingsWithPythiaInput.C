void simpleScalingsWithPythiaInput(){
  // instead of using A as scaling we use more realistic 0-10% central Nch, Npart and Ncoll factors from Pythia and Glauber
  // see note from David: https://www.dropbox.com/s/ze2j3txlyqprloj/ALICE_3_System_comparison.pdf?dl=0
  

  const Int_t nSystems = 4;
  TString Syst[nSystems]   = {"Ar","Kr","Xe","Pb"};
  Double_t A[nSystems]     = {40.,78.,129.,208.};
  Double_t lumi[nSystems]  = {1080.,123.,28.9,4.92};// from yellow report (p=1.5, in older version of this macro 1.9 was used)
  Double_t nCh[nSystems]   = {315.05,723.35,1114.45,1830.30};//from David's report (Pythia)
  Double_t nColl[nSystems] = {162.75,456.10,856.70,1598.70};//from David's report (Glauber)

  cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": A = "<<A[i]<<", Lumi(month) = "<<lumi[i]<<" nb-1 ("<<lumi[i]/lumi[nSystems-1]<<")"<<endl;
  }
  cout<<"================================================================"<<endl;
  cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": Signal factor(A) = "<<TMath::Power(A[i]/A[nSystems-1],1.4)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
    cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": Signal factor(Nch) = "<<TMath::Power(nCh[i]/nCh[nSystems-1],1.4)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
  cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": Cocktail factor(A) = "<<TMath::Power(A[i]/A[nSystems-1],1.)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
  cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": Cocktail factor(Nch) = "<<TMath::Power(nCh[i]/nCh[nSystems-1],1.)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
  cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": Ncoll factor(A) = "<<TMath::Power(A[i]/A[nSystems-1],4./3.)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
  cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": Ncoll factor(Ncoll) = "<<TMath::Power(nColl[i]/nColl[nSystems-1],1.)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
  cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": BG factor(A) = "<<TMath::Power(A[i]/A[nSystems-1],2)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
cout<<"================================================================"<<endl;
  for(Int_t i = 0; i < nSystems; i++){
    cout<<Syst[i].Data()<<": BG factor(Nch) = "<<TMath::Power(nCh[i]/nCh[nSystems-1],2)<<" "<<endl;
  }
  cout<<"================================================================"<<endl;
}
