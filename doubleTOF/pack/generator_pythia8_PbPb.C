#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"
#include "Generators/GeneratorPythia8.h"
#include "TRandom3.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include <map>
#include <unordered_set>
//#include <utility>	// for std::pair

using namespace Pythia8;

class GeneratorPythia8GunPbPb : public o2::eventgen::GeneratorPythia8{
public:
  /// default constructor
  GeneratorPythia8GunPbPb() = default;
  
  /// constructor
  GeneratorPythia8GunPbPb(int input_pdg){
    genMinPt=0.0;
    genMaxPt=20.0;
    genminY=-1.5;
    genmaxY=1.5;
    genminEta=-1.5;
    genmaxEta=1.5;
    
    UEOverSampling = 20;
    genEventCountUse = 2000; //start at large number: regen
    
    pdg = input_pdg;
    E=0;
    px=0;
    py=0;
    pz=0;
    p=0;
    y=0;
    eta=0;
    xProd=0;
    yProd=0;
    zProd=0;
    xProd=0.; yProd=0.; zProd=0.;
    //addFurtherPion=false;
    
    randomizePDGsign=false;
    
    fLVHelper = new TLorentzVector();
    
    if( input_pdg!=0 ) m = getMass(input_pdg);
    if( input_pdg==0 ) m = 1.0;
    
    //fSpectra = new TF1("fPtDist",myLevyPt,0.0,10,3);
    fSpectra = new TF1("fSpectra",this ,&GeneratorPythia8GunPbPb::myLevyPt, 0.0,20,3, "GeneratorPythia8GunPbPb","myLevyPt");
    fSpectra = STAR_BlastWave("fSpectra",  m, 20);
    fSpectra ->SetNpx( 1000 );
    
    /*
     Double_t pt = x[0];
     Double_t mass = p[0];
     Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
     Double_t beta_max = p[1];
     Double_t temp = p[2];
     Double_t n = p[3];
     Double_t norm = p[4];
     */
    
    // source: https://arxiv.org/pdf/1910.07678.pdf
    
    fSpectra->SetParameter(0, m);   //mass: automatic!
    fSpectra->SetParameter(1,0.6615);   //beta-max
    fSpectra->SetParameter(2,0.0905);   //T
    fSpectra->SetParameter(3,0.7355);   //n
    fSpectra->SetParameter(4,1000);   //norm (not relevant)
    
    furtherPrim={};
    keys_furtherPrim={};
    cout<<"Initalizing extra PYTHIA object"<<endl;
    // Read settings from external file.
    pythiaObject.readFile("pythia8_hi.cmnd");
    //pythiaObject.readFile("pp13.cmnd");
    pythiaObject.init();
    cout << "Done." << endl;

  }
  
  ///  Destructor
  ~GeneratorPythia8GunPbPb() = default;
  
  /// set PDG code
  void setPDG(int input_pdg){pdg=input_pdg;}
  
  /// randomize the PDG code sign of core particle
  void setRandomizePDGsign(){randomizePDGsign=true;}
  
  Double_t myLevyPt(const Double_t *pt, const Double_t *par)
  {
    //Levy Fit Function
    Double_t lMass  = 4.797; //pion Mass
    Double_t ldNdy  = par[0];
    Double_t lTemp = par[1];
    Double_t lPower = par[2];

    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);

    return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
  }
  
  //___________________________________________________________________

  Double_t STAR_BlastWave_Func(const Double_t *x, const Double_t *p) {
      /* dN/dpt */
      
      Double_t pt = x[0];
      Double_t mass = p[0];
      Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
      Double_t beta_max = p[1];
      Double_t temp = p[2];
      Double_t n = p[3];
      Double_t norm = p[4];
      
      if (!fIntegrand)
        fIntegrand = new TF1("fIntegrand", this, &GeneratorPythia8GunPbPb::STAR_BlastWave_Integrand_Improved, 0., 1., 5, "GeneratorPythia8GunPbPb", "STAR_BlastWave_Integrand_Improved");
      fIntegrand->SetParameters(mt, pt, beta_max, temp, n);

      Double_t integral = fIntegrand->Integral(0., 1.);
      return norm * pt * integral;
  }

  //___________________________________________________________________
  
  Double_t STAR_BlastWave_Integrand_Improved(const Double_t *x, const Double_t *p) {
      
      /*
       x[0] -> r (radius)
       p[0] -> mT (transverse mass)
       p[1] -> pT (transverse momentum)
       p[2] -> beta_max (surface velocity)
       p[3] -> T (freezout temperature)
       p[4] -> n (velocity profile)
       */
      
      Double_t r = x[0];
      Double_t mt = p[0];
      Double_t pt = p[1];
      Double_t beta_max = p[2];
      Double_t temp_1 = 1. / p[3];
      Double_t n = p[4];
      
      Double_t beta = beta_max * TMath::Power(r, n);
      Double_t rho = TMath::ATanH(beta);
      Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
      Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
      //  if (argI0 > 100 || argI0 < -100)
      //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
      return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
      
  }
  
  //___________________________________________________________________

  TF1 *STAR_BlastWave(const Char_t *name, Double_t mass,Float_t upperlim, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6) {
    
      //new TF1("fSpectra",this ,&GeneratorPythia8GunPbPb::myLevyPt, 0.0,20,3, "GeneratorPythia8GunPbPb","myLevyPt");
      TF1 *fBlastWave = new TF1(name, this, &GeneratorPythia8GunPbPb::STAR_BlastWave_Func, 0., upperlim, 5, "GeneratorPythia8GunPbPb", "STAR_BlastWave_Func");
      fBlastWave->SetParameters(mass, beta_max, temp, n, norm);
      fBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
      fBlastWave->FixParameter(0, mass);
      fBlastWave->SetParLimits(1, 0.1, 0.9); // don't touch :) adding some 99 youu get floating point exception
      fBlastWave->SetParLimits(2, 0.03,1.);//0.05, 1.);  // no negative values!! for the following as well
      fBlastWave->SetParLimits(3, 0.25,4.5); // was 2.5  // omega-->at limit even moving it to 4.5 but yield same
      return fBlastWave;
  }

  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  
  Double_t y2eta(Double_t pt, Double_t mass, Double_t y){
    Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
    return TMath::ASinH(mt / pt * TMath::SinH(y));
  }
  
  /// set mass
  void setMass(int input_m){m=input_m;}
  
  /// set 4-momentum
  void set4momentum(double input_px, double input_py, double input_pz){
    px = input_px;
    py = input_py;
    pz = input_pz;
    E  = sqrt( m*m+px*px+py*py+pz*pz );
    fourMomentum.px(px);
    fourMomentum.py(py);
    fourMomentum.pz(pz);
    fourMomentum.e(E);
    p   = sqrt( px*px+py*py+pz*pz );
    y   = 0.5*log( (E+pz)/(E-pz) );
    eta = 0.5*log( (p+pz)/(p-pz) );
    
    ////std::cout << "##### Particle #####" << std::endl;
    ////std::cout << " - PDG code: " << pdg << std::endl;
    ////std::cout << " - mass: "     << m   << std::endl;
    ////std::cout << " - (px,py,pz): (" << px << "," << py << "," << pz << ")" << std::endl;
    ////std::cout << " - momentum: " << p << std::endl;
    ////std::cout << " - energy: " << E << std::endl;
    ////std::cout << " - rapidity: " << y << std::endl;
    ////std::cout << " - pseudorapidity: " << eta << std::endl;
    ////std::cout << " - production vertex: (" << xProd << "," << yProd << "," << zProd << ")" << std::endl;
  }
  
  /// set 3-momentum
  void setMomentum(double input_p){p=input_p;}
  
  /// set x,y,z of production vertex
  void setProdVtx(double input_xProd, double input_yProd, double input_zProd){xProd=input_xProd; yProd=input_xProd; zProd=input_zProd;}
  
  /// setter to add further primary particles to the event
  void setAddFurtherPrimaries(const int pdgCode, const int howMany){
    /// check if this species has been already added
    const int map_counts = furtherPrim.count(pdgCode);
    if(map_counts==1){	// species already present
      const int howMany_already = furtherPrim[pdgCode];
      std::cout << "BEWARE: " << howMany_already << " particles of species " << pdgCode << " already required.";
      std::cout << " Ignoring the command setAddFurtherPrimaries(" << pdgCode << "," << howMany << ")" << std::endl;
      return;
    }
    /// add particles, if not yet present
    furtherPrim[pdgCode] = howMany;
    keys_furtherPrim.insert(pdgCode);
  }
  
  /// set add a further primary pion
  //void setAddFurtherPion(){addFurtherPion=true;}
  
  /// get mass from TParticlePDG
  double getMass(int input_pdg){
    double mass = 0;
    if(TDatabasePDG::Instance()){
      TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(input_pdg);
      if(particle)	mass = particle->Mass();
      else			std::cout << "===> particle mass equal to 0" << std::endl;
    }
    //Override database if needed, please
    if( input_pdg==4444 ) mass = 4.797; //OmegaCCC
    if( input_pdg==4432 ) mass = 3.746; //OmegaCC
    if( input_pdg==4422) mass = 3.621; //XiCC
    return mass;
  }
  
  //_________________________________________________________________________________
  /// generate uniform eta and uniform momentum
  void genUniformMomentumEta(double minP, double maxP, double minY, double maxY){
    // random generator
    std::unique_ptr<TRandom3> ranGenerator { new TRandom3() };
    ranGenerator->SetSeed(0);
    
    // momentum
    const double gen_p = ranGenerator->Uniform(minP,maxP);
    // eta
    const double gen_eta = ranGenerator->Uniform(minY,maxY);
    // z-component momentum from eta
    const double cosTheta = ( exp(2*gen_eta)-1 ) / ( exp(2*gen_eta)+1 );	// starting from eta = -ln(tan(theta/2)) = 1/2*ln( (1+cos(theta))/(1-cos(theta)) ) ---> NB: valid for cos(theta)!=1
    const double gen_pz = gen_p*cosTheta;
    // y-component: random uniform
    const double maxVal = sqrt( gen_p*gen_p-gen_pz*gen_pz );
    double sign_py = ranGenerator->Uniform(0,1);
    sign_py = (sign_py>0.5)?1.:-1.;
    const double gen_py = ranGenerator->Uniform(0.,maxVal)*sign_py;
    // x-component momentum
    double sign_px = ranGenerator->Uniform(0,1);
    sign_px = (sign_px>0.5)?1.:-1.;
    const double gen_px = sqrt( gen_p*gen_p-gen_pz*gen_pz-gen_py*gen_py )*sign_px;
    
    set4momentum(gen_px,gen_py,gen_pz);
  }
  
  //_________________________________________________________________________________
  /// generate uniform eta and uniform momentum
  void genSpectraMomentumEta(double minP, double maxP, double minY, double maxY){
    // random generator
    std::unique_ptr<TRandom3> ranGenerator { new TRandom3() };
    ranGenerator->SetSeed(0);
    
    // generate transverse momentum
    const double gen_pT = fSpectra->GetRandom(minP,maxP);
    
    //Actually could be something else without loss of generality but okay
    const double gen_phi = ranGenerator->Uniform(0,2*TMath::Pi());
    
    // sample flat in rapidity, calculate eta
    Double_t gen_Y=10, gen_eta=10;
    
    while( gen_eta>genmaxEta || gen_eta<genminEta ){
      gen_Y = ranGenerator->Uniform(minY,maxY);
      //(Double_t pt, Double_t mass, Double_t y)
      gen_eta = y2eta(gen_pT, m, gen_Y);
    }
    
    fLVHelper->SetPtEtaPhiM(gen_pT, gen_eta, gen_phi, m);
    set4momentum(fLVHelper->Px(),fLVHelper->Py(),fLVHelper->Pz());
  }
  
protected:
  
  //__________________________________________________________________
  Particle createParticle(){
    //std::cout << "createParticle() mass " << m << " pdgCode " << pdg << std::endl;
    Particle myparticle;
    myparticle.id(pdg);
    myparticle.status(11);
    myparticle.px(px);
    myparticle.py(py);
    myparticle.pz(pz);
    myparticle.e(E);
    myparticle.m(m);
    myparticle.xProd(xProd);
    myparticle.yProd(yProd);
    myparticle.zProd(zProd);
    
    return myparticle;
  }
  
  //__________________________________________________________________
  int randomizeSign(){
    
    std::unique_ptr<TRandom3> gen_random {new TRandom3(0)};
    const float n = gen_random->Uniform(-1,1);
    
    return n/abs(n);
  }
  
  //__________________________________________________________________
  Bool_t generateEvent() override {
    
    double original_m = m;
    int original_pdg  = pdg;
    
    /// reset event
    mPythia.event.reset();
    
    if(original_pdg!=211){
      for(Int_t ii=0; ii<15; ii++){
        xProd=0.0;
        yProd=0.0;
        zProd=0.0;
        genSpectraMomentumEta(genMinPt,genMaxPt,genminY,genmaxY);
        Particle lAddedParticle = createParticle();
        mPythia.event.append(lAddedParticle);
        lAddedParticles++;
      }
    }
    
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    Bool_t lPythiaOK = kFALSE;
    while (!lPythiaOK){
      lPythiaOK = pythiaObject.next();
      //Select rough central events, please, disregard
      //if( pythiaObject.info.hiInfo->b() > 6) lPythiaOK = kFALSE; //regenerate, please
    }
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // use event
    Long_t nParticles = pythiaObject.event.size();
    lAddedParticles = 0;
    for ( Long_t j=0; j < nParticles; j++ ) {
      Int_t pypid = pythiaObject.event[j].id();
      Float_t pypx = pythiaObject.event[j].px();
      Float_t pypy = pythiaObject.event[j].py();
      Float_t pypz = pythiaObject.event[j].pz();
      Float_t pypt = pythiaObject.event[j].pT();
      Float_t pyrap = pythiaObject.event[j].y();
      Float_t pyeta = pythiaObject.event[j].eta();
      Int_t pystate = pythiaObject.event[j].status();
      //if(TMath::Abs(state > 89)) {continue;}
      Float_t pyenergy = pythiaObject.event[j].e();
      Int_t pycharge = pythiaObject.event[j].charge();
           
      //Per-species loop: skip outside of mid-rapidity, please
      if ( TMath::Abs(pyeta) > 2.0 ) continue;
      
      //final only
      if (!pythiaObject.event[j].isFinal()) continue;
      
      pdg = pypid;
      px = pypx;
      py = pypy;
      pz = pypz;
      E = pyenergy;
      m = pythiaObject.event[j].m();
      xProd = pythiaObject.event[j].xProd();
      yProd = pythiaObject.event[j].yProd();
      zProd = pythiaObject.event[j].zProd();
      
      Particle lAddedParticle = createParticle();
      mPythia.event.append(lAddedParticle);
      lAddedParticles++;
    }
    //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
    //Revert back or else there'll be trouble next time!
    m = original_m;
    pdg = original_pdg;
    
    /// go to next Pythia event
    mPythia.next();
    
    return true;
  }
  
private:
  
  double genMinPt;			/// minimum 3-momentum for generated particles
  double genMaxPt;			/// maximum 3-momentum for generated particles
  double genminY;		/// minimum pseudorapidity for generated particles
  double genmaxY;		/// maximum pseudorapidity for generated particles
  double genminEta;
  double genmaxEta;
  int UEOverSampling; //number of times to repeat underlying event
  int genEventCountUse;
  
  Vec4   fourMomentum;	/// four-momentum (px,py,pz,E)
  double E;				/// energy: sqrt( m*m+px*px+py*py+pz*pz ) [GeV/c]
  double m;				/// particle mass [GeV/c^2]
  int    pdg;				/// particle pdg code
  double px;				/// x-component momentum [GeV/c]
  double py;				/// y-component momentum [GeV/c]
  double pz;				/// z-component momentum [GeV/c]
  double p;				/// momentum
  double y;				/// rapidity
  double eta;				/// pseudorapidity
  double xProd;			/// x-coordinate position production vertex [cm]
  double yProd;			/// y-coordinate position production vertex [cm]
  double zProd;			/// z-coordinate position production vertex [cm]
  
  //Max number: max number of particles to be added
  long lAddedParticles;
  
  bool randomizePDGsign;	/// bool to randomize the PDG code of the core particle
  
  TF1 *fSpectra; /// TF1 to store more realistic shape of spectrum
  
  //BW integrand
  TF1 *fIntegrand = NULL;
  
  TLorentzVector *fLVHelper;
  
  Pythia pythiaObject; ///Generate a full event if requested to do  so
  
  //bool   addFurtherPion;	/// bool to attach an additional primary pion
  std::map<int,int> furtherPrim;				/// key: PDG code; value: how many further primaries of this species to be added
  std::unordered_set<int> keys_furtherPrim;	/// keys of the above map (NB: only unique elements allowed!)
};

//The Omega Family
FairGenerator* generateNativeOmegaCCC(){
  auto myGen = new GeneratorPythia8GunPbPb(4444);
  return myGen;
}

FairGenerator* generateNativeOmegaCC(){
  auto myGen = new GeneratorPythia8GunPbPb(4432);
  return myGen;
}

FairGenerator* generateNativeOmegaC(){
  auto myGen = new GeneratorPythia8GunPbPb(4332);
  return myGen;
}

FairGenerator* generateNativeOmega(){
  auto myGen = new GeneratorPythia8GunPbPb(3334);
  return myGen;
}

//The Xi Family
FairGenerator* generateNativeXiCC(){
  auto myGen = new GeneratorPythia8GunPbPb(4422);
  return myGen;
}

FairGenerator* generateNativeXiC(){
  auto myGen = new GeneratorPythia8GunPbPb(4232);
  return myGen;
}

FairGenerator* generateNativeXi(){
  auto myGen = new GeneratorPythia8GunPbPb(3312);
  return myGen;
}

//Use this for minimum-bias
FairGenerator* generatePYTHIA(){
  auto myGen = new GeneratorPythia8GunPbPb(211);
  return myGen;
}
