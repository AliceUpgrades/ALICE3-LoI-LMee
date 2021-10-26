# ifndef anaEEstudy_H
# define anaEEstudy_H

R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libDelphesO2)

//#################################################################
//#                                                               #
//#      Analysis task for the dilepton studys for ALICE 3        #
//#                                                               #
//#  Authors:                                                     #
//#   Florian Eisenhut, Uni Frankfurt / florian.eisenhut@cern.ch  #
//#                                                               #
//#################################################################


class TH1F;
class TH2F;
class TH3F;


//#####################################
// initialize TClonesArray particles, read from Delphes file
//#####################################
TClonesArray* particles;



//#####################################
// initialize detector setup
//#####################################

o2::delphes::TOFLayer innertoflayer;
o2::delphes::TOFLayer toflayer;
o2::delphes::RICHdetector richdetector;
o2::delphes::PreShower preshower;


//#####################################
// initialize further variables
//#####################################

TString smearingfile;



//#####################################
//    initialize global Histograms
//#####################################

// initialize histograms used for weigths
TH1F* hweightHFtoe_lowerPt;
TH1F* hweightHFtoe_higherPt;
TH1F* hweightLFtoe;

// Pair histograms
TH3F* hMPtDCA_ULS_gen[6];
TH3F* hMPtDCA_ULS_rec[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle_charmTOe[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle_beautyTOe[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle_hfTOe[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle_lfTOee[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle_ccTOee[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle_bbTOee[6];
TH3F* hMPtDCA_ULS_rec_MCpidEle_hfTOee[6];
TH3F* hMPtDCA_LS_gen[6];
TH3F* hMPtDCA_LS_rec[6];
TH3F* hMPtDCA_LS_rec_MCpidEle[6];
// TH3F* hMPtDCA_LS_rec_misIDoneLeg[6];
// TH3F* hMPtDCA_LS_rec_misIDtwoLeg[6];
// TH3F* hMPtDCA_LS_rec_misIDPion[6];
// TH3F* hMPtDCA_LS_rec_misIDhf[6];

TH3F* hMPtDCA_LS_rec_charmTOe[6];
TH3F* hMPtDCA_LS_rec_beautyTOe[6];
TH3F* hMPtDCA_LS_rec_hfTOe[6];
TH3F* hMPtDCA_LS_rec_lfTOee[6];
TH3F* hMPtDCA_LS_rec_lfTOee_selectPDG[6];
TH3F* hMPtDCA_LS_rec_ccTOee[6];
TH3F* hMPtDCA_LS_rec_bbTOee[6];
TH3F* hMPtDCA_LS_rec_hfTOee[6];




//#####################################
//    initialize defined functions
//#####################################


void dileptonPairingRec(std::vector<Track *> vec_track_neg,std::vector<Track *> vec_track_pos, bool pairULS, bool MCpidEle, Int_t iSce);
void dileptonPairingGen(std::vector<GenParticle *> vec_track_neg,std::vector<GenParticle *> vec_track_pos, bool pairULS, Int_t iSce);
void PreFilter(std::vector<Track *>* vec_track_neg,std::vector<Track *>* vec_track_pos, Int_t iPID_scenario, std::vector<TH2F*> vecPFhistos);
void PreFilter(std::vector<GenParticle *>* vec_particle,std::vector<Track *>* vec_track, std::vector<Int_t> vecMotherID, std::vector<std::tuple<GenParticle*, GenParticle*>> vecTuples, Int_t iPID_scenario, std::vector<TH2F*> vecPFhistos);
void ApplyPID(std::vector<Track *>* vecTracks_pfPID, std::vector<Track *>* vecTracks, Int_t iPID_scenario);






# endif
