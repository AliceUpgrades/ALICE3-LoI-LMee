#include <iostream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "DetectorsCommonDataFormats/DetID.h"
#include "SimulationDataFormat/MCTrack.h"

bool isStable(int pdg);
bool isInFITacc(double eta);
bool isHF(int pdg);

enum isCharm { kIsNoCharm, kIsCharm, kIsCharmFromBeauty };

void ana(TString generator = "hijing")
{
	TChain mcTree("o2sim");
	mcTree.AddFile(Form("../../run/%s/tmp/hijing_PbPb_b45_Kine.root",generator.Data()));
	mcTree.SetBranchStatus("*", 0);
	mcTree.SetBranchStatus("MCTrack*", 1);

	std::vector<o2::MCTrack> *mcTracks = nullptr;
	mcTree.SetBranchAddress("MCTrack", &mcTracks);

	const int nEvents = mcTree.GetEntries();
	TH1D hEvents {"hEvents","nEvents",1,0,1};
	hEvents.SetBinContent(1,nEvents);
	const float r = 10.;
	TH1D hMult {"hMult", "Multiplicity",100,0,3500};
	TH1D hfePt {"hHfePt", "p_{T} spectrum of electrons from charm; p_{T} (GeV/c)",200,0,10};
	TH2D hVertex {"hVertex", "prod. vertices of e^{+}/e^{-} with photon mother;x (cm);y (cm)", 1000, -r, r, 1000, -r, r};
	TH2D hVertexR {"hVertexR", "prod. vertices of e^{+}/e^{-} with photon mother;z (cm);r (cm)", 1000, -20., 20., 1000, 0., 5.};
	TH1D hInvMass {"hInvMass", "invariant mass of e^{+}/e^{-} with photon mother;m (GeV/c);N / N_{ev}", 200, 0., 2.};
	TH1D hInvMassPrim {"hInvMassPrim", "invariant mass of primary e^{+}/e^{-};m (GeV/c);N / N_{ev}", 200, 0., 2.};

	TH1D hLS1 {"hLS1","Like Sign spectrum of -- pairs from all particles;m_{ee} (GeV/c^{2})", 400,0,4};
	TH1D hLS2 {"hLS2","Like Sign spectrum of ++ pairs from all particles;m_{ee} (GeV/c^{2})", 400,0,4};
	TH1D hLS1conv {"hLS1conv","Like Sign spectrum of -- pairs including at least one conversion electron;m_{ee} (GeV/c^{2})", 400,0,4};
	TH1D hLS2conv {"hLS2conv","Like Sign spectrum of ++ pairs including at least one conversion electron;m_{ee} (GeV/c^{2})", 400,0,4};
	TH1D hLS1prim {"hLS1prim","Like Sign spectrum of -- pairs from primary particles;m_{ee} (GeV/c^{2})", 400,0,4};
	TH1D hLS2prim {"hLS2prim","Like Sign spectrum of ++ pairs from primary particles;m_{ee} (GeV/c^{2})", 400,0,4};
	TH1D hLS1primHF {"hLS1primHF","Like Sign spectrum of -- pairs from primary electrons including at least one HFe;m_{ee} (GeV/c^{2})", 400,0,4};
	TH1D hLS2primHF {"hLS2primHF","Like Sign spectrum of ++ pairs from primary electrons including at least one HFe;m_{ee} (GeV/c^{2})", 400,0,4};


	double eMass = 0.000511;
	double ptCut = 0.2;
	double etaCut = 0.8;

	bool runPrefilter = false;
	double prfltrM_cut = 0.050;
	double prfltrPhi_cut = 0.050;

	// lambda to define prefilter criterion
	// this can be anything that returns a boolean
	auto prfltr = [&prfltrM_cut,&prfltrPhi_cut](o2::MCTrack &t1,o2::MCTrack &t2) {
		TLorentzVector lv1,lv2;
		t1.Get4Momentum(lv1);
		t2.Get4Momentum(lv2);
		if ( (lv1 + lv2).M() < prfltrM_cut && (fabs(t1.GetPhi() - t2.GetPhi()) < prfltrPhi_cut) ) return false;
		else return true;
	};


	std::vector<o2::MCTrack> ep, em, ep_prim, em_prim;
	for(int iEvent = 0; iEvent < nEvents; ++iEvent) {
		ep.clear();
		em.clear();
		ep_prim.clear();
		em_prim.clear();
		mcTree.GetEntry(iEvent);
		int nConv = 0;
		int Ntracks = 0;
		for (const auto track : *mcTracks) {
			// if (!track.isPrimary()) continue; //only tracks from generator, not geant
			const auto r_vtx = std::sqrt(std::pow(track.GetStartVertexCoordinatesX(), 2) + std::pow(track.GetStartVertexCoordinatesY(), 2));
			if (r_vtx > 0.1) continue; // additional 'primary' selection
			auto pdg = track.GetPdgCode();
			if (!pdg) continue; // if pdg code is 0 something is wrong
			TParticlePDG *pPDG = TDatabasePDG::Instance()->GetParticle(pdg);
			if (!pPDG) continue; // if no particle, something is wrong
			if (fabs(pPDG->Charge()) < 3.) continue; // charge of quarks is 1/2 hadrons is 3
			if (!isStable(pdg)) continue; // check list of 'stable particles'
			int mpdg = 0;
			auto motherId = track.getMotherTrackId();
			if (motherId < 0) mpdg = 1;
			else {auto mTrack = (*mcTracks)[motherId];
			      mpdg = mTrack.GetPdgCode();}
			if (isStable(pdg) && isStable(mpdg) ) continue; // if both are mother and daugther are stable, just fill the mother (other iteration)
			if (abs(track.GetEta()) > 0.5)  continue;
			// if (abs(track.GetPt()) < ptCut)  continue;
			Ntracks++; // check that we are in the right eta range
		}
		hMult.Fill(Ntracks);

		for (const auto track : *mcTracks) {
			if (abs(track.GetPdgCode()) != 11) continue;
			if (TMath::Abs(track.GetEta()) > etaCut) continue;
			if (track.GetPt() < ptCut) continue;

			const auto r_vtx = std::sqrt(std::pow(track.GetStartVertexCoordinatesX(), 2) + std::pow(track.GetStartVertexCoordinatesY(), 2));
			if (r_vtx > 0.55) continue;
			if (r_vtx < 0.05) // check for primaries... TODO: implement proper is physical primary based on particle species etc.
			{
				// fill primaries in both vectors
				if (track.GetPdgCode() > 0)
				{
					ep.emplace_back(track);
					ep_prim.emplace_back(track);
				}
				else
				{
					em.emplace_back(track);
					em_prim.emplace_back(track);
				}
			}
			auto motherId = track.getMotherTrackId();
			if (motherId < 0) continue;
			auto mTrack = (*mcTracks)[motherId];
			// our track is an electron and we have its mother!
			// check if its a charm mother:
			bool isCharmMother = false;
			bool isFromBeauty = false;
			int mpdg = abs(mTrack.GetPdgCode());
			if (((mpdg > 400) && (mpdg < 439)) || ((mpdg > 4000) && (mpdg < 4399))) isCharmMother = true;
			while (true){ // recursiv loop until beauty is found, or the stack is done.
				auto gmId = mTrack.getMotherTrackId(); // get the mother id
				if(gmId < 0) break; // break condition if done here
				mTrack =  (*mcTracks)[gmId]; // get track from id
				auto gmpdg = abs(mTrack.GetPdgCode()); // get the pdg code
				if (((gmpdg > 500) && (gmpdg < 549)) || ((gmpdg > 5000) && (gmpdg < 5499))) {isFromBeauty = true; break;} // if there is beauty, we can stop.
			}
			if(isCharmMother || isFromBeauty) hfePt.Fill(track.GetPt()); // we are looking for HFE, so both are fine
			if (mpdg != 22)  continue; // conversions will only be taken into account if they are produced in the inner most layer/ Be-foil and come from a photon... no weak decays
			hVertex.Fill(track.GetStartVertexCoordinatesX(), track.GetStartVertexCoordinatesY());
			hVertexR.Fill(track.GetStartVertexCoordinatesZ(), r_vtx);
			++nConv;
			if (track.GetPdgCode() == 11) ep.emplace_back(track);
			else em.emplace_back(track);
		}


		// vectors to store indicies of elements we want to get rid of
	  std::vector<int> vDelE,vDelP;
		for (unsigned long int e = 0; e<em.size(); ++e) {
	    for (unsigned long int p = 0; p<ep.size(); ++p) {
				if(!prfltr(em[e],ep[p])) {vDelE.push_back(e);vDelP.push_back(p);}
			}
		}
		// sort the indizies
	  // not sure if it is really needed, but better be safe
	  std::sort(vDelE.begin(),vDelE.end());
	  // erase entries that are not unique
	  vDelE.erase(std::unique(vDelE.begin(), vDelE.end()), vDelE.end());
	  // for second vector do as above
	  std::sort(vDelP.begin(),vDelP.end());
	  vDelP.erase(std::unique(vDelP.begin(), vDelP.end()), vDelP.end());
		// only delete if we run with a prefilter
		if (runPrefilter){
			// loop from end to begin of vector to not screw up ordering when deleting a elemt of the vector
		  for(int i = vDelE.size()-1; i >= 0; i--)
		  // erase element at position stored in index vector vDelE
		  { em.erase(em.begin()+vDelE.at(i));}
		  for(int i = vDelP.size()-1; i >= 0; i--){ ep.erase(ep.begin()+vDelP.at(i));}
		}
		// printf("nConv / nTracks = %i / %zu\n", nConv, mcTracks->size());
		for (auto p : ep) {
			TLorentzVector vp;
			p.Get4Momentum(vp);
			for (auto e : em) {
				TLorentzVector ve;
				e.Get4Momentum(ve);
				const float mass = (ve + vp).M();
				hInvMass.Fill(mass);
			}
		}
		for (auto p : ep_prim) {
			TLorentzVector vp;
			p.Get4Momentum(vp);
			for (auto e : em_prim) {
				TLorentzVector ve;
				e.Get4Momentum(ve);
				const float mass = (ve + vp).M();
				hInvMassPrim.Fill(mass);
			}
		}
		for (auto track1 = em.begin(); track1!=em.end(); ++track1)
		{
			TLorentzVector LV1;
			track1->Get4Momentum(LV1);
			auto motherId1 = track1->getMotherTrackId();
			if (motherId1 < 0) continue;
			auto mTrack1 = (*mcTracks)[motherId1];
			int mpdg1 = abs(mTrack1.GetPdgCode());

			for (auto track2 = track1+1; track2!=em.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS1.Fill((LV1+LV2).M());
				auto motherId2 = track2->getMotherTrackId();
				if (motherId2 < 0) continue;
				auto mTrack2 = (*mcTracks)[motherId2];
				int mpdg2 = abs(mTrack2.GetPdgCode());
				if((mpdg1 == 22) || (mpdg2 == 22)) hLS1conv.Fill((LV1+LV2).M());
			}
		}
		for (auto track1 = ep.begin(); track1!=ep.end(); ++track1)
		{
			TLorentzVector LV1;
			track1->Get4Momentum(LV1);
			auto motherId1 = track1->getMotherTrackId();
			if (motherId1 < 0) continue;
			auto mTrack1 = (*mcTracks)[motherId1];
			int mpdg1 = abs(mTrack1.GetPdgCode());
			for (auto track2 = track1+1; track2!=ep.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS2.Fill((LV1+LV2).M());
				auto motherId2 = track2->getMotherTrackId();
				if (motherId2 < 0) continue;
				auto mTrack2 = (*mcTracks)[motherId2];
				int mpdg2 = abs(mTrack2.GetPdgCode());
				if((mpdg1 == 22) || (mpdg2 == 22)) hLS2conv.Fill((LV1+LV2).M());
			}
		}

		for (auto track1 = em_prim.begin(); track1!=em_prim.end(); ++track1)
		{
			TLorentzVector LV1;
			track1->Get4Momentum(LV1);
			auto motherId1 = track1->getMotherTrackId();
			if (motherId1 < 0) continue;
			auto mTrack1 = (*mcTracks)[motherId1];
			int mpdg1 = abs(mTrack1.GetPdgCode());
			for (auto track2 = track1+1; track2!=em_prim.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS1prim.Fill((LV1+LV2).M());
				auto motherId2 = track2->getMotherTrackId();
				if (motherId2 < 0) continue;
				auto mTrack2 = (*mcTracks)[motherId2];
				int mpdg2 = abs(mTrack2.GetPdgCode());
				if (isHF(mpdg1) || isHF(mpdg2)) hLS1primHF.Fill((LV1+LV2).M());
			}
		}
		for (auto track1 = ep_prim.begin(); track1!=ep_prim.end(); ++track1)
		{
			TLorentzVector LV1;
			track1->Get4Momentum(LV1);
			auto motherId1 = track1->getMotherTrackId();
			if (motherId1 < 0) continue;
			auto mTrack1 = (*mcTracks)[motherId1];
			int mpdg1 = abs(mTrack1.GetPdgCode());
			for (auto track2 = track1+1; track2!=ep_prim.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS2prim.Fill((LV1+LV2).M());
				auto motherId2 = track2->getMotherTrackId();
				if (motherId2 < 0) continue;
				auto mTrack2 = (*mcTracks)[motherId2];
				int mpdg2 = abs(mTrack2.GetPdgCode());
				if (isHF(mpdg1) || isHF(mpdg2)) hLS2primHF.Fill((LV1+LV2).M());
			}
		}


	}

	// hLS1prim.Add(&hLS1prim);
	// FillLsHist(&hLS1,em);
	//
	// TCanvas c("c", "c", 1600, 1600);
	// hVertex.Draw("colz");
	// hVertex.SetStats(false);
	// c.SaveAs("conv_xy.png");
	// hVertexR.Draw("colz");
	// hVertexR.SetStats(false);
	// c.SaveAs("conv_rz.png");
	// // hLS1.Scale(1./nEvents);
	// hLS1.SetMarkerStyle(20);
	// hLS1.SetMarkerColor(kBlack);
	// hLS1.SetLineColor(kBlack);
	// // hLS1prim.Scale(1./nEvents);
	// hLS1prim.SetMarkerStyle(20);
	// hLS1prim.SetMarkerColor(kRed);
	// hLS1prim.SetLineColor(kRed);
	// hLS1.SetMaximum(10);
	// hLS1.SetMinimum(0.001);
	// c.SetLogy();
	// hLS1.SetStats(0);
	// hLS1.Draw();
	// hLS1prim.Draw("same");
	// TLegend l1(0.6, 0.7, .9, .9);
	// l1.AddEntry(&hLS1,"All -- pairs (w/ conversions)");
	// l1.AddEntry(&hLS1prim,"All primary -- pairs");
	// l1.Draw();
	// c.SaveAs("../output/massSpectrum.png");
	// c.SetLogy(kFALSE);
	//
	//
	// hInvMass.Scale(1./nEvents);
	// hInvMass.SetStats(kFALSE);
	// hInvMass.Draw();
	// hInvMassPrim.Scale(1./nEvents);
	// hInvMassPrim.Draw("same");
	// hInvMassPrim.SetLineColor(kRed);
	// TLegend l(0.6, 0.7, .9, .9);
	// l.AddEntry(&hInvMass, "conversions", "l");
	// l.AddEntry(&hInvMassPrim, "primary", "l");
	// l.Draw();
	// c.SaveAs("../output/invMass.png");

	std::unique_ptr<TFile> f {TFile::Open(Form("../output/ana_%1.2f_%1.1f_%s_prf_%d.root",ptCut,etaCut,generator.Data(), (int)runPrefilter), "RECREATE")};
	f->WriteTObject(&hEvents);
	f->WriteTObject(&hMult);
	f->WriteTObject(&hVertex);
	f->WriteTObject(&hfePt);
	f->WriteTObject(&hLS1);
	f->WriteTObject(&hLS2);
	f->WriteTObject(&hLS1prim);
	f->WriteTObject(&hLS2prim);
	f->WriteTObject(&hLS1conv);
	f->WriteTObject(&hLS2conv);
	f->WriteTObject(&hLS1primHF);
	f->WriteTObject(&hLS2primHF);
}

bool isStable(int pdg)
{
	if (abs(pdg) > 1000000000) return true;
	switch (abs(pdg)) {
		case 22: return true;      // Photon
    case 11: return true;      // Electron
    case 13: return true;      // Muon
    case 211: return true;     // Pion
    case 321: return true;     // Kaon
    case 310: return true;     // K0s
    case 130: return true;     // K0l
    case 2212: return true;    // Proton
    case 2112: return true;    // Neutron
    case 3122: return true;    // Lambda_0
    case 3212: return true;    // Sigma0
    case 3112: return true;    // Sigma Minus
    case 3222: return true;    // Sigma Plus
    case 3312: return true;    // Xsi Minus
    case 3322: return true;    // Xsi
    case 3334: return true;    // Omega
    case 12: return true;      // Electron Neutrino
    case 14: return true;      // Muon Neutrino
    case 16: return true;
		default: return false;
	}
}



bool isInFITacc(double eta)
{
	if ((2.2 < eta) && (eta< 5.0)) return true;
	if ((-3.4 < eta) && (eta< -2.3)) return true;
	else return false;
}

bool isHF(int pdg)
{
	return (((pdg > 400) && pdg < (600)) || ((pdg > 4000) && pdg < (6000)));
}
