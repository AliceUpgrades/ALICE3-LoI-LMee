#include <iostream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"
#include "TLegend.h"

#include "DetectorsCommonDataFormats/DetID.h"
#include "SimulationDataFormat/MCTrack.h"


void ana()
{
	TChain mcTree("o2sim");
	mcTree.AddFile("../input/hijing/o2sim_Kine.root");
	mcTree.SetBranchStatus("*", 0);
	mcTree.SetBranchStatus("MCTrack*", 1);

	std::vector<o2::MCTrack> *mcTracks = nullptr;
	mcTree.SetBranchAddress("MCTrack", &mcTracks);

	const int nEvents = mcTree.GetEntries();

	const float r = 10.;
	TH2D hVertex {"hVertex", "prod. vertices of e^{+}/e^{-} with photon mother;x (cm);y (cm)", 1000, -r, r, 1000, -r, r};
	TH2D hVertexR {"hVertexR", "prod. vertices of e^{+}/e^{-} with photon mother;z (cm);r (cm)", 1000, -20., 20., 1000, 0., 5.};
	TH1D hInvMass {"hInvMass", "invariant mass of e^{+}/e^{-} with photon mother;m (GeV/c);N / N_{ev}", 200, 0., 2.};
	TH1D hInvMassPrim {"hInvMassPrim", "invariant mass of primary e^{+}/e^{-};m (GeV/c);N / N_{ev}", 200, 0., 2.};

	TH1D hLS1 {"hLS1","Like Sign spectrum of -- pairs from all particles;m_{ee} (GeV/c^{2})", 100,0,1};
	TH1D hLS2 {"hLS2","Like Sign spectrum of ++ pairs from all particles;m_{ee} (GeV/c^{2})", 100,0,1};
	TH1D hLS1prim {"hLS1prim","Like Sign spectrum of -- pairs from primary particles;m_{ee} (GeV/c^{2})", 100,0,1};
	TH1D hLS2prim {"hLS2prim","Like Sign spectrum of ++ pairs from primary particles;m_{ee} (GeV/c^{2})", 100,0,1};


	double eMass = 0.000511;

	std::vector<o2::MCTrack> ep, em, ep_prim, em_prim;
	for(int iEvent = 0; iEvent < nEvents; ++iEvent) {
		ep.clear();
		em.clear();
		ep_prim.clear();
		em_prim.clear();
		mcTree.GetEntry(iEvent);
		int nConv = 0;
		for (const auto track : *mcTracks) {
			if (TMath::Abs(track.GetPdgCode()) != 11) continue;
			if (TMath::Abs(track.GetEta()) > 1.1) continue;
			if (track.GetPt() < 0.04) continue;
			if (sqrt(track.GetStartVertexCoordinatesX()*track.GetStartVertexCoordinatesX() + track.GetStartVertexCoordinatesY()*track.GetStartVertexCoordinatesY())> 0.55) continue;
			if (track.isPrimary()) {
				if (track.GetPdgCode() > 0) ep_prim.emplace_back(track);
				else em_prim.emplace_back(track);
				// continue;
			}
			auto motherId = track.getMotherTrackId();
			if (motherId < 0) continue;
			// auto mTrack = (*mcTracks)[motherId];
			// if (mTrack.GetPdgCode() != 22) continue;
			hVertex.Fill(track.GetStartVertexCoordinatesX(), track.GetStartVertexCoordinatesY());
			const auto r_vtx = std::sqrt(std::pow(track.GetStartVertexCoordinatesX(), 2) + std::pow(track.GetStartVertexCoordinatesY(), 2));
			hVertexR.Fill(track.GetStartVertexCoordinatesZ(), r_vtx);
			++nConv;
			if (track.GetPdgCode() == 11) ep.emplace_back(track);
			else em.emplace_back(track);
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
			for (auto track2 = track1+1; track2!=em.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS1.Fill((LV1+LV2).M());
			}
		}
		for (auto track1 = ep.begin(); track1!=ep.end(); ++track1)
		{
			TLorentzVector LV1;
			track1->Get4Momentum(LV1);
			for (auto track2 = track1+1; track2!=ep.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS2.Fill((LV1+LV2).M());
			}
		}

		for (auto track1 = em_prim.begin(); track1!=em_prim.end(); ++track1)
		{
			TLorentzVector LV1;
			track1->Get4Momentum(LV1);
			for (auto track2 = track1+1; track2!=em_prim.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS1prim.Fill((LV1+LV2).M());
			}
		}
		for (auto track1 = ep_prim.begin(); track1!=ep_prim.end(); ++track1)
		{
			TLorentzVector LV1;
			track1->Get4Momentum(LV1);
			for (auto track2 = track1+1; track2!=ep_prim.end(); ++track2)
			{
				TLorentzVector LV2;
				track2->Get4Momentum(LV2);
				hLS2prim.Fill((LV1+LV2).M());
			}
		}


	}

	// hLS1prim.Add(&hLS1prim);
	// FillLsHist(&hLS1,em);

	TCanvas c("c", "c", 1600, 1600);
	hVertex.Draw("colz");
	hVertex.SetStats(false);
	c.SaveAs("conv_xy.png");
	hVertexR.Draw("colz");
	hVertexR.SetStats(false);
	c.SaveAs("conv_rz.png");
	hLS1.Scale(1./nEvents);
	hLS1.SetMarkerStyle(20);
	hLS1.SetMarkerColor(kBlack);
	hLS1.SetLineColor(kBlack);
	hLS1prim.Scale(1./nEvents);
	hLS1prim.SetMarkerStyle(20);
	hLS1prim.SetMarkerColor(kRed);
	hLS1prim.SetLineColor(kRed);
	hLS1.SetMaximum(10);
	hLS1.SetMinimum(0.001);
	c.SetLogy();
	hLS1.SetStats(0);
	hLS1.Draw();
	hLS1prim.Draw("same");
	TLegend l1(0.6, 0.7, .9, .9);
	l1.AddEntry(&hLS1,"All -- pairs (w/ conversions)");
	l1.AddEntry(&hLS1prim,"All primary -- pairs");
	l1.Draw();
	c.SaveAs("massSpectrum.png");
	c.SetLogy(kFALSE);


	hInvMass.Scale(1./nEvents);
	hInvMass.SetStats(kFALSE);
	hInvMass.Draw();
	hInvMassPrim.Scale(1./nEvents);
	hInvMassPrim.Draw("same");
	hInvMassPrim.SetLineColor(kRed);
	TLegend l(0.6, 0.7, .9, .9);
	l.AddEntry(&hInvMass, "conversions", "l");
	l.AddEntry(&hInvMassPrim, "primary", "l");
	l.Draw();
	c.SaveAs("invMass.png");

	std::unique_ptr<TFile> f {TFile::Open("ana.root", "RECREATE")};
	f->WriteTObject(&hVertex);
}
