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
	mcTree.AddFile("o2sim_Kine.root");
	mcTree.SetBranchStatus("*", 0);
	mcTree.SetBranchStatus("MCTrack*", 1);

	std::vector<o2::MCTrack> *mcTracks = nullptr;
	mcTree.SetBranchAddress("MCTrack", &mcTracks);

	const int nEvents = mcTree.GetEntries();

	const float r = 10.;
	TH2D hVertex {"hVertex", "prod. vertices of e^{+}/e^{-} with photon mother;x (cm);y (cm)", 1000, -r, r, 1000, -r, r};
	TH2D hVertexR {"hVertexR", "prod. vertices of e^{+}/e^{-} with photon mother;z (cm);r (cm)", 1000, -200., 200., 1000, 0., 5.};
	TH1D hInvMass {"hInvMass", "invariant mass of e^{+}/e^{-} with photon mother;m (GeV/c);N / N_{ev}", 200, 0., 2.};
	TH1D hInvMassPrim {"hInvMassPrim", "invariant mass of primary e^{+}/e^{-};m (GeV/c);N / N_{ev}", 200, 0., 2.};

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
			if (track.isPrimary()) {
				if (track.GetPdgCode() > 0) ep_prim.emplace_back(track);
				else em_prim.emplace_back(track);
				continue;
			}
			auto motherId = track.getMotherTrackId();
			if (motherId < 0) continue;
			auto mTrack = (*mcTracks)[motherId];
			if (mTrack.GetPdgCode() != 22) continue;
			hVertex.Fill(track.GetStartVertexCoordinatesX(), track.GetStartVertexCoordinatesY());
			const auto r_vtx = std::sqrt(std::pow(track.GetStartVertexCoordinatesX(), 2) + std::pow(track.GetStartVertexCoordinatesY(), 2));
			hVertexR.Fill(track.GetStartVertexCoordinatesZ(), r_vtx);
			++nConv;
			if (TMath::Abs(track.GetEta()) > 1.1) continue;
			if (track.GetPt() < 0.04) continue;
			if (track.GetPdgCode() == 11) ep.emplace_back(track);
			else em.emplace_back(track);
		}
		printf("nConv / nTracks = %i / %zu\n", nConv, mcTracks->size());
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
	}

	TCanvas c("c", "c", 1600, 1600);
	hVertex.Draw("colz");
	hVertex.SetStats(false);
	c.SaveAs("conv_xy.png");
	hVertexR.Draw("colz");
	hVertexR.SetStats(false);
	c.SaveAs("conv_rz.png");

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
