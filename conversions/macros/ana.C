#include <iostream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"

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

	for(int iEvent = 0; iEvent < nEvents; ++iEvent) {
		mcTree.GetEntry(iEvent);
		int nConv = 0;
		for (const auto track : *mcTracks) {
			if (TMath::Abs(track.GetPdgCode()) != 11) continue;
			auto motherId = track.getMotherTrackId();
			if (motherId < 0) continue;
			auto mTrack = (*mcTracks)[motherId];
			if (mTrack.GetPdgCode() != 22) continue;
			hVertex.Fill(track.GetStartVertexCoordinatesX(), track.GetStartVertexCoordinatesY());
			++nConv;
		}
		printf("nConv / nTracks = %i / %zu\n", nConv, mcTracks->size());
	}
	TCanvas c("c", "c", 1600, 1600);
	hVertex.Draw("colz");
	c.SaveAs("conv.png");
	std::unique_ptr<TFile> f {TFile::Open("ana.root", "RECREATE")};
	f->WriteTObject(&hVertex);
}
