#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <vector>

int pitchToIndex(int pitch);
int indexToPitch(int index);

void clusterSizeMacro() {
    // access track processing for data and MC
    TFile* infile_tracks_data = TFile::Open("/home/stephan/MID_cluster/reco_2/track_processing.root", "READ");
    TFile* infile_tracks_MC = TFile::Open("/home/stephan/sims/mu_boxgen_100k_558801_anchor/track_processing.root", "READ");

    // access cluster size files for data and MC
    TFile* infile_clusters_data = TFile::Open("/home/stephan/MID_cluster/reco_2/cluster_size.root", "READ");
    TFile* infile_clusters_MC = TFile::Open("/home/stephan/sims/mu_boxgen_100k_558801_anchor/cluster_size.root", "READ");

    // access track processing distributions
    TH1D* MID_tracks_ROF_data = (TH1D*) infile_tracks_data->Get("mid_track_count_ROF");
    TH1D* MID_tracks_ROF_MC = (TH1D*) infile_tracks_MC->Get("mid_track_count_ROF");

    TH1D* muon_tracks_ROF_data = (TH1D*) infile_tracks_data->Get("muon_track_count_ROF");
    TH1D* muon_tracks_ROF_MC = (TH1D*) infile_tracks_MC->Get("muon_track_count_ROF");

    // access cluster size distributions
    const int nPitches = 3;
    std::vector<TH1D*> nStripsClusterBending_data(nPitches);
    std::vector<TH1D*> nStripsClusterNonBending_data(nPitches);

    std::vector<TH1D*> nStripsClusterBending_MC(nPitches);
    std::vector<TH1D*> nStripsClusterNonBending_MC(nPitches);

    for (int index = 0; index < nPitches; index++) {
        // data
        nStripsClusterBending_data[index] = (TH1D*) infile_clusters_data->Get(Form("cluster_strip_size_bending_%i", indexToPitch(index)));
        nStripsClusterNonBending_data[index] = (TH1D*) infile_clusters_data->Get(Form("cluster_strip_size_nonbending_%i", indexToPitch(index)));

        // MC
        nStripsClusterBending_MC[index] = (TH1D*) infile_clusters_MC->Get(Form("cluster_strip_size_bending_%i", indexToPitch(index)));
        nStripsClusterNonBending_MC[index] = (TH1D*) infile_clusters_MC->Get(Form("cluster_strip_size_nonbending_%i", indexToPitch(index)));
    }

    // plotting
    TFile* outFile = new TFile("cluster_size_data_MC_joined.root", "RECREATE");

    // MID tracks
    TCanvas *canvasMIDtracks = new TCanvas("mid_tracks_ROF","mid_tracks_ROF");
    gPad->SetLogy();

    MID_tracks_ROF_data->Scale(1/MID_tracks_ROF_data->Integral());
    MID_tracks_ROF_data->SetLineColor(1);
    MID_tracks_ROF_data->SetMinimum(0.000001);
    MID_tracks_ROF_data->Draw("HIST ][ SAME");

    MID_tracks_ROF_MC->Scale(1/MID_tracks_ROF_MC->Integral());
    MID_tracks_ROF_MC->SetLineColor(2);
    MID_tracks_ROF_MC->SetMinimum(0.000001);
    MID_tracks_ROF_MC->Draw("HIST ][ SAME");

    auto legendMIDtracks = new TLegend();
    legendMIDtracks->AddEntry(MID_tracks_ROF_data, "Data", "l");
    legendMIDtracks->AddEntry(MID_tracks_ROF_MC, "O2sim", "l");
    legendMIDtracks->Draw("SAME");

    canvasMIDtracks->Write();

    // muon tracks
    TCanvas *canvasMuonTracks = new TCanvas("muon_tracks_ROF","muon_tracks_ROF");
    gPad->SetLogy();

    muon_tracks_ROF_data->Scale(1/muon_tracks_ROF_data->Integral());
    muon_tracks_ROF_data->SetLineColor(1);
    //muon_tracks_ROF_data->SetMinimum(0.000001);
    muon_tracks_ROF_data->Draw("HIST ][ SAME");

    muon_tracks_ROF_MC->Scale(1/muon_tracks_ROF_MC->Integral());
    muon_tracks_ROF_MC->SetLineColor(2);
    //muon_tracks_ROF_MC->SetMinimum(0.000001);
    muon_tracks_ROF_MC->Draw("HIST ][ SAME");

    auto legendMuonTracks = new TLegend();
    legendMuonTracks->AddEntry(muon_tracks_ROF_data, "Data", "l");
    legendMuonTracks->AddEntry(muon_tracks_ROF_MC, "O2sim", "l");
    legendMuonTracks->Draw("SAME");

    canvasMuonTracks->Write();

    // clusters
    for (int index = 0; index < nPitches; index++) {
        // BP
        TCanvas *canvasBP = new TCanvas(Form("cluster_strip_size_BP_pitch_%i", indexToPitch(index)),Form("cluster_strip_size_BP_pitch_%i", indexToPitch(index)));
        gPad->SetLogy();

        nStripsClusterBending_data[index]->Scale(1/nStripsClusterBending_data[index]->Integral());
        nStripsClusterBending_data[index]->SetLineColor(1);
        nStripsClusterBending_data[index]->SetMinimum(0.000001);
        nStripsClusterBending_data[index]->Draw("HIST ][ SAME");

        nStripsClusterBending_MC[index]->Scale(1/nStripsClusterBending_MC[index]->Integral());
        nStripsClusterBending_MC[index]->SetLineColor(2);
        nStripsClusterBending_MC[index]->SetMinimum(0.000001);
        nStripsClusterBending_MC[index]->Draw("HIST ][ SAME");

        auto legendBP = new TLegend();
        legendBP->AddEntry(nStripsClusterBending_data[index], "Data", "l");
        legendBP->AddEntry(nStripsClusterBending_MC[index], "O2sim", "l");
        legendBP->Draw("SAME");

        canvasBP->Write();

        // NBP
        TCanvas *canvasNBP = new TCanvas(Form("cluster_strip_size_NBP_pitch_%i", indexToPitch(index)),Form("cluster_strip_size_NBP_pitch_%i", indexToPitch(index)));
        gPad->SetLogy();

        nStripsClusterNonBending_data[index]->Scale(1/nStripsClusterNonBending_data[index]->Integral());
        nStripsClusterNonBending_data[index]->SetLineColor(1);
        nStripsClusterNonBending_data[index]->SetMinimum(0.000001);
        nStripsClusterNonBending_data[index]->Draw("HIST ][ SAME");

        nStripsClusterNonBending_MC[index]->Scale(1/nStripsClusterNonBending_MC[index]->Integral());
        nStripsClusterNonBending_MC[index]->SetLineColor(2);
        nStripsClusterNonBending_MC[index]->SetMinimum(0.000001);
        nStripsClusterNonBending_MC[index]->Draw("HIST ][ SAME");

        auto legendNBP = new TLegend();
        legendNBP->AddEntry(nStripsClusterNonBending_data[index], "Data", "l");
        legendNBP->AddEntry(nStripsClusterNonBending_MC[index], "O2sim", "l");
        legendNBP->Draw("SAME");

        canvasNBP->Write();
    }

    delete outFile;
}

int pitchToIndex(int pitch) {
    int index = -1;
    if (pitch == 4) index = 2;
    else index = pitch - 1;

    return index;
}

int indexToPitch(int index) {
    int pitch = -1;
    if (index == 2) pitch = 4;
    else pitch = index + 1;

    return pitch;
}