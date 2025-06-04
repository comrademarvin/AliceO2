#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <tuple>
#include <gsl/span>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TGeoManager.h"
#include "TPolyLine.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsMID/ROFRecord.h"
#include "DataFormatsMID/ColumnData.h"
#include "DataFormatsMID/Cluster.h"
#include "DataFormatsMID/Track.h"
#include "MIDBase/DetectorParameters.h"
#include "MIDBase/Mapping.h"
#include "MIDBase/GeometryTransformer.h"
#include "MIDClustering/PreClusterizer.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"

std::tuple<TFile *, TTreeReader *> loadData(const char *fileName, const char *treeName);

int main() {
    // read in the MID/Muon track and digit infomation
    auto [digitFile, digitReader] = loadData("mid-digits-decoded.root", "middigits");
    auto [recoFileMID, recoReaderMID] = loadData("mid-reco.root", "midreco");
    auto [recoFileMCH, recoReaderMCH] = loadData("mchtracks.root", "mchtracks");
    auto [recoFileMuon, recoReaderMuon] = loadData("muontracks.root", "muontracks");

    // MID
    TTreeReaderValue<std::vector<o2::mid::ColumnData>> digits{*digitReader, "MIDDigit"};
    TTreeReaderValue<std::vector<o2::mid::ROFRecord>> digitRofs{*digitReader, "MIDROFRecords"};
    // TTreeReaderValue<std::vector<o2::mid::Cluster>> clusters{*recoReaderMID, "MIDTrackCluster"};
    // TTreeReaderValue<std::vector<o2::mid::ROFRecord>> clusterRofs{*recoReaderMID, "MIDTrackClusterROF"};
    TTreeReaderValue<std::vector<o2::mid::Track>> tracksMID{*recoReaderMID, "MIDTrack"};
    TTreeReaderValue<std::vector<o2::mid::ROFRecord>> trackRofsMID{*recoReaderMID, "MIDTrackROF"};

    // MCH
    TTreeReaderValue<std::vector<o2::mch::TrackMCH>> tracksMCH{*recoReaderMCH, "tracks"};
    TTreeReaderValue<std::vector<o2::mch::ROFRecord>> trackRofsMCH{*recoReaderMCH, "trackrofs"};

    if ((digitReader->GetEntries() != recoReaderMID->GetEntries()) && (digitReader->GetEntries() != recoReaderMCH->GetEntries())) // same number of TFs
    {
        std::cout << "Error: the digit and cluster readers do not contain the same number of TFs";
        exit(-1);
    }

    // output histograms
    TH1D* nTracksROF_MID = new TH1D("ROF_track_count", "Number of MID Tracks per ROF;nTracks/ROF;count", 50, 0, 50);
    TH1D* nClusterStrips = new TH1D("cluster_strips_size", "Number of Strips in MID Clusters;nStrips/Cluster;count", 30, 0, 30);

    // re-run preclusterizer on event digits
    o2::mid::PreClusterizer preClusterizer;

    // itterate over entries, where each entry is a TF
    int entriesCount = 0;
    while (digitReader->Next() && recoReaderMID->Next() && recoReaderMCH->Next()) {
        //auto clusterRofIt = (*clusterRofs).begin(); // itterator over cluster ROFs for one TF
        auto trackRofItMID = (*trackRofsMID).begin(); // itterator over MID track ROFs for one TF
        auto trackRofItMCH = (*trackRofsMCH).begin(); // itterator over MCH track ROFs for one TF

        gsl::span<o2::mid::ColumnData> sdigits(*digits); // all digits for the TF

        // itterate over digit ROFs for one TF
        for (auto digitRofIt = (*digitRofs).begin(), digitEnd = (*digitRofs).end(); digitRofIt != digitEnd; ++digitRofIt)
        {
            auto nTracksMID = trackRofItMID->nEntries; // number of MID tracks for the ROF
            auto nTracksMCH = trackRofItMCH->nEntries; // number of MCH tracks for the ROF
            if ((nTracksMID > 0) && (nTracksMCH > 0)) { // check whether there are any tracks for the ROF
                nTracksROF_MID->Fill(nTracksMID);

                // subspan of digits for the ROF where there are tracks
                auto eventDigits = sdigits.subspan(digitRofIt->firstEntry, digitRofIt->nEntries);

                // run pre-clusterizer on the event digits
                preClusterizer.process(eventDigits);
                auto preClusters = preClusterizer.getPreClusters();

                // itterate over pre-clusters
                for (auto& pc : preClusters) {
                    // number of strips in preCluster (from PreClusterHelper)
                    int nStripsInBetween = pc.lastStrip - pc.firstStrip + 16 * (pc.lastLine - pc.firstLine);

                    nClusterStrips->Fill(nStripsInBetween);
                }
            }

            //++clusterRofIt;
            ++trackRofItMID;
        }
        
        entriesCount++;
    }

    std::cout << "Number of entries (TFs): " << entriesCount << std::endl;

    // read out histograms
    auto outFile = new TFile("cluster_size.root", "RECREATE");

    TCanvas* nTracksCanvas = new TCanvas("ROF_tracks_count", "ROF_tracks_count");
    nTracksROF_MID->Draw();
    nTracksCanvas->Write();

    TCanvas* nStripsCanvas = new TCanvas("cluster_strip_size", "cluster_strip_size");
    nClusterStrips->Draw();
    nStripsCanvas->Write();

    delete outFile;

    return 0;
}

std::tuple<TFile *, TTreeReader *> loadData(const char *fileName, const char *treeName)
{
    /// open the input file and get the intput tree

    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cout << "opening file " << fileName << " failed";
        exit(-1);
    }

    TTreeReader *tr = new TTreeReader(treeName, file);
    if (tr->IsZombie())
    {
        std::cout << "tree " << treeName << " not found";
        exit(-1);
    }

    return std::make_tuple(file, tr);
}