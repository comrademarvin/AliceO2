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

std::tuple<TFile *, TTreeReader *> loadData(const char *fileName, const char *treeName);

int main() {
    std::cout << "Hello World, Cluster Size Task Exists!" << std::endl;

    // read in the track and digit infomation
    auto [digitFile, digitReader] = loadData("mid-digits-decoded.root", "middigits");
    auto [recoFile, recoReader] = loadData("mid-reco.root", "midreco");

    TTreeReaderValue<std::vector<o2::mid::ColumnData>> digits{*digitReader, "MIDDigit"};
    TTreeReaderValue<std::vector<o2::mid::ROFRecord>> digitRofs{*digitReader, "MIDROFRecords"};
    // TTreeReaderValue<std::vector<o2::mid::Cluster>> clusters{*recoReader, "MIDTrackCluster"};
    // TTreeReaderValue<std::vector<o2::mid::ROFRecord>> clusterRofs{*recoReader, "MIDTrackClusterROF"};
    TTreeReaderValue<std::vector<o2::mid::Track>> tracks{*recoReader, "MIDTrack"};
    TTreeReaderValue<std::vector<o2::mid::ROFRecord>> trackRofs{*recoReader, "MIDTrackROF"};

    if (digitReader->GetEntries() != recoReader->GetEntries()) // same number of TFs
    {
        std::cout << "Error: the digit and cluster readers do not contain the same number of TFs";
        exit(-1);
    }

    // output histograms
    TH1D* nTracksROF = new TH1D("ROF_track_count", "Number of Tracks per ROF;nTracks/ROF;count", 50, 0, 50);
    TH1D* nClusterStrips = new TH1D("cluster_strips_size", "Number of Strips in MID Clusters;nStrips/Cluster;count", 30, 0, 30);

    // re-run preclusterizer on event digits
    o2::mid::PreClusterizer preClusterizer;

    // itterate over entries, where each entry is a TF
    int entriesCount = 0;
    while (digitReader->Next() && recoReader->Next()) {
        //auto clusterRofIt = (*clusterRofs).begin(); // itterator over cluster ROFs for one TF
        auto trackRofIt = (*trackRofs).begin(); // itterator over track ROFs for one TF

        gsl::span<o2::mid::ColumnData> sdigits(*digits); // all digits for the TF

        // itterate over digit ROFs for one TF
        for (auto digitRofIt = (*digitRofs).begin(), digitEnd = (*digitRofs).end(); digitRofIt != digitEnd; ++digitRofIt)
        {
            auto nTracks = trackRofIt->nEntries; // number of tracks for the ROF
            if (nTracks > 0) { // check whether there are any tracks for the ROF
                //std::cout << "nTracks of ROF: " << nTracks << std::endl;
                nTracksROF->Fill(nTracks);

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
            ++trackRofIt;
        }
        
        entriesCount++;
    }

    std::cout << "Number of entries (TFs): " << entriesCount << std::endl;

    // read out histograms
    auto outFile = new TFile("cluster_size.root", "RECREATE");

    TCanvas* nTracksCanvas = new TCanvas("ROF_tracks_count", "ROF_tracks_count");
    nTracksROF->Draw();
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