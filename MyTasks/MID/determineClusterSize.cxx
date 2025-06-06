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
#include "MIDClustering/PreCluster.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

std::tuple<TFile *, TTreeReader *> loadData(const char *fileName, const char *treeName);
std::vector<o2::InteractionRecord> processMuonTracks(const char *fileName);
std::vector<o2::mid::PreCluster> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> tracksIR);
void processPreClusters(std::vector<o2::mid::PreCluster> preClustersMID);

int main() {
    auto muonTracksIR = processMuonTracks("muontracks.root");

    auto preClustersMID = processMIDdigits("mid-digits-decoded.root", "mid-reco.root", muonTracksIR);

    processPreClusters(preClustersMID);

    return 0;
}

void processPreClusters(std::vector<o2::mid::PreCluster> preClustersMID) {
    TH1D* nStripsCluster = new TH1D("cluster_strip_size", "Number of Strips in MID PreClusters;nStrips/PreCluster;count", 30, 0, 30);

    for (auto pc : preClustersMID) {
        int nStripsInBetween;
        if (pc.cathode == 0) { // is the cluster in the bending plane?
            nStripsInBetween = pc.lastStrip - pc.firstStrip + 16 * (pc.lastLine - pc.firstLine);
        } else {
            nStripsInBetween = pc.lastStrip - pc.firstStrip;
        }

        nStripsCluster->Fill(nStripsInBetween);
    }

    // read out histograms
    auto outFile = new TFile("cluster_size.root", "RECREATE");

    TCanvas* nStripsCanvas = new TCanvas("cluster_strip_size", "cluster_strip_size");
    nStripsCluster->Draw();
    nStripsCanvas->Write();

    delete outFile;
}

std::vector<o2::mid::PreCluster> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> tracksIR) {
    // read in the MID track and digit infomation
    auto [digitFile, digitReader] = loadData(fileMIDdigits, "middigits");
    auto [recoFileMID, recoReaderMID] = loadData(fileMIDtracks, "midreco");

    // branches of interest
    TTreeReaderValue<std::vector<o2::mid::ColumnData>> digits{*digitReader, "MIDDigit"};
    TTreeReaderValue<std::vector<o2::mid::ROFRecord>> digitRofs{*digitReader, "MIDROFRecords"};
    TTreeReaderValue<std::vector<o2::mid::Track>> tracksMID{*recoReaderMID, "MIDTrack"};
    TTreeReaderValue<std::vector<o2::mid::ROFRecord>> trackRofsMID{*recoReaderMID, "MIDTrackROF"};

    // check if MID digits and tracks have the same number of TFs
    if (digitReader->GetEntries() != recoReaderMID->GetEntries()) // same number of TFs
    {
        std::cout << "Error: the digit and cluster readers do not contain the same number of TFs";
        exit(-1);
    }

    // preclusterizer for MID event digits
    o2::mid::PreClusterizer preClusterizer;

    // vector of MID preClusters
    std::vector<o2::mid::PreCluster> midPreClusters;

    // itterate over entries, where each entry is a TF
    int timeframeCounter = 0;
    int selectedROFcounter = 0;
    while (digitReader->Next() && recoReaderMID->Next()) {
        auto trackRofItMID = (*trackRofsMID).begin(); // iterator over MID track ROFs for one TF

        gsl::span<o2::mid::ColumnData> sdigits(*digits); // all MID digits for the TF

        // itterate over MID digit ROFs for one TF
        for (auto digitRofIt = (*digitRofs).begin(), digitEnd = (*digitRofs).end(); digitRofIt != digitEnd; ++digitRofIt) {
            bool foundMuonTrack = false;
            auto nTracksMID = trackRofItMID->nEntries; // number of MID tracks for the ROF
            if (nTracksMID > 0) { // first check whether there are any MID tracks for the ROF
                // secondly check whether there are any matched muon tracks for the ROF
                auto rofIR = digitRofIt->interactionRecord;
                for (auto trackIR : tracksIR) {
                    if (rofIR == trackIR) {
                        foundMuonTrack = true;
                        break;
                    }
                }
            }

            if (foundMuonTrack) {
                // subspan of digits for the ROF where there are tracks
                auto eventDigits = sdigits.subspan(digitRofIt->firstEntry, digitRofIt->nEntries);

                // run pre-clusterizer on the event digits
                preClusterizer.process(eventDigits);
                auto preClusters = preClusterizer.getPreClusters();
                for (auto& pc : preClusters) midPreClusters.push_back(pc);

                selectedROFcounter++;
            }

            ++trackRofItMID;
        }

        timeframeCounter++;
    }

    std::cout << "Number of TF in MID track/digits file: " << timeframeCounter << std::endl;
    std::cout << "Number of selected digit ROFs: " << selectedROFcounter << std::endl;

    return midPreClusters;
}

std::vector<o2::InteractionRecord> processMuonTracks(const char *fileName) {
    // read in muon tracks from file
    auto [recoFileMuon, recoReaderMuon] = loadData(fileName, "o2sim");
    TTreeReaderValue<std::vector<o2::dataformats::TrackMCHMID>> muonTracks{*recoReaderMuon, "tracks"};

    // declare vector to store interaction records of muon tracks (MCH+MID match)
    std::vector<o2::InteractionRecord> muonTracksIR; 

    // iterate timeframe entries
    int timeframeCounter = 0;
    while (recoReaderMuon->Next()) {
        // iterate over muon tracks for one TF
        for (auto muonTrackIt = (*muonTracks).begin(), muonTrackEnd = (*muonTracks).end(); muonTrackIt != muonTrackEnd; ++muonTrackIt) {
            muonTracksIR.push_back(muonTrackIt->getIR());
        }

        timeframeCounter++;
    }

    std::cout << "Number of TF in muon track file: " << timeframeCounter << std::endl;
    std::cout << "Number of muon tracks: " << muonTracksIR.size() << std::endl;

    return muonTracksIR;
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