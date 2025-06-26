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
#include "TLegend.h"
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
std::vector<o2::mid::PreCluster> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> muonTracksIR);
void processPreClusters(std::vector<o2::mid::PreCluster> preClustersMID);
int pitchToIndex(int pitch);
int indexToPitch(int index);

int main() {
    auto muonTracksIR = processMuonTracks("muontracks.root");

    auto preClustersMID = processMIDdigits("mid-digits-decoded.root", "mid-reco.root", muonTracksIR);
    //auto preClustersMID = processMIDdigits("middigits.root", "mid-reco.root", muonTracksIR); // for MC

    processPreClusters(preClustersMID);

    return 0;
}

void processPreClusters(std::vector<o2::mid::PreCluster> preClustersMID) {
    // mapping object to extract strip size per column
    o2::mid::Mapping* mapper = new o2::mid::Mapping();

    // strip size distribution for each chamber
    const int nPitches = 3;
    std::vector<TH1D*> nStripsClusterBending(nPitches);
    std::vector<TH1D*> nStripsClusterNonBending(nPitches);

    for (int index = 0; index < nPitches; index++) {
        nStripsClusterBending[index] = new TH1D(Form("cluster_strip_size_bending_%i", indexToPitch(index)), Form("Number of Strips in MID PreClusters for Bending Plane (pitch = %i);nStrips/PreCluster;count", indexToPitch(index)), 30, 0, 30);
        nStripsClusterNonBending[index] = new TH1D(Form("cluster_strip_size_nonbending_%i", indexToPitch(index)), Form("Number of Strips in MID PreClusters for Non-Bending Plane (pitch = %i);nStrips/PreCluster;count", indexToPitch(index)), 30, 0, 30);
    }

    for (auto pc : preClustersMID) {
        // check the chamber of the clusters
        // int clusterChamber = o2::mid::detparams::getChamber(pc.deId);

        // find the strip pitch 
        o2::mid::MpArea mpArea_first = mapper->stripByLocation(pc.firstStrip, pc.cathode, pc.firstLine, pc.firstColumn, pc.deId);
        o2::mid::MpArea mpArea_last = mapper->stripByLocation(pc.lastStrip, pc.cathode, pc.lastLine, pc.lastColumn, pc.deId);
        int strip_pitch_first, strip_pitch_last;

        // determine strip size 
        int nStrips;
        if (pc.cathode == 0) { // is the cluster in the bending plane?
            strip_pitch_first = static_cast<int>(2*mpArea_first.getHalfSizeY());
            strip_pitch_last = static_cast<int>(2*mpArea_last.getHalfSizeY());

            nStrips = (int)(pc.lastStrip - pc.firstStrip + 16 * (pc.lastLine - pc.firstLine)) + 1;

            if (strip_pitch_first == strip_pitch_last) nStripsClusterBending[pitchToIndex(strip_pitch_first)]->Fill(nStrips);
        } 
        else {
            strip_pitch_first = static_cast<int>(2*mpArea_first.getHalfSizeX());
            strip_pitch_last = static_cast<int>(2*mpArea_last.getHalfSizeX());

            nStrips = (pc.lastStrip - pc.firstStrip) + 1;
            for (int column = pc.firstColumn; column < pc.lastColumn; column++) {
                nStrips += mapper->getNStripsNBP(column, pc.deId);
            }

            if (strip_pitch_first == strip_pitch_last) nStripsClusterNonBending[pitchToIndex(strip_pitch_first)]->Fill(nStrips);
        }
    }

    // read out histograms
    auto outFile = new TFile("cluster_size.root", "RECREATE");

    for (int index = 0; index < nPitches; index++) {
        nStripsClusterBending[index]->Write();
        nStripsClusterNonBending[index]->Write();
    }

    delete outFile;
}

std::vector<o2::mid::PreCluster> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> muonTracksIR) {
    // output histograms
    TH1D* nMuonTracksROF = new TH1D("muon_track_count_ROF", "Number of MCH+MID matched tracks per ROF;nTracks/ROF;count", 10, 0, 10);

    // read in the MID track and digit infomation
    auto [digitFile, digitReader] = loadData(fileMIDdigits, "middigits");
    //auto [digitFile, digitReader] = loadData(fileMIDdigits, "o2sim"); // for MC

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
            int muonTrackCount = 0;
            auto nTracksMID = trackRofItMID->nEntries; // number of MID tracks for the ROF
            if (nTracksMID > 0) { // first check whether there are any MID tracks for the ROF
                // secondly check whether there are any matched muon tracks for the ROF
                auto rofIR = digitRofIt->interactionRecord;
                for (auto trackIR : muonTracksIR) {
                    if (rofIR == trackIR) muonTrackCount++;
                }
            }

            if (muonTrackCount > 0) {
                nMuonTracksROF->Fill(muonTrackCount);

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

    // read out histograms
    auto outFile = new TFile("track_processing.root", "RECREATE");

    TCanvas* nMuonTracksCanvas = new TCanvas("muon_tracks_ROF", "muon_tracks_ROF");
    nMuonTracksROF->Draw();
    nMuonTracksCanvas->Write();

    delete outFile;

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