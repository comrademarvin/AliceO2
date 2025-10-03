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
#include "TF1.h"
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
#include "MIDBase/GeometryParameters.h"
#include "MIDClustering/PreClusterizer.h"
#include "MIDClustering/PreCluster.h"
#include "MIDSimulation/ChamberResponse.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

const int nPitches = 3;

struct clusterSizeHist {
    uint8_t deId;        ///< Detection element ID
    uint8_t cathode;     ///< Cathode
    int pitch;       ///< Strip pitch
    int chamber;     ///< Chamber number
    TH1F* clusterPosition; ///< Histogram of strip positions (x)
    int nClusters;   ///< Number of clusters
};
std::vector<clusterSizeHist*> initializeClusterSizeHist();

std::tuple<TFile *, TTreeReader *> loadData(const char *fileName, const char *treeName);
std::vector<o2::InteractionRecord> processMuonTracks(const char *fileName);
std::vector<o2::mid::PreCluster> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> muonTracksIR);
std::vector<clusterSizeHist*> processPreClusters(std::vector<o2::mid::PreCluster> preClustersMID);
void processClusterHist(std::vector<clusterSizeHist*> clusterSizeHistograms);
void fillStripPositions(clusterSizeHist* hist, int nStrips, int pitch);

int pitchToIndex(int pitch);
int indexToPitch(int index);

int main() {
    auto muonTracksIR = processMuonTracks("muontracks.root");

    //auto preClustersMID = processMIDdigits("mid-digits-decoded.root", "mid-reco.root", muonTracksIR);
    auto preClustersMID = processMIDdigits("middigits.root", "mid-reco.root", muonTracksIR); // for MC

    auto clusterSizeHist = processPreClusters(preClustersMID);

    processClusterHist(clusterSizeHist);

    return 0;
}

void fillStripPositions(clusterSizeHist* hist, int nStrips, int pitch) {
    double scaleFactor = o2::mid::geoparams::getStripUnitPitchSize(hist->chamber);

    double clusterPos; // position in mm
    div_t nStripsDiv2 = div(nStrips, 2);
    // number of strips even or odd?
    if (nStripsDiv2.rem == 0) { // even
        // assume symmetric
        for (int strip = 0; strip <= (nStripsDiv2.quot - 1); strip++) {
            clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * (static_cast<double>(strip) + 0.5);
            hist->clusterPosition->Fill(clusterPos);
            hist->clusterPosition->Fill(clusterPos); // fill twice for assumed symmetry of cluster position distribution
        }
        // assumme shifted asymmetry
        // //hist->clusterPosition->Fill(0.0); // middle strip
        // for (int strip = 0; strip <= nStripsDiv2.quot; strip++) {
        //     clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * static_cast<double>(strip);
        //     hist->clusterPosition->Fill(clusterPos);
        //     if (strip != nStripsDiv2.quot) hist->clusterPosition->Fill(clusterPos); // fill once for shifted asymmetry of cluster position distribution
        // }
    } else { // odd
        // assume symmetric
        // //hist->clusterPosition->Fill(0.0); // middle strip
        // for (int strip = 0; strip <= nStripsDiv2.quot; strip++) {
        //     clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * static_cast<double>(strip);
        //     hist->clusterPosition->Fill(clusterPos);
        //     hist->clusterPosition->Fill(clusterPos); // fill twice for assumed symmetry of cluster position distribution
        // }
        // assume shifted asymmetry
        for (int strip = 0; strip <= nStripsDiv2.quot; strip++) {
            clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * (static_cast<double>(strip) + 0.5);
            hist->clusterPosition->Fill(clusterPos);
            if (strip != nStripsDiv2.quot) hist->clusterPosition->Fill(clusterPos); // fill once for shifted asymmetry of cluster position distribution
        }
    }

    hist->nClusters++;
}

Double_t pdfFunc(Double_t *x, Double_t *par) // x = position (mm); par = {b, a0, a1, c0, c1, hv}
{
    Float_t xx = x[0];
    double hv = par[5];

    double a = par[2]*hv + par[1]; // a = a1*hv + a0
    double c = par[4]*hv + par[3]; // c = c1*hv + c0
    double b = par[0];

    return (1/(1+c))*((a/(a+pow(xx,b)))+c);
}

void processClusterHist(std::vector<clusterSizeHist*> clusterSizeHistograms) {
    if (clusterSizeHistograms.empty()) {
        std::cout << "No cluster size histograms available." << std::endl;
        return;
    }

    // add histograms with same cathode and deId together
    std::map<std::tuple<int, int>, std::pair<TH1F*, int>> combinedHists; // key: (cathode, deId), value: pair of combined histogram and cluster count
    for (auto& hist : clusterSizeHistograms) {
        auto key = std::make_tuple(hist->cathode, hist->deId);
        if (combinedHists.find(key) == combinedHists.end()) {
            combinedHists[key] = std::make_pair(
                (TH1F*)hist->clusterPosition->Clone(Form("combined_cathode%i_deId%i", hist->cathode, hist->deId)),
                hist->nClusters
            );
            combinedHists[key].first->SetTitle(Form("Combined Strip Position (mm) for Cathode %i, DE %i", hist->cathode, hist->deId));
        } else {
            combinedHists[key].first->Add(hist->clusterPosition);
            combinedHists[key].second += hist->nClusters;
        }
    }

    auto testKey = std::make_tuple(0, 3); // using cathode = 0, deId = 0 as a test
    (combinedHists[testKey].first)->Scale(1.0 / (2*static_cast<float>(combinedHists[testKey].second))); // normalize by number of clusters

    // pick test histogram for now (to fit parameters A,B,C)
    auto baseClusterHist = clusterSizeHistograms[343]; // index 240: deId = 48, cathode = 0, pitch = 1 | index 343: deId = 68, cathode = 1, pitch = 2
    //baseClusterHist->clusterPosition->Scale(1.0 / (2*static_cast<float>(baseClusterHist->nClusters))); // normalize by number of clusters
    baseClusterHist->clusterPosition->Scale(1 / (2*baseClusterHist->clusterPosition->Integral())); // normalize by area + for overcounting due to symmetry assumption


    // create current O2 PDF for the chosen hist (using default HV for now)
    o2::mid::ChamberResponse chamberResp = o2::mid::createDefaultChamberResponse();

    const int nPoints = 100; // Number of points in the grid
    double distances[nPoints];
    double firedProbabilities[nPoints];

    for (int i = 0; i < nPoints; ++i) {
        distances[i] = i; // Distance values from 0 to 100 in steps of 1 in mm
        firedProbabilities[i] = chamberResp.getFiredProbability(distances[i]/10, baseClusterHist->cathode, baseClusterHist->deId, 0.0); // need to pass in cm
    }

    TGraph* graph = new TGraph(nPoints, distances, firedProbabilities);

    // access HV values from the CCDB
    //auto ccdb = o2::ccdb::BasicCCDBManager::getInstance();

    // define my own fit function
    auto clusterPDF = new TF1("clusterPDF", pdfFunc, 0, 100, 6);
    clusterPDF->SetParNames("b", "a0", "a1", "c0", "c1", "hv");
    clusterPDF->SetParameters(1.97, -52.70, 6.089, -0.5e-3, 8.3e-4, 9.6); // initial parameters
    clusterPDF->FixParameter(5, 9.6); // Fix the high-voltage parameter

    //baseClusterHist->clusterPosition->Fit("clusterPDF", "R", "", 0, 100); // fit in range of 0-100 mm

    // read out histograms
    auto outFile = new TFile("cluster_hist_fitting.root", "RECREATE");

    // plot first hist and current PDF together as a sanity check
    TCanvas* canvasCheckFirst = new TCanvas("FiredProbabilityCanvas", "Fired Probability vs Distance", 800, 600);
    gPad->SetLogy();
    baseClusterHist->clusterPosition->SetLineColor(kBlack);
    baseClusterHist->clusterPosition->SetMaximum(1.0);
    baseClusterHist->clusterPosition->SetMinimum(0.0001);
    baseClusterHist->clusterPosition->Draw("SAME");

    graph->Draw("SAME");
    graph->SetLineColor(kRed);
    graph->SetLineWidth(2);

    canvasCheckFirst->Write();

    // plot test combined histogram with fit function for different b-parameter values
    TCanvas* canvasCheckCombine = new TCanvas("FiredProbabilityCanvasCombinedB", "Fired Probability vs Distance for cathode = 0 and pitch = 1", 800, 600);
    gPad->SetLogy();
    (combinedHists[testKey].first)->SetLineColor(kBlack);
    (combinedHists[testKey].first)->SetMaximum(1.0);
    (combinedHists[testKey].first)->SetMinimum(0.0001);
    (combinedHists[testKey].first)->Draw("SAME");

    double b_parameter_values[] = {1.97, 2.22, 2.47, 2.97}; // plot the PDF for different b-parameter values
    int colors[] = {kRed, kBlue, kGreen, kMagenta};
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetHeader("b-Parameter Values", "C");

    for (int i = 0; i < 4; ++i) {
        auto tempPDF = new TF1("tempPDF", pdfFunc, 0, 100, 6);
        tempPDF->SetParNames("b", "a0", "a1", "c0", "c1", "hv");
        tempPDF->SetParameters(b_parameter_values[i], -52.70, 6.089, -0.5e-3, 8.3e-4, 9.6); // initial parameters
        tempPDF->SetLineColor(colors[i]);               
        tempPDF->SetLineWidth(2);                      
        tempPDF->Draw("SAME"); // Overlay subsequent functions
        legend->AddEntry(tempPDF, Form("b = %.2f", b_parameter_values[i]), "l");
    }
    legend->Draw("SAME");

    canvasCheckCombine->Write();

    // plot test combined histogram with fit function for different HV values
    TCanvas* canvasCheckCombineHV = new TCanvas("FiredProbabilityCanvasCombinedHV", "Fired Probability vs Distance for cathode = 0 and pitch = 1", 800, 600);
    gPad->SetLogy();
    (combinedHists[testKey].first)->SetLineColor(kBlack);
    (combinedHists[testKey].first)->SetMaximum(1.0);
    (combinedHists[testKey].first)->SetMinimum(0.0001);
    (combinedHists[testKey].first)->Draw("SAME");

    double hv_values[] = {9.2, 9.4, 9.6, 9.8}; // plot the PDF for different hv values
    TLegend* legendHV = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendHV->SetHeader("HV Values (mV)", "C");

    for (int i = 0; i < 4; ++i) {
        auto tempPDF = new TF1("tempPDF", pdfFunc, 0, 100, 6);
        tempPDF->SetParNames("b", "a0", "a1", "c0", "c1", "hv");
        tempPDF->SetParameters(1.97, -52.70, 6.089, -0.5e-3, 8.3e-4, hv_values[i]); // initial parameters
        tempPDF->SetLineColor(colors[i]);               
        tempPDF->SetLineWidth(2);                      
        tempPDF->Draw("SAME"); // Overlay subsequent functions
        legendHV->AddEntry(tempPDF, Form("HV = %.1f", hv_values[i]), "l");
    }
    legendHV->Draw("SAME");

    canvasCheckCombineHV->Write();

    int histIndex = 0;
    for (auto& hist : clusterSizeHistograms) {
        hist->clusterPosition->Write(Form("strip_position_de%i_cathode%i_pitch%i_index%i", hist->deId, hist->cathode, hist->pitch, histIndex));
        histIndex++;
    }

    // for (const auto& [key, value] : combinedHists) {
    //     int cathode = std::get<0>(key);
    //     int deId = std::get<1>(key);
    //     value.first->Write(Form("combined_strip_position_de%i_cathode%i", deId, cathode));
    // }

    delete outFile;
}

std::vector<clusterSizeHist*> processPreClusters(std::vector<o2::mid::PreCluster> preClustersMID) {
    // mapping object to extract strip size per column
    o2::mid::Mapping* mapper = new o2::mid::Mapping();

    // strip size distribution for each chamber
    std::vector<TH1D*> nStripsClusterBending(nPitches);
    std::vector<TH1D*> nStripsClusterNonBending(nPitches);

    for (int index = 0; index < nPitches; index++) {
        nStripsClusterBending[index] = new TH1D(Form("cluster_strip_size_bending_%i", indexToPitch(index)), Form("Number of Strips in MID PreClusters for Bending Plane (pitch = %i);nStrips/PreCluster;count", indexToPitch(index)), 30, 0, 30);
        nStripsClusterNonBending[index] = new TH1D(Form("cluster_strip_size_nonbending_%i", indexToPitch(index)), Form("Number of Strips in MID PreClusters for Non-Bending Plane (pitch = %i);nStrips/PreCluster;count", indexToPitch(index)), 30, 0, 30);
    }

    // vector of cluster size histograms for fitting
    std::vector<clusterSizeHist*> clusterSizeHistograms = initializeClusterSizeHist();

    // loop over preClusters
    for (auto pc : preClustersMID) {
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

        // fill cluster size histograms
        if (strip_pitch_first == strip_pitch_last) { // only when the strip pitch is the same for the full cluster
            for (auto& hist : clusterSizeHistograms) {
                if (hist->deId == pc.deId && hist->cathode == pc.cathode && hist->pitch == strip_pitch_first) {
                    fillStripPositions(hist, nStrips, strip_pitch_first);
                    break;
                }
            }
        }
    }

    // read out histograms
    auto outFile = new TFile("cluster_size.root", "RECREATE");

    for (int index = 0; index < nPitches; index++) {
        nStripsClusterBending[index]->Write();
        nStripsClusterNonBending[index]->Write();
    }

    delete outFile;

    return clusterSizeHistograms; // to be used for fitting
}

std::vector<o2::mid::PreCluster> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> muonTracksIR) {
    // output histograms
    TH1D* nMuonTracksROF = new TH1D("muon_track_count_ROF", "Number of MCH+MID matched tracks per ROF;nTracks/ROF;count", 5, 0, 5);
    TH1D* nMIDTracksROF = new TH1D("mid_track_count_ROF", "Number of MID tracks per ROF;nTracks/ROF;count", 70, 0, 70);

    // read in the MID track and digit infomation
    //auto [digitFile, digitReader] = loadData(fileMIDdigits, "middigits");
    auto [digitFile, digitReader] = loadData(fileMIDdigits, "o2sim"); // for MC

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
                nMIDTracksROF->Fill(nTracksMID);
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

    nMIDTracksROF->Write();
    nMuonTracksROF->Write();

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

std::vector<clusterSizeHist*> initializeClusterSizeHist() {
    // initialize cluster size histograms
    std::vector<clusterSizeHist*> clusterSizeHistograms;

    // loop over all MID deIds
    for (int deId = 0; deId < o2::mid::detparams::NDetectionElements; deId++) {
        for (int cathode = 0; cathode < 2; cathode++) { // 0: bending, 1: non-bending
            for (int pitch = 0; pitch < nPitches; pitch++) { // strip pitches
                if (pitch == 0 && cathode == 1) continue; // skip non-bending plane for pitch 1
                auto hist = new clusterSizeHist();
                hist->deId = deId;
                hist->cathode = cathode;
                hist->pitch = indexToPitch(pitch);
                hist->chamber = o2::mid::detparams::getChamber(deId);
                hist->nClusters = 0;
                hist->clusterPosition = new TH1F(Form("strip_position_de%i_cathode%i_pitch%i", deId, cathode, indexToPitch(pitch)), 
                                                Form("Strip Position (mm) for DE %i, Cathode %i, Pitch %i", deId, cathode, indexToPitch(pitch)), 50, 0, 100);
                clusterSizeHistograms.push_back(hist);
            }
        }
    }

    return clusterSizeHistograms;
}