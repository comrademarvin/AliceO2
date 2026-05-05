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
#include "TLine.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TGeoManager.h"
#include "TPolyLine.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
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
#include "MIDSimulation/ChamberResponseParams.h"
#include "MIDConditions/DCSNamer.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

// CCDB HV data formats
using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;
using DPMAP = std::unordered_map<DPID, std::vector<DPVAL>>;

const int nPitches = 3;
int pitchToIndex(int pitch);
int indexToPitch(int index);

struct cluster {
    int rofID;
    int cathode;
    int deId;
    int pitch;
    int nStrips;
};

struct clusterSizeHist {
    uint8_t deId;     ///< Detection element ID
    uint8_t cathode;     ///< Cathode (Bending = 0; Non-bending = 1)
    int pitch;       ///< Strip pitch
    TH1D* clusterSize; ///< Histogram of strip positions (x)
};
std::vector<clusterSizeHist*> initializeClusterSizeHist();

struct clusterPosHist {
    uint8_t deId;        ///< Detection element ID
    uint8_t cathode;     ///< Cathode
    int pitch;       ///< Strip pitch
    int chamber;     ///< Chamber number
    TH1F* clusterPosition; ///< Histogram of strip positions (x)
    int nClusters;   ///< Number of clusters
};
std::vector<clusterPosHist*> initializeClusterPosHist();

struct fitParams {
    int deId;
    int cathode;
    int pitch;
    int nBins;
    int nEntries;
    std::array<double, 5> params; // b, a0, a1, c0, c1
};

std::tuple<TFile *, TTreeReader *> loadData(const char *fileName, const char *treeName);
std::vector<o2::InteractionRecord> processMuonTracks(const char *fileName);
std::vector<cluster*> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> muonTracksIR, bool isMC, bool removeMultColumns = true);
cluster* processPreCluster(int rofID, o2::mid::PreCluster pc, o2::mid::Mapping* mapper);
std::vector<clusterPosHist*> processClustersPosition(std::vector<cluster*> midClusters);
void processClustersSize(std::vector<cluster*> midClusters);
void processClusterPosHist(std::vector<clusterPosHist*> clusterPosHistograms, bool isMC);
void processFitParams(std::vector<fitParams*> fittedParams, const std::vector<double>& currentBParams, const std::array<double, 2>& currentAParams, const std::array<double, 2>& currentCParams);
void fillStripPositions(clusterPosHist* hist, int nStrips, int pitch);
DPMAP* accessHVObjectCCDB(int runNumber);
float accessHVperDE(DPMAP* HV_map, int deId);
float sum(float s, o2::dcs::DataPointValue v);

int main(int argc, char* argv[]) {
    // Check if the 'isMC' parameter is passed
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <isMC: 0 or 1>" << std::endl;
        return -1;
    }
    bool isMC = std::stoi(argv[1]) != 0;

    // run on multiple batches of data and combine the results for cluster size and position analysis
    const int numBatches = 1;
    //std::string recoFilesDirectory[numBatches] = {"/home/stephan/MID_cluster/run_558801_batch_1", "/home/stephan/MID_cluster/run_558801_batch_2"};
    std::string recoFilesDirectory[numBatches] = {"/home/stephan/sims/mu_boxgen_fitted_100k"};

    std::vector<cluster*> combinedClusters;

    for (int i = 0; i < numBatches; ++i) {
        std::string muonTracksFile = recoFilesDirectory[i] + "/muontracks.root";
        std::string midDigitsFile = recoFilesDirectory[i] + (isMC ? "/middigits.root" : "/mid-digits-decoded.root");
        std::string midRecoFile = recoFilesDirectory[i] + "/mid-reco.root";

        auto muonTracksIR = processMuonTracks(muonTracksFile.c_str());
        auto midClusters = processMIDdigits(midDigitsFile.c_str(), midRecoFile.c_str(), muonTracksIR, isMC, true);

        std::cout << "Number of MID clusters in batch " << i+1 << ": " << midClusters.size() << std::endl;

        combinedClusters.insert(combinedClusters.end(), midClusters.begin(), midClusters.end());
    }

    std::cout << "Total number of MID clusters across all batches: " << combinedClusters.size() << std::endl;

    processClustersSize(combinedClusters);

    auto clusterPosHist = processClustersPosition(combinedClusters);

    processClusterPosHist(clusterPosHist, isMC);

    return 0;
}

void fillStripPositions(clusterPosHist* hist, int nStrips, int pitch) {
    double scaleFactor = o2::mid::geoparams::getStripUnitPitchSize(hist->chamber);

    double clusterPos; // position in mm
    // second approach
    for (int strip = 0; strip < nStrips; strip++) {
        int bin = strip; // Bins are 1-indexed in ROOT histograms
        double binCenter = hist->clusterPosition->GetBinCenter(bin);
        hist->clusterPosition->Fill(binCenter);
    }

    hist->nClusters++;

    // first approach    
    // div_t nStripsDiv2 = div(nStrips, 2);
    // // number of strips even or odd?
    // if (nStripsDiv2.rem == 0) { // even
    //     // assume symmetric
    //     for (int strip = 0; strip <= (nStripsDiv2.quot - 1); strip++) {
    //         clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * (static_cast<double>(strip) + 0.5);
    //         hist->clusterPosition->Fill(clusterPos);
    //         hist->clusterPosition->Fill(clusterPos); // fill twice for assumed symmetry of cluster position distribution
    //     }
    //     // assumme shifted asymmetry
    //     //hist->clusterPosition->Fill(0.0); // middle strip
    //     // for (int strip = 0; strip <= nStripsDiv2.quot; strip++) {
    //     //     clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * static_cast<double>(strip);
    //     //     hist->clusterPosition->Fill(clusterPos);
    //     //     if (strip != nStripsDiv2.quot) hist->clusterPosition->Fill(clusterPos); // fill once for shifted asymmetry of cluster position distribution
    //     // }
    // } else { // odd
    //     // assume symmetric
    //     //hist->clusterPosition->Fill(0.0); // middle strip
    //     for (int strip = 0; strip <= nStripsDiv2.quot; strip++) {
    //         clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * static_cast<double>(strip);
    //         hist->clusterPosition->Fill(clusterPos);
    //         hist->clusterPosition->Fill(clusterPos); // fill twice for assumed symmetry of cluster position distribution
    //     }
    //     // // assume shifted asymmetry
    //     // for (int strip = 0; strip <= nStripsDiv2.quot; strip++) {
    //     //     clusterPos = (static_cast<double>(pitch) * scaleFactor * 10) * (static_cast<double>(strip) + 0.5);
    //     //     hist->clusterPosition->Fill(clusterPos);
    //     //     if (strip != nStripsDiv2.quot) hist->clusterPosition->Fill(clusterPos); // fill once for shifted asymmetry of cluster position distribution
    //     // }
    // }
}

Double_t pdfFunc(Double_t *x, Double_t *par) // x = position (mm); par = {b, a0, a1, c0, c1, hv}
{
    Float_t xx = x[0];
    double hv = par[5];
    double costheta = std::cos(par[6]);

    double a = par[2]*hv + par[1]; // a = a1*hv + a0
    double c = par[4]*hv + par[3]; // c = c1*hv + c0
    double b = par[0];

    return (1/(1+c)) * ((a / (a + (pow(xx,b)*costheta))) + c);
}

void processClusterPosHist(std::vector<clusterPosHist*> clusterPosHistograms, bool isMC) {
    if (clusterPosHistograms.empty()) {
        std::cout << "No cluster size histograms available." << std::endl;
        return;
    }

    // MID HV information
    DPMAP* HV_map = accessHVObjectCCDB(558801);

    // read out histograms
    const char* outputFileName = isMC ? "cluster_hist_position_MC.root" : "cluster_hist_position_data.root";
    auto outFile = new TFile(outputFileName, "RECREATE");

    const char* legendLabel = isMC ? "O2 Sim" : "Data (Run 3)";

    // store cumulative fitted parameters for each deId/cathode/pitch
    std::vector<fitParams*> fittedParams;

    // create current O2 PDF and parameters for comparison and fitting
    o2::mid::ChamberResponseParams chamberRespParam = o2::mid::createDefaultChamberResponseParams();

    // Arrays to store original parameter values
    std::vector<double> currentBParams;
    std::array<double, 2> currentAParams; // a0, a1
    currentAParams[0] = chamberRespParam.getParametersA()[0];
    currentAParams[1] = chamberRespParam.getParametersA()[1] * 1000.; // defined as [1/V] in ChamberResponseParams, but needs to be multiplied by 1000 for the PDF (HV in kV)
    std::array<double, 2> currentCParams; // c0, c1
    currentCParams[0] = chamberRespParam.getParametersC()[0];
    currentCParams[1] = chamberRespParam.getParametersC()[1] * 1000.;
    
    for (auto clusterHist : clusterPosHistograms) {
        if (clusterHist->nClusters == 0) continue; // skip empty histograms

        clusterHist->clusterPosition->Scale(1.0 / (static_cast<float>(clusterHist->nClusters))); // normalize by number of clusters

        // HV for the deID
        float HV_value = accessHVperDE(HV_map, clusterHist->deId);

        // Get the range of the current histogram
        double xMin = clusterHist->clusterPosition->GetXaxis()->GetXmin();
        double xMax = clusterHist->clusterPosition->GetXaxis()->GetXmax();

        // original parameters for comparison
        auto currenBparam = chamberRespParam.getParB(clusterHist->cathode, clusterHist->deId);
        currentBParams.push_back(currenBparam);

        // define current PDF for comparison
        auto clusterPDF_current = new TF1(Form("clusterPDF_current_de%i_cathode%i_pitch%i", clusterHist->deId, clusterHist->cathode, clusterHist->pitch), pdfFunc, 0, xMax, 7);
        clusterPDF_current->SetParNames("b", "a0", "a1", "c0", "c1", "hv", "theta");
        clusterPDF_current->SetParameters(currenBparam, currentAParams[0], currentAParams[1], currentCParams[0], currentCParams[1], HV_value, 0.0); // current parameters
        for (int i = 0; i < 7; ++i) clusterPDF_current->FixParameter(i, clusterPDF_current->GetParameter(i)); // fix all parameters (for comparison only)

        // define my own function for fitting
        auto clusterPDF_fit = new TF1(Form("clusterPDF_fit_de%i_cathode%i_pitch%i", clusterHist->deId, clusterHist->cathode, clusterHist->pitch), pdfFunc, 0, xMax, 7);
        clusterPDF_fit->SetParNames("b", "a0", "a1", "c0", "c1", "hv", "theta");
        clusterPDF_fit->SetParameters(currenBparam, currentAParams[0], currentAParams[1], currentCParams[0], currentCParams[1], HV_value, 0.0); // initial parameters

        // Set parameter limits for fitting (a,c limits were determined by previous fits)
        double percentageChangeB = 0.9; // 'b' parameter needs less constraining
        //clusterPDF_fit->SetParLimits(0, currenBparam * (1.0 - percentageChangeB), currenBparam * (1.0 + percentageChangeB)); // b
        //clusterPDF_fit->SetParLimits(1, -60.0, -15.0); // a0 (negative)
        clusterPDF_fit->FixParameter(1, currentAParams[0]); // a0 (negative)
        //clusterPDF_fit->SetParLimits(2, 6.0, 15.0); // a1
        clusterPDF_fit->FixParameter(2, currentAParams[1]); //a1
        //clusterPDF_fit->SetParLimits(3, -0.015, 0.005); // c0 (negative)
        clusterPDF_fit->FixParameter(3, -0.00475); // c0 (negative)
        //clusterPDF_fit->SetParLimits(4, -0.0015, 0.0015); // c1
        clusterPDF_fit->FixParameter(4, 0.0004); // c1
        clusterPDF_fit->FixParameter(5, HV_value); // fix HV parameter
        clusterPDF_fit->FixParameter(6, 0.0); // fix theta parameter

        // fit the histogram with the PDF
        double fitMin = clusterHist->clusterPosition->GetBinLowEdge(1); // Lower edge of the first bin
        int lastNonEmptyBin = clusterHist->clusterPosition->FindLastBinAbove(0); // Find the last non-empty bin
        if (lastNonEmptyBin > 8) {lastNonEmptyBin = 8;} // limit the fit range to the first n bins
        double fitMax = clusterHist->clusterPosition->GetBinLowEdge(lastNonEmptyBin + 1); // Upper edge of the last non-empty bin
        clusterHist->clusterPosition->Fit(Form("clusterPDF_fit_de%i_cathode%i_pitch%i", clusterHist->deId, clusterHist->cathode, clusterHist->pitch), "R", "", fitMin, fitMax);

        // add parameters for later analysis (only add when it has 4 or more bins with entries to ensure fit reliability)
        fitParams* params = new fitParams();
        params->deId = clusterHist->deId;
        params->cathode = clusterHist->cathode;
        params->pitch = clusterHist->pitch;
        params->nBins = clusterHist->clusterPosition->FindLastBinAbove(0);
        params->nEntries = clusterHist->clusterPosition->GetEntries();
        params->params[0] = clusterPDF_fit->GetParameter(0); // b
        params->params[1] = clusterPDF_fit->GetParameter(1); // a0
        params->params[2] = clusterPDF_fit->GetParameter(2); // a1
        params->params[3] = clusterPDF_fit->GetParameter(3); // c0
        params->params[4] = clusterPDF_fit->GetParameter(4); // c1
        fittedParams.push_back(params);

        // plot first hist, current PDF, and fitted function together
        TCanvas* canvasCheckFirst = new TCanvas(Form("strip_position_de%i_cathode%i_pitch%i", clusterHist->deId, clusterHist->cathode, clusterHist->pitch), "Fired Probability vs Distance", 800, 600);
        gPad->SetLogy();
        clusterHist->clusterPosition->SetLineColor(kBlack);
        clusterHist->clusterPosition->SetMaximum(1.0);
        clusterHist->clusterPosition->SetMinimum(0.0001);
        clusterHist->clusterPosition->SetStats(0);
        clusterHist->clusterPosition->Draw("SAME");

        clusterPDF_current->SetLineColor(kRed);
        clusterPDF_current->SetLineWidth(2);
        clusterPDF_current->Draw("SAME");

        clusterPDF_fit->SetLineColor(kBlue);
        clusterPDF_fit->SetLineWidth(2);
        clusterPDF_fit->DrawClone("SAME");

        // Add a legend
        TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend->AddEntry(clusterHist->clusterPosition, legendLabel, "l");
        legend->AddEntry(clusterPDF_current, "Current O2 PDF", "l");
        legend->AddEntry(clusterPDF_fit, "Fitted PDF", "l");
        legend->Draw("SAME");

        // Add a label with fit parameters and errors
        TPaveText* fitParamsLabel = new TPaveText(0.15, 0.6, 0.45, 0.9, "NDC");
        fitParamsLabel->SetBorderSize(1);
        fitParamsLabel->SetFillColor(0);
        fitParamsLabel->SetTextAlign(12);

        for (int i = 0; i < 7; ++i) {
            TString paramText = TString::Format("%s = %.3g #pm %.3g", clusterPDF_fit->GetParName(i), clusterPDF_fit->GetParameter(i), clusterPDF_fit->GetParError(i));
            fitParamsLabel->AddText(paramText);
        }
        fitParamsLabel->Draw("SAME");

        canvasCheckFirst->Write();
    }

    delete outFile;

    processFitParams(fittedParams, currentBParams, currentAParams, currentCParams);
}

void processFitParams(std::vector<fitParams*> fittedParams, const std::vector<double>& currentBParams, const std::array<double, 2>& currentAParams, const std::array<double, 2>& currentCParams) {
    // Check if there are any fitted parameters to process 
    if (fittedParams.empty()) {
        std::cout << "No fitted parameters to process." << std::endl;
        return;
    }

    // Histograms for each parameter
    std::map<std::string, TH1D*> histograms;
    histograms["b"] = new TH1D("b_distribution", "Distribution of b-parameter; b; Count", 40, 1.0, 4.0);
    histograms["b_original"] = new TH1D("b_original_distribution", "Original Distribution of parameter b; b; Count", 40, 1.0, 4.0);
    histograms["a0"] = new TH1D("a0_distribution", "Distribution of fitted a0-parameters; a0; Count", 40, -80.0, 120.0);
    histograms["a1"] = new TH1D("a1_distribution", "Distribution of fitted a1-parameters; a1; Count", 40, 0.0, 50.0);
    histograms["c0"] = new TH1D("c0_distribution", "Distribution of fitted c0-parameters; c0; Count", 40, -0.015, 0.005);
    histograms["c1"] = new TH1D("c1_distribution", "Distribution of fitted c1-parameters; c1; Count", 40, -0.0015, 0.0015);

    // Fill original b parameter histogram
    for (const auto& b : currentBParams) {
        histograms["b_original"]->Fill(b);
    }

    // Filter fitParams to select only one object per deId/cathode
    std::map<std::pair<int, int>, fitParams*> filteredParams; // Key: {deId, cathode}

    for (const auto& params : fittedParams) {
        auto key = std::make_pair(params->deId, params->cathode);

        // Check if the key already exists in the map
        if (filteredParams.find(key) == filteredParams.end()) {
            // If not, add the current params
            filteredParams[key] = params;
        } else {
            // If it exists, compare nBins and nEntries
            auto existingParams = filteredParams[key];
            if (params->nBins > existingParams->nBins || 
                (params->nBins == existingParams->nBins && params->nEntries > existingParams->nEntries)) {
                filteredParams[key] = params; // Replace with the current params
            }
        }
    }

    // output textfile with fit parameters for each deId/cathode
    std::ofstream outFileTxt("fitted_parameters_O2.txt");
    outFileTxt << "params.setParB(cathode, deId, b);"; // header

    // Fill fitted parameter histograms using the filtered parameters
    for (const auto& [key, params] : filteredParams) {
        histograms["b"]->Fill(params->params[0]); // b
        histograms["a0"]->Fill(params->params[1]); // a0
        histograms["a1"]->Fill(params->params[2]); // a1
        histograms["c0"]->Fill(params->params[3]); // c0
        histograms["c1"]->Fill(params->params[4]); // c1

        // Write b-parameter to text file for O2
        outFileTxt << "params.setParB(" << params->cathode << ", " << params->deId << ", " << std::fixed << std::setprecision(2) << params->params[0] << ");" << std::endl;
    }

    // Create a ROOT file to store histograms
    auto outFile = new TFile("fit_parameters_distributions.root", "RECREATE");

    // Plot b parameter distributions on the same canvas
    TCanvas* canvasB = new TCanvas("b_distributions_compare", "Distributions of original and fitted b-parameters", 800, 600);
    histograms["b_original"]->SetLineColor(kRed);
    histograms["b_original"]->SetStats(0);
    histograms["b_original"]->Draw();
    histograms["b"]->SetLineColor(kBlue);
    histograms["b"]->SetStats(0);
    histograms["b"]->Draw("SAME");

    // Add a legend to distinguish between the two distributions
    TLegend* legendB = new TLegend(0.6, 0.7, 0.9, 0.9);
    legendB->AddEntry(histograms["b_original"], "Original b-parameters", "l");
    legendB->AddEntry(histograms["b"], "Fitted b-parameters", "l");
    legendB->Draw("SAME");

    canvasB->Write();

    // Create a canvas for each histogram and draw the line on the same canvas
    for (auto& [name, hist] : histograms) {
        if (name == "b" || name == "b_original") {
            continue; // Skip b-parameter histograms
        }

        TCanvas* canvas = new TCanvas(Form("%s_distribution", name.c_str()), Form("Distribution of fitted %s-parameters", name.c_str()), 800, 600);
        hist->SetStats(0);
        hist->Draw();

        TLine* line = nullptr;
        if (name == "a0") {
            line = new TLine(currentAParams[0], 0, currentAParams[0], hist->GetMaximum());
        } else if (name == "a1") {
            line = new TLine(currentAParams[1], 0, currentAParams[1], hist->GetMaximum());
        } else if (name == "c0") {
            line = new TLine(currentCParams[0], 0, currentCParams[0], hist->GetMaximum());
        } else if (name == "c1") {
            line = new TLine(currentCParams[1], 0, currentCParams[1], hist->GetMaximum());
        }

        if (line) {
            line->SetLineColor(kRed);
            line->SetLineStyle(2);
            line->Draw("SAME");
        }

        // Add a legend to distinguish between the histogram and the original parameter line
        TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend->AddEntry(hist, "Fitted Parameters", "l");
        if (line) {
            legend->AddEntry(line, "Original Parameter", "l");
        }
        legend->Draw("SAME");

        canvas->Write(); // Save the canvas to the ROOT file
    }

    delete outFile;
}

float accessHVperDE(DPMAP* HV_map, int deId) {
    std::string deAlias = o2::mid::dcs::detElemId2DCSAlias(deId, o2::mid::dcs::MeasurementType::HV_V);
    for (auto const& [key, values] : *HV_map) {
        std::string entry_alias = key.get_alias();
        if (entry_alias == deAlias) {
            // std::cout << "HV values for DE " << deId << ": ";
            // for (auto val : values) {
            //     double HV_value = sum(0.0, val);
            //     std::cout << HV_value << ", ";
            // }
            // std::cout << std::endl;
            float HV_value = sum(0.0, values.back())/1000; // return the last HV value (in kV) - for now!
            //std::cout << "Accessed HV for DE " << deId << ": " << HV_value << " kV" << std::endl;
            return HV_value;
        }
    }
}

DPMAP* accessHVObjectCCDB(int runNumber) {
    // Access CCDB API to retrieve HV values
    o2::ccdb::CcdbApi api;
    std::string ccdbUrl = "http://alice-ccdb.cern.ch";
    api.init(ccdbUrl);

    // access info about the relevant run from the CCDB
    std::pair<uint64_t, uint64_t> runBoundaries = o2::ccdb::CCDBManagerInstance::getRunDuration(api, runNumber);

    // access HV values for run from CCDB
    std::map<std::string, std::string> metadata;
    DPMAP* HV_map = api.retrieveFromTFileAny<DPMAP>("MID/Calib/HV", metadata, runBoundaries.first);

    // for (auto const& [key, values] : *HV_map) {
    //     std::string entry_alias = key.get_alias();
    //     std::cout << "HV entry alias: " << entry_alias << std::endl;
    // }

    return HV_map;
}

void processClustersSize(std::vector<cluster*> midClusters) {
    // Initialize histograms segmented by deId/cathode/pitch
    std::vector<clusterSizeHist*> clusterSizeHistograms = initializeClusterSizeHist();

    // Initialize histograms aggregated by cathode/pitch
    std::map<std::pair<int, int>, TH1D*> aggregatedHistograms;

    std::cout << "Number of processed MID clusters: " << midClusters.size() << std::endl;

    // loop over preClusters
    for (auto cluster : midClusters) {
        // Fill histograms segmented by deId/cathode/pitch
        for (auto& hist : clusterSizeHistograms) {
            if (hist->deId == cluster->deId && hist->cathode == cluster->cathode && hist->pitch == cluster->pitch) {
                hist->clusterSize->Fill(cluster->nStrips);
                break;
            }
        }

        // Fill aggregated histograms by cathode/pitch
        auto key = std::make_pair(cluster->cathode, cluster->pitch);
        if (aggregatedHistograms.find(key) == aggregatedHistograms.end()) {
            Int_t nBins = 64 / cluster->pitch;
            aggregatedHistograms[key] = new TH1D(Form("cluster_size_cathode%i_pitch%i", cluster->cathode, cluster->pitch),
                                                 Form("Cluster Size for Cathode %i, Pitch %i; nStrips; Count", cluster->cathode, cluster->pitch),
                                                 nBins, 1, nBins + 1);
            aggregatedHistograms[key]->Sumw2();
        }
        aggregatedHistograms[key]->Fill(cluster->nStrips);
    }

    // read out histograms
    auto outFile = new TFile("cluster_size_histograms.root", "RECREATE");

    // Write histograms segmented by deId/cathode/pitch
    for (auto hist : clusterSizeHistograms) {
        if (hist->clusterSize->GetEntries() == 0) continue; // skip empty histograms
        TCanvas* canvas = new TCanvas(Form("cluster_size_de%i_cathode%i_pitch%i", hist->deId, hist->cathode, hist->pitch), "", 800, 600);
        canvas->SetLogy(); // Set y-axis to log scale
        hist->clusterSize->Draw("E1"); // Draw histogram with error bars and horizontal lines
        canvas->Write(); // Write the canvas to the output file
    }

    // Write aggregated histograms by cathode/pitch
    for (auto& [key, hist] : aggregatedHistograms) {
        hist->Write();
    }

    delete outFile;
}

std::vector<clusterPosHist*> processClustersPosition(std::vector<cluster*> midClusters) {
    // vector of cluster size histograms for fitting
    std::vector<clusterPosHist*> clusterPosHistograms = initializeClusterPosHist();

    // loop over preClusters
    for (auto cluster : midClusters) {
        // fill cluster size histograms
        for (auto& hist : clusterPosHistograms) {
            if (hist->deId == cluster->deId && hist->cathode == cluster->cathode && hist->pitch == cluster->pitch) {
                fillStripPositions(hist, cluster->nStrips, cluster->pitch);
                break;
            }
        }
    }

    return clusterPosHistograms; // to be used for fitting
}

cluster* processPreCluster(int rofID, o2::mid::PreCluster pc, o2::mid::Mapping* mapper) {
    // find the strip pitch 
    o2::mid::MpArea mpArea_first = mapper->stripByLocation(pc.firstStrip, pc.cathode, pc.firstLine, pc.firstColumn, pc.deId);
    o2::mid::MpArea mpArea_last = mapper->stripByLocation(pc.lastStrip, pc.cathode, pc.lastLine, pc.lastColumn, pc.deId);
    int strip_pitch_first, strip_pitch_last;

    // determine strip size 
    int nStrips;
    if (pc.cathode == 0) { // is the cluster in the bending plane?
        if (pc.firstColumn != pc.lastColumn) return NULL; // should not happen (different cluster per column)

        strip_pitch_first = static_cast<int>(2*mpArea_first.getHalfSizeY());
        strip_pitch_last = static_cast<int>(2*mpArea_last.getHalfSizeY());

        if (strip_pitch_first != strip_pitch_last) return NULL;

        nStrips = ((int)pc.lastStrip - (int)pc.firstStrip) + (16 * ((int)pc.lastLine - (int)pc.firstLine)) + 1;
    } 
    else {
        strip_pitch_first = static_cast<int>(2*mpArea_first.getHalfSizeX());
        strip_pitch_last = static_cast<int>(2*mpArea_last.getHalfSizeX());
        
        if (strip_pitch_first != strip_pitch_last) return NULL;

        nStrips = (pc.lastStrip - pc.firstStrip) + 1;
        for (int column = pc.firstColumn; column < pc.lastColumn; column++) {
            nStrips += mapper->getNStripsNBP(column, pc.deId);
        }
    }

    cluster* processedCluster = new cluster();
    processedCluster->rofID = rofID;
    processedCluster->cathode = pc.cathode;
    processedCluster->deId = pc.deId;
    processedCluster->pitch = strip_pitch_first;
    processedCluster->nStrips = nStrips;

    return processedCluster;
}

std::vector<cluster*> processMIDdigits(const char *fileMIDdigits, const char *fileMIDtracks, std::vector<o2::InteractionRecord> muonTracksIR, bool isMC, bool removeMultColumns) {
    // output histograms
    TH1D* nMuonTracksROF = new TH1D("muon_track_count_ROF", "Number of MCH+MID matched tracks per ROF;nTracks/ROF;count", 5, 0, 5);
    TH1D* nMIDTracksROF = new TH1D("mid_track_count_ROF", "Number of MID tracks per ROF;nTracks/ROF;count", 70, 0, 70);

    // read in the MID track and digit infomation
    const char* treeNameDigits = isMC ? "o2sim" : "middigits";
    auto [digitFile, digitReader] = loadData(fileMIDdigits, treeNameDigits);

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

    // mapping object to extract strip size per column
    o2::mid::Mapping* mapper = new o2::mid::Mapping();

    // vector of MID preClusters
    std::vector<cluster*> midClusters;

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

            if (muonTrackCount == 1) { // number of muon tracks per ROF to consider
                nMuonTracksROF->Fill(muonTrackCount);

                // subspan of digits for the ROF where there are tracks
                auto eventDigits = sdigits.subspan(digitRofIt->firstEntry, digitRofIt->nEntries);

                // run pre-clusterizer on the event digits
                preClusterizer.process(eventDigits);
                auto preClusters = preClusterizer.getPreClusters();

                bool bendingClusterPicker[nPitches][o2::mid::detparams::NDetectionElements] = {{false}}; // to pick only one bending plane cluster per pitch and DE

                for (auto& pc : preClusters)  { // process pre-clusters
                    cluster* processedCluster = processPreCluster(selectedROFcounter, pc, mapper);
                    if (processedCluster != NULL) {
                        if (processedCluster->cathode == 0 && removeMultColumns) { // bending plane
                            int pitchIndex = pitchToIndex(processedCluster->pitch);
                            if (!bendingClusterPicker[pitchIndex][processedCluster->deId]) {
                                midClusters.push_back(processedCluster);
                                bendingClusterPicker[pitchIndex][processedCluster->deId] = true;
                            } else {
                                delete processedCluster; // already have a cluster for this pitch and DE in this ROF
                            }
                        } else { // non-bending plane or no removal of multiple columns
                            midClusters.push_back(processedCluster);
                        }
                    }
                }

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

    return midClusters;
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
    std::vector<clusterSizeHist*> clusterSizeHistograms;

    // Loop over all MID deIds
    for (int deId = 0; deId < o2::mid::detparams::NDetectionElements; deId++) {
        for (int cathode = 0; cathode < 2; cathode++) { // 0: bending, 1: non-bending
            for (int pitch = 0; pitch < nPitches; pitch++) { // strip pitches
                if (pitch == 0 && cathode == 1) continue; // Skip non-bending plane for pitch 1
                auto hist = new clusterSizeHist();
                hist->deId = deId;
                hist->cathode = cathode;
                hist->pitch = indexToPitch(pitch);
                const Int_t nBins = 64 / hist->pitch; // Number of bins for the histogram
                hist->clusterSize = new TH1D(Form("cluster_size_de%i_cathode%i_pitch%i", deId, cathode, indexToPitch(pitch)),
                                             Form("Cluster Size for DE %i, Cathode %i, Pitch %i; nStrips; Count", deId, cathode, indexToPitch(pitch)),
                                             nBins, 1, nBins + 1);
                hist->clusterSize->Sumw2();
                clusterSizeHistograms.push_back(hist);
            }
        }
    }

    return clusterSizeHistograms;
}

std::vector<clusterPosHist*> initializeClusterPosHist() {
    // initialize cluster size histograms
    std::vector<clusterPosHist*> clusterPosHistograms;

    // loop over all MID deIds
    for (int deId = 0; deId < o2::mid::detparams::NDetectionElements; deId++) {
        for (int cathode = 0; cathode < 2; cathode++) { // 0: bending, 1: non-bending
            for (int pitch = 0; pitch < nPitches; pitch++) { // strip pitches
                if (pitch == 0 && cathode == 1) continue; // skip non-bending plane for pitch 1
                auto hist = new clusterPosHist();
                hist->deId = deId;
                hist->cathode = cathode;
                hist->pitch = indexToPitch(pitch);
                hist->chamber = o2::mid::detparams::getChamber(deId);
                hist->nClusters = 0;

                Double_t stripWidth = o2::mid::geoparams::getStripUnitPitchSize(hist->chamber) * static_cast<Double_t>(hist->pitch) * 10; // in mm
                const Int_t nBins = 2 * 32 / hist->pitch; // number of bins for the histogram
                hist->clusterPosition = new TH1F(Form("strip_position_de%i_cathode%i_pitch%i", deId, cathode, indexToPitch(pitch)), 
                                                Form("Strip Position (mm) for DE %i, Cathode %i, Pitch %i; x (mm); f(x)", deId, cathode, indexToPitch(pitch)), nBins, 0.0, (nBins/2) * stripWidth);
                hist->clusterPosition->Sumw2();
                clusterPosHistograms.push_back(hist);
            }
        }
    }

    return clusterPosHistograms;
}

float sum(float s, o2::dcs::DataPointValue v) {
    union Converter {
        uint64_t raw_data;
        double value;
    } converter;
    converter.raw_data = v.payload_pt1;
    return s + converter.value;
};