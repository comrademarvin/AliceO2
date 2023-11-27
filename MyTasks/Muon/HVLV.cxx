#include "CCDB/CcdbApi.h"
#include "DetectorsDCS/DataPointIdentifier.h"
#include "DetectorsDCS/DataPointValue.h"
#include "DetectorsDCS/DataPointCreator.h"
#include "DetectorsDCS/DataPointCompositeObject.h"
#include "MCHGlobalMapping/Mapper.h"
#include "MCHGlobalMapping/DsIndex.h"
#include "MCHGlobalMapping/ChannelCode.h"
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <numeric>
#include <vector>
#include <set>
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"

using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;
using DPMAP = std::unordered_map<DPID, std::vector<DPVAL>>;

float sum(float s, o2::dcs::DataPointValue v)
{
    union Converter {
        uint64_t raw_data;
        double value;
    } converter;
    converter.raw_data = v.payload_pt1;
    return s + converter.value;
};

int main() {
    TGraph* HV_entries = new TGraph();

    o2::ccdb::CcdbApi api;

    string ccdbUrl = "http://alice-ccdb.cern.ch";
    api.init(ccdbUrl);

    std::map<std::string, std::string> metadata;
    uint64_t timestamp = 1698299830000;
    auto* HV_map = api.retrieveFromTFileAny<DPMAP>("MCH/Calib/HV", metadata, timestamp);

    std::cout << "size of map = " << HV_map->size() << std::endl;

    auto outFile = new TFile("HVLV_ouput.root", "RECREATE");

    // as a naive start, use a HV threshold
    double HV_threshold = 1000.0;

    for (auto& entry : *HV_map) {
        std::string entry_alias = (entry.first).get_alias();

        if (entry_alias.find("vMon") != std::string::npos) { // only for voltage channels
            auto entry_values = entry.second;
            Int_t entries_size = entry_values.size();
            bool HV_threshold_trigger = false;

            if (entries_size > 1) {
                // std::cout << "Key: " << entry_alias;
                // std::cout << " |-| Values: ";
                int entry_count = 0;
                Int_t voltage_array[entries_size];
                Int_t timestamp_array[entries_size];

                for (auto value : entry_values) {
                    double HV_entry_value = sum(0.0, value);
                    // HV entry below threshold
                    if (HV_entry_value < HV_threshold) {
                        HV_threshold_trigger = true;
                    }

                    //std::cout << (Int_t)(value.get_epoch_time()*(0.001)) << std::endl;
                    timestamp_array[entry_count] = (Int_t)(value.get_epoch_time()*(0.001));
                    voltage_array[entry_count] = (Int_t)(HV_entry_value);
                    entry_count++;
                    //std::cout << sum(0.0, value) << " (" << entry_timestamp << ");";
                }

                auto canvasHV = new TCanvas("HV_timestamps","HV_timestamps");
                auto HV_entries = new TGraph(entries_size, timestamp_array, voltage_array);
                const std::string title_axis = entry_alias + ";Timestamp (s);High-Voltage (V)";
                HV_entries->SetTitle(title_axis.c_str());
                // HV_entries->GetXaxis()->SetTimeDisplay(1);
                // HV_entries->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
                // HV_entries->GetXaxis()->SetTimeOffset(0,"gmt");
                // HV_entries->GetXaxis()->SetNdivisions(505);
                HV_entries->Draw("AL*");

                canvasHV->Write();

                // auto mean = std::accumulate(entry_values.begin(), entry_values.end(), 0.0, sum);
                // if (entry_values.size()) {
                //     mean /= entry_values.size();
                // }
                // std::cout << " |-| Mean: " << mean << std::endl;
            }

            if (HV_threshold_trigger) { // flagged as bad channel
                std::cout << "Alias of channel: " << entry_alias << std::endl;
                std::set<int> channel_DsIndices = o2::mch::dcs::aliasToDsIndices("MchHvLvRight/Chamber00Right/Quad0Sect0.actual.vMon");
                std::cout << "Made it out!" << std::endl;
                //std::cout << "Number of channel DS indices: " << channel_DsIndices << std::endl;
            }
        }
    }

    delete outFile;
    
    return 0;
}