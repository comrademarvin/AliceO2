#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsDCS/DataPointIdentifier.h"
#include "DetectorsDCS/DataPointValue.h"
#include "DetectorsDCS/DataPointCreator.h"
#include "DetectorsDCS/DataPointCompositeObject.h"
#include "MCHGlobalMapping/Mapper.h"
#include "MCHGlobalMapping/DsIndex.h"
#include "MCHGlobalMapping/ChannelCode.h"
#include "MCHConditions/DCSAliases.h"
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <numeric>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"

using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;
using DPMAP = std::unordered_map<DPID, std::vector<DPVAL>>;
using DPMAP2 = std::map<std::string, std::map<uint64_t, double>>;
using RBMAP = std::map<int, std::pair<uint64_t, uint64_t>>;
using HVBMAP = std::map<uint64_t, uint64_t>;

float sum(float s, o2::dcs::DataPointValue v);
std::set<int> GetRuns(std::string runList);
RBMAP GetRunBoundaries(o2::ccdb::CcdbApi const& api, std::string runList);
void PrintRunBoundaries(const RBMAP& runBoundaries);
std::string GetTime(uint64_t ts);
uint64_t MSToS(uint64_t ts);
HVBMAP GetHVBoundaries(o2::ccdb::CcdbApi const& api, uint64_t tStart, uint64_t tStop);
void PrintHVBoundaries(const HVBMAP& hvBoundaries);
void PrintDataPoints(const DPMAP2 dpsMapsPerCh[10], bool all);

int main() {
    //TGraph* HV_entries = new TGraph();

    // read in runlist from textfile
    std::ifstream runfile("runlist.txt");
    std::stringstream buffer;
    buffer << runfile.rdbuf();

    std::string runlist = buffer.str();

    // Access CCDB API
    o2::ccdb::CcdbApi api;

    string ccdbUrl = "http://alice-ccdb.cern.ch";
    api.init(ccdbUrl);

    // look at run boundaries
    auto runBoundaries = GetRunBoundaries(api, runlist);

    if (runBoundaries.empty()) {
        printf("no run found from the list\n");
        return 0;
    }

    PrintRunBoundaries(runBoundaries);

    // extract the timestamp boundaries for each HV file in the full time range
    auto hvBoundaries = GetHVBoundaries(api, runBoundaries.begin()->second.first, runBoundaries.rbegin()->second.second);
    PrintHVBoundaries(hvBoundaries);


    // fetch HV object in timestamp range
    DPMAP2 dpsMapsPerCh[10];
    std::map<std::string, std::string> metadata;

    for (auto boundaries : hvBoundaries) {
      auto* HV_map = api.retrieveFromTFileAny<DPMAP>("MCH/Calib/HV", metadata, boundaries.first);

      for (auto& entry : *HV_map) {
        std::string entry_alias = (entry.first).get_alias();

        if (entry_alias.find("vMon") != std::string::npos) { // only for voltage channels
          auto entry_values = entry.second;
          Int_t entries_size = entry_values.size();
          int chamber = o2::mch::dcs::toInt(o2::mch::dcs::aliasToChamber(entry_alias));

          auto& dps2 = dpsMapsPerCh[chamber][entry_alias];
          for (const auto& value : entry_values) {
            double HV_value = sum(0.0, value);
            dps2.emplace(value.get_epoch_time(), HV_value);
          }
        }
      }
    }

    PrintDataPoints(dpsMapsPerCh, true);

    /*
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

    delete outFile; */
    
    return 0;
}

float sum(float s, o2::dcs::DataPointValue v)
{
    union Converter {
        uint64_t raw_data;
        double value;
    } converter;
    converter.raw_data = v.payload_pt1;
    return s + converter.value;
};

RBMAP GetRunBoundaries(o2::ccdb::CcdbApi const& api, std::string runList)
{
  /// return the SOR / EOR time stamps for every runs in the list

  RBMAP runBoundaries{};

  auto runs = GetRuns(runList);

  for (auto run : runs) {
    auto boundaries = o2::ccdb::CCDBManagerInstance::getRunDuration(api, run);
    runBoundaries.emplace(run, boundaries);
  }

  return runBoundaries;
}

std::set<int> GetRuns(std::string runList)
{
  /// read the runList from an ASCII file, or a comma separated run list, or a single run

  std::set<int> runs{};

  auto isNumber = [](std::string val) { return !val.empty() && val.find_first_not_of("0123456789") == val.npos; };

  if (isNumber(runList)) {

    runs.insert(std::stoi(runList));

  } else if (runList.find(",") != runList.npos) {

    std::istringstream input(runList);
    for (std::string run; std::getline(input, run, ',');) {
      if (isNumber(run)) {
        runs.insert(std::stoi(run));
      }
    }

  } else {

    std::ifstream input(runList);
    if (input.is_open()) {
      for (std::string run; std::getline(input, run);) {
        if (isNumber(run)) {
          runs.insert(std::stoi(run));
        }
      }
    }
  }

  return runs;
}

void PrintRunBoundaries(const RBMAP& runBoundaries)
{
  /// print the list of runs with their time boundaries

  printf("\nlist of runs with their boundaries:\n");
  printf("------------------------------------\n");

  for (const auto& [run, boundaries] : runBoundaries) {
    printf("%d: %lld - %lld (%s - %s)\n", run, boundaries.first, boundaries.second,
           GetTime(boundaries.first).c_str(), GetTime(boundaries.second).c_str());
  }

  printf("------------------------------------\n");
}

std::string GetTime(uint64_t ts)
{
  /// convert the time stamp (ms) to local time

  time_t t = MSToS(ts);

  std::string time = std::ctime(&t);
  time.pop_back(); // remove trailing \n

  return time;
}

uint64_t MSToS(uint64_t ts)
{
  /// convert the time stamp from ms to s

  return (ts + 500) / 1000;
}

HVBMAP GetHVBoundaries(o2::ccdb::CcdbApi const& api, uint64_t tStart, uint64_t tStop)
{
  /// get the time boundaries of every HV files found in the time range

  // add extra margin (ms) of ± 1 min to the creation time, which occurs every 30 min
  static const uint64_t timeMarging[2] = {60000, 1860000};

  HVBMAP hvBoundaries{};

  std::istringstream fileInfo(api.list("MCH/Calib/HV", false, "text/plain",
                                       tStop + timeMarging[1], tStart - timeMarging[0]));

  for (std::string line; std::getline(fileInfo, line);) {
    if (line.find("Validity:") == 0) {
      hvBoundaries.emplace(std::stoull(line.substr(10, 13)), std::stoull(line.substr(26, 13)));
    }
  }

  return hvBoundaries;
}

void PrintHVBoundaries(const HVBMAP& hvBoundaries)
{
  /// print the time boundaries of every HV files found in the full time range

  printf("\nlist of HV file time boundaries:\n");
  printf("------------------------------------\n");

  for (auto [tStart, tStop] : hvBoundaries) {
    printf("%lld - %lld (%s - %s)\n", tStart, tStop, GetTime(tStart).c_str(), GetTime(tStop).c_str());
  }

  printf("------------------------------------\n");
}

void PrintDataPoints(const DPMAP2 dpsMapsPerCh[10], bool all)
{
  /// print all the registered data points

  for (int ch = 0; ch < 10; ++ch) {

    printf("\n------------ chamber %d ------------\n", ch + 1);

    for (const auto& [alias, dps] : dpsMapsPerCh[ch]) {

      printf("- %s: %lu values", alias.c_str(), dps.size());

      if (all) {

        printf("\n");
        for (const auto& [ts, hv] : dps) {
          printf("  %lld (%s): %7.2f V\n", ts, GetTime(ts).c_str(), hv);
        }

      } else if (!dps.empty()) {

        const auto firstdt = dps.begin();
        const auto lastdt = dps.rbegin();
        printf(": %lld (%s): %7.2f V -- %lld (%s): %7.2f V\n",
               firstdt->first, GetTime(firstdt->first).c_str(), firstdt->second,
               lastdt->first, GetTime(lastdt->first).c_str(), lastdt->second);

      } else {
        printf("\n");
      }
    }
  }
}