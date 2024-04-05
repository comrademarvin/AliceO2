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
#include "MCHStatus/HVStatusCreator.h"
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <numeric>
#include <vector>
#include <set>
#include <tuple>
#include <fstream>
#include <sstream>
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"

using DPID = o2::dcs::DataPointIdentifier;
using DPVAL = o2::dcs::DataPointValue;
using DPMAP = std::unordered_map<DPID, std::vector<DPVAL>>;
using HVVALUES = std::map<uint64_t, double>;
using DPMAP2 = std::map<std::string, HVVALUES>;
using RBMAP = std::map<int, std::pair<uint64_t, uint64_t>>;
using HVBMAP = std::map<uint64_t, uint64_t>;
using BADHVLIST = std::vector<std::tuple<uint64_t, uint64_t, double, double, uint32_t, std::string>>;
using BADHVMAP = std::map<std::string, BADHVLIST>;

double yRange[2] = {-1., 1700.};
double hvLimits[10] = {1550., 1550., 1590., 1590., 1590., 1590., 1590., 1590., 1590., 1590.};
uint64_t minDuration = 3000; // Tune this for fluctuations

float sum(float s, o2::dcs::DataPointValue v);
std::set<int> GetRuns(std::string runList);
RBMAP GetRunBoundaries(o2::ccdb::CcdbApi const& api, std::string runList);
void PrintRunBoundaries(const RBMAP& runBoundaries);
std::string GetTime(uint64_t ts);
uint64_t MSToS(uint64_t ts);
std::string GetDE(std::string alias);
HVBMAP GetHVBoundaries(o2::ccdb::CcdbApi const& api, uint64_t tStart, uint64_t tStop);
void SelectHVBoundaries(HVBMAP& hvBoundaries, const RBMAP& runBoundaries);
void PrintHVBoundaries(const HVBMAP& hvBoundaries);
void PrintDataPoints(const DPMAP2 dpsMapsPerCh[10], bool all);
void SelectDataPoints(DPMAP2 dpsMapsPerCh[10], uint64_t tStart, uint64_t tStop);
TGraph* MapToGraph(std::string alias, const std::map<uint64_t, double>& dps);
TCanvas* DrawDataPoints(TMultiGraph* mg);
void DrawRunBoudaries(const RBMAP& runBoundaries, TCanvas* c);
void FindHVIssues(const HVVALUES& HV_values, double hvLimit, BADHVLIST& hvIssuesList);
void SelectHVIssues(BADHVMAP& hvIssuesList, const RBMAP& runBoundaries, uint64_t minDuration);
void PrintHVIssues(const BADHVMAP hvIssuesPerCh[10]);
void FindHVIssuesMacro(o2::ccdb::CcdbApi const& api, HVBMAP& hvBoundaries, const RBMAP& runBoundaries);
void FindHVIssuesPerObjectMacro(o2::ccdb::CcdbApi const& api, HVBMAP& hvBoundaries, const RBMAP& runBoundaries, bool plotHVObjects);
void FindHVIssuesWithClassMacro(o2::ccdb::CcdbApi const& api, HVBMAP& hvBoundaries, const RBMAP& runBoundaries);
void SelectPrintHVIssuesFromClass(const o2::mch::HVStatusCreator::BADHVMAP& hvIssues, const RBMAP& runBoundaries);

int main() {
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

    // Compare the macros from per object manual -> per object by class
    SelectHVBoundaries(hvBoundaries, runBoundaries);
    PrintHVBoundaries(hvBoundaries);

    std::cout << "==== Issues found from per object MACRO ====\n";
    FindHVIssuesPerObjectMacro(api, hvBoundaries, runBoundaries, true);

    std::cout << "==== Issues found from per object through CLASS ==== \n";
    FindHVIssuesWithClassMacro(api, hvBoundaries, runBoundaries);
    
    return 0;
}

void FindHVIssuesWithClassMacro(o2::ccdb::CcdbApi const& api, HVBMAP& hvBoundaries, const RBMAP& runBoundaries) {
  std::map<std::string, std::string> metadata;

  for (auto boundaries : hvBoundaries) {
    auto* HV_map = api.retrieveFromTFileAny<DPMAP>("MCH/Calib/HV", metadata, boundaries.first);

    // Create HVStatusCreator object
    auto* HVStatusCreator = new o2::mch::HVStatusCreator();
    HVStatusCreator->findAllIssues(*HV_map);

    // remove issues that are not within a run range (to get rid of ramp up/down)
    auto hvIssues = HVStatusCreator->getHVIssuesList();

    printf("Issues of HV Object: %lld - %lld\n", boundaries.first, boundaries.second);
    SelectPrintHVIssuesFromClass(hvIssues, runBoundaries);
    std::cout << std::endl;
  };
};

void FindHVIssuesPerObjectMacro(o2::ccdb::CcdbApi const& api, HVBMAP& hvBoundaries, const RBMAP& runBoundaries, bool plotHVObjects) {
  // Case of iterating through one HV object at a time and spot HV issue
  std::map<std::string, std::string> metadata;

  for (auto boundaries : hvBoundaries) {
    auto* HV_map = api.retrieveFromTFileAny<DPMAP>("MCH/Calib/HV", metadata, boundaries.first);

    DPMAP2 dpsMapsPerCh[10];

    // Fill HV values for single object
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
    };

    SelectDataPoints(dpsMapsPerCh, boundaries.first, boundaries.second);

    // printf("%lld - %lld (%s - %s)\n", boundaries.first, boundaries.second, GetTime(boundaries.first).c_str(), GetTime(boundaries.second).c_str());
    // PrintDataPoints(dpsMapsPerCh, true);
    // std::cout << std::endl;

    // // Find HV Issues in one object
    BADHVMAP hvIssuesPerCh[10];

    for (int ch = 0; ch < 10; ++ch) {
      // iterate through each alias of the channel
      for (const auto& [alias, HV_values] : dpsMapsPerCh[ch]) {
        FindHVIssues(HV_values, hvLimits[ch], hvIssuesPerCh[ch][alias]);
      }
      SelectHVIssues(hvIssuesPerCh[ch], runBoundaries, minDuration); // do a selection without run boundaries
    }

    printf("Issues of HV Object: %lld - %lld\n", boundaries.first, boundaries.second);
    PrintHVIssues(hvIssuesPerCh);
    std::cout << std::endl;

    // Plotting
    if (plotHVObjects) {
      auto outFile = new TFile(fmt::format("Plots/HV_ouput_{}.root", boundaries.first).c_str(), "RECREATE");

      TMultiGraph* mg[10];
      for (int ch = 0; ch < 10; ++ch) {
        mg[ch] = new TMultiGraph;
        mg[ch]->SetNameTitle(fmt::format("ch{}", ch).c_str(), fmt::format("Chamber {};time;HV (V)", ch + 1).c_str());
        for (const auto& [alias, dps] : dpsMapsPerCh[ch]) {
          mg[ch]->Add(MapToGraph(alias, dps), "lp");
        }
      }

      // Draw canvas for each channel
      for (int ch = 0; ch < 10; ++ch) {
        TCanvas* c = DrawDataPoints(mg[ch]);
        //DrawRunBoudaries(runBoundaries, c);
        c->Write();
      }

      delete outFile;
    }
  };
}

void FindHVIssuesMacro(o2::ccdb::CcdbApi const& api, HVBMAP& hvBoundaries, const RBMAP& runBoundaries) {
  // Case of iterating through all HV points continuously to spot issues

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

  //PrintDataPoints(dpsMapsPerCh, true);

  SelectDataPoints(dpsMapsPerCh, runBoundaries.begin()->second.first, runBoundaries.rbegin()->second.second);
  //PrintDataPoints(dpsMapsPerCh, true);

  // Look for HV issues for each channel
  BADHVMAP hvIssuesPerCh[10];
  
  for (int ch = 0; ch < 10; ++ch) {
    // iterate through each alias of the channel
    for (const auto& [alias, HV_values] : dpsMapsPerCh[ch]) {
      FindHVIssues(HV_values, hvLimits[ch], hvIssuesPerCh[ch][alias]); // How does it know alias entries already? Need to initialize?
    }
    SelectHVIssues(hvIssuesPerCh[ch], runBoundaries, 2);
  }

  PrintHVIssues(hvIssuesPerCh);

  // Plotting
  auto outFile = new TFile("HVLV_ouput.root", "RECREATE");

  TMultiGraph* mg[10];
  for (int ch = 0; ch < 10; ++ch) {
    mg[ch] = new TMultiGraph;
    mg[ch]->SetNameTitle(fmt::format("ch{}", ch).c_str(), fmt::format("chamber {};time;HV (V)", ch + 1).c_str());
    for (const auto& [alias, dps] : dpsMapsPerCh[ch]) {
      mg[ch]->Add(MapToGraph(alias, dps), "lp");
    }
  }

  // Draw canvas for each channel
  for (int ch = 0; ch < 10; ++ch) {
    TCanvas* c = DrawDataPoints(mg[ch]);
    DrawRunBoudaries(runBoundaries, c);
    c->Write();
  }

  delete outFile;
}

void FindHVIssues(const HVVALUES& HV_values, double hvLimit, BADHVLIST& hvIssuesList) 
{
  uint64_t start_stamp = 0;
  uint64_t end_stamp = 0;
  bool ongoing_issue = false;
  double hv_mean = 0.0;
  double hv_min = 0.0;
  uint32_t hv_count = 0;
  uint64_t boundaryExtension = 10000;
  bool isFirst = true;

  for (auto& [timestamp, HV_value] : HV_values) {
    if (HV_value < hvLimit) {
      if (!ongoing_issue) {
        start_stamp = timestamp;
        // if it is the first point in the object, extend it backwards in time?
        if (isFirst) {
          //start_stamp = start_stamp - boundaryExtension;
          start_stamp = 0;
        }
        end_stamp = start_stamp;
        hv_mean = HV_value;
        hv_min = HV_value;
        hv_count = 1;
        ongoing_issue = true;
      } else {
        end_stamp = timestamp;
        hv_mean += HV_value;
        if (HV_value < hv_min) {
          hv_min = HV_value;
        }
        hv_count += 1;
      }
    } else {
      if (ongoing_issue) {
        end_stamp = timestamp;
        hv_mean = hv_mean/hv_count;
        // Possibly do the min duration selection here?
        hvIssuesList.emplace_back(start_stamp, end_stamp, hv_mean, hv_min, hv_count, "");
        ongoing_issue = false;
        hv_mean = 0.0;
        hv_min = 0.0;
        hv_count = 0;
      }
    }
    isFirst = false;
  };

  // The issue is still ongoing at the end of the object boundary, extend it forwards in time?
  if (ongoing_issue && (start_stamp != end_stamp)) {
    hv_mean = hv_mean/hv_count;
    //end_stamp = end_stamp + boundaryExtension;
    end_stamp = 9999999999999;
    hvIssuesList.emplace_back(start_stamp, end_stamp, hv_mean, hv_min, hv_count, "");
  };
};

void SelectHVIssues(BADHVMAP& hvIssuesList, const RBMAP& runBoundaries, uint64_t minDuration) 
{
  // itterate through all identified HV issues to select only those within a run range 
  for (auto& hvIssues: hvIssuesList) {
    if (!hvIssues.second.empty()) {
      for (auto itIssue = hvIssues.second.begin(); itIssue != hvIssues.second.end();) {
        auto tStart = std::get<0>(*itIssue);
        auto tStop = std::get<1>(*itIssue);

        if ((tStop - tStart) < minDuration) {
          itIssue = hvIssues.second.erase(itIssue);
        } else {
          // itterate through run boundaries
          bool inRunRange = false;
          std::string runs = "";
          for (const auto& [run, boundaries] : runBoundaries) {
            if (tStart >= boundaries.second) {
              continue;
            } else if (tStop <= boundaries.first) {
              break;
            }
            inRunRange = true;
            runs += fmt::format("{},", run);
          }

          if (!inRunRange) {
            itIssue = hvIssues.second.erase(itIssue);
          } else {
            std::get<5>(*itIssue) = runs;
            ++itIssue;
          }
        }
      }
    }
  }
}

void SelectPrintHVIssuesFromClass(const o2::mch::HVStatusCreator::BADHVMAP& hvIssues, const RBMAP& runBoundaries) {
  for (const auto& [alias, timeRanges] : hvIssues) {
    bool foundIssue = false;

    for (const auto& timeRange : timeRanges) {
      auto tStart = timeRange.begin;
      auto tStop = timeRange.end;

      // itterate through run boundaries
      bool inRunRange = false;
      for (const auto& [run, boundaries] : runBoundaries) {
        if (tStart >= boundaries.second) {
          continue;
        } else if (tStop <= boundaries.first) {
          break;
        }
        inRunRange = true;
      }

      if (inRunRange) {
        if (!foundIssue) {
          std::cout << "Issues for the alias: " << alias << std::endl;
          foundIssue = true;
        }
        std::cout << tStart << " -> " << tStop << std::endl;
      }
    }
  }
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

void SelectHVBoundaries(HVBMAP& hvBoundaries, const RBMAP& runBoundaries) 
{
  // Only select HV objects that fall within a run range
  for (auto itBound = hvBoundaries.begin(), itNext = itBound; itBound != hvBoundaries.end(); itBound = itNext) {
    auto tStart = itBound->first;
    auto tStop = itBound->second;
    bool inRunRange = false;

    for (const auto& [run, boundaries] : runBoundaries) {
      if (tStart >= boundaries.second) {
        continue;
      } else if (tStop <= boundaries.first) {
        break;
      }
      inRunRange = true;
    }

    ++itNext;
    if (!inRunRange) {
      hvBoundaries.erase(itBound);
    }
  }
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

void SelectDataPoints(DPMAP2 dpsMapsPerCh[10], uint64_t tStart, uint64_t tStop)
{
  /// remove the data points outside of the given time range and, if needed,
  /// add a data point at the boundaries with HV equal to the preceding value

  for (int ch = 0; ch < 10; ++ch) {
    for (auto& [alias, dps] : dpsMapsPerCh[ch]) {

      // get the first data point in the time range, remove the previous ones
      // and add a data point with HV equal to the preceding value if needed
      auto itFirst = dps.lower_bound(tStart);
      if (itFirst != dps.begin()) {
        double previousHV = std::prev(itFirst)->second;
        for (auto it = dps.begin(); it != itFirst;) {
          it = dps.erase(it);
        }
        dps.emplace(tStart, previousHV);
      } else if (itFirst->first != tStart) {
        printf("error (%s) at object (%lld): first data point is after the start of the time range\n", alias.c_str(), tStart);
        // The first point is after start of object - don't know what was before.
      }

      // get the first data point exceeding the time range, remove it and the next ones
      // and add a data point with HV equal to the preceding value if needed
      auto itLast = dps.upper_bound(tStop);
      if (itLast != dps.begin()) {
        double previousHV = std::prev(itLast)->second;
        for (auto it = itLast; it != dps.end();) {
          it = dps.erase(it);
        }
        dps.emplace(tStop, previousHV);
      } else {
        printf("error (%s) at object (%lld): all data points are after the end of the time range\n", alias.c_str(), tStart);
        //dps.clear();
      }
    }
  }
}

TGraph* MapToGraph(std::string alias, const std::map<uint64_t, double>& dps)
{
  /// create a graph for the DCS channel and add the data points

  TGraph* g = new TGraph(dps.size());

  auto shortAlias = alias.substr(0, alias.size() - 12);
  auto title = fmt::format("{} ({})", GetDE(alias).c_str(), shortAlias.c_str());
  g->SetNameTitle(alias.c_str(), title.c_str());

  int i(0);
  for (auto [ts, hv] : dps) {
    g->SetPoint(i, MSToS(ts), hv);
    ++i;
  }

  g->SetMarkerSize(1.5);
  g->SetMarkerStyle(2);
  g->SetLineStyle(2);

  return g;
}

TCanvas* DrawDataPoints(TMultiGraph* mg)
{
  /// display the data points of the given chamber

  TCanvas* c = new TCanvas(mg->GetName(), mg->GetHistogram()->GetTitle(), 1500, 900);

  mg->Draw("A plc pmc");
  mg->SetMinimum(yRange[0]);
  mg->SetMaximum(yRange[1]);
  mg->GetXaxis()->SetTimeDisplay(1);
  mg->GetXaxis()->SetTimeFormat("%d/%m %H:%M");
  mg->GetXaxis()->SetTimeOffset(0, "local");
  mg->GetXaxis()->SetNdivisions(21010);

  c->BuildLegend();
  c->Update();

  return c;
}

std::string GetDE(std::string alias)
{
  /// get the DE (and sector) corresponding to the DCS alias

  auto de = o2::mch::dcs::aliasToDetElemId(alias);

  return (o2::mch::dcs::isQuadrant(o2::mch::dcs::aliasToChamber(alias)))
           ? fmt::format("DE{}-{}", *de, o2::mch::dcs::aliasToNumber(alias) % 10)
           : fmt::format("DE{}", *de);
}

void DrawRunBoudaries(const RBMAP& runBoundaries, TCanvas* c)
{
  /// draw the run time boundaries

  c->cd();

  for (const auto& [run, boundaries] : runBoundaries) {

    TLine* startRunLine = new TLine(MSToS(boundaries.first), yRange[0], MSToS(boundaries.first), yRange[1]);
    startRunLine->SetUniqueID(run);
    startRunLine->SetLineColor(4);
    startRunLine->SetLineWidth(1);
    startRunLine->Draw();

    TLine* endRunLine = new TLine(MSToS(boundaries.second), yRange[0], MSToS(boundaries.second), yRange[1]);
    endRunLine->SetUniqueID(run);
    endRunLine->SetLineColor(2);
    endRunLine->SetLineWidth(1);
    endRunLine->Draw();
  }
}

void PrintHVIssues(const BADHVMAP hvIssuesPerCh[10])
{
  for (int ch = 0; ch < 10; ++ch) {
    for (const auto& [alias, hv_issues]: hvIssuesPerCh[ch]) {
      if (!hv_issues.empty()) {
        printf("Issues for the alias: %s (Chamber %d)\n", alias.c_str(), (ch+1));
        for (const auto& [start, stop, mean, min, count, runs] : hv_issues) {
          printf("%lld -> %lld \n", start, stop);
          //printf("Start: %s; End: %s; Mean: %f; Min: %f; Count: %d; Runs: %s \n", GetTime(start).c_str(), GetTime(stop).c_str(), mean, min, count, (runs).c_str());
        }
      }
    }
  }
}