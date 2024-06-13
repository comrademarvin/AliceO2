// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <algorithm>

#include "MCHStatus/HVStatusCreator.h"

#include "MCHConditions/DetectionElement.h"
#include "MCHGlobalMapping/Mapper.h"
#include "MCHStatus/StatusMap.h"

float sum(float s, o2::dcs::DataPointValue v);

namespace o2::mch
{

void HVStatusCreator::findBadHVs(const DPMAP& dpMap)
{
  // decode the DCS DPMAP
  DPMAP2 dpsMapPerAlias = decodeDPMAP(dpMap);
  
  // Find list of HV issues per alias 
  for (const auto& [alias, dpsHV] : dpsMapPerAlias) {
    int chamber = o2::mch::dcs::toInt(o2::mch::dcs::aliasToChamber(alias));
    std::vector<TimeRange> hvIssuesList;

    uint64_t tStart = 0;
    uint64_t tStop = 0;
    bool ongoingIssue = false;

    for (auto& [timestamp, valueHV] : dpsHV) {
      if (valueHV < StatusMapCreatorParam::Instance().hvLimits[chamber]) { // check whether HV point is below set threshold for chamber
        if (!ongoingIssue) {
          tStart = timestamp;
          tStop = tStart;
          ongoingIssue = true;
        } else {
          tStop = timestamp;
        }
      } else {
        if (ongoingIssue) {
          tStop = timestamp;
          TimeRange newIssue;
          newIssue.begin = tStart;
          newIssue.end = tStop;
          hvIssuesList.push_back(newIssue);
          ongoingIssue = false;
        }
      }
    }
    // ongoing issue at the end of the object
    if (ongoingIssue && (tStart != tStop)) {
      TimeRange newIssue;
      newIssue.begin = tStart;
      newIssue.end = tStop;
      hvIssuesList.push_back(newIssue);
    }

    // add issues for the alias if non-empty
    if (!hvIssuesList.empty()) {
      this->mBadHVTimeRanges.emplace(alias, hvIssuesList);
    }
  }
}

bool HVStatusCreator::findCurrentBadHVs(uint64_t timestamp)
{
  // list issues at the given time stamp
  std::set<std::string> currentBadHVs{};
  for (const auto& [alias, timeRanges] : mBadHVTimeRanges) {
    auto it = std::find_if(timeRanges.begin(), timeRanges.end(),
                           [timestamp](const TimeRange& timeRange) { return timeRange.contains(timestamp); });
    if (it != timeRanges.end()) {
      currentBadHVs.emplace(alias);
    }
  }

  // check if the list of issues has changed and update it in this case
  if (currentBadHVs != mCurrentBadHVs) {
    mCurrentBadHVs.swap(currentBadHVs);
    return true;
  }

  return false;
}

void HVStatusCreator::updateStatusMap(StatusMap& statusMap)
{
  for (const auto& alias : mCurrentBadHVs) {
    int deId = dcs::aliasToDetElemId(alias).value();
    if (deId < 500) {
      for (auto dsIndex : dcs::aliasToDsIndices(alias)) {
        statusMap.addDS(dsIndex, StatusMap::kBadHV);
      }
    } else {
      statusMap.addDE(deId, StatusMap::kBadHV);
    }
  }
}

HVStatusCreator::DPMAP2 HVStatusCreator::decodeDPMAP(const DPMAP& dpMap) {
  DPMAP2 dpsMapPerAlias;

  for (auto& entry : dpMap) {
    std::string alias = (entry.first).get_alias();

    if (alias.find("vMon") != std::string::npos) { // only consider voltage channels
      auto dpsHV = entry.second;
      auto& dps2 = dpsMapPerAlias[alias];

      // copy first point to the beginning of time
      auto firstPoint = dpsHV.front();
      dps2.emplace(0, sum(0.0, firstPoint));

      for (const auto& value : dpsHV) {
        double valueConverted = sum(0.0, value);
        dps2.emplace(value.get_epoch_time(), valueConverted);
      }

      // copy last point to the end of time
      auto lastPoint = dpsHV.back();
      dps2.emplace(9999999999999, sum(0.0, lastPoint));
    }
  }

  return dpsMapPerAlias;
}

} // namespace o2::mch

// converts DCS data point value to double HV value
float sum(float s, o2::dcs::DataPointValue v) 
{
  union Converter {
    uint64_t raw_data;
    double value;
  } converter;
  converter.raw_data = v.payload_pt1;
  return s + converter.value;
};