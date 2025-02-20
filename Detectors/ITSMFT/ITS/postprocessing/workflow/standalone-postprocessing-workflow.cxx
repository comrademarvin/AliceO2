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

#include "Framework/CallbacksPolicy.h"
#include "CommonUtils/ConfigurableParam.h"
#include "Framework/CompletionPolicy.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CompletionPolicyHelpers.h"
#include "GlobalTrackingWorkflowHelpers/InputHelper.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"

// Include studies hereafter
#include "ITSStudies/ImpactParameter.h"
#include "ITSStudies/AvgClusSize.h"
#include "ITSStudies/PIDStudy.h"
#include "ITSStudies/AnomalyStudy.h"
#include "ITSStudies/Efficiency.h"
#include "ITSStudies/TrackCheck.h"
#include "ITSStudies/TrackExtension.h"
#include "Steer/MCKinematicsReader.h"

using namespace o2::framework;
using GID = o2::dataformats::GlobalTrackID;
using DetID = o2::detectors::DetID;

void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  // option allowing to set parameters
  std::vector<o2::framework::ConfigParamSpec> options{
    {"input-from-upstream", VariantType::Bool, false, {"read clusters from the clusterer"}},
    {"track-sources", VariantType::String, std::string{"ITS,ITS-TPC-TRD-TOF,ITS-TPC-TOF,ITS-TPC,ITS-TPC-TRD"}, {"comma-separated list of track sources to use"}},
    {"cluster-sources", VariantType::String, std::string{"ITS"}, {"comma-separated list of cluster sources to use"}},
    {"disable-root-input", VariantType::Bool, false, {"disable root-files input reader"}},
    {"disable-mc", VariantType::Bool, false, {"disable MC propagation even if available"}},
    {"cluster-size-study", VariantType::Bool, false, {"Perform the average cluster size study"}},
    {"pid-study", VariantType::Bool, false, {"Perform the PID study"}},
    {"track-study", VariantType::Bool, false, {"Perform the track study"}},
    {"impact-parameter-study", VariantType::Bool, false, {"Perform the impact parameter study"}},
    {"anomaly-study", VariantType::Bool, false, {"Perform the anomaly study"}},
    {"track-extension-study", VariantType::Bool, false, {"Perform the track extension study"}},
    {"efficiency-study", VariantType::Bool, false, {"Perform the efficiency study"}},
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings ..."}}};
  // o2::raw::HBFUtilsInitializer::addConfigOption(options, "o2_tfidinfo.root");
  std::swap(workflowOptions, options);
}

#include <Framework/runDataProcessing.h>

WorkflowSpec defineDataProcessing(ConfigContext const& configcontext)
{
  WorkflowSpec specs;
  GID::mask_t srcTrc, srcCls;

  o2::conf::ConfigurableParam::updateFromString(configcontext.options().get<std::string>("configKeyValues"));
  auto useMC = !configcontext.options().get<bool>("disable-mc");

  std::shared_ptr<o2::steer::MCKinematicsReader> mcKinematicsReader;
  if (useMC) {
    mcKinematicsReader = std::make_shared<o2::steer::MCKinematicsReader>("collisioncontext.root");
  }
  bool anyStudy{false};
  // Declare specs related to studies hereafter
  if (configcontext.options().get<bool>("impact-parameter-study")) {
    anyStudy = true;
    srcTrc = GID::getSourcesMask(configcontext.options().get<std::string>("track-sources"));
    srcCls = GID::getSourcesMask(configcontext.options().get<std::string>("cluster-sources"));
    o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, srcCls, srcTrc, srcTrc, useMC, srcCls, srcTrc);
    specs.emplace_back(o2::its::study::getImpactParameterStudy(srcTrc, srcCls, useMC));
  }
  if (configcontext.options().get<bool>("cluster-size-study")) {
    anyStudy = true;
    srcTrc = GID::getSourcesMask(configcontext.options().get<std::string>("track-sources"));
    srcCls = GID::getSourcesMask(configcontext.options().get<std::string>("cluster-sources"));
    o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, srcCls, srcTrc, srcTrc, useMC, srcCls, srcTrc);
    specs.emplace_back(o2::its::study::getAvgClusSizeStudy(srcTrc, srcCls, useMC, mcKinematicsReader));
  }
  if (configcontext.options().get<bool>("pid-study")) {
    anyStudy = true;
    srcTrc = GID::getSourcesMask(configcontext.options().get<std::string>("track-sources"));
    srcCls = GID::getSourcesMask(configcontext.options().get<std::string>("cluster-sources"));
    o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, srcCls, srcTrc, srcTrc, useMC, srcCls, srcTrc);
    specs.emplace_back(o2::its::study::getPIDStudy(srcTrc, srcCls, useMC, mcKinematicsReader));
  }
  if (configcontext.options().get<bool>("track-study")) {
    anyStudy = true;
    srcTrc = GID::getSourcesMask(configcontext.options().get<std::string>("track-sources"));
    srcCls = GID::getSourcesMask(configcontext.options().get<std::string>("cluster-sources"));
    if (!configcontext.options().get<bool>("input-from-upstream")) {
      o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, srcCls, srcTrc, srcTrc, useMC, srcCls, srcTrc);
    }
    specs.emplace_back(o2::its::study::getTrackCheckStudy(GID::getSourcesMask("ITS"), GID::getSourcesMask("ITS"), useMC, mcKinematicsReader));
  }
  if (configcontext.options().get<bool>("anomaly-study")) {
    anyStudy = true;
    srcCls = GID::getSourcesMask(configcontext.options().get<std::string>("cluster-sources"));
    if (!configcontext.options().get<bool>("input-from-upstream")) {
      o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, srcCls, srcTrc, srcTrc, useMC, srcCls, srcTrc);
    }
    specs.emplace_back(o2::its::study::getAnomalyStudy(srcCls, useMC));
  }
  if (configcontext.options().get<bool>("track-extension-study")) {
    if (!useMC) {
      LOGP(fatal, "Track Extension Study needs MC!");
    }
    anyStudy = true;
    srcTrc = GID::getSourcesMask(configcontext.options().get<std::string>("track-sources"));
    srcCls = GID::getSourcesMask("ITS");
    o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, srcCls, srcTrc, srcTrc, true, srcCls, srcTrc);
    specs.emplace_back(o2::its::study::getTrackExtensionStudy(srcTrc, srcCls, mcKinematicsReader));
  }
  if (configcontext.options().get<bool>("efficiency-study")) {
    anyStudy = true;
    srcTrc = GID::getSourcesMask(configcontext.options().get<std::string>("track-sources"));
    srcCls = GID::getSourcesMask(configcontext.options().get<std::string>("cluster-sources"));
    if (!configcontext.options().get<bool>("input-from-upstream")) {
      o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, srcCls, srcTrc, srcTrc, useMC, srcCls, srcTrc);
    }
    specs.emplace_back(o2::its::study::getEfficiencyStudy(GID::getSourcesMask("ITS"), GID::getSourcesMask("ITS"), useMC, mcKinematicsReader));
  }
  if (!anyStudy) {
    LOGP(info, "No study selected, dryrunning");
  }

  o2::raw::HBFUtilsInitializer hbfIni(configcontext, specs);
  // write the configuration used for the studies workflow
  o2::conf::ConfigurableParam::writeINI("o2_its_standalone_configuration.ini");

  return std::move(specs);
}
