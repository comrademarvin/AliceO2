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

/// \file GPUTPCRow.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCRow.h"
using namespace o2::gpu;

#if !defined(GPUCA_GPUCODE)
GPUTPCRow::GPUTPCRow() : mNHits(0), mX(0), mMaxY(0), mGrid(), mHy0(0), mHz0(0), mHstepY(0), mHstepZ(0), mHstepYi(0), mHstepZi(0), mHitNumberOffset(0), mFirstHitInBinOffset(0)
{
  // dummy constructor
}

#endif
