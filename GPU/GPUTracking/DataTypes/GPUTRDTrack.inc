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

/// \file GPUTRDTrack.inc
/// \author Ole Schmidt, Sergey Gorbunov

#ifndef GPUTRDTRACK_INC_H
#define GPUTRDTRACK_INC_H

#include "GPUTRDTrack.h"

template <typename T>
GPUd() GPUTRDTrack_t<T>::GPUTRDTrack_t()
{
  // default constructor
  initialize();
}

template <typename T>
GPUd() void GPUTRDTrack_t<T>::initialize()
{
  // TODO: Use in-class initialization instead of setting all members to their default values
  mChi2 = 0.f;
  mSignal = -1.f;
  mRefGlobalTrackId = 0;
  mCollisionId = -1;
  mFlags = 0;
  mIsCrossingNeighbor = 0;
  for (int32_t i = 0; i < kNLayers; ++i) {
    mAttachedTracklets[i] = -1;
  }
}

#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "DataFormatsTPC/TrackTPC.h"

template <typename T>
GPUd() GPUTRDTrack_t<T>::GPUTRDTrack_t(const o2::dataformats::TrackTPCITS& t) : T(t)
{
  initialize();
}

template <typename T>
GPUd() GPUTRDTrack_t<T>::GPUTRDTrack_t(const o2::tpc::TrackTPC& t) : T(t)
{
  initialize();
}

template <typename T>
GPUd() GPUTRDTrack_t<T>::GPUTRDTrack_t(const GPUTRDTrack_t<T>& t)
  : T(t), mChi2(t.mChi2), mSignal(t.mSignal), mRefGlobalTrackId(t.mRefGlobalTrackId), mCollisionId(t.mCollisionId), mFlags(t.mFlags), mIsCrossingNeighbor(t.mIsCrossingNeighbor)
{
  // copy constructor
  for (int32_t i = 0; i < kNLayers; ++i) {
    mAttachedTracklets[i] = t.mAttachedTracklets[i];
  }
}

template <typename T>
GPUd() GPUTRDTrack_t<T>::GPUTRDTrack_t(const T& t) : T(t)
{
  // copy constructor from anything
  initialize();
}

template <typename T>
GPUd() GPUTRDTrack_t<T>& GPUTRDTrack_t<T>::operator=(const GPUTRDTrack_t<T>& t)
{
  // assignment operator
  if (&t == this) {
    return *this;
  }
  *(T*)this = t;
  mChi2 = t.mChi2;
  mSignal = t.mSignal;
  mRefGlobalTrackId = t.mRefGlobalTrackId;
  mCollisionId = t.mCollisionId;
  mFlags = t.mFlags;
  mIsCrossingNeighbor = t.mIsCrossingNeighbor;
  for (int32_t i = 0; i < kNLayers; ++i) {
    mAttachedTracklets[i] = t.mAttachedTracklets[i];
  }
  return *this;
}

template <typename T>
GPUd() int32_t GPUTRDTrack_t<T>::getNlayersFindable() const
{
  // returns number of layers in which the track is in active area of TRD
  int32_t retVal = 0;
  for (int32_t iLy = 0; iLy < kNLayers; iLy++) {
    if ((mFlags >> iLy) & 0x1) {
      ++retVal;
    }
  }
  return retVal;
}

template <typename T>
GPUd() int32_t GPUTRDTrack_t<T>::getNmissingConsecLayers(int32_t iLayer) const
{
  // returns number of consecutive layers in which the track was
  // inside the deadzone up to (and including) the given layer
  int32_t retVal = 0;
  while (!getIsFindable(iLayer)) {
    ++retVal;
    --iLayer;
    if (iLayer < 0) {
      break;
    }
  }
  return retVal;
}

#endif // GPUTRDTRACK_INC_H
