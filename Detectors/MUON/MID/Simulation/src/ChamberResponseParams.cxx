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

/// \file   MID/Simulation/src/ChamberResponseParams.cxx
/// \brief  Implementation of the parameters for MID RPC response
/// \author Diego Stocco <Diego.Stocco at cern.ch>
/// \date   26 April 2018

/// This class implements the parameters for the parameterization of the RPC spatial resolution.
/// The parameters were tuned by Massimiliano Marchisone in his PhD thesis:
/// http://www.theses.fr/2013CLF22406
/// See ChamberResponse for further details

#include "MIDSimulation/ChamberResponseParams.h"

namespace o2
{
namespace mid
{

//______________________________________________________________________________
double ChamberResponseParams::getParA(int cathode, int deId, double hv) const
{
  /// Gets first parameter
  /// \par cathode Cathode
  /// \par deId Detection element ID
  /// \par hv RPC HV in volts
  return mParA[72 * cathode + deId].second * hv + mParA[72 * cathode + deId].first;
}

//______________________________________________________________________________
double ChamberResponseParams::getParC(int cathode, int deId, double hv) const
{
  /// Get third parameter
  /// \par cathode Cathode
  /// \par deId Detection element ID
  /// \par hv RPC HV in volts
  return mParC[72 * cathode + deId].second * hv + mParC[72 * cathode + deId].first;
}

//______________________________________________________________________________
double ChamberResponseParams::getParB(int cathode, int deId) const
{
  /// Gets the second parameter
  /// \par cathode Cathode
  /// \par deId Detection element ID
  return mParB[72 * cathode + deId];
}

//______________________________________________________________________________
void ChamberResponseParams::setParA(int cathode, int deId, double a0, double a1)
{
  /// Sets parameter A
  mParA[72 * cathode + deId].first = a0;
  mParA[72 * cathode + deId].second = a1;
}

//______________________________________________________________________________
void ChamberResponseParams::setParC(int cathode, int deId, double c0, double c1)
{
  /// Sets parameter C
  mParC[72 * cathode + deId].first = c0;
  mParC[72 * cathode + deId].second = c1;
}

//______________________________________________________________________________
void ChamberResponseParams::setParB(int cathode, int deId, double val)
{
  /// Sets parameter B
  mParB[72 * cathode + deId] = val;
}

void ChamberResponseParams::setParAll(int cathode, int deId, double b, double a0, double a1, double c0, double c1)
{
  /// Sets all parameters
  setParB(cathode, deId, b);
  setParA(cathode, deId, a0, a1);
  setParC(cathode, deId, c0, c1);
}

const std::pair<double,double> ChamberResponseParams::getParametersA(int cathode, int deId) const
{
  /// Gets the parameters to compute A
  /// \par cathode Cathode
  /// \par deId Detection element ID
  return mParA[72 * cathode + deId];
}


const std::pair<double,double> ChamberResponseParams::getParametersC(int cathode, int deId) const
{
  /// Gets the parameters to compute C
  /// \par cathode Cathode
  /// \par deId Detection element ID
  return mParC[72 * cathode + deId];
}

ChamberResponseParams createDefaultChamberResponseParams()
{
  /// Creates the default parameters
  ChamberResponseParams params;
  
  params.setParAll(0,0,2.35,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,0,2.08,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,1,2.47,-26.35,9.133,-0.000545,0.0004582);
  params.setParAll(1,1,2.39,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,2,2.48,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,2,2.25,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,3,2.53,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,3,2.08,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,4,2.53,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,4,2.06,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,5,2.72,-29.25,8.736,-0.00075,0.000415);
  params.setParAll(1,5,2.23,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,6,2.41,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,6,2.34,-26.52,9.049,-0.00075,0.000415);
  params.setParAll(0,7,2.56,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,7,2.5,-26.35,9.134,-0.0005037,0.0005275);
  params.setParAll(0,8,2.28,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,8,2.37,-26.35,9.133,-0.0005667,0.0007745);
  params.setParAll(0,9,2.49,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,9,1.79,-38.35,7.733,-0.00075,0.000415);
  params.setParAll(0,10,2.54,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,10,1.64,-39.64,7.6,-0.000674,0.0007972);
  params.setParAll(0,11,2.27,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,11,2.22,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,12,2.47,-26.35,9.134,-0.0006655,0.0005404);
  params.setParAll(1,12,2.11,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,13,2.46,-27.47,8.958,-0.0006267,0.0005422);
  params.setParAll(1,13,2.1,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,14,2.83,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,14,2.33,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,15,2.43,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,15,2.17,-34.42,8.138,-0.00075,0.000415);
  params.setParAll(0,16,2.46,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,16,2.16,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,17,2.5,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,17,2.12,-37.17,7.85,-0.0005,0.00083);
  params.setParAll(0,18,2.31,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,18,1.76,-39.02,7.664,-0.0005,0.00083);
  params.setParAll(0,19,2.24,-31.77,8.428,-0.0006793,0.0005029);
  params.setParAll(1,19,2.19,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,20,2.41,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,20,2.11,-28.34,8.841,-0.00075,0.000415);
  params.setParAll(0,21,2.63,-29.53,8.715,-0.000648,0.0004326);
  params.setParAll(1,21,2.12,-33.83,8.212,-0.00075,0.000415);
  params.setParAll(0,22,2.36,-36.73,7.868,-0.00075,0.000415);
  params.setParAll(1,22,1.91,-40.85,7.367,-0.000576,0.0004216);
  params.setParAll(0,23,2.57,-27.42,8.957,-0.0006155,0.0005939);
  params.setParAll(1,23,2.13,-33.49,8.246,-0.00075,0.000415);
  params.setParAll(0,24,2.44,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,24,2.28,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,25,2.53,-26.35,9.133,-0.0007499,0.0004623);
  params.setParAll(1,25,1.94,-39.2,7.646,-0.00075,0.000415);
  params.setParAll(0,26,2.44,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,26,2.26,-31.13,8.498,-0.00075,0.000415);
  params.setParAll(0,27,2.37,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,27,2.03,-37.08,7.869,-0.0005216,0.0004332);
  params.setParAll(0,28,2.32,-34.29,8.162,-0.00075,0.000415);
  params.setParAll(1,28,2.18,-27.19,8.97,-0.00075,0.000415);
  params.setParAll(0,29,2.58,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,29,2.54,-28.43,8.783,-0.00075,0.000415);
  params.setParAll(0,30,2.36,-37.53,7.773,-0.0005833,0.0006917);
  params.setParAll(1,30,2.18,-31.4,8.473,-0.00075,0.000415);
  params.setParAll(0,31,2.24,-38.07,7.722,-0.00075,0.000415);
  params.setParAll(1,31,1.89,-38.91,7.675,-0.00075,0.000415);
  params.setParAll(0,32,2.55,-43.92,7.104,-0.0007311,0.0004315);
  params.setParAll(1,32,2.13,-35.35,8.048,-0.00075,0.000415);
  params.setParAll(0,33,2.48,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,33,2.34,-29.46,8.688,-0.00075,0.000415);
  params.setParAll(0,34,2.29,-34.71,8.116,-0.00075,0.000415);
  params.setParAll(1,34,2.15,-34.87,8.101,-0.00075,0.000415);
  params.setParAll(0,35,2.33,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,35,2.21,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,36,2.33,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,36,2.09,-26.35,9.134,-0.0005,0.00083);
  params.setParAll(0,37,2.34,-33.36,8.255,-0.00075,0.000415);
  params.setParAll(1,37,1.8,-38.81,7.685,-0.0005002,0.0008277);
  params.setParAll(0,38,2.33,-26.35,9.133,-0.000629,0.0004404);
  params.setParAll(1,38,2.16,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,39,2.72,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,39,2.37,-27.78,8.859,-0.00075,0.000415);
  params.setParAll(0,40,2.56,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,40,1.89,-33.84,8.203,-0.00075,0.000415);
  params.setParAll(0,41,2.17,-27.05,9.003,-0.0007011,0.0006288);
  params.setParAll(1,41,1.85,-26.35,9.133,-0.0007499,0.0004572);
  params.setParAll(0,42,2.43,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,42,2.22,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,43,2.36,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,43,2.11,-31.59,8.448,-0.00075,0.000415);
  params.setParAll(0,44,2.31,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,44,2.05,-36,7.975,-0.0005001,0.00083);
  params.setParAll(0,45,2.43,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,45,1.93,-38.83,7.683,-0.0007498,0.000415);
  params.setParAll(0,46,2.29,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,46,2.05,-35.16,8.062,-0.0007498,0.000415);
  params.setParAll(0,47,2.31,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,47,2.23,-31.62,8.438,-0.00075,0.000415);
  params.setParAll(0,48,2.46,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,48,2.2,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,49,2.52,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,49,1.82,-35.25,8.051,-0.00075,0.000415);
  params.setParAll(0,50,2.37,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,50,2.03,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,51,2.46,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,51,2.3,-30.99,8.514,-0.00075,0.000415);
  params.setParAll(0,52,2.26,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,52,2.14,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,53,2.39,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,53,2.04,-38.1,7.756,-0.00075,0.000415);
  params.setParAll(0,54,2.35,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,54,1.93,-39.02,7.664,-0.00075,0.000415);
  params.setParAll(0,55,2.3,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(1,55,2.14,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,56,2.31,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,56,2.07,-31.52,8.459,-0.00075,0.000415);
  params.setParAll(0,57,2.48,-30.63,8.649,-0.0005833,0.0006917);
  params.setParAll(1,57,2.02,-35.12,8.073,-0.00075,0.000415);
  params.setParAll(0,58,2.53,-28.95,8.772,-0.0007262,0.0005205);
  params.setParAll(1,58,2.37,-30.49,8.564,-0.00075,0.000415);
  params.setParAll(0,59,2.66,-27.52,8.998,-0.00075,0.000415);
  params.setParAll(1,59,2.24,-34.25,8.176,-0.00075,0.000415);
  params.setParAll(0,60,2.32,-32.05,8.399,-0.0005019,0.0004276);
  params.setParAll(1,60,1.98,-36.41,7.935,-0.00075,0.000415);
  params.setParAll(0,61,2.33,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,61,2.13,-26.37,9.11,-0.0007244,0.0007135);
  params.setParAll(0,62,2.25,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,62,1.92,-38.08,7.76,-0.00075,0.000415);
  params.setParAll(0,63,2.31,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,63,1.71,-41.36,7.339,-0.0005,0.00083);
  params.setParAll(0,64,2.31,-26.35,9.134,-0.0006399,0.0004901);
  params.setParAll(1,64,2.14,-29.83,8.648,-0.00075,0.000415);
  params.setParAll(0,65,2.29,-28.42,8.814,-0.00075,0.000415);
  params.setParAll(1,65,2.36,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(0,66,2.13,-36.55,7.901,-0.0005833,0.0006917);
  params.setParAll(1,66,2,-33.82,8.225,-0.00075,0.000415);
  params.setParAll(0,67,2.47,-34.45,8.096,-0.0007404,0.0004244);
  params.setParAll(1,67,2.17,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,68,2.36,-35.06,8.022,-0.000744,0.0006501);
  params.setParAll(1,68,2.09,-30.66,8.496,-0.00075,0.000415);
  params.setParAll(0,69,2.29,-26.35,9.134,-0.0005405,0.0004364);
  params.setParAll(1,69,2.05,-28.02,8.847,-0.00075,0.000415);
  params.setParAll(0,70,2.18,-34.07,8.18,-0.0006373,0.0004252);
  params.setParAll(1,70,2.3,-26.35,9.133,-0.00075,0.000415);
  params.setParAll(0,71,2.22,-26.35,9.134,-0.00075,0.000415);
  params.setParAll(1,71,1.96,-35.16,8.059,-0.00075,0.000415);

  // old parameters for fitting, to be removed when the new ones will be fully validated
  //   // BP
  // // MT11R
  // params.setParAll(0, 0, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 1, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 2, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 3, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 4, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 5, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 6, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 7, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 8, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT12R
  // params.setParAll(0, 9, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 10, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 11, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 12, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 13, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 14, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 15, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 16, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 17, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT21R
  // params.setParAll(0, 18, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 19, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 20, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 21, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 22, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 23, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 24, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 25, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 26, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT22R
  // params.setParAll(0, 27, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 28, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 29, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 30, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 31, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 32, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 33, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 34, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 35, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT11L
  // params.setParAll(0, 36, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 37, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 38, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 39, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 40, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 41, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 42, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 43, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 44, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT12L
  // params.setParAll(0, 45, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 46, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 47, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 48, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 49, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 50, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 51, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 52, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 53, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT21L
  // params.setParAll(0, 54, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 55, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 56, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 57, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 58, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 59, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 60, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 61, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 62, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT22L
  // params.setParAll(0, 63, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 64, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 65, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 66, 1.72, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 67, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 68, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 69, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 70, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(0, 71, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);

  // // NBP
  // // MT11R
  // params.setParAll(1, 0, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 1, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 2, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 3, 1.72, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 4, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 5, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 6, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 7, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 8, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT12R
  // params.setParAll(1, 9, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 10, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 11, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 12, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 13, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 14, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 15, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 16, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 17, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT21R
  // params.setParAll(1, 18, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 19, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 20, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 21, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 22, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 23, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 24, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 25, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 26, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT22R
  // params.setParAll(1, 27, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 28, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 29, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 30, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 31, 1.72, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 32, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 33, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 34, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 35, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT11L
  // params.setParAll(1, 36, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 37, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 38, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 39, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 40, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 41, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 42, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 43, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 44, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT12L
  // params.setParAll(1, 45, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 46, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 47, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 48, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 49, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 50, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 51, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 52, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 53, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT21L
  // params.setParAll(1, 54, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 55, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 56, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 57, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 58, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 59, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 60, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 61, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 62, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // // MT22L
  // params.setParAll(1, 63, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 64, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 65, 2.47, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 66, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 67, 2.22, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 68, 1.72, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 69, 1.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 70, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);
  // params.setParAll(1, 71, 2.97, -52.70, 6.089 / 1000., -0.5e-3, 8.3e-4 / 1000.);

  return std::move(params);
}

} // namespace mid
} // namespace o2
