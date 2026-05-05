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
double ChamberResponseParams::getParA(double hv) const
{
  /// Gets first parameter
  /// \par hv RPC HV in volts
  return mParA[1] * hv + mParA[0];
}

//______________________________________________________________________________
double ChamberResponseParams::getParC(double hv) const
{
  /// Get third parameter
  /// \par hv RPC HV in volts
  return mParC[1] * hv + mParC[0];
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
void ChamberResponseParams::setParA(double a0, double a1)
{
  /// Sets parameter A
  mParA[0] = a0;
  mParA[1] = a1;
}

//______________________________________________________________________________
void ChamberResponseParams::setParC(double c0, double c1)
{
  /// Sets parameter C
  mParC[0] = c0;
  mParC[1] = c1;
}

//______________________________________________________________________________
void ChamberResponseParams::setParB(int cathode, int deId, double val)
{
  /// Sets parameter B
  mParB[72 * cathode + deId] = val;
}

ChamberResponseParams createDefaultChamberResponseParams()
{
  /// Creates the default parameters
  ChamberResponseParams params;
  params.setParA(-52.70, 6.089 / 1000.);   // par1 in 1/V
  params.setParC(-4.75e-3, 4.0e-4 / 1000.); // par1 in 1/V

  // if (isStreamer) {
  //   mParB.fill(2.966);
  //   return;
  // }

  params.setParB(1, 0, 1.52);
  params.setParB(0, 1, 1.87);
  params.setParB(1, 1, 1.52);
  params.setParB(0, 2, 1.75);
  params.setParB(1, 2, 1.48);
  params.setParB(0, 3, 1.88);
  params.setParB(1, 3, 1.41);
  params.setParB(0, 4, 1.83);
  params.setParB(1, 4, 1.30);
  params.setParB(0, 5, 2.18);
  params.setParB(1, 5, 1.51);
  params.setParB(0, 6, 1.65);
  params.setParB(1, 6, 1.52);
  params.setParB(0, 7, 1.95);
  params.setParB(1, 7, 1.54);
  params.setParB(0, 8, 1.56);
  params.setParB(1, 8, 1.50);
  params.setParB(0, 9, 1.73);
  params.setParB(1, 9, 1.55);
  params.setParB(0, 10, 1.90);
  params.setParB(1, 10, 1.44);
  params.setParB(0, 11, 1.69);
  params.setParB(1, 11, 1.45);
  params.setParB(0, 12, 1.91);
  params.setParB(1, 12, 1.40);
  params.setParB(0, 13, 1.64);
  params.setParB(1, 13, 1.36);
  params.setParB(0, 14, 2.09);
  params.setParB(1, 14, 1.56);
  params.setParB(0, 15, 1.75);
  params.setParB(1, 15, 1.55);
  params.setParB(0, 16, 1.91);
  params.setParB(1, 16, 1.36);
  params.setParB(0, 17, 1.78);
  params.setParB(1, 17, 1.60);
  params.setParB(0, 18, 1.55);
  params.setParB(1, 18, 1.50);
  params.setParB(0, 19, 1.80);
  params.setParB(1, 19, 1.41);
  params.setParB(0, 20, 1.78);
  params.setParB(1, 20, 1.43);
  params.setParB(0, 21, 2.00);
  params.setParB(1, 21, 1.49);
  params.setParB(0, 22, 1.64);
  params.setParB(1, 22, 1.47);
  params.setParB(0, 23, 1.87);
  params.setParB(1, 23, 1.45);
  params.setParB(0, 24, 1.66);
  params.setParB(1, 24, 1.53);
  params.setParB(0, 25, 1.88);
  params.setParB(1, 25, 1.52);
  params.setParB(0, 26, 1.66);
  params.setParB(1, 26, 1.55);
  params.setParB(0, 27, 1.48);
  params.setParB(1, 27, 1.41);
  params.setParB(0, 28, 1.77);
  params.setParB(1, 28, 1.53);
  params.setParB(0, 29, 1.81);
  params.setParB(1, 29, 1.64);
  params.setParB(0, 30, 1.89);
  params.setParB(1, 30, 1.44);
  params.setParB(0, 31, 1.62);
  params.setParB(1, 31, 1.46);
  params.setParB(0, 32, 1.73);
  params.setParB(1, 32, 1.50);
  params.setParB(0, 33, 1.80);
  params.setParB(1, 33, 1.60);
  params.setParB(0, 34, 1.82);
  params.setParB(1, 34, 1.45);
  params.setParB(0, 35, 1.56);
  params.setParB(1, 35, 1.45);
  params.setParB(0, 36, 1.60);
  params.setParB(1, 36, 1.53);
  params.setParB(0, 37, 1.87);
  params.setParB(1, 37, 1.49);
  params.setParB(0, 38, 1.73);
  params.setParB(1, 38, 1.46);
  params.setParB(0, 39, 2.05);
  params.setParB(1, 39, 1.58);
  params.setParB(0, 40, 1.94);
  params.setParB(1, 40, 1.34);
  params.setParB(0, 41, 1.71);
  params.setParB(1, 41, 1.25);
  params.setParB(0, 42, 1.77);
  params.setParB(1, 42, 1.48);
  params.setParB(0, 43, 1.86);
  params.setParB(1, 43, 1.42);
  params.setParB(0, 44, 1.58);
  params.setParB(1, 44, 1.53);
  params.setParB(0, 45, 1.66);
  params.setParB(1, 45, 1.62);
  params.setParB(0, 46, 1.39);
  params.setParB(1, 46, 1.47);
  params.setParB(0, 47, 1.65);
  params.setParB(1, 47, 1.54);
  params.setParB(0, 48, 1.89);
  params.setParB(1, 48, 1.43);
  params.setParB(0, 49, 1.95);
  params.setParB(1, 49, 1.32);
  params.setParB(0, 50, 1.69);
  params.setParB(1, 50, 1.36);
  params.setParB(0, 51, 1.73);
  params.setParB(1, 51, 1.60);
  params.setParB(0, 52, 1.31);
  params.setParB(1, 52, 1.44);
  params.setParB(0, 53, 1.66);
  params.setParB(1, 53, 1.73);
  params.setParB(0, 54, 1.59);
  params.setParB(1, 54, 1.67);
  params.setParB(0, 55, 1.74);
  params.setParB(1, 55, 1.45);
  params.setParB(0, 56, 1.71);
  params.setParB(1, 56, 1.45);
  params.setParB(0, 57, 1.82);
  params.setParB(1, 57, 1.43);
  params.setParB(0, 58, 1.70);
  params.setParB(1, 58, 1.54);
  params.setParB(0, 59, 1.86);
  params.setParB(1, 59, 1.46);
  params.setParB(0, 60, 1.81);
  params.setParB(1, 60, 1.56);
  params.setParB(0, 61, 1.75);
  params.setParB(1, 61, 1.48);
  params.setParB(0, 62, 1.50);
  params.setParB(1, 62, 1.44);
  params.setParB(0, 63, 1.52);
  params.setParB(1, 63, 1.40);
  params.setParB(0, 64, 1.70);
  params.setParB(1, 64, 1.43);
  params.setParB(0, 65, 1.76);
  params.setParB(1, 65, 1.55);
  params.setParB(0, 66, 1.89);
  params.setParB(1, 66, 1.42);
  params.setParB(0, 67, 1.68);
  params.setParB(1, 67, 1.46);
  params.setParB(0, 68, 1.92);
  params.setParB(1, 68, 1.33);
  params.setParB(0, 69, 1.71);
  params.setParB(1, 69, 1.46);
  params.setParB(0, 70, 1.84);
  params.setParB(1, 70, 1.50);
  params.setParB(0, 71, 1.50);
  params.setParB(1, 71, 1.53);

  return std::move(params);
}

} // namespace mid
} // namespace o2