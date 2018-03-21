/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ObsBiasCovariance.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdlib.h>     /* srand, rand */
#include <string>
#include <vector>

#include "util/Logger.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "eckit/config/LocalConfiguration.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace soca {
// -----------------------------------------------------------------------------
ObsBiasCovariance::ObsBiasCovariance(const eckit::Configuration & conf)
  : conf_(conf), variance_(ObsBias::ntypes, 0.0)
{
//  if (!conf.empty()) {
    std::vector<double> zz(4, 0.0);
    if (conf.has("fraction")) zz[0] = conf.getDouble("fraction");
    if (conf.has("freeboard"))  zz[1] = conf.getDouble("freeboard");
    if (conf.has("temp"))  zz[2] = conf.getDouble("temp");
    if (conf.has("salt")) zz[3] = conf.getDouble("salt");

    std::string strn = "";
    for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
      if (jj > 0) strn += ", ";
      if (std::abs(zz[jj]) > 1.0e-8) {
        variance_[jj] = zz[jj] * zz[jj];
        std::ostringstream strs;
        strs << variance_[jj];
        strn += strs.str();
      } else {
        variance_[jj] = 0.0;
        strn += "0.0";
      }
    }
    Log::info() << "ObsBiasCovariance created, variances = " << strn << std::endl;
//  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::multiply(const ObsBiasIncrement & dxin,
                                 ObsBiasIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] * variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & dxin,
                                        ObsBiasIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] / variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::randomize(ObsBiasIncrement & dx) const {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      double zz = static_cast<double>(std::rand()) / RAND_MAX;
      dx[jj] = zz * std::sqrt(variance_[jj]);
    } else {
      dx[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsBiasCovariance::print(std::ostream & os) const {
  os << "ObsBiasCovariance::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace soca
