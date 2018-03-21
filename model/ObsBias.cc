/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ObsBias.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "util/Logger.h"
#include "model/ObsBiasIncrement.h"
#include "eckit/config/Configuration.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  ObsBias::ObsBias(const eckit::Configuration & conf) : bias_(ntypes, 0.0), active_(false) {
    Log::info() << "ObsBias: conf = " << conf << std::endl;
    active_ = conf.has("Fraction");// || conf.has("freeboard") || conf.has("temp") || conf.has("salt");
    if (active_) {
      if (conf.has("fraction")) bias_[0] = conf.getDouble("fraction");
      //if (conf.has("freeboard"))  bias_[1] = conf.getDouble("freeboard");
      //if (conf.has("temp"))  bias_[2] = conf.getDouble("temp");
      //if (conf.has("salt")) bias_[3] = conf.getDouble("salt");
      std::string strn = "";
      for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
	if (jj > 0) strn += ", ";
	std::ostringstream strs;
	strs << bias_[jj];
	strn += strs.str();
      }
      Log::info() << "ObsBias::ObsBias created, bias = " << strn << std::endl;
    }
  }
  // -----------------------------------------------------------------------------
  ObsBias::ObsBias(const ObsBias & other, const bool copy)
    : bias_(ntypes, 0.0), active_(other.active_)
  {
    if (active_ && copy) {
      for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] = other.bias_[jj];
    }
  }
  // -----------------------------------------------------------------------------
  ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
    if (active_) {
      for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] += dx[jj];
    }
    return *this;
  }
  // -----------------------------------------------------------------------------
  double ObsBias::norm() const {
    double zz = 0.0;
    if (active_) {
      for (unsigned int jj = 0; jj < ntypes; ++jj) zz += bias_[jj]*bias_[jj];
      zz = std::sqrt(zz/3.0);
    }
    return zz;
  }
  // -----------------------------------------------------------------------------
  void ObsBias::print(std::ostream & os) const {
    if (active_) {
      std::string strn = "";
      for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
	if (jj > 0) strn += ", ";
	std::ostringstream strs;
	strs << bias_[jj];
	strn += strs.str();
      }
      os << std::endl << "ObsBias = " << strn;
    }
  }
  // -----------------------------------------------------------------------------
}  // namespace soca

