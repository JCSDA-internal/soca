
#include "model/ObsBiasIncrement.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "util/Logger.h"
#include "model/ObsBias.h"
#include "model/ObsBiasCovariance.h"
#include "eckit/config/Configuration.h"

using oops::Log;


// -----------------------------------------------------------------------------
namespace mom5cice5 {
// -----------------------------------------------------------------------------
ObsBiasIncrement::ObsBiasIncrement(const eckit::Configuration & conf)
  : bias_(ObsBias::ntypes, 0.0), active_(ObsBias::ntypes, false)
{
  active_[0] = conf.has("fraction");
  active_[1] = conf.has("freeboard");
  active_[2] = conf.has("temp");
  active_[3] = conf.has("salt");
  bool on = false;
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (active_[jj]) {
      strn += "on";
      on = true;
    } else {
      strn += "off";
    }
  }
  if (on) {Log::trace() << "ObsBiasIncrement created : " << strn << std::endl;}
}
// -----------------------------------------------------------------------------
ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const bool copy)
  : bias_(ObsBias::ntypes, 0.0), active_(other.active_)
{
  if (copy) {
    for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = other.bias_[jj];
  }
  this->makePassive();
}
// -----------------------------------------------------------------------------
ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const eckit::Configuration &)
  : bias_(ObsBias::ntypes, 0.0), active_(other.active_)
{
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = other.bias_[jj];
  this->makePassive();
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::makePassive() {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (!active_[jj]) bias_[jj] = 0.0;
  }
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    bias_[jj] = b1[jj] - b2[jj];
  }
  this->makePassive();
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::zero() {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = 0.0;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] = rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] += rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] -= rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] *= fact;
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) bias_[jj] += fact * rhs.bias_[jj];
  this->makePassive();
}
// -----------------------------------------------------------------------------
double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (active_[jj]) zz += bias_[jj] * rhs.bias_[jj];
  }
  return zz;
}
// -----------------------------------------------------------------------------
double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  int ii = 0;
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (active_[jj]) {
      zz += bias_[jj] * bias_[jj];
      ++ii;
    }
  }
  if (ii>0) zz = std::sqrt(zz/ii);
  return zz;
}
// -----------------------------------------------------------------------------
void ObsBiasIncrement::print(std::ostream & os) const {
  bool on = false;
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (active_[jj]) {
      on = true;
      std::ostringstream strs;
      strs << bias_[jj];
      strn += strs.str();
    } else {
      strn += "0.0";
    }
  }
  if (on) os << std::endl << "ObsBiasIncrement = " << strn;
}
// -----------------------------------------------------------------------------
}  // namespace mom5cice5
