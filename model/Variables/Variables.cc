/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "Variables.h"

#include<vector>

#include "oops/base/Variables.h"
#include "eckit/config/Configuration.h"
#include "util/Logger.h"

namespace soca {

// -----------------------------------------------------------------------------

Variables::Variables(const oops::Variables & oopsvars) {
  oops::Log::debug() << "Variables oopsvar:" << oopsvars.variables() << std::endl;
  this->setF90(oopsvars.variables());
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & config) {
  oops::Log::debug() << "Variables config:" << config << std::endl;
  std::vector<std::string> vars;
  config.get("variables", vars);
  this->setF90(vars);
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

void Variables::setF90(const std::vector<std::string> vars) {
  size_t nv = vars.size();
  oops::Log::debug() << "setF90 " << nv << " vars = " << vars << std::endl;
  fvars_.resize(nv + 2);
  fvars_[0] = nv;
  for (size_t jj = 0; jj < nv; ++jj) {
     int ii = 0;
     if (vars[jj]=="cicen") ii = 1;
     if (vars[jj]=="hicen") ii = 2;
     if (vars[jj]=="hsnon") ii = 3;
     if (vars[jj]=="tsfcn") ii = 4;
     if (vars[jj]=="qsnon") ii = 5;
     if (vars[jj]=="sicnk") ii = 6;
     if (vars[jj]=="qicnk") ii = 7;
     if (vars[jj]=="socn") ii = 8;     
     if (vars[jj]=="tocn") ii = 9;
     if (vars[jj]=="ssh") ii = 10;
     //std::cout << vars[jj] << std::endl;     
     ASSERT(ii > 0);
     fvars_[jj+1] = ii;
  }
  fvars_[nv+1] = 999;  // just for checking
  oops::Log::debug() << "setF90 " << nv << " fvars = " << fvars_ << std::endl;
}

// -----------------------------------------------------------------------------

Variables::~Variables() {}

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other): fvars_(other.fvars_) {}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  os << "soca::Variables: vars = " << fvars_;
}

// -----------------------------------------------------------------------------

}  // namespace soca


