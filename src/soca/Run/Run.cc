/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Run/Run.h"

#include "oops/util/Logger.h"
#include "oops/runs/Run.h"
#include "eckit/config/Configuration.h"

namespace soca {

// -----------------------------------------------------------------------------

Run::Run(int argc, char ** argv): oops::Run(argc, argv) {
  oops::Log::trace() << "Creating SOCA Run" << std::endl;
  const eckit::Configuration * conf = &config();
  oops::Log::trace() << "SOCA Run created" << std::endl;
}

// -----------------------------------------------------------------------------

Run::~Run() {
  oops::Log::trace() << "Destroying SOCA Run" << std::endl;
  oops::Log::trace() << "SOCA Run destroyed" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace soca
