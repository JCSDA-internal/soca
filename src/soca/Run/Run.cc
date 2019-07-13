/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Run/Run.h"

#include "oops/util/Logger.h"
#include "oops/runs/Run.h"
#include "eckit/config/Configuration.h"

#include "soca/Fortran.h"

namespace soca {

// -----------------------------------------------------------------------------

Run::Run(int argc, char ** argv): oops::Run(argc, argv) {
  oops::Log::trace() << "Creating Run" << std::endl;
  const eckit::Configuration * conf = &config();
  soca_setup_f(&conf);
  oops::Log::trace() << "Run created" << std::endl;
}

// -----------------------------------------------------------------------------

Run::~Run() {
  oops::Log::trace() << "Destructing Run" << std::endl;
  soca_finalize_f();
  oops::Log::trace() << "MPI finalized, Run destructed" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace soca
