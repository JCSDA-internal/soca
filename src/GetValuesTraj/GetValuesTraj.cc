/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>
#include "src/GetValuesTraj/GetValuesTraj.h"
#include "oops/util/Logger.h"
#include "Fortran.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace soca {
// -----------------------------------------------------------------------------
GetValuesTraj::GetValuesTraj() {
  oops::Log::trace() << "GetValuesTraj constructor starting"
                     << std::endl;
  soca_getvaltraj_setup_f90(keyGetValuesTraj_);
  oops::Log::trace() << "GetValuesTraj constructor done"
                     << keyGetValuesTraj_ << std::endl;
}
// -----------------------------------------------------------------------------
GetValuesTraj::~GetValuesTraj() {
  oops::Log::trace() << "GetValuesTraj destructor starting"
                     << std::endl;
  soca_getvaltraj_delete_f90(keyGetValuesTraj_);
  oops::Log::trace() << "GetValuesTraj destructor done" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
