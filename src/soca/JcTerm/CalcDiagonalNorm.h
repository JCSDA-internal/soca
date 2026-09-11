/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

namespace soca {

class CalcDiagonalNorm {
 public:
  CalcDiagonalNorm(const State &,
                   const Geometry &,
                   Increment &,
                   Increment &,
                   const eckit::Configuration &) {
    ABORT("CalcDiagonalNorm::CalcDiagonalNorm not implemented.");
  }
};

}  // namespace soca

