/*
 * (C) Crown Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Norm/Norm.h"

#include "oops/util/abor1_cpp.h"

namespace soca {

Norm::Norm(const oops::Variables &,
           const eckit::Configuration &) {
  ABORT("Norm::Norm not implemented.");
}

void Norm::calculate(const State &) {
  ABORT("Norm::calculate not implemented.");
}

void Norm::apply(Increment &) const {
  ABORT("Norm::apply not implemented.");
}

void Norm::applyInverse(Increment &) const {
  ABORT("Norm::applyInverse not implemented.");
}

void Norm::print(std::ostream &) const {
  ABORT("Norm::print not implemented.");
}

}  // namespace soca
