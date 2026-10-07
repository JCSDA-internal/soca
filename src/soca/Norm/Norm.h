/*
 * (C) Crown Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

namespace soca {

class Norm : public util::Printable,
             private util::ObjectCounter<Norm> {
 public:
  static const std::string classname() {return "soca::Norm";}

// Constructor, destructor
  Norm(const oops::Variables &,
       const eckit::Configuration &);
  ~Norm() = default;

// Compute values for norm
  void calculate(const State &);

// Apply norm to increment
  void apply(Increment &) const;
  void applyInverse(Increment &) const;

 private:
  void print(std::ostream &) const override;
};

}  // namespace soca
