/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

// #include "atlas/field.h"
// #include "atlas/functionspace.h"

// #include "eckit/config/Configuration.h"

// #include "oops/base/Variables.h"
#include "oops/generic/InterpolatorUnstructured.h"
// #include "oops/util/Logger.h"
#include "oops/util/Printable.h"


namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace soca {
  class Geometry;
  class Increment;
  class State;
}


namespace soca {

// -----------------------------------------------------------------------------

class LocalUnstructuredInterpolator : public util::Printable {
 public:
  LocalUnstructuredInterpolator(const eckit::Configuration &, const Geometry &,
                                const std::vector<double> &);
  ~LocalUnstructuredInterpolator() {}

  void apply(const oops::Variables &, const State &, std::vector<double> &) const;
  void apply(const oops::Variables &, const Increment &, std::vector<double> &) const;
  void applyAD(const oops::Variables &, Increment &, const std::vector<double> &) const;

 private:
  void print(std::ostream &) const;

  std::shared_ptr<const Geometry> geom_;
  std::unique_ptr<oops::InterpolatorUnstructured> interp_[6];
  //const atlas::FunctionSpace * fspace_;
};


// }

// // -----------------------------------------------------------------------------

// template<typename GEOMETRY, typename STATE, typename INCREMENT>
// void LocalUnstructuredInterpolator<GEOMETRY, STATE, INCREMENT>::
//   apply(const Variables & vars, const INCREMENT & dx, std::vector<double> & locvals) const
// {
//   Log::trace() << "LocalUnstructuredInterpolator::apply INCREMENT start" << std::endl;

//   atlas::FieldSet fset;
//   dx.getFieldSet(vars, fset);

//   interp_->apply(fset, locvals);

//   Log::trace() << "LocalUnstructuredInterpolator::apply INCREMENT done" << std::endl;
// }

// // -----------------------------------------------------------------------------

// template<typename GEOMETRY, typename STATE, typename INCREMENT>
// void LocalUnstructuredInterpolator<GEOMETRY, STATE, INCREMENT>::
//   applyAD(const Variables & vars, INCREMENT & dx, const std::vector<double> & locvals) const
// {
//   Log::trace() << "LocalUnstructuredInterpolator::applyAD starting" << std::endl;

//   atlas::FieldSet fset;
//   std::vector<size_t> levs(dx.geometry()->variableSizes(vars));
//   for (size_t jf = 0; jf < vars.size(); ++jf) {
//     const std::string fname = vars[jf];
//     if (levs[jf] > 1) {
//       atlas::Field incfld = fspace_->createField<double>(atlas::option::name(fname) |
//                                                          atlas::option::levels(levs[jf]));
//       fset.add(incfld);
//     } else {
//       atlas::Field incfld = fspace_->createField<double>(atlas::option::name(fname));
//       fset.add(incfld);
//     }
//   }

//   interp_->applyAD(fset, locvals);  // set fields to zero inside

//   dx.getFieldSetAD(vars, fset);  // should add to dx, not over-write it

//   Log::trace() << "LocalUnstructuredInterpolator::applyAD done" << std::endl;
// }

// // -----------------------------------------------------------------------------

// template<typename GEOMETRY, typename STATE, typename INCREMENT>
// void LocalUnstructuredInterpolator<GEOMETRY, STATE, INCREMENT>::print(std::ostream & os) const
// {
//   os << "LocalUnstructuredInterpolator<"
//      << GEOMETRY::classname() << STATE::classname() << INCREMENT::classname()
//      << ">";
// }

// -----------------------------------------------------------------------------

}  // namespace oops
