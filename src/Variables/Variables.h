/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_SOCAVARIABLES_H_
#define SOCA_MODEL_SOCAVARIABLES_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace soca {

// -----------------------------------------------------------------------------

class Variables : public util::Printable,
                    private util::ObjectCounter<Variables> {
 public:
  static const std::string classname() {return "soca::Variables";}

  explicit Variables(const oops::Variables &);
  explicit Variables(const eckit::Configuration &);

  ~Variables();

  Variables(const Variables &);

  const int * toFortran() const {return &fvars_[0];}

 private:
  void print(std::ostream &) const;
  void setF90(const std::vector<std::string>);
  std::vector<int> fvars_;
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif    // SOCA_MODEL_SOCAVARIABLES_H_

