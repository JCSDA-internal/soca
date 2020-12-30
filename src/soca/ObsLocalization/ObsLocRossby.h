/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
#define SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_

#include <string>
#include <vector>

#include "oops/util/Printable.h"

// Forward declarations

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsVector;
  class ObsSpace;
}

namespace soca {


  class ObsLocRossby: public util::Printable {
   public:
    static const std::string classname() {return "soca::ObsLocRossby";}

    ObsLocRossby(const eckit::Configuration &, const ioda::ObsSpace &);
    ~ObsLocRossby() {}
    void multiply(ioda::ObsVector &) const;

   private:
    void print(std::ostream &) const override {}
    const ioda::ObsSpace & obsdb_;
    const double rscale_;

    const double locMult_;
    const double locMin_;
  };

}  // namespace soca

#endif  // SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
