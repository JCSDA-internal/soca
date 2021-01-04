/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_OBSLOCALIZATION_OBSLOCFROMMETADATA_H_
#define SOCA_OBSLOCALIZATION_OBSLOCFROMMETADATA_H_

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


  class ObsLocFromMetadata: public util::Printable {
   public:
    static const std::string classname() {return "soca::ObsLocFromMetadata";}

    ObsLocFromMetadata(const eckit::Configuration &, const ioda::ObsSpace &);
    void multiply(ioda::ObsVector &) const;

   private:
    void print(std::ostream &) const override {}
    const ioda::ObsSpace & obsdb_;
  };

}  // namespace soca

#endif  // SOCA_OBSLOCALIZATION_OBSLOCFROMMETADATA_H_
