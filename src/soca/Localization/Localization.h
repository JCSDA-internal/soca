/*
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_LOCALIZATION_LOCALIZATION_H_
#define SOCA_LOCALIZATION_LOCALIZATION_H_

#include <ostream>
#include <string>

#include "soca/Fortran.h"
#include "soca/Traits.h"

#include "oops/interface/LocalizationBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace soca {
  class Geometry;
  class Increment;
}

// -----------------------------------------------------------------------------

namespace soca {

  // Localization for SOCA model.
  class Localization : public oops::interface::LocalizationBase<soca::Traits> {
   public:
      static const std::string classname() {return "soca::Localization";}

      Localization(const Geometry &, const eckit::Configuration &);
      ~Localization();
      void randomize(Increment &) const override;
      void multiply(Increment &) const override;

   private:
      void print(std::ostream &) const override;
      int keyFtnConfig_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_LOCALIZATION_LOCALIZATION_H_
