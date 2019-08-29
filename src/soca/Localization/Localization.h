/*
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_LOCALIZATION_LOCALIZATION_H_
#define SOCA_LOCALIZATION_LOCALIZATION_H_

#include <ostream>
#include <string>
#include <boost/scoped_ptr.hpp>

#include "soca/Geometry/Geometry.h"
#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

#include "soca/Fortran.h"

// Forward declarations
namespace soca {
  class Geometry;
  class Increment;

  // Localization for SOCA model.

  // -----------------------------------------------------------------------------
  class Localization : public util::Printable,
     private boost::noncopyable,
     private util::ObjectCounter<Localization>{
   public:
      static const std::string classname() {return "soca::Localization";}

      Localization(const Geometry &, const eckit::Configuration &);
      ~Localization();
      void multiply(Increment &) const;

   private:
      void print(std::ostream &) const;
      int keyFtnConfig_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_LOCALIZATION_LOCALIZATION_H_
