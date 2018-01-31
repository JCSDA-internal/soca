/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SOCA_MODEL_LOCALIZATIONMATRIXSOCA_H_
#define SOCA_MODEL_LOCALIZATIONMATRIXSOCA_H_

#include <ostream>
#include <string>
#include <boost/scoped_ptr.hpp>

#include "model/Geometry/Geometry.h"
#include "eckit/config/Configuration.h"
//#include "oops/interface/LocalizationBase.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"

#include "model/Fortran.h"
//#include "model/Traits.h"

// Forward declarations
namespace soca {
  class Geometry;  
  class Increment;

  /// Localization matrix for SOCA model.

  // -----------------------------------------------------------------------------
  class LocalizationMatrix : public util::Printable,
                            private boost::noncopyable,
                            private util::ObjectCounter<LocalizationMatrix>{
  public:
      static const std::string classname() {return "soca::LocalizationMatrix";}

      LocalizationMatrix(const Geometry &, const eckit::Configuration &);
      ~LocalizationMatrix();
      void multiply(Increment &) const;

  private:
      void print(std::ostream &) const;
      int keyFtnConfig_;
    };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_LOCALIZATIONMATRIXSOCA_H_
