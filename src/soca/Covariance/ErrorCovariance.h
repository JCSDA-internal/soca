/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_COVARIANCE_ERRORCOVARIANCE_H_
#define SOCA_COVARIANCE_ERRORCOVARIANCE_H_

#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace soca {
  class Geometry;
  class Increment;
  class State;
}

// ----------------------------------------------------------------------------

namespace soca {

  // Background error covariance matrix for SOCA model.
  class ErrorCovariance : public util::Printable,
    private boost::noncopyable,
    private util::ObjectCounter<ErrorCovariance> {
   public:
      static const std::string classname() {return "soca::ErrorCovariance";}

      ErrorCovariance(const Geometry &, const oops::Variables &,
                      const eckit::Configuration &,
                      const State &, const State &);
      ~ErrorCovariance();

      void linearize(const State &, const Geometry &);
      void multiply(const Increment &, Increment &) const;
      void inverseMultiply(const Increment &, Increment &) const;
      void randomize(Increment &) const;

   private:
      void print(std::ostream &) const;
      int keyFtnConfig_;
      boost::scoped_ptr<const Geometry> geom_;
      boost::scoped_ptr<const State> traj_;
      util::DateTime time_;
      oops::Variables vars_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_COVARIANCE_ERRORCOVARIANCE_H_
