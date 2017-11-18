
#ifndef SOCA_MODEL_SOCAERRORCOVARIANCE_H_
#define SOCA_MODEL_SOCAERRORCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "model/Fortran.h"
#include "model/Geometry/Geometry.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace soca {
  class Increment;
  class State;
  class Variables;

  // -----------------------------------------------------------------------------
  /// Background error covariance matrix for SOCA model.

  class ErrorCovariance : public util::Printable,
    private boost::noncopyable,
    private util::ObjectCounter<ErrorCovariance> {
  public:
      static const std::string classname() {return "mom5cic5::ErrorCovariance";}

      ErrorCovariance(const Geometry &, const Variables &, const eckit::Configuration &, const State &);
      ~ErrorCovariance();

      void linearize(const State &, const Geometry &);
      void multiply(const Increment &, Increment &) const;
      void inverseMultiply(const Increment &, Increment &) const;
      void randomize(Increment &) const;

  private:
      void print(std::ostream &) const;
      int keyFtnConfig_;
      boost::scoped_ptr<const Geometry> geom_;
      util::DateTime time_;
    };
  // -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_MODEL_SOCAERRORCOVARIANCE_H_
