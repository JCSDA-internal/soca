
#include "model/ErrorCovariance.h"

#include <cmath>

#include "util/Logger.h"
#include "model/Fields.h"
#include "model/Fortran.h"
#include "model/Geometry.h"
#include "model/Increment.h"
#include "model/State.h"
#include "eckit/config/Configuration.h"


using oops::Log;

// -----------------------------------------------------------------------------
namespace mom5cice5 {
  // -----------------------------------------------------------------------------

  ErrorCovariance::ErrorCovariance(const Geometry & resol, const Variables &,
				   const eckit::Configuration & conf, const State &) {
    time_ = util::DateTime(conf.getString("date"));
    const eckit::Configuration * configc = &conf;
    //mom5cice5_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
    Log::trace() << "ErrorCovariance created" << std::endl;
  }

  // -----------------------------------------------------------------------------

  ErrorCovariance::~ErrorCovariance() {
    //mom5cice5_b_delete_f90(keyFtnConfig_);
    Log::trace() << "ErrorCovariance destructed" << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::linearize(const State &, const Geometry & resol) {
    geom_.reset(new Geometry(resol));
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::multiply(const Increment & dxin, Increment & dxout) const {
    //mom5cice5_b_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
    //			 dxout.fields().toFortran());
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::inverseMultiply(const Increment & dxin, Increment & dxout) const {
    //mom5cice5_b_invmult_f90(keyFtnConfig_, dxin.fields().toFortran(),
    //			    dxout.fields().toFortran());
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::randomize(Increment & dx) const {
    //mom5cice5_b_randomize_f90(keyFtnConfig_, dx.fields().toFortran());
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::print(std::ostream & os) const {
    os << "ErrorCovariance::print not implemented";
  }

  // -----------------------------------------------------------------------------

}  // namespace mom5cice5
