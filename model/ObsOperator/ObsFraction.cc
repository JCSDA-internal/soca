
#include "model/ObsOperator/ObsFraction.h"

#include "util/Logger.h"
#include "model/ModelAtLocations/Gom.h"
#include "model/ObsBias.h"
#include "model/ObsSpace/ObsSpace.h"
#include "model/ObsVector/ObsVec.h"
#include "model/Variables/Variables.h"
#include "eckit/config/Configuration.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace mom5cice5 {
  // -----------------------------------------------------------------------------

  ObsFraction::ObsFraction(ObsSpace & odb, const eckit::Configuration & config)
    : obsdb_(odb), obsname_("Fraction"), varin_()
  {
    const eckit::Configuration * configc = &config;
    mom5cice5_fraction_setup_f90(keyOperStrm_, &configc);
    int keyVarin;
    mom5cice5_obsoper_inputs_f90(keyOperStrm_, keyVarin);
    varin_.reset(new Variables(keyVarin));
    Log::trace() << "ObsFraction created " << obsname_ << std::endl;
  }

  // -----------------------------------------------------------------------------

  ObsFraction::~ObsFraction() {
    mom5cice5_fraction_delete_f90(keyOperStrm_);
  }

  // -----------------------------------------------------------------------------

  void ObsFraction::obsEquiv(const Gom & gom, ObsVec & ovec,
			     const ObsBias & bias) const {
    mom5cice5_fraction_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.fraction());
  }

  // -----------------------------------------------------------------------------

  void ObsFraction::generateObsError(const eckit::Configuration & conf) {
    const double err = conf.getDouble("obs_error");
    mom5cice5_obsdb_seterr_f90(obsdb_.toFortran(), keyOperStrm_, err);
  }

  // -----------------------------------------------------------------------------

  void ObsFraction::print(std::ostream & os) const {
    os << "ObsFraction::print not implemented";
  }

  // -----------------------------------------------------------------------------

}  // namespace mom5cice5
