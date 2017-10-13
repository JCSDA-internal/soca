
#include "model/ObsFractionMOM5CICE5.h"

#include "util/Logger.h"
#include "model/Gom.h"
//#include "model/ObsBias.h"
#include "model/ObsSpace.h"
#include "model/ObsVec.h"
#include "model/Variables.h"
#include "eckit/config/Configuration.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace mom5cice5 {
// -----------------------------------------------------------------------------

ObsFractionMOM5CICE5::ObsFractionMOM5CICE5(ObsSpace & odb, const eckit::Configuration & config)
  : obsdb_(odb), obsname_("Fraction"), varin_()
{
  const eckit::Configuration * configc = &config;
  //mom5cice5_stream_setup_f90(keyOperStrm_, &configc);
  int keyVarin;
  //mom5cice5_obsoper_inputs_f90(keyOperStrm_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  Log::trace() << "ObsFractionMOM5CICE5 created " << obsname_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsFractionMOM5CICE5::~ObsFractionMOM5CICE5() {
  //mom5cice5_stream_delete_f90(keyOperStrm_);
}

// -----------------------------------------------------------------------------

void ObsFractionMOM5CICE5::obsEquiv(const Gom & gom, ObsVec & ovec,
                           const ObsBias & bias) const {
  //mom5cice5_stream_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsFractionMOM5CICE5::generateObsError(const eckit::Configuration & conf) {
  const double err = conf.getDouble("obs_error");
  mom5cice5_obsdb_seterr_f90(obsdb_.toFortran(), keyOperStrm_, err);
}

// -----------------------------------------------------------------------------

void ObsFractionMOM5CICE5::print(std::ostream & os) const {
  os << "ObsFractionMOM5CICE5::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace mom5cice5
