
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
namespace soca {
  // -----------------------------------------------------------------------------

  static oops::ObsOperatorMaker<Traits, ObsFraction> makerFraction_("Fraction");
  
  ObsFraction::ObsFraction(const ObsSpace & odb, const eckit::Configuration & config)
    : obsname_("Fraction"), varin_()
  {
    Log::trace() << "ObsFraction start create " << obsname_ << std::endl;    
    const eckit::Configuration * configc = &config;
    soca_fraction_setup_f90(keyOperFraction_, &configc);
    int keyVarin;
    soca_obsoper_inputs_f90(keyOperFraction_, keyVarin);
    varin_.reset(new Variables(keyVarin));
    Log::trace() << "ObsFraction created " << obsname_ << std::endl;
  }

  // -----------------------------------------------------------------------------

  ObsFraction::~ObsFraction() {
    soca_fraction_delete_f90(keyOperFraction_);
  }

  // -----------------------------------------------------------------------------

  void ObsFraction::obsEquiv(const Gom & gom, ObsVec & ovec,
			     const ObsBias & bias) const {
    Log::trace() << "Starting ObsEquiv ... " << obsname_ << std::endl;        
    soca_fraction_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.fraction());
    Log::trace() << "Out of ObsEquiv ... " << obsname_ << std::endl;            
  }

  // -----------------------------------------------------------------------------

  //void ObsFraction::generateObsError(const eckit::Configuration & conf) {
  //  const double err = conf.getDouble("obs_error");
  //  soca_obsdb_seterr_f90(obsdb_.toFortran(), keyOperFraction_, err);
  //}

  // -----------------------------------------------------------------------------

  void ObsFraction::print(std::ostream & os) const {
    os << "ObsFraction::print not implemented";
  }

  // -----------------------------------------------------------------------------

}  // namespace soca
