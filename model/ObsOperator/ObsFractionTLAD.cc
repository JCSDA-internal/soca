
#include "model/ObsOperator/ObsFractionTLAD.h"

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "model/ModelAtLocations/Gom.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsSpace/ObsSpace.h"
#include "model/ObsVector/ObsVec.h"
#include "model/Variables/Variables.h"
#include "model/Fortran.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace mom5cice5 {
// -----------------------------------------------------------------------------
static oops::LinearObsOpMaker<Traits, ObsFractionTLAD> makerFractionTL_("Fraction");
// -----------------------------------------------------------------------------
//ObsFractionTLAD::ObsFractionTLAD(const ObsSpace &, const int & keyOperStrm)
//  : keyOperFraction_(keyOperStrm), varin_()
ObsFractionTLAD::ObsFractionTLAD(const ObsSpace &, const eckit::Configuration & config)
  : keyOperFraction_(0), varin_()    
{
  const eckit::Configuration * configc = &config;    
  mom5cice5_fraction_setup_f90(keyOperFraction_, &configc);
  int keyVarin; 
  mom5cice5_obsoper_inputs_f90(keyOperFraction_, keyVarin);
  varin_.reset(new Variables(keyVarin));
  Log::trace() << "ObsFractionTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsFractionTLAD::~ObsFractionTLAD() {
  Log::trace() << "ObsFractionTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsFractionTLAD::setTrajectory(const Gom &, const ObsBias &) {}

// -----------------------------------------------------------------------------

void ObsFractionTLAD::obsEquivTL(const Gom & gom, ObsVec & ovec,
                               const ObsBiasIncrement & bias) const {
  // NOT IMPLEMENTED YET
  mom5cice5_fraction_equiv_tl_f90(gom.toFortran(), ovec.toFortran(), bias.fraction());
}

// -----------------------------------------------------------------------------

void ObsFractionTLAD::obsEquivAD(Gom & gom, const ObsVec & ovec,
                               ObsBiasIncrement & bias) const {
  mom5cice5_fraction_equiv_ad_f90(gom.toFortran(), ovec.toFortran(), bias.fraction());
}

// -----------------------------------------------------------------------------
void ObsFractionTLAD::print(std::ostream & os) const {
  os << "ObsFractionTLAD::print not implemented" << std::endl;  
}
}  // namespace mom5cice5
