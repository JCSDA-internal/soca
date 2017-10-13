
#include "model/ObsFractionTLAD.h"

#include "util/Logger.h"
#include "model/Gom.h"
//#include "model/ObsBias.h"
//#include "model/ObsBiasIncrement.h"
#include "model/ObsSpace.h"
#include "model/ObsVec.h"
#include "model/Variables.h"


using oops::Log;

// -----------------------------------------------------------------------------
namespace mom5cice5 {
// -----------------------------------------------------------------------------

ObsFractionTLAD::ObsFractionTLAD(const ObsSpace &, const int & keyOperStrm)
  : keyOperStrm_(keyOperStrm), varin_()
{
  int keyVarin;
  //mom5cice5_obsoper_inputs_f90(keyOperStrm_, keyVarin);
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
  //mom5cice5_stream_equiv_tl_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsFractionTLAD::obsEquivAD(Gom & gom, const ObsVec & ovec,
                               ObsBiasIncrement & bias) const {
  //mom5cice5_stream_equiv_ad_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

}  // namespace mom5cice5
