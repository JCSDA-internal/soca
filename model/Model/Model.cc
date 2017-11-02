
#include "model/Model/Model.h"

#include "util/Logger.h"
#include "model/ModelBias.h"
#include "model/Fields/Fields.h"
#include "model/Fortran.h"
#include "model/Geometry/Geometry.h"
#include "model/State/State.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"

using oops::Log;

namespace mom5cice5 {
  // -----------------------------------------------------------------------------
  Model::Model(const Geometry & resol, const eckit::Configuration & model)
    : keyConfig_(0), tstep_(0), geom_(resol)
  {
    Log::trace() << "Model::Model" << std::endl;
    tstep_ = util::Duration(model.getString("tstep"));
    const eckit::Configuration * configc = &model;
    mom5cice5_setup_f90(&configc, geom_.toFortran(), keyConfig_);
    Log::trace() << "Model created" << std::endl;
  }
  // -----------------------------------------------------------------------------
  Model::~Model() {
    mom5cice5_delete_f90(keyConfig_);
    Log::trace() << "Model destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::initialize(State & xx) const {
    xx.activateModel();
    ASSERT(xx.fields().isForModel(true));
    mom5cice5_prepare_integration_f90(keyConfig_, xx.fields().toFortran());
    Log::debug() << "Model::initialize" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::step(State & xx, const ModelBias &) const {
    ASSERT(xx.fields().isForModel(true));
    Log::debug() << "Model::step fields in" << xx.fields() << std::endl;
    mom5cice5_propagate_f90(keyConfig_, xx.fields().toFortran());
    xx.validTime() += tstep_;
    Log::debug() << "Model::step fields out" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::finalize(State & xx) const {
    xx.deactivateModel();
    Log::debug() << "Model::finalize" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  int Model::saveTrajectory(State & xx, const ModelBias &) const {
    ASSERT(xx.fields().isForModel(true));
    int ftraj = 0;
    Log::debug() << "Model::saveTrajectory fields in" << xx.fields() << std::endl;
    //mom5cice5_prop_traj_f90(keyConfig_, xx.fields().toFortran(), ftraj);
    xx.validTime() += tstep_;
    ASSERT(ftraj != 0);
    Log::debug() << "Model::saveTrajectory fields out" << xx.fields() << std::endl;
    return ftraj;
  }
  // -----------------------------------------------------------------------------
  void Model::print(std::ostream & os) const {
    os << "Model::print not implemented";
  }
  // -----------------------------------------------------------------------------
}  // namespace mom5cice5
