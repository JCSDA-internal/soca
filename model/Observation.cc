
#include "model/Observation.h"

#include <map>
#include <string>

#include "util/Logger.h"
#include "util/abor1_cpp.h"
#include "eckit/config/Configuration.h"


using oops::Log;

namespace mom5cice5 {

  // -----------------------------------------------------------------------------

  // We have to use a pointer here because we don't know the order in
  // which objects are created and this causes problems with xlc.

  std::map < std::string, ObsFactory * > * ObsFactory::makers_ = 0;

  // -----------------------------------------------------------------------------

  ObsFactory::ObsFactory(const std::string & name) {
    if (!makers_) makers_=new std::map < std::string, ObsFactory * >();

    if (makers_->find(name) != makers_->end()) {
      Log::error() << name << " already registered in MOM5CICE5 observation factory." << std::endl;
      ABORT("Element already registered in ObsFactory.");
    }
    (*makers_)[name] = this;
  }

  // -----------------------------------------------------------------------------

  Observation * ObsFactory::create(ObsSpace & odb, const eckit::Configuration & conf) {
    if (!makers_) makers_=new std::map < std::string, ObsFactory * >();

    std::string id = conf.getString("ObsType");
    Log::trace() << "MOM5CICE5 ObsFactory ObsType =" << id << std::endl;
    std::map<std::string, ObsFactory*>::iterator j = makers_->find(id);
    if (j == makers_->end()) {
      Log::error() << id << " does not exist in MOM5CICE5 observation factory." << std::endl;
      ABORT("Element does not exist in ObsFactory.");
    }
    return (*j).second->make(odb, conf);
  }

  // -----------------------------------------------------------------------------
}  // namespace mom5cice5
