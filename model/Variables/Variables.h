
#ifndef SOCA_MODEL_SOCAVARIABLES_H_
#define SOCA_MODEL_SOCAVARIABLES_H_

#include <ostream>
#include <string>

#include "util/Logger.h"
#include "model/Fortran.h"
#include "eckit/config/Configuration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace soca {

  // -----------------------------------------------------------------------------
  /// Variables class to handle variables for SOCA model.

  class Variables : public util::Printable,
    private util::ObjectCounter<Variables> {
  public:
      static const std::string classname() {return "soca::Variables";}

      explicit Variables(const eckit::Configuration & config) {
	using oops::Log;
	Log::debug() << "Variables config:" << config << std::endl;
	const eckit::Configuration * conf = &config;
	soca_var_create_f90(keyVar_, &conf);
      }
      explicit Variables(const int keyVar): keyVar_(keyVar) {}

      ~Variables() {soca_var_delete_f90(keyVar_);}

      Variables(const Variables & other) {soca_var_clone_f90(other.keyVar_, keyVar_);}

      int& toFortran() {return keyVar_;}
      const int& toFortran() const {return keyVar_;}

  private:
      void print(std::ostream &) const;
      int keyVar_;
    };

  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_SOCAVARIABLES_H_
