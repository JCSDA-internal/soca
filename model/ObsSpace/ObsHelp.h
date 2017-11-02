
#ifndef MOM5CICE5_MODEL_OBSHELP_H_
#define MOM5CICE5_MODEL_OBSHELP_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "model/Fortran.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mom5cice5 {

/// Observation Data Handler for MOM5CICE5 Model

class ObsHelp : private boost::noncopyable,
                  private util::ObjectCounter<ObsHelp> {
 public:
  static const std::string classname() {return "mom5cice5::ObsHelp";}

  explicit ObsHelp(const eckit::Configuration &);
  ~ObsHelp();

  void getdb(const std::string &, const std::string &, int & keyOvec) const;
  void putdb(const std::string &, const std::string &, const int & keyOvec);

  int locations(const std::string &, const util::DateTime &, const util::DateTime &) const;
  void generateDistribution(const eckit::Configuration &, const std::string &,
                            const util::DateTime &, const util::DateTime &, unsigned int &);
  int nobs(const std::string &) const;

  int & toFortran() {return keyHelp_;}
  const int & toFortran() const {return keyHelp_;}

 private:
  int keyHelp_;
};

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_OBSHELP_H_
