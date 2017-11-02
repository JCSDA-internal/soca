#ifndef MOM5CICE5_MODEL_OBSSPACE_H_
#define MOM5CICE5_MODEL_OBSSPACE_H_

#include <map>
#include <ostream>
#include <string>

#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/Printable.h"

#include "model/ObsSpace/ObsHelp.h"
#include "model/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace mom5cice5 {
  class ObsVec;
  class Observation;

  /// Wrapper around ObsHelp, mostly to hide the factory

  class ObsSpace : public util::Printable {
  public:
    ObsSpace(const eckit::Configuration &, const util::DateTime &, const util::DateTime &);
    ObsSpace(const ObsSpace &);
    ~ObsSpace();

    void getdb(const std::string & col, int & keyData) const {
      helper_->getdb(obsname_, col, keyData);
    }
    void putdb(const std::string & col, const int & keyData) const {
      helper_->putdb(obsname_, col, keyData);
    }

    int locations(const util::DateTime & t1, const util::DateTime & t2) const {
      int key_locs;
      key_locs = helper_->locations(obsname_, t1, t2);
      return key_locs;
    }

    void generateDistribution(const eckit::Configuration & conf) {
      helper_->generateDistribution(conf, obsname_, winbgn_, winend_, nobs_);
    }
 
    void printJo(const ObsVec &, const ObsVec &);

    int nobs() const {return nobs_;}
    int nvin() const {return nvin_;}
    int nout() const {return nout_;}
    const std::string & obsname() const {return obsname_;}
    const util::DateTime & windowStart() const {return winbgn_;}
    const util::DateTime & windowEnd() const {return winend_;}

    int & toFortran() {return helper_->toFortran();}
    const int & toFortran() const {return helper_->toFortran();}

  private:
    void print(std::ostream &) const;
    ObsSpace & operator= (const ObsSpace &);
    std::string ref_;
    mutable ObsHelp * helper_;
    std::string obsname_;
    unsigned int nobs_;
    unsigned int nvin_;
    unsigned int nout_;
    const util::DateTime winbgn_;
    const util::DateTime winend_;

    static std::map < std::string, int > theObsFileCount_;
  };

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_OBSSPACE_H_
