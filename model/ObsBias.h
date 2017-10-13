
#ifndef MOM5CICE5_MODEL_OBSBIAS_H_
#define MOM5CICE5_MODEL_OBSBIAS_H_

#include <iostream>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mom5cice5 {
  class ObsBiasIncrement;

/// Class to handle observation bias parameters.

// -----------------------------------------------------------------------------

class ObsBias : public util::Printable,
                private boost::noncopyable,
                private util::ObjectCounter<ObsBias> {
 public:
  static const unsigned int ntypes = 4;
  static const std::string classname() {return "mom5cice5::ObsBias";}

  explicit ObsBias(const eckit::Configuration &);
  ObsBias(const ObsBias &, const bool);
  ~ObsBias() {}

  ObsBias & operator+=(const ObsBiasIncrement &);

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const;

  const double & operator[](const unsigned int ii) const {return bias_[ii];}

  const double & fraction() const {return bias_[0];}
  const double & temp() const {return bias_[1];}
  const double & salt() const {return bias_[3];}

 private:
  void print(std::ostream &) const;
  std::vector<double> bias_;
  bool active_;
};

// -----------------------------------------------------------------------------

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_OBSBIAS_H_
