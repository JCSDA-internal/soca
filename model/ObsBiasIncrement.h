
#ifndef MOM5CICE5_MODEL_OBSBIASINCREMENT_H_
#define MOM5CICE5_MODEL_OBSBIASINCREMENT_H_

#include <iostream>
#include <vector>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mom5cice5 {
  class ObsBias;

  // -----------------------------------------------------------------------------

  class ObsBiasIncrement : public util::Printable {
  public:
    /// Constructor, destructor
    explicit ObsBiasIncrement();
    explicit ObsBiasIncrement(const eckit::Configuration &);
    ObsBiasIncrement(const ObsBiasIncrement &, const bool copy = true);
    ObsBiasIncrement(const ObsBiasIncrement &, const eckit::Configuration &);
    ~ObsBiasIncrement() {}

    /// Linear algebra operators
    void diff(const ObsBias &, const ObsBias &);
    void zero();
    ObsBiasIncrement & operator=(const ObsBiasIncrement &);
    ObsBiasIncrement & operator+=(const ObsBiasIncrement &);
    ObsBiasIncrement & operator-=(const ObsBiasIncrement &);
    ObsBiasIncrement & operator*=(const double);
    void axpy(const double, const ObsBiasIncrement &);
    double dot_product_with(const ObsBiasIncrement &) const;

    /// I/O and diagnostics
    void read(const eckit::Configuration &) {}
    void write(const eckit::Configuration &) const {}
    double norm() const;

    double & operator[](const unsigned int ii) {return bias_[ii];}
    const double & operator[](const unsigned int ii) const {return bias_[ii];}

    double & stream() {return bias_[0];}
    double & temp() {return bias_[1];}
    double & salt() {return bias_[3];}
    const double & fraction() const {return bias_[0];}
    const double & wind() const {return bias_[1];}
    const double & wspd() const {return bias_[3];}

  private:
    void print(std::ostream &) const;
    void makePassive();

    std::vector<double> bias_;
    std::vector<bool> active_;
  };

  // -----------------------------------------------------------------------------

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_OBSBIASINCREMENT_H_
