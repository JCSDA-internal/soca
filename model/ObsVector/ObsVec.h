
#ifndef SOCA_MODEL_OBSVEC_H_
#define SOCA_MODEL_OBSVEC_H_

#include <ostream>
#include <string>

#include "model/Fortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace soca {
  class ObsSpace;

  // -----------------------------------------------------------------------------
  /// ObsVec class to handle vectors in observation space for SOCA model.

  class ObsVec : public util::Printable,
    private util::ObjectCounter<ObsVec> {
  public:
      static const std::string classname() {return "soca::ObsVec";}

      explicit ObsVec(const ObsSpace &);
      ObsVec(const ObsVec &, const bool copy = true);
      ~ObsVec();

      ObsVec & operator = (const ObsVec &);
      ObsVec & operator*= (const double &);
      ObsVec & operator+= (const ObsVec &);
      ObsVec & operator-= (const ObsVec &);
      ObsVec & operator*= (const ObsVec &);
      ObsVec & operator/= (const ObsVec &);

      void zero();
      void axpy(const double &, const ObsVec &);
      void invert();
      void random();
      double dot_product_with(const ObsVec &) const;
      double rms() const;

      unsigned int size() const;

      int & toFortran() {return keyOvec_;}
      const int & toFortran() const {return keyOvec_;}

      // I/O
      void read(const std::string &);
      void save(const std::string &) const;

  private:
      void print(std::ostream &) const;

      const ObsSpace & obsdb_;
      int keyOvec_;
    };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_OBSVEC_H_
