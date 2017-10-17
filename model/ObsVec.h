
#ifndef MOM5CICE5_MODEL_OBSVEC_H_
#define MOM5CICE5_MODEL_OBSVEC_H_

#include <ostream>
#include <string>

#include "model/Fortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace mom5cice5 {
  class ObsSpace;

  // -----------------------------------------------------------------------------
  /// ObsVec class to handle vectors in observation space for MOM5CICE5 model.

  class ObsVec : public util::Printable,
    private util::ObjectCounter<ObsVec> {
  public:
      static const std::string classname() {return "mom5cice5::ObsVec";}

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

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_OBSVEC_H_
