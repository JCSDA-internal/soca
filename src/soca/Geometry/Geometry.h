/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GEOMETRY_GEOMETRY_H_
#define SOCA_GEOMETRY_GEOMETRY_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "soca/Fortran.h"
#include "soca/Geometry/GeometryFortran.h"
#include "soca/GeometryIterator/GeometryIterator.h"


namespace eckit {
  class Configuration;
}

namespace soca {

  class GeometryIterator;

  // -----------------------------------------------------------------------------
  /// Geometry handles geometry for SOCA model.

  class Geometry : public util::Printable,
    private util::ObjectCounter<Geometry> {
   public:
      static const std::string classname() {return "soca::Geometry";}

      explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
      Geometry(const Geometry &);
      ~Geometry();

      GeometryIterator begin() const;
      GeometryIterator end() const;


      int& toFortran() {return keyGeom_;}
      const int& toFortran() const {return keyGeom_;}
      void gridgen(const eckit::Configuration &) const;
      const eckit::mpi::Comm & getComm() const {return comm_;}

   private:
      Geometry & operator=(const Geometry &);
      void print(std::ostream &) const;
      int keyGeom_;
      const eckit::mpi::Comm & comm_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_GEOMETRY_GEOMETRY_H_
