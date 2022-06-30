/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GEOMETRY_GEOMETRY_H_
#define SOCA_GEOMETRY_GEOMETRY_H_

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/KDTree.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "soca/Fortran.h"
#include "soca/Geometry/FmsInput.h"
#include "soca/Geometry/GeometryFortran.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/GeometryIterator/GeometryIteratorFortran.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace atlas {
  class FieldSet;
  class FunctionSpace;
  namespace functionspace {
    class PointCloud;
  }
}
namespace oops {
  class Variables;
}

// -----------------------------------------------------------------------------

namespace soca {

  /// Geometry handles geometry for SOCA model.
  class Geometry : public util::Printable,
    private util::ObjectCounter<Geometry> {
   public:
      static const std::string classname() {return "soca::Geometry";}

      explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &, const bool = false);
      Geometry(const Geometry &);
      ~Geometry();

      GeometryIterator begin() const;
      GeometryIterator end() const;
      int IteratorDimension() const;
      std::vector<size_t> variableSizes(const oops::Variables & vars) const;
      std::vector<double> verticalCoord(std::string &) const {return {};}

      int& toFortran() {return keyGeom_;}
      const int& toFortran() const {return keyGeom_;}
      const eckit::mpi::Comm & getComm() const {return comm_;}

      const atlas::FunctionSpace & functionSpace() const {return functionSpaceIncHalo_;}
      atlas::FunctionSpace & functionSpace() {return functionSpaceIncHalo_;}
      const atlas::FieldSet & extraFields() const {return extraFields_;}
      atlas::FieldSet & extraFields() {return extraFields_;}

      void latlon(std::vector<double> &, std::vector<double> &, const bool) const;
      void latlon(std::vector<double> &, std::vector<double> &, const bool,
                  const char, const bool) const;
      atlas::util::KDTree<size_t>::ValueList closestPoints(
        const double, const double, const int,
        const char, const bool) const;

      void getVarGrid(const std::string &, char &, bool &) const;

   private:
      void mask(std::vector<double> &, const bool, const char) const;
      Geometry & operator=(const Geometry &);
      void print(std::ostream &) const;

      int keyGeom_;
      const eckit::mpi::Comm & comm_;
      FmsInput fmsinput_;
      atlas::FunctionSpace functionSpace_;
      atlas::FunctionSpace functionSpaceIncHalo_;
      atlas::FieldSet extraFields_;

      atlas::util::IndexKDTree localTree_[6];
  };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_GEOMETRY_GEOMETRY_H_
