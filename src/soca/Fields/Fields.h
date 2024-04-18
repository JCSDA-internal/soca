/*
 * (C) Copyright 2024-2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Serializable.h"

namespace soca {

// --------------------------------------------------------------------------------------

class Geometry;

// --------------------------------------------------------------------------------------

class Fields : public util::Serializable {

 public:
  explicit Fields(const Geometry &, const oops::Variables &,
                  const util::DateTime &);

  virtual void toFieldSet(atlas::FieldSet &) const = 0;
  virtual void fromFieldSet(const atlas::FieldSet &) = 0;

  // serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  // accessors
  const Geometry & geometry() const {return geom_;}
  const util::DateTime & time() const {return time_;}

 protected:
  util::DateTime time_;
  oops::Variables vars_;
  const Geometry & geom_;
};

// --------------------------------------------------------------------------------------

}