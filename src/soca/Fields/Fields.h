/*
 * (C) Copyright 2024-2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

 */

#pragma once

#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Serializable.h"

namespace soca {

// --------------------------------------------------------------------------------------

class Geometry;

// --------------------------------------------------------------------------------------

/// @brief The `Fields` class contains the components that are common to both
/// the `State` and `Increment` classes.
class Fields : public util::Serializable,
               public util::Printable {
 public:
  explicit Fields(const Geometry &, const oops::Variables &, const util::DateTime &);

  // These create copies of FieldSets
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);
  // TODO(travis) remove this direct accessor when I figure out why sharing in
  // toFieldSet causes a change of answers.
  const atlas::FieldSet & fieldSet() const {return fieldSet_;}
  atlas::FieldSet & fieldSet() {return fieldSet_;}

  // math operators
  void zero();
  double norm() const;
  void updateTime(const util::Duration & dt) {time_ += dt;}
  void accumul(const double &, const Fields &);

  // serialization
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

  // accessors
  const Geometry & geometry() const {return geom_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  const oops::Variables & variables() const {return vars_;}

 protected:
  void print(std::ostream &) const override;

  atlas::FieldSet fieldSet_;
  util::DateTime time_;
  oops::Variables vars_;
  const Geometry & geom_;
};

// --------------------------------------------------------------------------------------

}  // namespace soca
