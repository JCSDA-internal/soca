/*
 * (C) Copyright 2024-2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/field.h"

#include "soca/Fields/Fields.h"

namespace soca
{

// -----------------------------------------------------------------------------

Fields::Fields(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & vt)
 : geom_(geom), vars_(vars), time_(vt)  {
}

// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t Fields::serialSize() const {
  size_t nn = 1;  // plus magic factor
  atlas::FieldSet fs; toFieldSet(fs);
  for(const auto & field : fs) {
      nn += field.size();
  }

  // Date and time
  nn += time_.serialSize();
  return nn;
}

// -----------------------------------------------------------------------------

constexpr double SerializeCheckValue = -54321.98765;
void Fields::serialize(std::vector<double> & vect) const {
  atlas::FieldSet fs; toFieldSet(fs);

  // Serialize the field
  size_t n = 0;
  vect.reserve(serialSize());
  for(const auto & field : fs) {
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t i=0; i < field.shape(0); i++) {
      for (size_t j = 0; j < field.shape(1); j++) {
        vect.push_back(view(i,j));
        n++;
      }
    }
  }

  // Magic value placed in serialization; used to validate deserialization
  vect.push_back(SerializeCheckValue);

  // Serialize the date and time
  time_.serialize(vect);
}

// -----------------------------------------------------------------------------

void Fields::deserialize(const std::vector<double> & vect, size_t & index) {
  // Deserialize the field
  atlas::FieldSet fs; toFieldSet(fs);
  for (auto & field : fs) {
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t i=0; i < view.shape(0); i++) {
      for (size_t j = 0; j < view.shape(1); j++) {
        view(i,j) = vect[index++];
      }
    }
  }
  // Use magic value to validate deserialization
  ASSERT(vect.at(index) == SerializeCheckValue);
  ++index;

  // Deserialize the date and time
  time_.deserialize(vect, index);

  fromFieldSet(fs);  // TODO temp
}

}