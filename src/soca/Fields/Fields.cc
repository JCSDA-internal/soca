/*
 * (C) Copyright 2024-2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <algorithm>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <string>

#include "atlas/field.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Fields/Fields.h"

namespace soca
{

// Note, until we deal with U/V de-staggering, the print value is kinda
// slightly wrong. Oh well.
const char MASK_METADATA[] = "mask";

// -----------------------------------------------------------------------------

Fields::Fields(const Geometry & geom, const oops::Variables & vars, const util::DateTime & vt)
  : geom_(geom), vars_(vars), time_(vt)  {
}


// -----------------------------------------------------------------------------

void Fields::zero() {
  util::zeroFieldSet(fieldSet_);
}

// -----------------------------------------------------------------------------

void Fields::accumul(const double & zz, const Fields & xx) {
  atlas::FieldSet fs1, fs2; xx.toFieldSet(fs1);
  util::copyFieldSet(fs1, fs2);
  util::multiplyFieldSet(fs2, zz);
  util::addFieldSets(fieldSet_, fs2);
}

// -----------------------------------------------------------------------------

size_t Fields::serialSize() const {
  size_t nn = 1;  // plus magic factor
  for (const auto & field : fieldSet_) {
    nn += field.size();
  }

  // Date and time
  nn += time_.serialSize();
  return nn;
}

// -----------------------------------------------------------------------------

constexpr double SerializeCheckValue = -54321.98765;
void Fields::serialize(std::vector<double> & vect) const {
  size_t n = 0;
  vect.reserve(serialSize());
  for (const auto & field : fieldSet_) {
    const auto & view = atlas::array::make_view<double, 2>(field);
    for (size_t i=0; i < field.shape(0); i++) {
      for (size_t j = 0; j < field.shape(1); j++) {
        vect.push_back(view(i, j));
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
  for (auto & field : fieldSet_) {
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t i=0; i < view.shape(0); i++) {
      for (size_t j = 0; j < view.shape(1); j++) {
        view(i, j) = vect[index++];
      }
    }
  }
  // Use magic value to validate deserialization
  ASSERT(vect.at(index++) == SerializeCheckValue);

  // Deserialize the date and time
  time_.deserialize(vect, index);
}

// -----------------------------------------------------------------------------

double Fields::norm() const {
  double zz = 0.0;
  for (const auto & field : fieldSet_) {
    const auto & vGhost = atlas::array::make_view<int, 1>(field.functionspace().ghost());
    const auto & view = atlas::array::make_view<double, 2>(field);
    std::unique_ptr<atlas::array::ArrayView<double, 2> > mask;
    if (field.metadata().has(MASK_METADATA)) {
      // optionally get the mask field, if one is given
      const std::string & maskName = field.metadata().getString(MASK_METADATA);
      mask.reset(new atlas::array::ArrayView<double, 2>(
        atlas::array::make_view<double, 2>(geom_.fields().field(maskName))));
    }

    // traverse all points to calculate stats
    for (size_t i = 0; i < field.shape(0); i++) {
      if (vGhost(i)) continue;
      if (mask && (*mask)(i, 0) == 0.0) continue;
      for (size_t lvl = 0; lvl < field.shape(1); lvl++) {
        zz += view(i, lvl) * view(i, lvl);
      }
    }
  }

  geom_.getComm().allReduceInPlace(zz, eckit::mpi::sum());
  zz = sqrt(zz);
  return zz;
}

// -----------------------------------------------------------------------------

void Fields::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();

  // find the longest field name, for use in formatting of the print
  size_t maxNameLen = 0;
  for (const auto & field : fieldSet_) {
    maxNameLen = std::max(maxNameLen, field.name().size());
  }

  // for each field
  for (const auto & field : fieldSet_) {
    size_t count = 0;
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double sum = 0.0;

    const auto & vGhost = atlas::array::make_view<int, 1>(field.functionspace().ghost());
    const auto & view = atlas::array::make_view<double, 2>(field);
    std::unique_ptr<atlas::array::ArrayView<double, 2> > mask;
    if (field.metadata().has(MASK_METADATA)) {
      // optionally get the mask field, if one is given
      const std::string & maskName = field.metadata().getString(MASK_METADATA);
      mask.reset(new atlas::array::ArrayView<double, 2>(
        atlas::array::make_view<double, 2>(geom_.fields().field(maskName))));
    }

    // traverse all points to calculate stats
    for (size_t i = 0; i < field.shape(0); i++) {
      if (vGhost(i)) continue;
      if (mask && (*mask)(i, 0) == 0.0) continue;

      count++;
      for (size_t lvl = 0; lvl < field.shape(1); lvl++) {
        min = std::min(min, view(i, lvl));
        max = std::max(max, view(i, lvl));
        sum += view(i, lvl) / field.shape(1);
      }
    }

    // combine stats across PEs
    geom_.getComm().allReduceInPlace(count, eckit::mpi::sum());
    geom_.getComm().allReduceInPlace(min, eckit::mpi::min());
    geom_.getComm().allReduceInPlace(max, eckit::mpi::max());
    geom_.getComm().allReduceInPlace(sum, eckit::mpi::sum());
    if (count > 0) sum /= count;

    // done with this field, print information
    os << std::endl << std::right << std::setw(maxNameLen) << field.name()
        << "   min="  <<  std::fixed << std::setw(12) <<
                          std::right << min
        << "   max="  <<  std::fixed << std::setw(12) <<
                          std::right << max
        << "   mean=" <<  std::fixed << std::setw(12) <<
                          std::right << sum;
  }
}

// -----------------------------------------------------------------------------

void Fields::toFieldSet(atlas::FieldSet & fset) const {
  // copy seems like a waste, it would be nice if we could get away with "share"
  util::copyFieldSet(fieldSet_, fset);
}

// -----------------------------------------------------------------------------

void Fields::fromFieldSet(const atlas::FieldSet &fset) {
  // keep a copy of the metadata to copy back
  std::map<std::string, atlas::util::Metadata> metadata;
  for (const auto & f : fieldSet_) metadata[f.name()] = f.metadata();

  util::copyFieldSet(fset, fieldSet_);

  for (auto & f : fieldSet_) {
    if (metadata[f.name()].has("name")) {  // that's weird
      // THERE IS A BUG IN SOCA, i'm getting fieldsets with one or two fields
      // (hocn) that have have no metadata even though the rest of the field
      // do. Eh, i'll figure it out at a later date.
      f.metadata() = metadata[f.name()];
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace soca
