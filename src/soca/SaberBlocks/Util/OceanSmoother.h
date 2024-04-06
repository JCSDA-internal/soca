/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>

#include "soca/SaberBlocks/Util/OceanSmootherParameters.h"
#include "soca/Utils/Diffusion.h"

// forward declarations
namespace atlas {
  class Field;
  class FieldSet;
}

namespace oops {
  class GeometryData;
}

// --------------------------------------------------------------------------------------

namespace soca
{


class OceanSmoother {

 public:
  typedef OceanSmootherParameters Parameters_;

  explicit OceanSmoother(const oops::GeometryData &, const Parameters_ &);
  void multiply(atlas::FieldSet &) const;
  void multiply(atlas::Field &) const;

 private:
  const oops::GeometryData & geom_;
  Diffusion diffusion_;
};

}