/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

// Forward declarations
namespace atlas {
  class Field;
  class FieldSet;
}

namespace oops {
  class GeometryData;
  class Diffusion;
}

// -----------------------------------------------------------------------------

namespace soca {

// -----------------------------------------------------------------------------
/// @brief A class to smooth arbitrary ocean fields.
///
/// This class uses the oops::Diffusion class to do the actual smoothing. The
/// scales are calculated from our usual ocean related input parameters (e.g.
/// Rossby radius multiplier and minimum grid multiplier). Note that all scales
/// are defined as a Gaussian sigma.
class OceanSmoother {
 public:

  // -----------------------------------------------------------------------------
  // Yaml parameters for OceanSmoother
  class Parameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Parameters, oops::Parameters)
   public:

    // ---------------------------------------------------------------------------

    class Mask : public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Mask, oops::Parameters)
     public:
      oops::Parameter<std::string> maskVariable{"mask variable",
        "The name of the geometry variable to use as a 2D horizontal mask",
        "interp_mask", this};
      oops::Parameter<double> minThickness{"min thickness",
        "The minimum layer thickness. Anything smaller will be masked out",
         0.01, this, {oops::minConstraint(0.0)}};
      oops::Parameter<std::string> thicknessVariable{"thickness variable",
        "The state variable to use for the layer thickness",
        "hocn", this};
      oops::Parameter<bool> setMaskedToZero{"set masked to zero",
        "If true, set masked out values to zero. Otherwise, leave them alone",
        true, this};
    };

    // ---------------------------------------------------------------------------
    class Horizontal : public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Horizontal, oops::Parameters)
     public:
      oops::Parameter<std::string> rossbyVariable{"rossby radius variable", "rossby_radius", this};

      oops::Parameter<double> base{"base value", 0.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> rossbyMult{"rossby mult", 1.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> minGridMult{"min grid mult", 1.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> min{"min", 0.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> max{"max", std::numeric_limits<double>::max(), this,
                                  {oops::minConstraint(0.0)}};
    };

    // ---------------------------------------------------------------------------

    class Vertical : public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Vertical, oops::Parameters)
     public:
      oops::Parameter<double> base{"base value", 0.0, this, {oops::minConstraint(0.0)}};
    };

    // ---------------------------------------------------------------------------

    oops::Parameter<Mask> mask{"mask", Mask(), this};
    oops::OptionalParameter<Horizontal> horizontal{"horizontal", this};
    oops::OptionalParameter<Vertical> vertical{"vertical", this};
  };

  // -----------------------------------------------------------------------------
  explicit OceanSmoother(
    const oops::GeometryData & geom,
    const Parameters & params,
    const int levels,
    const atlas::FieldSet & bkg = atlas::FieldSet()
    );

  void multiply(atlas::FieldSet &);
  void multiply(atlas::Field &);

 private:
  const oops::GeometryData & geom_;
  std::shared_ptr<oops::Diffusion> diffusion_;
  std::shared_ptr<atlas::Field> mask_;
};

// -----------------------------------------------------------------------------

}  // namespace soca