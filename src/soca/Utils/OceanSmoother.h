/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

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
/// @brief A class to smooth arbitrary ocean fields using the commonly used scales.
///
/// This class uses the oops::Diffusion class to do the actual smoothing. The
/// scales are calculated from our usual ocean related input parameters (e.g.
/// Rossby radius multiplier and minimum grid multiplier for the horizontal, and
/// mixed layer depth interpolated to a constant value for the vertical). Note
/// that all scales are defined as a Gaussian sigma.
class OceanSmoother {
 public:
  // -----------------------------------------------------------------------------
  // Yaml parameters for OceanSmoother
  class Parameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Parameters, oops::Parameters)

   public:
    // ---------------------------------------------------------------------------
    /// @brief Parameters for generating the 3D mask
    ///
    /// In short, the mask is a combination of the 2D landmask, and
    /// a 3D mask based on the layer thickness. The 3D mask is generated
    /// by masking out any layers with thickness less than a certain value.
    class Mask : public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Mask, oops::Parameters)
     public:
      oops::Parameter<std::string> maskVariable{"mask variable",
        "The name of the geometry variable to use as a 2D horizontal mask",
        "interp_mask", this};
      oops::Parameter<double> minThickness{"min thickness",
        "The minimum layer thickness. Anything smaller will be masked out",
         0.01, this, {oops::minConstraint(0.0)}};
      oops::Parameter<bool> setMaskedToZero{"set masked to zero",
        "If true, set masked out values to zero during multiply(). Otherwise, leave them alone",
        true, this};
    };

    // ---------------------------------------------------------------------------
    /// @brief Parameters for the optional horizontal scales (meters)
    ///
    /// The horizontal scales are calculated as a Gaussian sigma based on the
    /// Rossby radius and the minimum grid size. The scales are calculated as
    /// follows:
    ///   - `base value` + `rossby Mult` * rossbyRadius
    ///   - `min grid mult` * gridSize
    ///
    /// The final scale is the maximum of the two, with extra clamping between
    /// the `min` and `max` values.
    class Horizontal : public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Horizontal, oops::Parameters)
     public:
      oops::Parameter<std::string> rossbyVariable{"rossby radius variable",
        "The name of the geometry variable to use as the Rossby radius",
        "rossby_radius", this};
      oops::Parameter<double> base{"base value", 0.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> rossbyMult{"rossby mult", 1.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> minGridMult{"min grid mult", 1.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> min{"min", 0.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> max{"max", 150e3, this, {oops::minConstraint(0.0)}};
    };

    // ---------------------------------------------------------------------------
    /// @brief Parameters for the optional vertical scales (number of levels)
    ///
    /// The vertical scales can be set at the surface by the mixed layer depth,
    /// interpolating down to the `min` value at the bottom of the mixed layer.
    /// The MLD can optionally be smoothed by the horizontal scales if desired.
    class Vertical : public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Vertical, oops::Parameters)
     public:
      oops::Parameter<double> min{"min", 1.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> max{"max", 10.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<bool> useMld{"use mld",
        "If true, use the mixed layer depth to set the vertical scales at the top of the ocean",
       false, this};
      oops::Parameter<std::string> mldVariable{"mld variable",
        "The name of the MLD variable in the fieldset passed to the constructor",
        "ocean_mixed_layer_thickness", this};
      oops::Parameter<bool> mldSmooth{"mld smooth",
        "If true, smooth the MLD by the horizontal scales before using it for the vertical scales",
        true, this};
      oops::Parameter<double> mldMax{"mld max",
        "The maximum MLD value to use for the vertical scales",
        250.0, this, {oops::minConstraint(0.0)}};
    };

    // ---------------------------------------------------------------------------

    oops::Parameter<Mask> mask{"mask", Mask(), this};
    oops::OptionalParameter<Horizontal> horizontal{"horizontal", this};
    oops::OptionalParameter<Vertical> vertical{"vertical", this};
    oops::OptionalParameter<eckit::LocalConfiguration> saveDiags{"save diagnostics",
      "If present, save the calculated scales and the mld parameters to the given file",
      this};
    oops::Parameter<std::string> thicknessVariable{"thickness variable",
      "The name of the layer thickness state variable in the fieldset passed to the constructor",
      "sea_water_cell_thickness", this};
  };

  // -----------------------------------------------------------------------------
  /// @brief Constructor for the ocean smoother
  /// @param geom geometry data for the model
  /// @param params parameters for the ocean smoother
  /// @param levels number of vertical levels for the fields IF vertical smoothing is requested
  /// @param bkg background fieldset that contains layer thickness and MLD fields (if needed)
  explicit OceanSmoother(
    const oops::GeometryData & geom,
    const Parameters & params,
    const int levels,
    const atlas::FieldSet & bkg = atlas::FieldSet());

  void multiply(atlas::FieldSet &);
  void multiply(atlas::Field &);

 private:
  const oops::GeometryData & geom_;
  std::shared_ptr<oops::Diffusion> diffusion_;
  std::shared_ptr<atlas::Field> mask_;
};

// -----------------------------------------------------------------------------

}  // namespace soca
