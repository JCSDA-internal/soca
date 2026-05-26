/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <string>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/Increment/Increment.h"
#include "soca/ObsLocalization/ObsLocRossby.h"
#include "soca/ObsLocalization/ObsLocRossbyParameters.h"
#include "soca/Traits.h"

namespace soca {

// -----------------------------------------------------------------------------
/// \brief Application that saves the Rossby-radius based observation
/// localization length scales at every grid point to a file.
///
/// This evaluates the same length-scale formula used by soca::ObsLocRossby
/// (base + rossby mult * rossby_radius, bounded by min/max and a minimum
/// grid-relative scale) for a given geometry and a "rossby localization"
/// configuration block, and writes the result as a single-field increment.
class SaveObsLocScales : public oops::Application {
 public:
  explicit SaveObsLocScales(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
  static const std::string classname() {return "soca::SaveObsLocScales";}

  int execute(const eckit::Configuration & fullConfig) const {
    // Setup geometry
    const eckit::LocalConfiguration geomConfig(fullConfig, "geometry");
    const Geometry geom(geomConfig, this->getComm());

    // Rossby-based localization specification
    ObsLocRossbyParameters options;
    options.deserialize(eckit::LocalConfiguration(fullConfig, "rossby localization"));

    // Whether to apply the gaussian -> Gaspari-Cohn width conversion, so the
    // saved value matches the width actually passed to the GC99 localization.
    // Defaults to true, if set to false, saves the raw gaussian length scale.
    const bool gcWidth = fullConfig.getBool("gaspari-cohn width", true);

    // Create a single-field 2D increment to hold the length scales
    const util::DateTime date("2000-01-01T00:00:00Z");  // dummy date
    const oops::Variables vars({"localization_lengthscale"});
    Increment scales(geom, vars, date);
    scales.zero();

    // Pull the increment into an atlas FieldSet so we can fill it point by point
    atlas::FieldSet fset;
    scales.toFieldSet(fset);

    // Compute the length scale at every (non-ghost) grid point, using the
    // shared formula in soca::ObsLocRossby (the single source of truth).
    for (auto & field : fset) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (GeometryIterator it = geom.begin(); it != geom.end(); ++it) {
        double lengthscale = ObsLocRossby::lengthScale(
            options, it.getFieldValue("rossby_radius"), it.getFieldValue("area"));
        if (gcWidth) lengthscale *= 2.0 / std::sqrt(0.3);

        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(it.i(), jlevel) = lengthscale;
        }
      }
    }

    // Push the filled FieldSet back into the increment and write it out
    scales.fromFieldSet(fset);

    const eckit::LocalConfiguration outputConfig(fullConfig, "output");
    scales.write(outputConfig);
    oops::Log::info() << "Output localization length scales: " << scales << std::endl;

    return 0;
  }
  // ---------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "soca::SaveObsLocScales";
  }
  // ---------------------------------------------------------------------------
};

}  // namespace soca
