/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "soca/Utils/readNcAndInterp.h"
#include "soca/SaberBlocks/ParametricOceanStdDev/ParametricOceanStdDev.h"

namespace soca {

// ----------------------------------------------------------------------------------------

static saber::SaberOuterBlockMaker<ParametricOceanStdDev>
  makerParametricOceanStdDev_("SOCAParametricOceanStdDev");

// ----------------------------------------------------------------------------------------

ParametricOceanStdDev::ParametricOceanStdDev(
    const oops::GeometryData & geom,
    const oops::Variables & outerVars,
    const eckit::Configuration & covarConfig,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberOuterBlockBase(params, xb.validTime()),
    geom_(geom),
    innerVars_(outerVars)
{
  util::Timer timer("soca::ParametricOceanStdDev", "ParametricOceanStdDev");

  // validate the configuration, using a hacky way
  {
    auto fullConfig = params.fullConfig.value();
    Parameters params2;
    params2.validate(fullConfig);
  }

  // things we'll need for later
  const auto & fs = geom_.functionSpace();
  const auto & depth = xb[params.depthVariable.value()];
  const auto & v_depth = atlas::array::make_view<double, 2>(depth);
  const auto & v_ghost = atlas::array::make_view<int, 1>(fs.ghost());
  const auto & v_mask = atlas::array::make_view<double, 2>(
    geom_.getField(params.maskVariable.value()));
  const int levels = depth.shape(1);

  // diagnostic fields that will be optionally saved to a file at the end
  atlas::FieldSet diags;
  diags.add(depth);

  // setup the smoother. If the user provides a smoother, use that, otherwise
  // create a smoother based on the default parameters for hz and vt smoothing.
  oops::Log::info() << "ParametricOceanStdDev: using user-provided smoother" << std::endl;
  eckit::LocalConfiguration smootherConfig;
  if (params.smoother.value()) {
    smootherConfig = params.smoother.value()->toConfiguration();
  } else {
    smootherConfig.set("horizontal", eckit::LocalConfiguration());
    smootherConfig.set("vertical", eckit::LocalConfiguration());
  }
  OceanSmoother::Parameters smootherParams;
  smootherParams.deserialize(smootherConfig);
  OceanSmoother smoother3D(geom_, smootherParams, levels, xb.fieldSet());

  // --------------------------------------------------------------------------------------------
  // create temperature error.
  // This is done by calculating the temperature gradient, and then using that
  // to calculate the temperature error as a function of the e-folding scale,
  // and min/max values. The min value at the surface is read in from a file
  // The temperature error is then smoothed.
  const auto & tocnParams = params.tocn.value();
  if (xb.has(tocnParams.varName.value())) {
    oops::Log::info() << "ParametricOceanStdDev: creating temperature error" << std::endl;

    const auto & tocn = xb[params.tocn.value().varName.value()];
    const auto & v_tocn = atlas::array::make_view<double, 2>(tocn);
    const double minVal = tocnParams.min.value();
    const double maxVal = tocnParams.max.value();
    const double efold = tocnParams.efold.value();
    const double minLayerThickness = tocnParams.minLayerThickness.value();

    // create empty error field
    atlas::Field tocnErr = tocn.clone();
    auto v_tocnErr = atlas::array::make_view<double, 2>(tocnErr);
    v_tocnErr.assign(0.0);
    bkgErr_.add(tocnErr);
    diags.add(tocnErr);

    // get the sst values
    atlas::Field sstErr;
    const auto & sstParams = tocnParams.sst.value();
    if (sstParams.has("fixed value") == sstParams.has("filepath")) {
      // sanity check to make sure only one option is specified
      throw eckit::BadParameter("ParametricOceanStdDev: tocn parameter must have either"
                                " 'fixed value' or 'filepath'", Here());
    } else if (sstParams.has("fixed value")) {
      // create a fixed value for the sst
      sstErr = fs.createField<double>(atlas::option::levels(1));
      auto v_sst = atlas::array::make_view<double, 2>(sstErr);
      v_sst.assign(sstParams.getDouble("fixed value"));
    } else if (sstParams.has("filepath")) {
      // read the sst field from a file and interpolate it to the destination function space
      const std::string & sstVar = sstParams.getString("variable");
      auto fset = readNcAndInterp(sstParams.getString("filepath"),
                                  {sstVar}, fs);
      sstErr = fset[sstVar];
    }
    sstErr.rename("min sst");
    diags.add(sstErr);
    const auto & v_sstErr = atlas::array::make_view<double, 2>(sstErr);

    // initialize empty dtdz and efold fields
    atlas::Field dtdz = fs.createField<double>(atlas::option::levels(levels) |
                                               atlas::option::name("dtdz"));
    auto v_dtdz = atlas::array::make_view<double, 2>(dtdz);
    v_dtdz.assign(0.0);
    diags.add(dtdz);

    // iterate over all points
    for (size_t i = 0; i < tocnErr.shape(0); i++) {
      if (v_ghost(i)) continue;  // skip ghost points
      if (v_mask(i, 0) == 0) continue;  // skip land points

      // calculate dt/dz
      for (size_t z = 1; z < levels-1; z++) {
        v_dtdz(i, z) = (v_tocn(i, z+1) - v_tocn(i, z-1)) /
                       (v_depth(i, z+1) - v_depth(i, z-1));

        // ignore dt/dz where layers are too thin
        if (v_depth(i, z+1) - v_depth(i, z) <= minLayerThickness) {
          v_dtdz(i, z) = 0.0;
        }
      }
      v_dtdz(i, 0) = v_dtdz(i, 1);
      v_dtdz(i, levels-1) = 0;

      // calculate the temperature error as a function of dt/dz, e-folding scale,
      // and min/max values
      double sstVal = v_sstErr(i, 0);
      for (size_t z = 0; z < levels; z++) {
        // step 1: calc value from dT/dz
        auto val = abs(tocnParams.dz.value() * v_dtdz(i, z));

        // step 2: calc a min from the e-folding scale, and min SST value
        auto localMin = minVal + (sstVal - minVal) * exp(-v_depth(i, z) / efold);

        // step 3: min/max
        v_tocnErr(i, z) = std::clamp(val, localMin, maxVal);;
      }
    }
    sstErr.set_dirty();

    // smooth the temperature error
    if (tocnParams.smooth.value()) {
      smoother3D.multiply(tocnErr);
    }
  }

  // --------------------------------------------------------------------------------------------
  // create unbalanced salinity error
  const auto & socnParams = params.socn.value();
  if (xb.has(socnParams.varName.value())) {
    oops::Log::info() << "ParametricOceanStdDev: creating salinity error" << std::endl;

    const double min = socnParams.min.value();
    const double max = socnParams.max.value();
    const double maxMld = socnParams.mldMax.value();

    // create empty error field
    atlas::Field socnErr = fs.createField<double>(atlas::option::levels(levels) |
                                                  atlas::option::name(socnParams.varName.value()));
    auto v_socnErr = atlas::array::make_view<double, 2>(socnErr);
    v_socnErr.assign(0.0);
    bkgErr_.add(socnErr);
    diags.add(socnErr);

    // find the MLD, and smooth it if requested
    atlas::Field mld = xb[socnParams.mldVariableName.value()].clone(); mld.rename("mld");
    const auto & v_mld = atlas::array::make_view<double, 2>(mld);
    diags.add(mld);

    // iterate over each grid point
    for (size_t i = 0; i < socnErr.shape(0); i++) {
      if (v_ghost(i)) continue;  // skip ghost points
      if (v_mask(i, 0) == 0) continue;  // skip land points

      // calculate the salinity error as a function of the mixed layer depth
      // and min/max values
      auto mld = std::min(v_mld(i, 0), maxMld);
      for (size_t z = 0; z < levels; z++) {
        if (v_depth(i, z) <= mld) {
          v_socnErr(i, z) = max;
        } else {
          // NOTE, this isn't quite right, it's not taking into account the min val
          const double r = 0.1 + 0.45 * (1.0-tanh(2.0 * log(v_depth(i, z) / mld)));
          v_socnErr(i, z) = r * max;
        }
      }
    }
    socnErr.set_dirty();

    // smooth the unbalanced salinity error
    if (socnParams.smooth.value()) {
      smoother3D.multiply(socnErr);
    }
  }

  // --------------------------------------------------------------------------------------------
  // create unbalanced SSH error
  // This is done by setting the SSH error to a maximum value in the
  // extra-tropics poles, and decreasing exponentially toward 0.0 the equator
  const auto & sshParams = params.ssh.value();
  if (xb.has(sshParams.varName.value())) {
    oops::Log::info() << "ParametricOceanStdDev: creating SSH error" << std::endl;

    const double min = sshParams.min.value();
    const double max = sshParams.max.value();
    const double phiEx = sshParams.phiEx.value();
    const auto & v_lonlat = atlas::array::make_view<double, 2>(fs.lonlat());
    const double pi = 2.0 * acos(0.0);

    // create empty error field
    atlas::Field sshErr = fs.createField<double>(atlas::option::levels(1) |
                                                 atlas::option::name(sshParams.varName.value()));
    auto v_sshErr = atlas::array::make_view<double, 2>(sshErr);
    v_sshErr.assign(0.0);
    bkgErr_.add(sshErr);
    diags.add(sshErr);

    // find the SSH error
    for (size_t i = 0 ; i < sshErr.shape(0); i++) {
      if (v_ghost(i)) continue;  // skip ghost points
      if (v_mask(i, 0) == 0) continue;  // skip land points

      // calculate the SSH error as a function of the exponential scale
      auto absLat = std::abs(v_lonlat(i, 1));
      if (absLat >= phiEx) {
        // in the extra-tropics, set to max value
        v_sshErr(i, 0) = max;
      } else {
        // otherwise, set to a value that increases exponentially toward equator
        v_sshErr(i, 0) = min + 0.5 * (max - min) * (1.0 - cos(pi * absLat / phiEx));
      }
    }
    sshErr.set_dirty();
  }

  // --------------------------------------------------------------------------------------------
  // deal with other generic variables
  // Which only have a min/max, and a fraction of the background
  if (params.otherVars.value()) {
    for (const auto & varParams : *params.otherVars.value()) {
      std::string varName = varParams.varName.value();
      if (!xb.has(varName)) continue;  // skip if the variable isn't
                                       // in the background (give warning?)

      oops::Log::info() << "ParametricOceanStdDev: creating "+varName+" error" << std::endl;
      const double min = varParams.min.value();
      const double max = varParams.max.value();
      const double fraction = varParams.fractionOfBkg.value();

      const auto & varBg = xb[varName];
      const auto & v_varBg = atlas::array::make_view<double, 2>(varBg);
      const auto levels = varBg.shape(1);

      auto varErr = varBg.clone();
      auto v_varErr = atlas::array::make_view<double, 2>(varErr);

      bkgErr_.add(varErr);
      diags.add(varErr);

      // iterate over each grid point
      for (size_t i = 0; i < varErr.shape(0); i++) {
        if (v_ghost(i)) continue;  // skip ghost points
        if (v_mask(i, 0) == 0) continue;  // skip land points

        // calculate the error as a fraction of the background
        for (size_t z = 0; z < levels; z++) {
          v_varErr(i, z) = std::clamp(std::abs(fraction * v_varBg(i, z)), min, max);
        }
      }
      varErr.set_dirty();

      // smooth the error
      if (varParams.smooth.value()) {
        smoother3D.multiply(varErr);
      }
    }
  }

  // save the diagnostics if requested
  if (params.saveDiags.value()) {
    oops::Log::info() << "ParametricOceanStdDev: saving diagnostics to file " << std::endl;
    util::writeFieldSet(geom_.comm(), *params.saveDiags.value(), diags);
  }
}

// ----------------------------------------------------------------------------------------
// The TL and AD multiply are the same, so we can use a common function for
// both. We could just have multiplyAD() call multiply(), but we want to keep
// the timings separate.
void ParametricOceanStdDev::commonMultiply(oops::FieldSet3D &fset) const {
  const auto & v_ghost = atlas::array::make_view<int, 1>(geom_.functionSpace().ghost());

  for (auto & field : fset) {
    if (!bkgErr_.has(field.name())) continue;

    const auto levels = field.shape(1);
    auto v_field = atlas::array::make_view<double, 2>(field);
    const auto & v_err  = atlas::array::make_view<double, 2>(bkgErr_[field.name()]);
    for (size_t i = 0; i < field.shape(0); i++) {
      if (v_ghost(i)) continue;  // skip ghost points

      for (size_t z = 0; z < levels; z++) {
        v_field(i, z) *= v_err(i, z);
      }
    }
    field.set_dirty();
  }
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiply(oops::FieldSet3D &fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::multiply starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "multiply");

  commonMultiply(fset);

  oops::Log::trace() << "ParametricOceanStdDev::multiply done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::multiplyAD starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "multiplyAD");

  commonMultiply(fset);

  oops::Log::trace() << "ParametricOceanStdDev::multiplyAD done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::leftInverseMultiply(oops::FieldSet3D &fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::leftInverseMultiply starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "leftInverseMultiply");

  const auto & v_ghost = atlas::array::make_view<int, 1>(geom_.functionSpace().ghost());

  for (auto & field : fset) {
    if (!bkgErr_.has(field.name())) continue;

    const auto levels = field.shape(1);
    auto v_field = atlas::array::make_view<double, 2>(field);
    const auto & v_err  = atlas::array::make_view<double, 2>(bkgErr_[field.name()]);
    for (size_t i = 0; i < field.shape(0); i++) {
      if (v_ghost(i)) continue;  // skip ghost points

      for (size_t z = 0; z < levels; z++) {
        if (v_err(i, z) == 0.0) {
          v_field(i, z) = 0.0;
        } else {
          v_field(i, z) /= v_err(i, z);
        }
      }
    }
    field.set_dirty();
  }

  oops::Log::trace() << "ParametricOceanStdDev::leftInverseMultiply done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::print(std::ostream & os) const {
  os << "soca::ParametricOceanStdDev";
}

// ----------------------------------------------------------------------------------------

}  // namespace soca
