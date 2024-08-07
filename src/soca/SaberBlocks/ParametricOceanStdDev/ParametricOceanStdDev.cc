/*
 * (C) Copyright 2017-2024 UCAR
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

static saber::SaberOuterBlockMaker<ParametricOceanStdDev> makerParametricOceanStdDev_("SOCAParametricOceanStdDev");

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

  // TODO get rid of hardcoded variable names

  // things we'll need for later
  const auto & fs = geom_.functionSpace();
  const auto & v_ghost = atlas::array::make_view<int, 1>(fs.ghost());
  const auto & v_mask = atlas::array::make_view<double, 2>(
    geom_.getField(params.maskVariable.value()));
  const auto & v_hocn = atlas::array::make_view<double, 2>(xb["hocn"]);
  const int levels = xb["hocn"].shape(1);
  const double minLayerThickness = 0.01;
  atlas::FieldSet diags;

  // setup the smoother
  std::unique_ptr<OceanSmoother> smoother3D;
  if (params.smoother.value()) {
    oops::Log::info() << "ParametricOceanStdDev: using user-provided smoother" << std::endl;
    smoother3D = std::make_unique<OceanSmoother>(geom_, *params.smoother.value(), levels, xb.fieldSet());
  }

  // --------------------------------------------------------------------------------------------
  // create temperature error (TODO only if temperature is active)
  {
    oops::Log::info() << "ParametricOceanStdDev: creating temperature error" << std::endl;
    const auto & tocnParams = params.tocn.value();
    const auto & tocn = xb["tocn"];
    const auto & v_tocn = atlas::array::make_view<double, 2>(tocn);
    const auto & v_depth = atlas::array::make_view<double, 2>(xb["layer_depth"]);
    const double minVal = tocnParams.min.value();
    const double maxVal = tocnParams.max.value();
    const double efold = tocnParams.efold.value();

    // create empty error field
    atlas::Field tocn_err = tocn.clone();
    tocn_err.rename("temperature");
    auto v_tocn_err = atlas::array::make_view<double, 2>(tocn_err);
    v_tocn_err.assign(0.0);
    bkgErr_.add(tocn_err);
    diags.add(tocn_err);

    // get the sst values
    atlas::Field sstErr;
    const auto & sstParams = tocnParams.sst.value();
    if(sstParams.has("fixed value") == sstParams.has("filepath")) {
      //sanity check to make sure only one option is specified
      throw eckit::BadParameter("ParametricOceanStdDev: tocn parameter must have either"
                                " 'fixed value' or 'filepath'", Here());
    } else if (sstParams.has("fixed value")) {
      // create a fixed value for the sst
      sstErr = fs.createField<double>(atlas::option::levels(1) |
                                   atlas::option::name("temperature_surface"));
      auto v_sst = atlas::array::make_view<double, 2>(sstErr);
      v_sst.assign(sstParams.getDouble("fixed value"));
    } else if (sstParams.has("filepath")) {
      // read the sst field from a file and interpolate it to the destination function space
      const std::string & sstVar = sstParams.getString("variable");
      auto fset = readNcAndInterp(sstParams.getString("filepath"),
                                  {sstVar}, fs);
      sstErr = fset[sstVar];
      sstErr.rename("temperature_surface");
    }
    diags.add(sstErr);
    const auto & v_sstErr = atlas::array::make_view<double, 2>(sstErr);

    // initialize empty dtdz and efold fields
    atlas::Field dtdz = fs.createField<double>(atlas::option::levels(tocn.shape(1)) |
                                               atlas::option::name("dtdz"));
    auto v_dtdz = atlas::array::make_view<double, 2>(dtdz);
    v_dtdz.assign(0.0);
    // atlas::Field efold = fs.createField<double>(atlas::option::levels(tocn.shape(1)) |
    //                                             atlas::option::name("efold"));
    // auto v_dtdz = atlas::array::make_view<double, 2>(dtdz);
    // v_dtdz.assign(0.0);

    diags.add(dtdz);

    // iterate over all points
    for (size_t i = 0; i < tocn_err.shape(0); i++) {
      if (v_ghost(i)) continue;  // skip ghost points
      if (v_mask(i, 0) == 0) continue;  // skip land points

      // calculate dt/dz
      v_dtdz(i, 0) = 2.0 * (v_tocn(i,1) - v_tocn(i,0)) / (v_hocn(i,0) + v_hocn(i,1));
      for (size_t z = 1; z < levels-1; z++) {
        v_dtdz(i, z) = (v_tocn(i, z+1) - v_tocn(i, z-1)) /
                       (v_hocn(i,z) + 0.5*(v_hocn(i, z+1) + v_hocn(i, z-1)));

        // ignore dt/dz where layers are too thin
        if(v_hocn(i,z) <= minLayerThickness) v_dtdz(1, z) = 0.0;
      }
      v_dtdz(1, levels-1) = 0;

      // calculate the temperature error as a function of dt/dz, e-folding scale,
      // and min/max values
      double sstVal = v_sstErr(i, 0);
      for (size_t z = 0; z < levels; z++) {

        // step 1: calc value from dT/dz
        auto val = abs(tocnParams.dz.value() * v_dtdz(i, z));

        // step 2: calc a min from the e-folding scale, and min SST value
        auto localMin = minVal + (sstVal - minVal) * exp(-v_depth(i,z) / efold );

        // step 3: min/max
        v_tocn_err(i, z) = std::clamp(val, localMin, maxVal);;
      }
    }

    // smooth the temperature error
    if (smoother3D) {
      smoother3D->multiply(tocn_err);
    }
  }

  // save the diagnostics if requested
  if (params.saveDiags.value()) {
    oops::Log::info() << "ParametricOceanStdDev: saving diagnostics to file " << std::endl;
    util::writeFieldSet(geom_.comm(), *params.saveDiags.value(), diags);
  }
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiply(oops::FieldSet3D &fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::multiply starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "multiply");

  // TODO implement this

  oops::Log::trace() << "ParametricOceanStdDev::multiply done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::multiplyAD starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "multiplyAD");

  // TODO implement this

  oops::Log::trace() << "ParametricOceanStdDev::multiplyAD done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::leftInverseMultiply(oops::FieldSet3D &fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::leftInverseMultiply starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "leftInverseMultiply");

  // TODO implement this

  oops::Log::trace() << "ParametricOceanStdDev::leftInverseMultiply done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::print(std::ostream & os) const {
  os << "soca::ParametricOceanStdDev";
}

// ----------------------------------------------------------------------------------------

}  // namespace soca
