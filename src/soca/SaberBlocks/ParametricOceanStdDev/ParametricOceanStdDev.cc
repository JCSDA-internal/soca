/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "soca/SaberBlocks/ParametricOceanStdDev/ParametricOceanStdDev.h"
#include "soca/SaberBlocks/Util/OceanSmoother.h"
#include "soca/Utils/readAndInterp.h"

#include "eckit/exception/Exceptions.h"

#include "oops/generic/AtlasInterpolator.h"

#include "saber/blocks/SaberOuterBlockBase.h"

namespace soca {

static saber::SaberOuterBlockMaker<ParametricOceanStdDev> makerParametricOceanStdDev_(
  "ParametricOceanStdDev");

// ------------------------------------------------------------------------------------------------

ParametricOceanStdDev::ParametricOceanStdDev(
    const oops::GeometryData & geom,
    const oops::Variables & outerVars,
    const eckit::Configuration & covarConfig,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(geom),
    innerVars_(outerVars),
    bkgErr_(xb.validTime(),xb.commGeom())
{
  const int levels = xb["hocn"].levels();
  const double MIN_LAYER_THICKNESS = 0.1;

  xb.fieldSet().haloExchange();
  auto v_lonlat = atlas::array::make_view<double, 2>(geom.functionSpace().lonlat());
  auto v_mask = atlas::array::make_view<double, 2>(geom.getField("interp_mask"));
  auto v_hocn = atlas::array::make_view<double, 2>(xb["hocn"]);
  auto v_tocn = atlas::array::make_view<double, 2>(xb["tocn"]);
  auto v_depth = atlas::array::make_view<double, 2>(xb["layer_depth"]);

  // initialize the smoother
  std::unique_ptr<OceanSmoother> smoother;
  if (params.smoother.value() != boost::none) {
    smoother.reset(new OceanSmoother(innerGeometryData_, *params.smoother.value(), levels));
  }

  //*************************************************************************************
  // calculate T background error
  //*************************************************************************************
  atlas::Field tocn_err = geom.functionSpace().createField<double>(
    atlas::option::levels(levels) | atlas::option::name("tocn"));
  bkgErr_.add(tocn_err);
  {
    oops::Log::info() << "ParametricOceanStdDev: Creating temperature background error:" << std::endl;

    // parameters from config
    auto minVal = params.tocn.value().min.value();
    auto maxVal = params.tocn.value().max.value();
    auto efold = params.tocn.value().efold.value();

    // output field
    auto v_tocn_err = atlas::array::make_view<double, 2>(tocn_err);
    v_tocn_err.assign(0.0);

    // Get SST values interpolated
    const std::string & sstFilepath = params.tocn.value().sst.value().getString("filepath") ;
    const std::string & sstVar = params.tocn.value().sst.value().getString("variable") ;
    auto sstFields = readNcAndInterp(sstFilepath, {sstVar}, geom.functionSpace());
    auto & sstErr = sstFields[sstVar];
    sstErr.rename("sst_surf");
    bkgErr_.add(sstErr);
    const auto & v_sstErr = atlas::array::make_view<double, 2>(sstErr);

    // loop over all points
    std::vector<double> dtdz(levels, 0.0);
    for (size_t i = 0; i < tocn_err.shape(0); i++) {
      // skip land
      if (v_mask(i, 0) == 0.0) continue;

      // calculate dt/dz
      dtdz[0] = 2.0 * (v_tocn(i,1) - v_tocn(i,0)) / (v_hocn(i,0) + v_hocn(i,1));
      for (size_t z = 1; z < levels-1; z++) {
        dtdz[z] = (v_tocn(i, z+1) - v_tocn(i, z-1)) /
                  (v_hocn(i,z) + 0.5*(v_hocn(i, z+1) + v_hocn(i, z-1)));

        // ignore dt/dz where layers are too thin
        if(v_hocn(i,z) <= MIN_LAYER_THICKNESS) dtdz[z] = 0.0;
      }
      dtdz[levels-1] = 0;

      // calculate value as function of dt/dz, efolding scale, and min/max
      double sstVal = v_sstErr(i, 0);
      for (size_t z = 0; z < levels; z++) {

        // step 1: calc value from dT/dz
        auto val = abs(params.tocn.value().dz * dtdz[z]);

        // step 2: calc a minimum from the efolding scale, and min SST value
        auto localMin = minVal + (sstVal - minVal) * exp(-v_depth(i,z) / efold );

        // step 3: min/max
        val = std::clamp(val, localMin,  maxVal);

        // done
        v_tocn_err(i, z) = val;
      }
    }
  }

  // smooth the resulting field
  if(smoother) smoother->multiply(tocn_err);

  //*************************************************************************************
  // calculate unbalanced S background error
  //*************************************************************************************
  auto findMLD=[this, &geom, &v_mask, &v_depth, &xb](){
    const atlas::Field &tocn =  xb["tocn"];
    auto v_tocn = atlas::array::make_view<double, 2>(tocn);

    atlas::Field mld = geom.functionSpace().createField<double>(
      atlas::option::name("mld_T02"));
    auto v_mld = atlas::array::make_view<double, 1>(mld);
    v_mld.assign(0.0);

    const double tDelta = 0.2;
    for (size_t i = 0; i < mld.shape(0); i++) {
      if (v_mask(i, 0) == 0.0) continue;

      v_mld(i) = v_depth(i,0);

      const double t_surf = v_tocn(i,0);
      for (size_t lvl = 1; lvl < tocn.shape(1); lvl++) {
        if (abs(t_surf - v_tocn(i, lvl)) > tDelta ) {
          const double t2 = v_tocn(i, lvl);
          const double t1 = v_tocn(i, lvl-1);
          const double d2 = v_depth(i, lvl);
          const double d1 = v_depth(i, lvl-1);

          const double td2 = t2 < t1 ? tDelta : -tDelta;
          //const double r = (t_surf-tDelta - t1) / (t2-t1);
          const double r = (t_surf-td2 - t1) / (t2-t1);
          v_mld(i) = v_depth(i, lvl-1)   +  r * (v_depth(i, lvl) - v_depth(i, lvl-1));
          break;
        }
      }
    }
    return mld;
  };

  std::cout << "DBG " << "creating S error"<<std::endl;
  atlas::Field socn_err = geom.functionSpace().createField<double>(
    atlas::option::levels(levels) | atlas::option::name("socn"));
  bkgErr_.add(socn_err);
  {
    auto mld = findMLD();
    auto v_mld = atlas::array::make_view<double, 1> (mld);
    auto v_socn = atlas::array::make_view<double, 2> (socn_err);
    v_socn.assign(0.0);

    const double maxS = 0.25;
    const double minS = 0.05;
    const double maxMld = 400;
    for (size_t i = 0; i < socn_err.shape(0); i++) {
      //v_socn(i, 0) = v_mld(i);
      if (v_mask(i, 0) == 0.0) continue;
       for (size_t lvl = 0; lvl < socn_err.shape(1); lvl++) {
        const double mldVal = std::min(v_mld(i), maxMld);
        if( v_depth(i, lvl) <= mldVal) {
          v_socn(i, lvl) = maxS;
        } else {
          const double r = 0.1 + 0.45 *(1-tanh(2 * log(v_depth(i, lvl) / mldVal)));
          v_socn(i, lvl) = r * maxS;
        }
      }
    }
  }
  // TODO only do hz smoothing
  if(smoother) smoother->multiply(socn_err);


  //*************************************************************************************
  // calculate unbalanced SSH background error
  //*************************************************************************************
  std::cout << "DBG " << "creating SSH error"<<std::endl;
  atlas::Field ssh_err = geom.functionSpace().createField<double>(
    atlas::option::levels(1) |
    atlas::option::name("ssh"));
  bkgErr_.add(ssh_err);
  {
    // parameter from input config
    auto minVal = params.ssh.value().min.value();
    auto maxVal = params.ssh.value().max.value();

    // output field
    auto v_ssh_err = atlas::array::make_view<double, 2>(ssh_err);
    v_ssh_err.assign(0.0);

    // loop over all points
    for (size_t i = 0; i < ssh_err.shape(0); i++) {
      // skip land
      if (v_mask(i, 0) == 0.0) continue;

      auto absLat = std::abs(v_lonlat(i, 1));
      if (absLat >= params.ssh.value().phi_ex.value()) {
        // In the extratropics, set to max value
        v_ssh_err(i, 0) = maxVal;
      } else {
        // otherwise, taper to min value towards the equator
        constexpr double pi = 3.14159265358979323846;
        v_ssh_err(i, 0) = minVal + 0.5 * (maxVal - minVal)*(1 - std::cos(pi * absLat / params.ssh.value().phi_ex.value()));
      }
    }
  }

  // **********************************************************************************************
  // Other fields. Which only have min/max and "fraction of background" for
  // defining what the parametric stddev should be set to
  // **********************************************************************************************
  if (params.others.value() != boost::none) {
    oops::Log::info() << "Inventing parametric standard deviation for other variables:" << std::endl;
    const auto & otherVarsConfigs = *params.others.value();
    for(auto varConfig : otherVarsConfigs) {
      oops::Log::info() << "  " << varConfig.first << ":"<< std::endl;
      double frac = varConfig.second.fractionOfBkg.value();
      double minVal = varConfig.second.min.value();
      double maxVal = varConfig.second.max.value();
      oops::Log::info() << "    fraction of background: " << frac << std::endl;
      oops::Log::info() << "    min: " << minVal << std::endl;
      oops::Log::info() << "    max: " << maxVal << std::endl;

      atlas::Field bkg = xb[varConfig.first];
      auto v_bkg = atlas::array::make_view<double, 2>(bkg);

      // create the error field
      atlas::Field err = geom.functionSpace().createField<double>(
        atlas::option::levels(bkg.shape(1)) | atlas::option::name(varConfig.first));
      auto v_err = atlas::array::make_view<double, 2>(err);
      bkgErr_.add(err);

      // assign values
      v_err.assign(0.0);
      for (size_t i = 0; i < err.shape(0); i++) {
        if (v_mask(i, 0) == 0.0) continue;
        for (size_t lvl = 0; lvl < err.shape(1); lvl++) {
          v_err(i, lvl) = std::clamp(frac * v_bkg(i, lvl), minVal, maxVal);
        }
      }

      // all done, optionally smooth the field
      if(smoother) smoother->multiply(err);
    }
  }

  // optionally write the output
  std::cout << "DBG " << "writing output"<<std::endl;
  if (params.diagnostics.value() != boost::none) {
    oops::Log::info() << "ParametricStdDev values have been generated. Writing final"
                          " values to diagnostic file." << std::endl;
    bkgErr_.write(*params.diagnostics.value());
  }

}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiply(oops::FieldSet3D & fset) const {

  for (atlas::Field field : fset) {
    if (! bkgErr_.has(field.name())) continue;

    auto v_field = atlas::array::make_view<double, 2>(field);
    auto v_mult = atlas::array::make_view<double, 2>(bkgErr_[field.name()]);

    for (size_t i = 0; i < field.shape(0); i++) {
      for (size_t z = 0; z < field.levels(); z++) {
        // TEMP, change this to *=
        v_field(i, z) = v_mult(i, z);
      }
    }
  }
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiplyAD(oops::FieldSet3D & fset) const
{
  multiply(fset);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::leftInverseMultiply(oops::FieldSet3D &) const
{
  ASSERT(1 == 3);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::print(std::ostream & os) const {
  os << classname();
}

// ------------------------------------------------------------------------------------------------
}  // namespace soca
