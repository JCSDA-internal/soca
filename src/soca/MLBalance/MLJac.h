/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "torch/torch.h"

#include "daml/IceEmul/IceEmul.h"

#include "soca/Geometry/Geometry.h"

namespace soca {
  class MLJac {
   private:
    daml::IceEmul iceEmulArctic_;
    daml::IceEmul iceEmulAntarctic_;

   public:
    MLJac(const eckit::Configuration & config,
          const oops::FieldSet3D & xb,
          atlas::FieldSet jacobian,
          const oops::GeometryData & GeometryData,
          const eckit::mpi::Comm & comm) :
      iceEmulArctic_(getConf(config, "arctic"), comm),
      iceEmulAntarctic_(getConf(config, "antarctic"), comm)
    {
      oops::Log::trace() << "In MLJack" << std::endl;
      // Geometry info
      const auto lonlat =
        atlas::array::make_view<double, 2>(GeometryData.functionSpace().lonlat());

      // Pointers to the background
      auto cicen = atlas::array::make_view<double, 2>(xb["cicen"]);
      auto hi = atlas::array::make_view<double, 2>(xb["hicen"]);
      auto hs = atlas::array::make_view<double, 2>(xb["hsnon"]);
      auto sst = atlas::array::make_view<double, 2>(xb["tocn"]);
      auto sss = atlas::array::make_view<double, 2>(xb["socn"]);
      auto sice = atlas::array::make_view<double, 2>(xb["sice"]);
      auto tair = atlas::array::make_view<double, 2>(xb["tair"]);
      auto tsfc = atlas::array::make_view<double, 2>(xb["tsfc"]);

      // Pointers to the Jacobian
      auto dcdsst = atlas::array::make_view<double, 2>(jacobian["dc/dsst"]);
      auto dcdsss = atlas::array::make_view<double, 2>(jacobian["dc/dsss"]);
      auto dcdhi = atlas::array::make_view<double, 2>(jacobian["dc/dhi"]);
      auto dcdhs = atlas::array::make_view<double, 2>(jacobian["dc/dhs"]);

      // Containerize the model's bkg in a torch tensor and compute Jacobian
      torch::Tensor pattern = torch::zeros({iceEmulArctic_.inputSize_});
      for (atlas::idx_t jnode = 0; jnode < xb["tocn"].shape(0); ++jnode) {
        // TODO(G): add tair, tsfc & ice salinity to the soca background
        pattern[0] = tair(jnode, 0);
        pattern[1] = tsfc(jnode, 0);
        pattern[2] = sst(jnode, 0);
        pattern[3] = sss(jnode, 0);
        pattern[4] = hs(jnode, 0);
        pattern[5] = hi(jnode, 0);
        pattern[6] = sice(jnode, 0);
        torch::Tensor dcdx = torch::zeros({iceEmulArctic_.inputSize_});
        if ( lonlat(jnode, 1) > 40.0 ) {
          dcdx = iceEmulArctic_.model_->jac(pattern);
        }
        if ( lonlat(jnode, 1) < -40.0 ) {
          dcdx = iceEmulAntarctic_.model_->jac(pattern);
        }
        dcdsst(jnode, 0) = dcdx[2].item<float>();
        dcdsss(jnode, 0) = dcdx[3].item<float>();
        dcdhs(jnode, 0) = dcdx[4].item<float>();
        dcdhi(jnode, 0) = dcdx[5].item<float>();
      }
      oops::Log::trace() << "In MLJack" << std::endl;
    }

    // Utility for initializer
    const eckit::LocalConfiguration getConf(const eckit::Configuration & conf,
                                            std::string str) {
      const eckit::LocalConfiguration localConf(conf, str);
      return localConf;
    }
  };
}  // namespace soca
