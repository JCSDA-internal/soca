/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "atlas/array.h"
#include "atlas/util/KDTree.h"
#include "eckit/config/Configuration.h"
#include "oops/util/Printable.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace soca {

class Geometry;
class State;

class PostProcessIce: public util::Printable {
 public:
  // ---------------------------------------------------------------------------
  /// @brief Parameters for rescaling CICE analysis in the ice pack
  class RescaleParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(RescaleParameters, oops::Parameters)
   public:
    oops::Parameter<bool> rescale{"rescale",
      "rescale analysis in the ice pack", false, this};
    oops::Parameter<float> min_hice{"min hice",
      "min ice thickness to trigger adjusting ice volume", 0.5, this};
    oops::Parameter<float> min_hsno{"min hsno",
      "min snow thickness to trigger adjusting snow volume", 0.1, this};
  };

  /// @brief Parameters for the redistribution of sea ice
  class RedistributionParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(RedistributionParameters, oops::Parameters)
   public:
    oops::Parameter<double> edge{"seaice edge",
      "Threshold for sea ice edge, used in shuffle and rescale prior options",
      0.15, this, {oops::minConstraint(0.0), oops::maxConstraint(1.0)}};
    oops::Parameter<bool> shuffle{"shuffle",
      "Option to shuffle sea ice in the marginal ice zone (where ice concentration"
      " < seaice edge)",
      true, this};
    oops::Parameter<RescaleParameters> rescale{"rescale prior",
      "Option to rescale sea ice in the ice pack zone (where ice concentration"
      " > seaice edge)", {}, this};
    oops::Parameter<bool> adjustSST{"update SST",
      "Option to update sea surface temperature", true, this};
    oops::Parameter<double> sstDiffMax{"max positive SST update",
      "Maximum update to sea surface temperature (K)", 1.0, this,
      {oops::minConstraint(0.0)}};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for adding soca increment to CICE restart files
  class Parameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Parameters, oops::Parameters)
   public:
    oops::RequiredParameter<int> ncat{"ncat", this, {oops::minConstraint(1)}};
    oops::RequiredParameter<int> ice_lev{"ice_lev", this, {oops::minConstraint(1)}};
    oops::RequiredParameter<int> sno_lev{"sno_lev", this, {oops::minConstraint(1)}};
    oops::RequiredParameter<RedistributionParameters> arctic{"arctic", this};
    oops::RequiredParameter<RedistributionParameters> antarctic{"antarctic", this};
    oops::Parameter<float> icepackTstep{"icepack time step",
            "icepack time step used for thickness categories rebinning", 300, this};
    oops::Parameter<int> shuffleStencilSize{"shuffle stencil depth",
            "stencil depth used in the shuffle search", 1, this,
            {oops::minConstraint(1)}};
    oops::Parameter<std::string> tFrzOption{"icepack tfrz_option",
            "icepack option to compute freezing temperature", "mushy", this};
    oops::Parameter<double> minAice{"min aice",
            "minimum allowable ice concentration", 0.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> minVice{"min vice",
            "minimum allowable ice volume", 0.00001, this, {oops::minConstraint(0.0)}};
  };

  const std::string classname() {return "soca::PostProcessIce";}

  PostProcessIce(const Geometry &, const eckit::Configuration &);

  void postProcess(State & pproc, const State & restart,
                   const State & analysis) const;

 private:
  const Geometry & geom_;
  Parameters params_;
  size_t ncat_;
  atlas::util::IndexKDTree kdTree_;

  atlas::array::ArrayView<double, 2> lonlat_;
  atlas::array::ArrayView<double, 2> mask_;

  // Helpers to compute totals and mean thicknesses from category fields
  double totalAice(const std::vector<atlas::array::ArrayView<double, 2>> & aiceCat,
                   size_t jnode) const;
  double meanHice(const std::vector<atlas::array::ArrayView<double, 2>> & viceCat,
                  double aice, size_t jnode) const;
  double meanHsno(const std::vector<atlas::array::ArrayView<double, 2>> & vsnoCat,
                  double aice, size_t jnode) const;
  void print(std::ostream &) const override;
};

}  // namespace soca
