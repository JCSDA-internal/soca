/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_OBSLOCALIZATION_ROSSBYLOCFILTER_H_
#define SOCA_OBSLOCALIZATION_ROSSBYLOCFILTER_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace soca {

class RossbyLocFilterParameters : public ufo::FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(RossbyLocFilterParameters, ufo::FilterParametersBase)

 public:
  oops::RequiredParameter<float> multiplier{"multiplier", this};
  oops::RequiredParameter<float> minvalue{"minvalue", this};
  oops::RequiredParameter<float> maxvalue{"maxvalue", this};
};

class RossbyLocFilter : public ufo::FilterBase {
 public:
  typedef RossbyLocFilterParameters Parameters_;
  static const std::string classname() {return "soca::RossbyLocFilter";}

  RossbyLocFilter(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~RossbyLocFilter() {}

 private:
  void print(std::ostream &) const override {}
  void applyFilter(const std::vector<bool> &, const ufo::Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return ufo::QCflags::preQC;}

  Parameters_ parameters_;
};

}  // namespace soca

#endif  // SOCA_OBSLOCALIZATION_ROSSBYLOCFILTER_H_
