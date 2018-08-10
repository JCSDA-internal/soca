/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SOCA_SRC_LINEARMODEL_TLMID_H_
#define SOCA_SRC_LINEARMODEL_TLMID_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "oops/interface/LinearModelBase.h"

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "src/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
// -----------------------------------------------------------------------------
///  linear identity model definition.
/*!
 *   linear identity model definition and configuration parameters.
 */

class TlmId: public oops::LinearModelBase<Traits>,
              private util::ObjectCounter<TlmId> {
 public:
  static const std::string classname() {return "soca::TlmId";}

  TlmId(const Geometry &, const eckit::Configuration &);
  ~TlmId();

/// Model trajectory computation
  void setTrajectory(const State &, State &, const ModelBias &) override;

/// Run TLM and its adjoint
  void initializeTL(Increment &) const override;
  void stepTL(Increment &, const ModelBiasIncrement &) const override;
  void finalizeTL(Increment &) const override;

  void initializeAD(Increment &) const override;
  void stepAD(Increment &, ModelBiasIncrement &) const override;
  void finalizeAD(Increment &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const Geometry & resolution() const {return resol_;}

 private:
  void print(std::ostream &) const override;

// Data
  int keyConfig_;
  util::Duration tstep_;
  const Geometry resol_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_SRC_LINEARMODEL_TLMID_H_
