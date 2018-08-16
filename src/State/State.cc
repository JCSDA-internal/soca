/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/State/State.h"

#include <algorithm>
#include <string>
#include <vector>
#include <utility>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "src/ModelBias.h"
#include "src/Fields/Fields.h"
#include "src/Geometry/Geometry.h"
#include "src/Increment/Increment.h"
#include "src/Model/Model.h"
#include "src/GetValuesTraj/GetValuesTraj.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"


using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  State::State(const Geometry & resol, const oops::Variables & vars,
               const util::DateTime & vt)
    : fields_(new Fields(resol, vars, vt)), stash_()
  {
    Log::trace() << "State::State created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & resol, const oops::Variables & vars,
	       const eckit::Configuration & file)
    : fields_(new Fields(resol, vars, util::DateTime())), stash_()
  {
    const std::vector<std::string> vv = {
          "cicen",
          "hicen",
          "hsnon",
          "tsfcn",
          "qsnon",
          "sicnk",
          "qicnk",
          "socn",
          "tocn",
          "ssh",
          "hocn"
    };
    //oops::Variables vars(vv);
    //fields_.reset(new Fields(resol, vars, util::DateTime()));
    fields_->read(file);
    ASSERT(fields_);
    Log::trace() << "State::State created and read in." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & resol, const State & other)
    : fields_(new Fields(*other.fields_, resol)), stash_()
  {
    ASSERT(fields_);
    Log::trace() << "State::State created by interpolation." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const State & other)
    : fields_(new Fields(*other.fields_)), stash_()
  {
    ASSERT(fields_);
    Log::trace() << "State::State copied." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::~State() {
    Log::trace() << "State::State destructed." << std::endl;
  }
  // -----------------------------------------------------------------------------
  void State::activateModel() {
    const std::vector<std::string> vv{
          "cicen",
          "hicen",
          "hsnon",
          "tsfcn",
          "qsnon",
          "sicnk",
          "qicnk",
          "socn",
          "tocn",
          "ssh",
          "hocn"
    };
    oops::Variables vars(vv);

    stash_.reset(new Fields(*fields_, vars));
    swap(fields_, stash_);
    ASSERT(fields_);
    ASSERT(stash_);
    Log::trace() << "State activated for Model" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void State::deactivateModel() {
    swap(fields_, stash_);
    *fields_ = *stash_;
    stash_.reset();
    ASSERT(fields_);
    ASSERT(!stash_);
    Log::trace() << "State deactivated for Model" << std::endl;
  }
  // -----------------------------------------------------------------------------
  /// Basic operators
  // -----------------------------------------------------------------------------
  State & State::operator=(const State & rhs) {
    ASSERT(fields_);
    *fields_ = *rhs.fields_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  /// Interpolate to observation location
  // -----------------------------------------------------------------------------
  void State::getValues(const ioda::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & cols) const {
    fields_->getValues(locs, vars, cols);
  }
  // -----------------------------------------------------------------------------
  void State::getValues(const ioda::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & cols,
                        GetValuesTraj &) const {
    fields_->getValues(locs, vars, cols);
  }
  // -----------------------------------------------------------------------------
  /// Interactions with Increments
  // -----------------------------------------------------------------------------
  State & State::operator+=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    ASSERT(fields_);
    fields_->add(dx.fields());
    return *this;
  }
  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------
  void State::read(const eckit::Configuration & files) {
    fields_->read(files);
  }
  // -----------------------------------------------------------------------------
  void State::write(const eckit::Configuration & files) const {
    fields_->write(files);
  }
  // -----------------------------------------------------------------------------
  void State::print(std::ostream & os) const {
    os << std::endl << "  Valid time: " << validTime();
    os << *fields_;
  }
  // -----------------------------------------------------------------------------
  /// For accumulator
  // -----------------------------------------------------------------------------
  void State::zero() {
    fields_->zero();
  }
  // -----------------------------------------------------------------------------
  void State::accumul(const double & zz, const State & xx) {
    fields_->axpy(zz, *xx.fields_);
  }
  // -----------------------------------------------------------------------------

}  // namespace soca
