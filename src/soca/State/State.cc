/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/State/State.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "soca/ModelBias.h"
#include "soca/Fields/Fields.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/Model/Model.h"
#include "soca/GetValuesTraj/GetValuesTraj.h"

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
  {
    fields_.reset(new Fields(resol, vars, util::DateTime()));
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
  void State::getValues(const ufo::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & cols) const {
    if (fields_->geometry()->getAtmInit())
      {
        // Get atm geovals
        // The variables in vars that are also defined in soca will be
        // over-written in the interpolation call bellow
        getValuesFromFile(locs, vars, cols);
        Log::trace() << cols << std::endl;
        }
    // Get ocean geovals
    fields_->getValues(locs, vars, cols);
    Log::trace() << cols << std::endl;
  }
  // -----------------------------------------------------------------------------
  void State::getValues(const ufo::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & cols,
                        GetValuesTraj & traj) const {
    fields_->getValues(locs, vars, cols, traj);
  }
  // -----------------------------------------------------------------------------
  /// Read Interpolated GeoVaLs from file
  // -----------------------------------------------------------------------------
  void State::getValuesFromFile(const ufo::Locations & locs,
                                const oops::Variables & vars,
                                ufo::GeoVaLs & atmgom) const {
    // Get Atm configuration
    eckit::LocalConfiguration conf(fields_->geometry()->getAtmConf());

    // Get Time Bounds
    util::DateTime bgn = util::DateTime(conf.getString("notocean.date_begin"));
    util::DateTime end = util::DateTime(conf.getString("notocean.date_end"));

    // Create the Atmospheric Geometry in Observation Space
    eckit::LocalConfiguration confatmobs(conf, "notocean.ObsSpace");
    ioda::ObsSpace atmobs(confatmobs, fields_->geometry()->getComm(),
                                bgn, end);

    // Get GeoVaLs from file
    eckit::LocalConfiguration confatm(conf, "notocean");
    atmgom.read(confatm, atmobs);
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
    Log::trace() << "State::State read started." << std::endl;
    fields_->read(files);
    Log::trace() << "State::State read done." << std::endl;
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
